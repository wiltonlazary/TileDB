/**
 * @file   global_order_writer.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2022 TileDB, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file implements class GlobalOrderWriter.
 */

#include "tiledb/sm/query/writers/global_order_writer.h"
#include "tiledb/common/common.h"
#include "tiledb/common/heap_memory.h"
#include "tiledb/common/logger.h"
#include "tiledb/sm/array/array.h"
#include "tiledb/sm/array_schema/array_schema.h"
#include "tiledb/sm/array_schema/dimension.h"
#include "tiledb/sm/fragment/fragment_metadata.h"
#include "tiledb/sm/misc/comparators.h"
#include "tiledb/sm/misc/hilbert.h"
#include "tiledb/sm/misc/parallel_functions.h"
#include "tiledb/sm/misc/tdb_math.h"
#include "tiledb/sm/misc/tdb_time.h"
#include "tiledb/sm/misc/utils.h"
#include "tiledb/sm/misc/uuid.h"
#include "tiledb/sm/query/hilbert_order.h"
#include "tiledb/sm/query/query_macros.h"
#include "tiledb/sm/stats/global_stats.h"
#include "tiledb/sm/storage_manager/storage_manager.h"
#include "tiledb/sm/tile/generic_tile_io.h"
#include "tiledb/sm/tile/tile_metadata_generator.h"
#include "tiledb/sm/tile/writer_tile.h"
#include "tiledb/storage_format/uri/generate_uri.h"

using namespace tiledb;
using namespace tiledb::common;
using namespace tiledb::sm::stats;

namespace tiledb {
namespace sm {

class GlobalOrderWriterStatusException : public StatusException {
 public:
  explicit GlobalOrderWriterStatusException(const std::string& message)
      : StatusException("GlobalOrderWriter", message) {
  }
};

/* ****************************** */
/*   CONSTRUCTORS & DESTRUCTORS   */
/* ****************************** */

GlobalOrderWriter::GlobalOrderWriter(
    stats::Stats* stats,
    shared_ptr<Logger> logger,
    StorageManager* storage_manager,
    Array* array,
    Config& config,
    std::unordered_map<std::string, QueryBuffer>& buffers,
    Subarray& subarray,
    Layout layout,
    uint64_t fragment_size,
    std::vector<WrittenFragmentInfo>& written_fragment_info,
    bool disable_checks_consolidation,
    std::vector<std::string>& processed_conditions,
    Query::CoordsInfo& coords_info,
    bool remote_query,
    optional<std::string> fragment_name,
    bool skip_checks_serialization)
    : WriterBase(
          stats,
          logger,
          storage_manager,
          array,
          config,
          buffers,
          subarray,
          layout,
          written_fragment_info,
          disable_checks_consolidation,
          coords_info,
          remote_query,
          fragment_name,
          skip_checks_serialization)
    , processed_conditions_(processed_conditions)
    , fragment_size_(fragment_size)
    , current_fragment_size_(0) {
  if (layout != Layout::GLOBAL_ORDER) {
    throw GlobalOrderWriterStatusException(
        "Failed to initialize global order writer. Layout " +
        layout_str(layout) + " is not global order.");
  }
}

GlobalOrderWriter::~GlobalOrderWriter() {
}

/* ****************************** */
/*               API              */
/* ****************************** */

Status GlobalOrderWriter::dowork() {
  get_dim_attr_stats();

  auto timer_se = stats_->start_timer("write");

  // In case the user has provided a coordinates buffer
  RETURN_NOT_OK(split_coords_buffer());

  if (check_coord_oob_) {
    RETURN_NOT_OK(check_coord_oob());
  }

  try {
    auto status = global_write();
    if (!status.ok()) {
      clean_up();
      return status;
    }
  } catch (...) {
    clean_up();
    std::throw_with_nested(std::runtime_error("[GlobalOrderWriter::dowork] "));
  }

  return Status::Ok();
}

Status GlobalOrderWriter::finalize() {
  auto timer_se = stats_->start_timer("finalize");

  if (global_write_state_ != nullptr) {
    try {
      throw_if_not_ok(finalize_global_write_state());
    } catch (...) {
      clean_up();
      std::throw_with_nested(
          std::runtime_error("[GlobalOrderWriter::dowork] "));
    }
  }

  return Status::Ok();
}

void GlobalOrderWriter::reset() {
  if (global_write_state_ != nullptr) {
    nuke_global_write_state();
  }
}

Status GlobalOrderWriter::alloc_global_write_state() {
  // Create global array state object
  if (global_write_state_ != nullptr)
    return logger_->status(
        Status_WriterError("Cannot initialize global write state; State not "
                           "properly finalized"));
  global_write_state_.reset(new GlobalWriteState);

  // Alloc FragmentMetadata object
  global_write_state_->frag_meta_ = make_shared<FragmentMetadata>(HERE());
  // Used in serialization when FragmentMetadata is built from ground up
  global_write_state_->frag_meta_->set_storage_manager(storage_manager_);

  return Status::Ok();
}

Status GlobalOrderWriter::init_global_write_state() {
  auto uri = global_write_state_->frag_meta_->fragment_uri();

  // Initialize global write state for attribute and coordinates
  for (const auto& it : buffers_) {
    // Initialize last tiles
    const auto& name = it.first;
    const auto var_size = array_schema_.var_size(name);
    const auto nullable = array_schema_.is_nullable(name);
    const auto cell_size = array_schema_.cell_size(name);
    const auto type = array_schema_.type(name);
    const auto& domain{array_schema_.domain()};
    const auto capacity = array_schema_.capacity();
    const auto cell_num_per_tile =
        coords_info_.has_coords_ ? capacity : domain.cell_num_per_tile();
    auto last_tile_vector =
        std::pair<std::string, WriterTileVector>(name, WriterTileVector());
    try {
      last_tile_vector.second.emplace_back(WriterTile(
          array_schema_,
          cell_num_per_tile,
          var_size,
          nullable,
          cell_size,
          type));
    } catch (const std::logic_error& le) {
      return Status_WriterError(le.what());
    }
    global_write_state_->last_tiles_.emplace(std::move(last_tile_vector));

    // Initialize cells written
    global_write_state_->cells_written_[name] = 0;

    // Initialize last var offsets
    global_write_state_->last_var_offsets_[name] = 0;
  }

  return Status::Ok();
}

GlobalOrderWriter::GlobalWriteState* GlobalOrderWriter::get_global_state() {
  return global_write_state_.get();
}

std::pair<Status, std::unordered_map<std::string, VFS::MultiPartUploadState>>
GlobalOrderWriter::multipart_upload_state(bool client) {
  if (client) {
    return {Status::Ok(), global_write_state_->multipart_upload_state_};
  }

  auto meta = global_write_state_->frag_meta_;
  std::unordered_map<std::string, VFS::MultiPartUploadState> result;

  // TODO: to be refactored, there are multiple places in writers where
  // we iterate over the internal fragment files manually
  const auto buf_names = buffer_names();
  for (const auto& name : buf_names) {
    auto&& [st1, uri] = meta->uri(name);
    RETURN_NOT_OK_TUPLE(st1, {});

    auto&& [st2, state] = storage_manager_->vfs()->multipart_upload_state(*uri);
    RETURN_NOT_OK_TUPLE(st2, {});
    // If there is no entry for this uri, probably multipart upload is disabled
    // or no write was issued so far
    if (!state.has_value()) {
      return {Status::Ok(), {}};
    }
    result[uri->remove_trailing_slash().last_path_part()] = std::move(*state);

    if (array_schema_.var_size(name)) {
      auto&& [status, var_uri] = meta->var_uri(name);
      RETURN_NOT_OK_TUPLE(status, {});

      auto&& [st, var_state] =
          storage_manager_->vfs()->multipart_upload_state(*var_uri);
      RETURN_NOT_OK_TUPLE(st, {});
      result[var_uri->remove_trailing_slash().last_path_part()] =
          std::move(*var_state);
    }
    if (array_schema_.is_nullable(name)) {
      auto&& [status, validity_uri] = meta->validity_uri(name);
      RETURN_NOT_OK_TUPLE(status, {});

      auto&& [st, val_state] =
          storage_manager_->vfs()->multipart_upload_state(*validity_uri);
      RETURN_NOT_OK_TUPLE(st, {});
      result[validity_uri->remove_trailing_slash().last_path_part()] =
          std::move(*val_state);
    }
  }

  return {Status::Ok(), result};
}

Status GlobalOrderWriter::set_multipart_upload_state(
    const std::string& uri,
    const VFS::MultiPartUploadState& state,
    bool client) {
  if (client) {
    global_write_state_->multipart_upload_state_[uri] = state;
    return Status::Ok();
  }

  // uri in this case holds only the buffer name
  auto absolute_uri =
      global_write_state_->frag_meta_->fragment_uri().join_path(uri);
  return storage_manager_->vfs()->set_multipart_upload_state(
      absolute_uri, state);
}

/* ****************************** */
/*        PRIVATE METHODS         */
/* ****************************** */

Status GlobalOrderWriter::check_coord_dups() const {
  auto timer_se = stats_->start_timer("check_coord_dups");

  // Check if applicable
  if (array_schema_.allows_dups() || !check_coord_dups_ || dedup_coords_)
    return Status::Ok();

  if (!coords_info_.has_coords_) {
    return logger_->status(
        Status_WriterError("Cannot check for coordinate duplicates; "
                           "Coordinates buffer not found"));
  }

  if (coords_info_.coords_num_ < 2)
    return Status::Ok();

  // Prepare auxiliary vectors for better performance
  auto dim_num = array_schema_.dim_num();
  std::vector<const unsigned char*> buffs(dim_num);
  std::vector<uint64_t> coord_sizes(dim_num);
  std::vector<const unsigned char*> buffs_var(dim_num);
  std::vector<uint64_t*> buffs_var_sizes(dim_num);
  for (unsigned d = 0; d < dim_num; ++d) {
    const auto& dim_name{array_schema_.dimension_ptr(d)->name()};
    buffs[d] = (const unsigned char*)buffers_.find(dim_name)->second.buffer_;
    coord_sizes[d] = array_schema_.cell_size(dim_name);
    buffs_var[d] =
        (const unsigned char*)buffers_.find(dim_name)->second.buffer_var_;
    buffs_var_sizes[d] = buffers_.find(dim_name)->second.buffer_var_size_;
  }

  auto status = parallel_for(
      storage_manager_->compute_tp(),
      1,
      coords_info_.coords_num_,
      [&](uint64_t i) {
        // Check for duplicate in adjacent cells
        bool found_dup = true;
        for (unsigned d = 0; d < dim_num; ++d) {
          auto dim{array_schema_.dimension_ptr(d)};
          if (!dim->var_size()) {  // Fixed-sized dimensions
            if (memcmp(
                    buffs[d] + i * coord_sizes[d],
                    buffs[d] + (i - 1) * coord_sizes[d],
                    coord_sizes[d]) != 0) {  // Not the same
              found_dup = false;
              break;
            }
          } else {
            auto offs = (uint64_t*)buffs[d];
            auto off_i_plus_1 = (i == coords_info_.coords_num_ - 1) ?
                                    *(buffs_var_sizes[d]) :
                                    offs[i + 1];
            auto size_i_minus_1 = offs[i] - offs[i - 1];
            auto size_i = off_i_plus_1 - offs[i];

            // Compare sizes
            if (size_i != size_i_minus_1) {  // Not same
              found_dup = false;
              break;
            }

            // Compare var values
            if (memcmp(
                    buffs_var[d] + offs[i - 1],
                    buffs_var[d] + offs[i],
                    size_i) != 0) {  // Not the same
              found_dup = false;
              break;
            }
          }
        }

        // Found duplicate
        if (found_dup) {
          std::stringstream ss;
          ss << "Duplicate coordinates " << coords_to_str(i);
          ss << " are not allowed";
          return Status_WriterError(ss.str());
        }

        return Status::Ok();
      });

  RETURN_NOT_OK(status);

  return Status::Ok();
}

Status GlobalOrderWriter::check_global_order() const {
  auto timer_se = stats_->start_timer("check_global_order");

  // Check if applicable
  if (!check_global_order_) {
    return Status::Ok();
  }

  // Applicable only to sparse writes - exit if coordinates do not exist
  if (!coords_info_.has_coords_ || coords_info_.coords_num_ == 0) {
    return Status::Ok();
  }

  // Special case for Hilbert
  if (array_schema_.cell_order() == Layout::HILBERT) {
    return check_global_order_hilbert();
  }

  // Make sure the last cell written by a previous write comes before the
  // first cell of this write in the global order
  auto& domain{array_schema_.domain()};
  DomainBuffersView domain_buffs{array_schema_, buffers_};
  if (global_write_state_->cells_written_.begin()->second > 0) {
    DomainBuffersView last_cell_buffs{array_schema_,
                                      *global_write_state_->last_cell_coords_};
    auto left{last_cell_buffs.domain_ref_at(domain, 0)};
    auto right{domain_buffs.domain_ref_at(domain, 0)};

    auto tile_cmp = domain.tile_order_cmp(left, right);
    auto fail = (tile_cmp > 0) ||
                ((tile_cmp == 0) && domain.cell_order_cmp(left, right) > 0);
    if (fail) {
      std::stringstream ss;
      ss << "Write failed; Coordinate " << coords_to_str(0);
      ss << " comes before last written coordinate in the global order";
      if (tile_cmp > 0)
        ss << " due to writes across tiles";
      return Status_WriterError(ss.str());
    }
  }

  // Check if all coordinates are in global order in parallel
  auto status = parallel_for(
      storage_manager_->compute_tp(),
      0,
      coords_info_.coords_num_ - 1,
      [&](uint64_t i) {
        auto left{domain_buffs.domain_ref_at(domain, i)};
        auto right{domain_buffs.domain_ref_at(domain, i + 1)};
        auto tile_cmp = domain.tile_order_cmp(left, right);
        auto fail = (tile_cmp > 0) ||
                    ((tile_cmp == 0) && domain.cell_order_cmp(left, right) > 0);
        if (fail) {
          std::stringstream ss;
          ss << "Write failed; Coordinates " << coords_to_str(i);
          ss << " succeed " << coords_to_str(i + 1);
          ss << " in the global order";
          if (tile_cmp > 0)
            ss << " due to writes across tiles";
          return Status_WriterError(ss.str());
        }
        return Status::Ok();
      });

  RETURN_NOT_OK(status);

  // Save the last cell's coordinates.
  auto last_cell_coords{
      domain_buffs.domain_ref_at(domain, coords_info_.coords_num_ - 1)};
  global_write_state_->last_cell_coords_ =
      SingleCoord(array_schema_, last_cell_coords);

  return Status::Ok();
}

Status GlobalOrderWriter::check_global_order_hilbert() const {
  // Compute hilbert values
  DomainBuffersView domain_buffs{array_schema_, buffers_};
  std::vector<uint64_t> hilbert_values(coords_info_.coords_num_);
  RETURN_NOT_OK(calculate_hilbert_values(domain_buffs, hilbert_values));

  // Make sure the last cell written by a previous write comes before the
  // first cell of this write in the hilbert order
  if (global_write_state_->cells_written_.begin()->second > 0) {
    if (global_write_state_->last_hilbert_value_ > hilbert_values[0]) {
      std::stringstream ss;
      ss << "Write failed; Coordinates " << coords_to_str(0);
      ss << " comes before last written coordinate in the hilbert order";
      return Status_WriterError(ss.str());
    }
  }

  // Check if all coordinates are in hilbert order in parallel
  auto status = parallel_for(
      storage_manager_->compute_tp(),
      0,
      coords_info_.coords_num_ - 1,
      [&](uint64_t i) {
        if (hilbert_values[i] > hilbert_values[i + 1]) {
          std::stringstream ss;
          ss << "Write failed; Coordinates " << coords_to_str(i);
          ss << " succeed " << coords_to_str(i + 1);
          ss << " in the hilbert order";
          return Status_WriterError(ss.str());
        }
        return Status::Ok();
      });

  RETURN_NOT_OK(status);

  // Save the last hilbert value
  global_write_state_->last_hilbert_value_ =
      hilbert_values[coords_info_.coords_num_ - 1];

  return Status::Ok();
}

void GlobalOrderWriter::clean_up() {
  if (global_write_state_ != nullptr) {
    auto meta = global_write_state_->frag_meta_;
    const auto& uri = meta->fragment_uri();
    throw_if_not_ok(storage_manager_->vfs()->remove_dir(uri));
    global_write_state_.reset(nullptr);

    // Cleanup all fragments pending commit.
    for (auto& uri : frag_uris_to_commit_) {
      throw_if_not_ok(storage_manager_->vfs()->remove_dir(uri));
    }
    frag_uris_to_commit_.clear();
  }
}

Status GlobalOrderWriter::filter_last_tiles(uint64_t cell_num) {
  // Adjust cell num
  for (auto& last_tiles : global_write_state_->last_tiles_) {
    last_tiles.second[0].set_final_size(cell_num);
  }

  // Compute coordinates metadata
  auto meta = global_write_state_->frag_meta_;
  auto mbrs = compute_mbrs(global_write_state_->last_tiles_);
  set_coords_metadata(0, 1, global_write_state_->last_tiles_, mbrs, meta);

  // Compute tile metadata.
  RETURN_NOT_OK(compute_tiles_metadata(1, global_write_state_->last_tiles_));

  // Gather stats
  stats_->add_counter(
      "cell_num",
      global_write_state_->last_tiles_.begin()->second[0].cell_num());
  stats_->add_counter("tile_num", 1);

  // Filter tiles
  RETURN_NOT_OK(filter_tiles(&global_write_state_->last_tiles_));

  return Status::Ok();
}

Status GlobalOrderWriter::compute_coord_dups(
    std::set<uint64_t>* coord_dups) const {
  auto timer_se = stats_->start_timer("compute_coord_dups");

  if (!coords_info_.has_coords_) {
    return logger_->status(
        Status_WriterError("Cannot check for coordinate duplicates; "
                           "Coordinates buffer not found"));
  }

  if (coords_info_.coords_num_ < 2)
    return Status::Ok();

  // Prepare auxiliary vectors for better performance
  auto dim_num = array_schema_.dim_num();
  std::vector<const unsigned char*> buffs(dim_num);
  std::vector<uint64_t> coord_sizes(dim_num);
  std::vector<const unsigned char*> buffs_var(dim_num);
  std::vector<uint64_t*> buffs_var_sizes(dim_num);
  for (unsigned d = 0; d < dim_num; ++d) {
    const auto& dim_name{array_schema_.dimension_ptr(d)->name()};
    buffs[d] = (const unsigned char*)buffers_.find(dim_name)->second.buffer_;
    coord_sizes[d] = array_schema_.cell_size(dim_name);
    buffs_var[d] =
        (const unsigned char*)buffers_.find(dim_name)->second.buffer_var_;
    buffs_var_sizes[d] = buffers_.find(dim_name)->second.buffer_var_size_;
  }

  std::mutex mtx;
  auto status = parallel_for(
      storage_manager_->compute_tp(),
      1,
      coords_info_.coords_num_,
      [&](uint64_t i) {
        // Check for duplicate in adjacent cells
        bool found_dup = true;
        for (unsigned d = 0; d < dim_num; ++d) {
          auto dim{array_schema_.dimension_ptr(d)};
          if (!dim->var_size()) {  // Fixed-sized dimensions
            if (memcmp(
                    buffs[d] + i * coord_sizes[d],
                    buffs[d] + (i - 1) * coord_sizes[d],
                    coord_sizes[d]) != 0) {  // Not the same
              found_dup = false;
              break;
            }
          } else {
            auto offs = (uint64_t*)buffs[d];
            auto off_i_plus_1 = (i == coords_info_.coords_num_ - 1) ?
                                    *(buffs_var_sizes[d]) :
                                    offs[i + 1];
            auto size_i_minus_1 = offs[i] - offs[i - 1];
            auto size_i = off_i_plus_1 - offs[i];

            // Compare sizes
            if (size_i != size_i_minus_1) {  // Not same
              found_dup = false;
              break;
            }

            // Compare var values
            if (memcmp(
                    buffs_var[d] + offs[i - 1],
                    buffs_var[d] + offs[i],
                    size_i) != 0) {  // Not the same
              found_dup = false;
              break;
            }
          }
        }

        // Found duplicate
        if (found_dup) {
          std::lock_guard<std::mutex> lock(mtx);
          coord_dups->insert(i);
        }

        return Status::Ok();
      });

  RETURN_NOT_OK(status);

  return Status::Ok();
}

Status GlobalOrderWriter::finalize_global_write_state() {
  assert(layout_ == Layout::GLOBAL_ORDER);
  auto meta = global_write_state_->frag_meta_;
  const auto& uri = meta->fragment_uri();

  // Handle last tile
  Status st = global_write_handle_last_tile();
  if (!st.ok()) {
    throw_if_not_ok(close_files(meta));
    return st;
  }

  // Close all files
  RETURN_NOT_OK(close_files(meta));

  // Check that the same number of cells was written across attributes
  // and dimensions
  auto cell_num = global_write_state_->cells_written_[buffers_.begin()->first];
  for (const auto& it : buffers_) {
    const auto& name = it.first;
    if (global_write_state_->cells_written_[name] != cell_num) {
      return logger_->status(Status_WriterError(
          "Failed to finalize global write state; Different "
          "number of cells written across attributes and coordinates"));
    }
  }

  // No cells written, clean up empty fragment.
  if (cell_num == 0) {
    return Status::Ok();
  }

  // Check if the total number of cells written is equal to the subarray size
  if (!coords_info_.has_coords_) {  // This implies a dense array
    auto expected_cell_num =
        array_schema_.domain().cell_num(subarray_.ndrange(0));
    if (cell_num != expected_cell_num) {
      std::stringstream ss;
      ss << "Failed to finalize global write state; Number "
         << "of cells written (" << cell_num
         << ") is different from the number of cells expected ("
         << expected_cell_num << ") for the query subarray";
      return logger_->status(Status_WriterError(ss.str()));
    }
  }

  // Set the processed conditions
  meta->set_processed_conditions(processed_conditions_);

  // Compute fragment min/max/sum/null count and flush fragment metadata to
  // storage
  meta->compute_fragment_min_max_sum_null_count();
  meta->store(array_->get_encryption_key());

  // Add written fragment infos
  for (auto& frag_uri : frag_uris_to_commit_) {
    RETURN_NOT_OK(add_written_fragment_info(frag_uri));
  }
  RETURN_NOT_OK(add_written_fragment_info(uri));

  // The following will make the fragment visible
  URI commit_uri = array_->array_directory().get_commit_uri(uri);

  // Write either one commit file or a consolidated commit file if multiple
  // fragments were written.
  if (frag_uris_to_commit_.size() == 0) {
    RETURN_NOT_OK(storage_manager_->vfs()->touch(commit_uri));
  } else {
    std::vector<URI> commit_uris;
    commit_uris.reserve(frag_uris_to_commit_.size() + 1);
    for (auto& uri : frag_uris_to_commit_) {
      commit_uris.emplace_back(array_->array_directory().get_commit_uri(uri));
    }
    commit_uris.emplace_back(commit_uri);

    auto write_version = array_->array_schema_latest().write_version();
    storage_manager_->write_consolidated_commits_file(
        write_version, array_->array_directory(), commit_uris);
  }

  // Delete global write state
  global_write_state_.reset(nullptr);

  return st;
}

Status GlobalOrderWriter::global_write() {
  // Applicable only to global write on dense/sparse arrays
  assert(layout_ == Layout::GLOBAL_ORDER);

  // Initialize the global write state if this is the first invocation
  if (!global_write_state_) {
    RETURN_CANCEL_OR_ERROR(alloc_global_write_state());
    RETURN_CANCEL_OR_ERROR(create_fragment(
        !coords_info_.has_coords_, global_write_state_->frag_meta_));
    RETURN_CANCEL_OR_ERROR(init_global_write_state());
  }

  // Check for coordinate duplicates
  if (coords_info_.has_coords_) {
    RETURN_CANCEL_OR_ERROR(check_coord_dups());
    RETURN_CANCEL_OR_ERROR(check_global_order());
  }

  // Retrieve coordinate duplicates
  std::set<uint64_t> coord_dups;
  if (dedup_coords_) {
    RETURN_CANCEL_OR_ERROR(compute_coord_dups(&coord_dups));
  }

  std::unordered_map<std::string, WriterTileVector> tiles;
  RETURN_CANCEL_OR_ERROR(prepare_full_tiles(coord_dups, &tiles));

  // Find number of tiles and gather stats
  uint64_t tile_num = 0;
  if (!tiles.empty()) {
    auto it = tiles.begin();
    tile_num = it->second.size();

    uint64_t cell_num = 0;
    for (size_t t = 0; t < tile_num; ++t) {
      cell_num += it->second[t].cell_num();
    }
    stats_->add_counter("cell_num", cell_num);
    stats_->add_counter("tile_num", tile_num);
  }

  // No cells to be written
  if (tile_num == 0) {
    return Status::Ok();
  }

  // Compute coordinate metadata (if coordinates are present)
  auto mbrs = compute_mbrs(tiles);

  // Compute tile metadata.
  RETURN_CANCEL_OR_ERROR(compute_tiles_metadata(tile_num, tiles));

  // Filter all tiles
  RETURN_CANCEL_OR_ERROR(filter_tiles(&tiles));

  uint64_t idx = 0;
  while (idx < tile_num) {
    auto frag_meta = global_write_state_->frag_meta_;

    // Compute the number of tiles that will fit in this fragment.
    auto num = num_tiles_to_write(idx, tile_num, tiles);

    // Set new number of tiles in the fragment metadata
    auto new_num_tiles = frag_meta->tile_index_base() + num;
    throw_if_not_ok(frag_meta->set_num_tiles(new_num_tiles));

    if (new_num_tiles == 0) {
      throw GlobalOrderWriterStatusException(
          "Fragment size is too small to write a single tile");
    }

    set_coords_metadata(idx, idx + num, tiles, mbrs, frag_meta);

    // Write tiles for all attributes
    RETURN_CANCEL_OR_ERROR(write_tiles(idx, idx + num, frag_meta, &tiles));
    idx += num;

    // If we didn't write all tiles, close this fragment and start another.
    if (idx != tile_num) {
      RETURN_CANCEL_OR_ERROR(start_new_fragment());
    }

    // Increment the tile index base for the next global order write.
    frag_meta->set_tile_index_base(new_num_tiles);
  }

  return Status::Ok();
}

Status GlobalOrderWriter::global_write_handle_last_tile() {
  auto capacity = array_schema_.capacity();
  auto& domain = array_schema_.domain();
  auto cell_num_per_tile =
      coords_info_.has_coords_ ? capacity : domain.cell_num_per_tile();
  auto cell_num_last_tiles =
      global_write_state_->cells_written_[buffers_.begin()->first] %
      cell_num_per_tile;
  if (cell_num_last_tiles == 0)
    return Status::Ok();

  // Reserve space for the last tile in the fragment metadata
  auto meta = global_write_state_->frag_meta_;
  throw_if_not_ok(meta->set_num_tiles(meta->tile_index_base() + 1));

  // Filter last tiles
  RETURN_CANCEL_OR_ERROR(filter_last_tiles(cell_num_last_tiles));

  // Write the last tiles
  RETURN_CANCEL_OR_ERROR(
      write_tiles(0, 1, meta, &global_write_state_->last_tiles_));

  // Increment the tile index base.
  meta->set_tile_index_base(meta->tile_index_base() + 1);

  return Status::Ok();
}

void GlobalOrderWriter::nuke_global_write_state() {
  auto meta = global_write_state_->frag_meta_;
  throw_if_not_ok(close_files(meta));
  throw_if_not_ok(storage_manager_->vfs()->remove_dir(meta->fragment_uri()));
  global_write_state_.reset(nullptr);
}

Status GlobalOrderWriter::prepare_full_tiles(
    const std::set<uint64_t>& coord_dups,
    std::unordered_map<std::string, WriterTileVector>* tiles) const {
  auto timer_se = stats_->start_timer("prepare_tiles");

  // Initialize attribute and coordinate tiles
  for (const auto& it : buffers_) {
    (*tiles)[it.first] = WriterTileVector();
  }

  auto num = buffers_.size();
  auto status =
      parallel_for(storage_manager_->compute_tp(), 0, num, [&](uint64_t i) {
        auto buff_it = buffers_.begin();
        std::advance(buff_it, i);
        const auto& name = buff_it->first;
        RETURN_CANCEL_OR_ERROR(
            prepare_full_tiles(name, coord_dups, &(*tiles)[name]));
        return Status::Ok();
      });

  RETURN_NOT_OK(status);

  return Status::Ok();
}

Status GlobalOrderWriter::prepare_full_tiles(
    const std::string& name,
    const std::set<uint64_t>& coord_dups,
    WriterTileVector* tiles) const {
  return array_schema_.var_size(name) ?
             prepare_full_tiles_var(name, coord_dups, tiles) :
             prepare_full_tiles_fixed(name, coord_dups, tiles);
}

Status GlobalOrderWriter::prepare_full_tiles_fixed(
    const std::string& name,
    const std::set<uint64_t>& coord_dups,
    WriterTileVector* tiles) const {
  // For easy reference
  auto nullable = array_schema_.is_nullable(name);
  auto type = array_schema_.type(name);
  auto it = buffers_.find(name);
  auto buffer = (unsigned char*)it->second.buffer_;
  auto buffer_validity = (unsigned char*)it->second.validity_vector_.buffer();
  auto buffer_size = it->second.buffer_size_;
  auto cell_size = array_schema_.cell_size(name);
  auto capacity = array_schema_.capacity();
  auto cell_num = *buffer_size / cell_size;
  auto& domain{array_schema_.domain()};
  auto cell_num_per_tile =
      coords_info_.has_coords_ ? capacity : domain.cell_num_per_tile();

  // Do nothing if there are no cells to write
  if (cell_num == 0) {
    return Status::Ok();
  }

  // First fill the last tile
  auto& last_tile = global_write_state_->last_tiles_[name][0];
  uint64_t cell_idx = 0;
  uint64_t last_tile_cell_idx =
      global_write_state_->cells_written_[name] % cell_num_per_tile;
  if (last_tile_cell_idx != 0) {
    if (coord_dups.empty()) {
      do {
        RETURN_NOT_OK(last_tile.fixed_tile().write(
            buffer + cell_idx * cell_size,
            last_tile_cell_idx * cell_size,
            cell_size));
        if (nullable) {
          RETURN_NOT_OK(last_tile.validity_tile().write(
              buffer_validity + cell_idx * constants::cell_validity_size,
              last_tile_cell_idx * constants::cell_validity_size,
              constants::cell_validity_size));
        }
        ++cell_idx;
        ++last_tile_cell_idx;
      } while (last_tile_cell_idx != cell_num_per_tile && cell_idx != cell_num);
    } else {
      do {
        if (coord_dups.find(cell_idx) == coord_dups.end()) {
          RETURN_NOT_OK(last_tile.fixed_tile().write(
              buffer + cell_idx * cell_size,
              last_tile_cell_idx * cell_size,
              cell_size));
          if (nullable) {
            RETURN_NOT_OK(last_tile.validity_tile().write(
                buffer_validity + cell_idx * constants::cell_validity_size,
                last_tile_cell_idx * constants::cell_validity_size,
                constants::cell_validity_size));
          }
          ++last_tile_cell_idx;
        }
        ++cell_idx;
      } while (last_tile_cell_idx != cell_num_per_tile && cell_idx != cell_num);
    }
  }

  // Initialize full tiles and set previous last tile as first tile
  auto full_tile_num = (cell_num - cell_idx) / cell_num_per_tile +
                       (int)(last_tile_cell_idx == cell_num_per_tile);
  auto cell_num_to_write =
      (full_tile_num - (last_tile_cell_idx == cell_num_per_tile)) *
      cell_num_per_tile;

  if (full_tile_num > 0) {
    tiles->reserve(full_tile_num);
    for (uint64_t i = 0; i < full_tile_num; i++) {
      tiles->emplace_back(WriterTile(
          array_schema_, cell_num_per_tile, false, nullable, cell_size, type));
    }

    // Handle last tile (it must be either full or empty)
    auto tile_it = tiles->begin();
    if (last_tile_cell_idx == cell_num_per_tile) {
      tile_it->fixed_tile().swap(last_tile.fixed_tile());
      if (nullable) {
        tile_it->validity_tile().swap(last_tile.validity_tile());
      }
      tile_it++;
    } else if (last_tile_cell_idx != 0) {
      return Status_WriterError(
          "Last tile was not empty when it should have been");
    }

    // Write all remaining cells one by one
    if (coord_dups.empty()) {
      for (uint64_t i = 0; i < cell_num_to_write;) {
        RETURN_NOT_OK(tile_it->fixed_tile().write(
            buffer + cell_idx * cell_size, 0, cell_size * cell_num_per_tile));

        if (nullable) {
          RETURN_NOT_OK(tile_it->validity_tile().write(
              buffer_validity + cell_idx * constants::cell_validity_size,
              0,
              constants::cell_validity_size * cell_num_per_tile));
        }

        cell_idx += cell_num_per_tile;
        i += cell_num_per_tile;
        tile_it++;
      }
    } else {
      uint64_t current_tile_cell_idx = 0;
      for (uint64_t i = 0; i < cell_num_to_write; ++cell_idx, ++i) {
        if (current_tile_cell_idx == cell_num_per_tile) {
          tile_it++;
          current_tile_cell_idx = 0;
        }

        if (coord_dups.find(cell_idx) == coord_dups.end()) {
          RETURN_NOT_OK(tile_it->fixed_tile().write(
              buffer + cell_idx * cell_size, cell_idx * cell_size, cell_size));

          if (nullable) {
            RETURN_NOT_OK(tile_it->validity_tile().write(
                buffer_validity + cell_idx * constants::cell_validity_size,
                cell_idx * constants::cell_validity_size,
                constants::cell_validity_size));
          }
        }
      }
    }
  }

  // Potentially fill the last tile
  last_tile_cell_idx = 0;
  if (coord_dups.empty()) {
    for (; cell_idx < cell_num; ++cell_idx, ++last_tile_cell_idx) {
      RETURN_NOT_OK(last_tile.fixed_tile().write(
          buffer + cell_idx * cell_size,
          last_tile_cell_idx * cell_size,
          cell_size));
      if (nullable) {
        RETURN_NOT_OK(last_tile.validity_tile().write(
            buffer_validity + cell_idx * constants::cell_validity_size,
            last_tile_cell_idx * constants::cell_validity_size,
            constants::cell_validity_size));
      }
    }
  } else {
    for (; cell_idx < cell_num; ++cell_idx) {
      if (coord_dups.find(cell_idx) == coord_dups.end()) {
        RETURN_NOT_OK(last_tile.fixed_tile().write(
            buffer + cell_idx * cell_size,
            last_tile_cell_idx * cell_size,
            cell_size));
        if (nullable) {
          RETURN_NOT_OK(last_tile.validity_tile().write(
              buffer_validity + cell_idx * constants::cell_validity_size,
              last_tile_cell_idx * constants::cell_validity_size,
              constants::cell_validity_size));
        }
        ++last_tile_cell_idx;
      }
    }
  }

  global_write_state_->cells_written_[name] += cell_num;

  return Status::Ok();
}

Status GlobalOrderWriter::prepare_full_tiles_var(
    const std::string& name,
    const std::set<uint64_t>& coord_dups,
    WriterTileVector* tiles) const {
  // For easy reference
  auto it = buffers_.find(name);
  auto nullable = array_schema_.is_nullable(name);
  auto cell_size = array_schema_.cell_size(name);
  auto type = array_schema_.type(name);
  auto buffer = (uint64_t*)it->second.buffer_;
  auto buffer_var = (unsigned char*)it->second.buffer_var_;
  auto buffer_validity = (uint8_t*)it->second.validity_vector_.buffer();
  auto buffer_size = get_offset_buffer_size(*it->second.buffer_size_);
  auto buffer_var_size = it->second.buffer_var_size_;
  auto capacity = array_schema_.capacity();
  auto cell_num = buffer_size / constants::cell_var_offset_size;
  auto& domain{array_schema_.domain()};
  auto cell_num_per_tile =
      coords_info_.has_coords_ ? capacity : domain.cell_num_per_tile();
  auto attr_datatype_size = datatype_size(array_schema_.type(name));

  // Do nothing if there are no cells to write
  if (cell_num == 0)
    return Status::Ok();

  // First fill the last tile
  auto& last_tile = global_write_state_->last_tiles_[name][0];
  auto& last_var_offset = global_write_state_->last_var_offsets_[name];
  uint64_t cell_idx = 0;
  uint64_t last_tile_cell_idx =
      global_write_state_->cells_written_[name] % cell_num_per_tile;
  if (last_tile_cell_idx != 0) {
    if (coord_dups.empty()) {
      do {
        // Write offset.
        RETURN_NOT_OK(last_tile.offset_tile().write(
            &last_var_offset,
            last_tile_cell_idx * sizeof(last_var_offset),
            sizeof(last_var_offset)));

        // Write var-sized value(s).
        auto buff_offset =
            prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
        uint64_t var_size = (cell_idx == cell_num - 1) ?
                                *buffer_var_size - buff_offset :
                                prepare_buffer_offset(
                                    buffer, cell_idx + 1, attr_datatype_size) -
                                    buff_offset;
        RETURN_NOT_OK(last_tile.var_tile().write_var(
            buffer_var + buff_offset, last_var_offset, var_size));
        last_var_offset += var_size;

        // Write validity value(s).
        if (nullable) {
          RETURN_NOT_OK(last_tile.validity_tile().write(
              buffer_validity + cell_idx,
              last_tile_cell_idx * constants::cell_validity_size,
              constants::cell_validity_size));
        }

        ++cell_idx;
        ++last_tile_cell_idx;
      } while (last_tile_cell_idx != cell_num_per_tile && cell_idx != cell_num);
    } else {
      do {
        if (coord_dups.find(cell_idx) == coord_dups.end()) {
          // Write offset.
          RETURN_NOT_OK(last_tile.offset_tile().write(
              &last_var_offset,
              last_tile_cell_idx * sizeof(last_var_offset),
              sizeof(last_var_offset)));

          // Write var-sized value(s).
          auto buff_offset =
              prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
          uint64_t var_size =
              (cell_idx == cell_num - 1) ?
                  *buffer_var_size - buff_offset :
                  prepare_buffer_offset(
                      buffer, cell_idx + 1, attr_datatype_size) -
                      buff_offset;
          RETURN_NOT_OK(last_tile.var_tile().write_var(
              buffer_var + buff_offset, last_var_offset, var_size));
          last_var_offset += var_size;

          // Write validity value(s).
          if (nullable) {
            RETURN_NOT_OK(last_tile.validity_tile().write(
                buffer_validity + cell_idx,
                last_tile_cell_idx * constants::cell_validity_size,
                constants::cell_validity_size));
          }
          ++last_tile_cell_idx;
        }

        ++cell_idx;
      } while (last_tile_cell_idx != cell_num_per_tile && cell_idx != cell_num);
    }

    last_tile.var_tile().set_size(last_var_offset);
  }

  // Initialize full tiles and set previous last tile as first tile
  auto full_tile_num = (cell_num - cell_idx) / cell_num_per_tile +
                       (last_tile_cell_idx == cell_num_per_tile);
  auto cell_num_to_write =
      (full_tile_num - (last_tile_cell_idx == cell_num_per_tile)) *
      cell_num_per_tile;

  if (full_tile_num > 0) {
    tiles->reserve(full_tile_num);
    for (uint64_t i = 0; i < full_tile_num; i++) {
      tiles->emplace_back(WriterTile(
          array_schema_, cell_num_per_tile, true, nullable, cell_size, type));
    }

    // Handle last tile (it must be either full or empty)
    auto tile_it = tiles->begin();
    if (last_tile_cell_idx == cell_num_per_tile) {
      last_var_offset = 0;
      tile_it->offset_tile().swap(last_tile.offset_tile());
      tile_it->var_tile().swap(last_tile.var_tile());
      if (nullable) {
        tile_it->validity_tile().swap(last_tile.validity_tile());
      }
      tile_it++;
    } else if (last_tile_cell_idx != 0) {
      return Status_WriterError(
          "Last tile was not empty when it should have been");
    }

    // Write all remaining cells one by one
    uint64_t current_tile_cell_idx = 0;
    if (cell_num_to_write != 0) {
      if (coord_dups.empty()) {
        for (uint64_t i = 0; i < cell_num_to_write;
             ++cell_idx, ++i, ++current_tile_cell_idx) {
          if (current_tile_cell_idx == cell_num_per_tile) {
            tile_it->var_tile().set_size(last_var_offset);
            current_tile_cell_idx = 0;
            last_var_offset = 0;
            tile_it++;
          }

          // Write offset.
          RETURN_NOT_OK(tile_it->offset_tile().write(
              &last_var_offset,
              current_tile_cell_idx * sizeof(last_var_offset),
              sizeof(last_var_offset)));

          // Write var-sized value(s).
          auto buff_offset =
              prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
          uint64_t var_size =
              (cell_idx == cell_num - 1) ?
                  *buffer_var_size - buff_offset :
                  prepare_buffer_offset(
                      buffer, cell_idx + 1, attr_datatype_size) -
                      buff_offset;
          RETURN_NOT_OK(tile_it->var_tile().write_var(
              buffer_var + buff_offset, last_var_offset, var_size));
          last_var_offset += var_size;

          // Write validity value(s).
          if (nullable) {
            RETURN_NOT_OK(tile_it->validity_tile().write(
                buffer_validity + cell_idx,
                current_tile_cell_idx * constants::cell_validity_size,
                constants::cell_validity_size));
          }
        }
      } else {
        for (uint64_t i = 0; i < cell_num_to_write; ++cell_idx, ++i) {
          if (coord_dups.find(cell_idx) == coord_dups.end()) {
            if (current_tile_cell_idx == cell_num_per_tile) {
              tile_it->var_tile().set_size(last_var_offset);
              current_tile_cell_idx = 0;
              last_var_offset = 0;
              tile_it++;
            }

            // Write offset.
            RETURN_NOT_OK(tile_it->offset_tile().write(
                &last_var_offset,
                current_tile_cell_idx * sizeof(last_var_offset),
                sizeof(last_var_offset)));

            // Write var-sized value(s).
            auto buff_offset =
                prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
            uint64_t var_size =
                (cell_idx == cell_num - 1) ?
                    *buffer_var_size - buff_offset :
                    prepare_buffer_offset(
                        buffer, cell_idx + 1, attr_datatype_size) -
                        buff_offset;
            RETURN_NOT_OK(tile_it->var_tile().write_var(
                buffer_var + buff_offset, last_var_offset, var_size));
            last_var_offset += var_size;

            // Write validity value(s).
            if (nullable) {
              RETURN_NOT_OK(tile_it->validity_tile().write(
                  buffer_validity + cell_idx,
                  current_tile_cell_idx * constants::cell_validity_size,
                  constants::cell_validity_size));
            }

            ++current_tile_cell_idx;
          }
        }
      }

      tile_it->var_tile().set_size(last_var_offset);
      last_var_offset = 0;
    }
  }

  // Potentially fill the last tile
  last_tile_cell_idx = 0;
  if (coord_dups.empty()) {
    for (; cell_idx < cell_num; ++cell_idx, ++last_tile_cell_idx) {
      // Write offset.
      RETURN_NOT_OK(last_tile.offset_tile().write(
          &last_var_offset,
          last_tile_cell_idx * sizeof(last_var_offset),
          sizeof(last_var_offset)));

      // Write var-sized value(s).
      auto buff_offset =
          prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
      uint64_t var_size =
          (cell_idx == cell_num - 1) ?
              *buffer_var_size - buff_offset :
              prepare_buffer_offset(buffer, cell_idx + 1, attr_datatype_size) -
                  buff_offset;
      RETURN_NOT_OK(last_tile.var_tile().write_var(
          buffer_var + buff_offset, last_var_offset, var_size));
      last_var_offset += var_size;

      // Write validity value(s).
      if (nullable) {
        RETURN_NOT_OK(last_tile.validity_tile().write(
            buffer_validity + cell_idx,
            last_tile_cell_idx * constants::cell_validity_size,
            constants::cell_validity_size));
      }
    }
  } else {
    for (; cell_idx < cell_num; ++cell_idx) {
      if (coord_dups.find(cell_idx) == coord_dups.end()) {
        // Write offset.
        RETURN_NOT_OK(last_tile.offset_tile().write(
            &last_var_offset,
            last_tile_cell_idx * sizeof(last_var_offset),
            sizeof(last_var_offset)));

        // Write var-sized value(s).
        auto buff_offset =
            prepare_buffer_offset(buffer, cell_idx, attr_datatype_size);
        uint64_t var_size = (cell_idx == cell_num - 1) ?
                                *buffer_var_size - buff_offset :
                                prepare_buffer_offset(
                                    buffer, cell_idx + 1, attr_datatype_size) -
                                    buff_offset;
        RETURN_NOT_OK(last_tile.var_tile().write_var(
            buffer_var + buff_offset, last_var_offset, var_size));
        last_var_offset += var_size;

        // Write validity value(s).
        if (nullable) {
          RETURN_NOT_OK(last_tile.validity_tile().write(
              buffer_validity + cell_idx,
              last_tile_cell_idx * constants::cell_validity_size,
              constants::cell_validity_size));
        }

        ++last_tile_cell_idx;
      }
    }
  }

  last_tile.var_tile().set_size(last_var_offset);

  global_write_state_->cells_written_[name] += cell_num;

  return Status::Ok();
}

uint64_t GlobalOrderWriter::num_tiles_to_write(
    uint64_t start,
    uint64_t tile_num,
    std::unordered_map<std::string, WriterTileVector>& tiles) {
  // Cache variables to prevent map lookups.
  const auto buf_names = buffer_names();
  std::vector<bool> var_size;
  std::vector<bool> nullable;
  std::vector<WriterTileVector*> writer_tile_vectors;
  var_size.reserve(buf_names.size());
  nullable.reserve(buf_names.size());
  writer_tile_vectors.reserve(buf_names.size());
  for (auto& name : buf_names) {
    var_size.emplace_back(array_schema_.var_size(name));
    nullable.emplace_back(array_schema_.is_nullable(name));
    writer_tile_vectors.emplace_back(&tiles[name]);
  }

  // Make sure we don't write more than the desired fragment size.
  for (uint64_t t = start; t < tile_num; t++) {
    uint64_t tile_size = 0;
    for (uint64_t a = 0; a < buf_names.size(); a++) {
      if (var_size[a]) {
        tile_size += writer_tile_vectors[a]
                         ->at(t)
                         .offset_tile()
                         .filtered_buffer()
                         .size();
        tile_size +=
            writer_tile_vectors[a]->at(t).var_tile().filtered_buffer().size();
      } else {
        tile_size +=
            writer_tile_vectors[a]->at(t).fixed_tile().filtered_buffer().size();
      }

      if (nullable[a]) {
        tile_size += writer_tile_vectors[a]
                         ->at(t)
                         .validity_tile()
                         .filtered_buffer()
                         .size();
      }
    }

    if (current_fragment_size_ + tile_size > fragment_size_) {
      return t - start;
    }

    current_fragment_size_ += tile_size;
  }

  return tile_num - start;
}

Status GlobalOrderWriter::start_new_fragment() {
  auto frag_meta = global_write_state_->frag_meta_;
  auto& uri = frag_meta->fragment_uri();

  // Close all files
  RETURN_NOT_OK(close_files(frag_meta));

  // Set the processed conditions
  frag_meta->set_processed_conditions(processed_conditions_);

  // Compute fragment min/max/sum/null count
  frag_meta->compute_fragment_min_max_sum_null_count();

  // Flush fragment metadata to storage
  frag_meta->store(array_->get_encryption_key());

  frag_uris_to_commit_.emplace_back(uri);

  // Make a new fragment URI.
  const auto write_version = array_->array_schema_latest().write_version();
  auto frag_dir_uri =
      array_->array_directory().get_fragments_dir(write_version);
  auto new_fragment_str = storage_format::generate_uri(
      fragment_timestamp_range_.first,
      fragment_timestamp_range_.second,
      write_version);
  fragment_uri_ = frag_dir_uri.join_path(new_fragment_str);

  // Create a new fragment.
  current_fragment_size_ = 0;
  RETURN_NOT_OK(create_fragment(
      !coords_info_.has_coords_, global_write_state_->frag_meta_));

  return Status::Ok();
}

}  // namespace sm
}  // namespace tiledb
