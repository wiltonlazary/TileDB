/**
 * @file   reading_sparse_global.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2018-2021 TileDB, Inc.
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
 * When run, this program will create a simple 2D sparse array, write some data
 * to it, and read a slice of the data back in the global layout.
 */

#include <chrono>
#include <cmath>
#include <iostream>
#include <tiledb/tiledb>
//#include "tiledb/common/status.h"

using namespace std::chrono;
using namespace tiledb;

// Name of array.
std::string array_name("sparse_global_order_reader_array");

template <typename T>
struct test_dim_t {
  test_dim_t(
      const std::string& name,
      const std::array<T, 2>& domain,
      const uint64_t tile_extent)
      : name_(name)
      , domain_(domain)
      , tile_extent_(tile_extent) {
  }

  std::string name_;
  std::array<T, 2> domain_;
  uint64_t tile_extent_;
};

template <typename T>
struct test_attr_t {
  test_attr_t(const std::string& name)
      : name_(name) {
  }

  std::string name_;
};

template <typename T>
struct test_query_buffer_t {
  test_query_buffer_t(const std::string& name, std::vector<T>* const data)
      : name_(name)
      , data_(data) {
  }

  std::string name_;
  std::vector<T>* data_;
};

template <typename DIM_T, typename ATTR_T>
void create_array(
    const std::vector<test_dim_t<DIM_T>>& test_dims,
    const std::vector<test_attr_t<ATTR_T>>& test_attrs) {
  // Create domain.
  Context ctx;
  Domain domain(ctx);

  // Create the dimensions.
  for (const auto& test_dim : test_dims) {
    domain.add_dimension(Dimension::create<DIM_T>(
        ctx, test_dim.name_, test_dim.domain_, test_dim.tile_extent_));
  }

  // The array will be sparse.
  ArraySchema schema(ctx, TILEDB_SPARSE);
  schema.set_domain(domain)
      .set_order({{TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR}})
      .set_allows_dups(true);

  // Create the attributes.
  for (const auto& test_attr : test_attrs) {
    Attribute attr = Attribute::create<ATTR_T>(ctx, test_attr.name_);
    schema.add_attribute(attr);
  }

  // Check the array schema.
  schema.check();

  // Create the array.
  Array::create(array_name, schema);
}

template <typename ATTR_T>
Query::Status write(
    std::vector<test_query_buffer_t<ATTR_T>>& test_query_buffers) {
  // Open the array for writing.
  Config config;
  config["sm.use_refactored_readers"] = true;
  Context ctx(config);
  Array array(ctx, array_name, TILEDB_WRITE);

  // Create the queries with unordered layouts.
  Query query(ctx, array);
  query.set_layout(TILEDB_UNORDERED);

  // Set the query buffers.
  for (auto& test_query_buffer : test_query_buffers) {
    query.set_data_buffer(test_query_buffer.name_, *test_query_buffer.data_);
  }

  // Submit and finalize the queries.
  auto s = query.submit();
  std::cerr << s << std::endl;
  query.finalize();

  for (auto& test_query_buffer : test_query_buffers) {
    test_query_buffer.data_->clear();
  }
  test_query_buffers.clear();

  // Close the array.
  array.close();
  return s;
}

template <typename ATTR_T>
Query::Status read_array(
    std::vector<test_query_buffer_t<ATTR_T>>& test_query_buffers) {  //,
  // uint64_t full_domain) {
  Config config;
  config["sm.use_refactored_readers"] = "true";
  Context ctx(config);

  // Open the array for reading and create the read query
  Array array(ctx, array_name, TILEDB_READ);
  Query query(ctx, array, TILEDB_READ);
  query.set_layout(TILEDB_GLOBAL_ORDER);

  // Set the query buffers.
  for (auto& test_query_buffer : test_query_buffers)
    query.set_data_buffer(test_query_buffer.name_, *test_query_buffer.data_);

  // Submit the query.
  auto status = query.submit();

  // std::string stats = query.stats();
  // std::cerr << stats << std::endl;

  // Validate
  /*auto result_elements = query.result_buffer_elements();
  auto num_a_elements = result_elements["a"].second;
  if (num_a_elements != full_domain) {
    std::cerr << "Reader Error: The number of fragments read ("
      << num_a_elements << ") does not match the expected value ("
      << full_domain << ")." << std::endl;
  }*/

  // Close the array.
  array.close();

  return status;
}

template <typename DIM_T>
void add_tile_coords(
    uint64_t i,
    const std::vector<test_dim_t<DIM_T>>& test_dims,
    uint64_t* coords) {
  const auto domain_extent =
      (test_dims[0].domain_[1] - test_dims[0].domain_[0] + 1);
  const int tile_extent = test_dims[0].tile_extent_;
  const int dim_num = test_dims.size();
  const int tiles_per_row_column = domain_extent / tile_extent;

  auto tile_pos = i / (tile_extent * tile_extent);

  auto div = (uint64_t)pow(tiles_per_row_column, dim_num - 1);
  for (int64_t dim = 0; dim < dim_num; dim++) {
    coords[dim] += ((tile_pos / div) % tiles_per_row_column) * tile_extent;
    div /= tiles_per_row_column;
  }
}

template <typename DIM_T>
void add_cell_coords(
    uint64_t i,
    const std::vector<test_dim_t<DIM_T>>& test_dims,
    uint64_t* coords) {
  const uint64_t domain_min = test_dims[0].domain_[0];
  const int tile_extent = test_dims[0].tile_extent_;
  const int dim_num = test_dims.size();

  auto cell_pos = i % (tile_extent * tile_extent);

  auto div = (uint64_t)pow(tile_extent, dim_num - 1);
  for (int64_t dim = 0; dim < dim_num; dim++) {
    coords[dim] += (cell_pos / div) % tile_extent + domain_min;
    div /= tile_extent;
  }
}

void create_write_query_buffer(
    int min_bound,
    int max_bound,
    std::vector<test_dim_t<uint64_t>> dims,
    std::vector<uint64_t>* a,
    std::vector<uint64_t>* rows,
    std::vector<uint64_t>* cols,
    std::vector<test_query_buffer_t<uint64_t>>& query_buffers) {
  for (int i = min_bound; i < max_bound; i++) {
    uint64_t coords[2] = {0, 0};
    add_tile_coords(i, dims, coords);
    add_cell_coords(i, dims, coords);

    rows->emplace_back(coords[0]);
    cols->emplace_back(coords[1]);
    a->emplace_back(i);
  }

  query_buffers.emplace_back("a", a);
  query_buffers.emplace_back("rows", rows);
  query_buffers.emplace_back("cols", cols);
}

void define_dimensions(
    const uint64_t domain_min,
    const uint64_t domain_max,
    const uint64_t tile_extent,
    std::vector<test_dim_t<uint64_t>>& dims) {
  std::array<uint64_t, 2> rows_domain = {domain_min, domain_max};
  dims.emplace_back("rows", rows_domain, tile_extent);
  std::array<uint64_t, 2> cols_domain = {domain_min, domain_max};
  dims.emplace_back("cols", cols_domain, tile_extent);
}

void print(
    // uint64_t full_domain,
    std::vector<uint64_t> data,
    std::vector<uint64_t> coords_rows,
    std::vector<uint64_t> coords_cols) {
  for (uint64_t i = 0; i < 7; i++)
    std::cerr << "{" << coords_rows[i] << "," << coords_cols[i] << "}"
              << " = " << data[i] << std::endl;
}

void validate_data(
    uint64_t validation_min,
    uint64_t validation_max,
    std::string layout,
    std::vector<uint64_t> data,
    std::vector<uint64_t> coords_rows,
    std::vector<uint64_t> coords_cols) {
  if (layout == "ordered" || layout == "interleaved") {
    for (uint64_t i = validation_min; i < validation_max; i++) {
      if (data[i - validation_min] != i) {
        std::cerr << "Error: Data " << data[i - validation_min]
                  << " starting at coordinate {"
                  << coords_rows[i - validation_min] << ","
                  << coords_cols[i - validation_min]
                  << "} is inconsistent with the anticipated value of " << i
                  << "." << std::endl;
        break;
      }
    }
  } else if (layout == "duplicated") {
    uint64_t count = 0;
    for (uint64_t i = validation_min; i < validation_max; i += 2) {
      if (data[i] != count && data[i + 1] != count) {
        std::cerr << "Error: Data at coordinate {" << coords_rows[i] << ","
                  << coords_cols[i]
                  << "} is inconsistent with the anticipated value."
                  << std::endl;
      }

      count++;
    }
  } else {
    std::cerr << "Error: Invalid fragment layout. "
              << "Must be \"ordered\", \"interleaved\", or \"duplicated\""
              << std::endl;
  }
}

// Assume there are always 3 fragments (for now)
// Thus full_domain should be some multiple of 3 and 6
void test(uint64_t full_domain, uint64_t num_fragments, std::string layout) {
  Context ctx;

  // Define the dimensions.
  uint64_t domain_max = ceil(sqrt(4 * full_domain));
  uint64_t tile_extent = ceil(0.2 * domain_max);
  std::vector<test_dim_t<uint64_t>> dims;
  define_dimensions(1, domain_max, tile_extent, dims);

  // Define the attributes.
  std::vector<test_attr_t<uint64_t>> attrs;
  attrs.emplace_back("a");

  // Create the array only if it does not exist.
  if (Object::object(ctx, array_name).type() != Object::Type::Array)
    create_array(dims, attrs);
  std::cerr << "array created" << std::endl;

  // Create buffers for the fragment write queries
  std::vector<test_query_buffer_t<uint64_t>> write_query_buffers;
  std::vector<uint64_t> a_write, row_write, col_write;
  uint64_t last = 0;
  uint64_t iterator = full_domain / num_fragments;
  uint64_t iterator_lower = iterator / 2;
  uint64_t iterator_size = iterator;
  Query::Status status;

  auto write_start = high_resolution_clock::now();
  if (layout == "ordered") {
    while (iterator <= full_domain) {
      create_write_query_buffer(
          last,
          iterator,
          dims,
          &a_write,
          &row_write,
          &col_write,
          write_query_buffers);
      status = write(write_query_buffers);
      last = iterator;
      iterator += iterator_size;
      if (status != Query::Status::COMPLETE)  // Error status
        return;
    }

  } else if (layout == "interleaved") {
    // We need 2 different query buffers
    // One for the lower half and one for the upper half of the domain
    std::vector<std::pair<uint64_t, uint64_t>> domains;
    uint64_t iterator = 0;
    uint64_t fragment_multiplier = 2;
    uint64_t iterator_size =
        full_domain / (fragment_multiplier * num_fragments);

    while (iterator <= full_domain) {
      domains.emplace_back(iterator, iterator + iterator_size);
      iterator += iterator_size;
    }

    while (!domains.empty()) {
      uint64_t rand_iterator_1 = rand() % domains.size();
      std::pair<uint64_t, uint64_t> domain_1 = domains[rand_iterator_1];
      domains.erase(domains.begin() + rand_iterator_1);
      uint64_t rand_iterator_2 = rand() % domains.size();
      std::pair<uint64_t, uint64_t> domain_2 = domains[rand_iterator_2];
      domains.erase(domains.begin() + rand_iterator_2);

      if (domain_1 > domain_2) {
        std::swap(domain_1, domain_2);
      }

      create_write_query_buffer(
          domain_1.first,
          domain_1.second,
          dims,
          &a_write,
          &row_write,
          &col_write,
          write_query_buffers);
      create_write_query_buffer(
          domain_2.first,
          domain_2.second,
          dims,
          &a_write,
          &row_write,
          &col_write,
          write_query_buffers);

      status = write(write_query_buffers);
      if (status != Query::Status::COMPLETE)  // Error status
        return;
    }

  } else if (layout == "duplicated") {
    // We need to write the same query buffer twice.
    // This will result in the first "half" of the data duplicated across
    // the entire domain.
    uint64_t iterator_size = iterator_lower;
    while (iterator_lower <= full_domain / 2) {
      create_write_query_buffer(
          last,
          iterator_lower,
          dims,
          &a_write,
          &row_write,
          &col_write,
          write_query_buffers);
      create_write_query_buffer(
          last,
          iterator_lower,
          dims,
          &a_write,
          &row_write,
          &col_write,
          write_query_buffers);
      status = write(write_query_buffers);
      last = iterator_lower;
      iterator_lower += iterator_size;
      if (status != Query::Status::COMPLETE)  // Error status
        return;
    }

  } else {
    std::cerr << "Error: Invalid fragment layout. "
              << "Must be \"ordered\", \"interleaved\", or \"duplicated\""
              << std::endl;
  }
  auto write_stop = high_resolution_clock::now();
  auto write_duration = duration_cast<milliseconds>(write_stop - write_start);
  std::cerr << "\n[Performance][Write]: " << write_duration.count()
            << " milliseconds." << std::endl;

  // Set up buffers
  std::vector<test_query_buffer_t<uint64_t>> read_query_buffers;
  std::vector<uint64_t> data(full_domain), coords_rows(full_domain),
      coords_cols(full_domain);
  read_query_buffers.emplace_back("a", &data);
  read_query_buffers.emplace_back("rows", &coords_rows);
  read_query_buffers.emplace_back("cols", &coords_cols);

  // Read array.
  auto start = high_resolution_clock::now();
  status = read_array(read_query_buffers);  //, full_domain);
  // if (status != Query::Status::COMPLETE)  // Error status
  //  return;

  // Performance analysis.
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  std::cerr << "\n[Performance][Read]: " << duration.count() << " milliseconds."
            << std::endl;

  // Validate.
  validate_data(0, full_domain, layout, data, coords_rows, coords_cols);
  // print(data, coords_rows, coords_cols);
  // Clear the query buffers
  for (auto& read_query_buffer : read_query_buffers) {
    read_query_buffer.data_->clear();
  }
  read_query_buffers.clear();

  return;
}

int main() {
  Context ctx;

  // Remove the array if it already exists.
  if (Object::object(ctx, array_name).type() == Object::Type::Array)
    Object::remove(ctx, array_name);

  // Note: num_fragments should be about 1% of full_domain
  // Note: full_domain must be divisible by num_fragments
  // Note: if using interleaved or duplicated order, full_domain must also be
  // divisible by num_fragments * 2
  test(1000, 10, "interleaved");
  // Object::remove(ctx, array_name);

  return 0;
}