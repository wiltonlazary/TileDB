/**
 * @file   query.h
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
 * This file defines class Query.
 */

#ifndef TILEDB_QUERY_H
#define TILEDB_QUERY_H

#include <atomic>
#include <functional>
#include <sstream>
#include <utility>
#include <vector>

#include "tiledb/common/logger_public.h"
#include "tiledb/common/status.h"
#include "tiledb/sm/array_schema/array_schema.h"
#include "tiledb/sm/array_schema/dimension.h"
#include "tiledb/sm/array_schema/domain.h"
#include "tiledb/sm/enums/query_status_details.h"
#include "tiledb/sm/fragment/written_fragment_info.h"
#include "tiledb/sm/query/iquery_strategy.h"
#include "tiledb/sm/query/query_buffer.h"
#include "tiledb/sm/query/query_condition.h"
#include "tiledb/sm/query/update_value.h"
#include "tiledb/sm/query/validity_vector.h"
#include "tiledb/sm/storage_manager/storage_manager_declaration.h"
#include "tiledb/sm/subarray/subarray.h"

using namespace tiledb::common;

namespace tiledb {
namespace sm {

class Array;
class ArrayDimensionLabelQueries;

enum class QueryStatus : uint8_t;
enum class QueryType : uint8_t;

/** Processes a (read/write) query. */
class Query {
 public:
  /* ********************************* */
  /*          PUBLIC DATATYPES         */
  /* ********************************* */

  /**
   * Contains any current state related to (de)serialization of this query.
   * Mostly this supports setting buffers on this query that were allocated
   * internally as a part of deserialization (as opposed to user-set buffers).
   */
  struct SerializationState {
    /** Serialization state for a single attribute. */
    struct AttrState {
      /**
       * Buffer holding (or wrapping) fixed-length data, either attribute or
       * offset data.
       */
      Buffer fixed_len_data;

      /** Buffer holding (or wrapping) variable-length data. */
      Buffer var_len_data;

      /** Buffer holding (or wrapping) validity vector data. */
      Buffer validity_len_data;

      /** Value holding the length of the fixed-length data. */
      uint64_t fixed_len_size = 0;

      /** Value holding the length of the variable-length data. */
      uint64_t var_len_size = 0;

      /** Value holding the length of the validity vector data. */
      uint64_t validity_len_size = 0;
    };

    /** Serialization state per attribute. */
    std::unordered_map<std::string, AttrState> attribute_states;
  };

  /**
   * Contains current state related to coords of this query.
   */
  struct CoordsInfo {
    /**
     * True if either zipped coordinates buffer or separate coordinate
     * buffers are set.
     */
    bool has_coords_;

    /** The zipped coordinates buffer potentially set by the user. */
    void* coords_buffer_;

    /** The zipped coordinates buffer size potentially set by the user. */
    uint64_t* coords_buffer_size_;

    /** Keeps track of the number of coordinates across coordinate buffers. */
    uint64_t coords_num_;
  };

  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  /**
   * Constructor. The query type is inherited from the query type of the
   * input array. An optional fragment URI is given as input in
   * case the query will be used as writes and the given URI should be used
   * for the name of the new fragment to be created.
   *
   * @note Array must be a properly opened array.
   *
   * @param array The array that is being queried.
   * @param fragment_uri The full URI for the new fragment. Only used for
   * writes.
   * @param fragment_base_uri Optional base name for new fragment. Only used for
   *     writes and only if fragment_uri is empty.
   */
  Query(
      StorageManager* storage_manager,
      shared_ptr<Array> array,
      optional<std::string> fragment_name = nullopt);

  /** Destructor. */
  virtual ~Query();

  DISABLE_COPY_AND_COPY_ASSIGN(Query);
  DISABLE_MOVE_AND_MOVE_ASSIGN(Query);

  /* ********************************* */
  /*                 API               */
  /* ********************************* */

  /**
   * Adds a range to the (read/write) query on the input dimension by index,
   * in the form of (start, end, stride).
   * The range components must be of the same type as the domain type of the
   * underlying array.
   */
  Status add_range(
      unsigned dim_idx, const void* start, const void* end, const void* stride);

  /**
   * Adds a variable-sized range to the (read/write) query on the input
   * dimension by index, in the form of (start, end).
   */
  Status add_range_var(
      unsigned dim_idx,
      const void* start,
      uint64_t start_size,
      const void* end,
      uint64_t end_size);

  /** Retrieves the number of ranges of the subarray for the given dimension
   * index. */
  Status get_range_num(unsigned dim_idx, uint64_t* range_num) const;

  /**
   * Retrieves a range from a dimension index in the form (start, end, stride).
   *
   * @param dim_idx The dimension to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start The range start to retrieve.
   * @param end The range end to retrieve.
   * @param stride The range stride to retrieve.
   * @return Status
   */
  Status get_range(
      unsigned dim_idx,
      uint64_t range_idx,
      const void** start,
      const void** end,
      const void** stride) const;

  /**
   * Retrieves a range's sizes for a variable-length dimension index
   *
   * @param dim_idx The dimension to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start_size range start size in bytes
   * @param end_size range end size in bytes
   * @return Status
   */
  Status get_range_var_size(
      unsigned dim_idx,
      uint64_t range_idx,
      uint64_t* start_size,
      uint64_t* end_size) const;

  /**
   * Retrieves a range from a variable-length dimension index in the form
   * (start, end).
   *
   * @param dim_idx The dimension to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start The range start to retrieve.
   * @param end The range end to retrieve.
   * @return Status
   */
  Status get_range_var(
      unsigned dim_idx, uint64_t range_idx, void* start, void* end) const;

  /**
   * Adds a range to the (read/write) query on the input dimension by name,
   * in the form of (start, end, stride).
   * The range components must be of the same type as the domain type of the
   * underlying array.
   */

  Status add_range_by_name(
      const std::string& dim_name,
      const void* start,
      const void* end,
      const void* stride);

  /**
   * Adds a variable-sized range to the (read/write) query on the input
   * dimension by name, in the form of (start, end).
   */
  Status add_range_var_by_name(
      const std::string& dim_name,
      const void* start,
      uint64_t start_size,
      const void* end,
      uint64_t end_size);

  /** Retrieves the number of ranges of the subarray for the given dimension
   * name. */
  Status get_range_num_from_name(
      const std::string& dim_name, uint64_t* range_num) const;

  /**
   * Retrieves a range from a dimension name in the form (start, end, stride).
   *
   * @param dim_name The dimension to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start The range start to retrieve.
   * @param end The range end to retrieve.
   * @param stride The range stride to retrieve.
   * @return Status
   */
  Status get_range_from_name(
      const std::string& dim_name,
      uint64_t range_idx,
      const void** start,
      const void** end,
      const void** stride) const;

  /**
   * Retrieves a range's sizes for a variable-length dimension name
   *
   * @param dim_name The dimension name to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start_size range start size in bytes
   * @param end_size range end size in bytes
   * @return Status
   */
  Status get_range_var_size_from_name(
      const std::string& dim_name,
      uint64_t range_idx,
      uint64_t* start_size,
      uint64_t* end_size) const;

  /**
   * Retrieves a range from a variable-length dimension name in the form (start,
   * end).
   *
   * @param dim_name The dimension name to retrieve the range from.
   * @param range_idx The id of the range to retrieve.
   * @param start The range start to retrieve.
   * @param end The range end to retrieve.
   * @return Status
   */
  Status get_range_var_from_name(
      const std::string& dim_name,
      uint64_t range_idx,
      void* start,
      void* end) const;

  /**
   * Gets the estimated result size (in bytes) for the input fixed-sized
   * attribute/dimension.
   */
  Status get_est_result_size(const char* name, uint64_t* size);

  /**
   * Gets the estimated result size (in bytes) for the input var-sized
   * attribute/dimension.
   */
  Status get_est_result_size(
      const char* name, uint64_t* size_off, uint64_t* size_val);

  /**
   * Gets the estimated result size (in bytes) for the input fixed-sized,
   * nullable attribute.
   */
  Status get_est_result_size_nullable(
      const char* name, uint64_t* size_val, uint64_t* size_validity);

  /**
   * Gets the estimated result size (in bytes) for the input var-sized,
   * nullable attribute.
   */
  Status get_est_result_size_nullable(
      const char* name,
      uint64_t* size_off,
      uint64_t* size_val,
      uint64_t* size_validity);

  /** Retrieves the number of written fragments. */
  Status get_written_fragment_num(uint32_t* num) const;

  /** Retrieves the URI of the written fragment with the input index. */
  Status get_written_fragment_uri(uint32_t idx, const char** uri) const;

  /**
   * Retrieves the timestamp range [t1,t2] of the written fragment with the
   * input index.
   */
  Status get_written_fragment_timestamp_range(
      uint32_t idx, uint64_t* t1, uint64_t* t2) const;

  /** Returns the array's smart pointer. */
  inline shared_ptr<Array> array_shared() {
    return array_shared_;
  }

  /** Returns the array. */
  const Array* array() const;

  /** Returns the array. */
  Array* array();

  /** Returns the array schema. */
  const ArraySchema& array_schema() const;

  /** Returns the array schema as a shared_ptr */
  const std::shared_ptr<const ArraySchema> array_schema_shared() const;

  /** Returns the names of the buffers set by the user for the query. */
  std::vector<std::string> buffer_names() const;

  /**
   * Gets the query buffer for the input attribute/dimension.
   * An empty string means the special default attribute.
   */
  QueryBuffer buffer(const std::string& name) const;

  /**
   * Marks a query that has not yet been started as failed. This should not be
   * called asynchronously to cancel an in-progress query; for that use the
   * parent StorageManager's cancellation mechanism.
   * @return Status
   */
  Status cancel();

  /**
   * Finalizes the query, flushing all internal state. Applicable only to global
   * layout writes. It has no effect for any other query type.
   */
  Status finalize();

  /**
   * Submits and finalizes a query, flushing all internal state. Applicable
   * only to global layout writes, returns an error otherwise.
   */
  Status submit_and_finalize();

  /**
   * This is a deprecated API.
   * Retrieves the buffer of a fixed-sized attribute/dimension.
   *
   * @param name The buffer attribute/dimension name. An empty string means
   *     the special default attribute/dimension.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size A pointer to the buffer size to be retrieved.
   * @return Status
   */
  Status get_buffer(
      const char* name, void** buffer, uint64_t** buffer_size) const;

  /**
   * This is a deprecated API.
   * Retrieves the offsets and values buffers of a var-sized
   * attribute/dimension.
   *
   * @param name The attribute/dimension name. An empty string means
   *     the special default attribute/dimension.
   * @param buffer_off The offsets buffer to be retrieved.
   * @param buffer_off_size A pointer to the offsets buffer size to be
   * retrieved.
   * @param buffer_val The values buffer to be retrieved.
   * @param buffer_val_size A pointer to the values buffer size to be retrieved.
   * @return Status
   */
  Status get_buffer(
      const char* name,
      uint64_t** buffer_off,
      uint64_t** buffer_off_size,
      void** buffer_val,
      uint64_t** buffer_val_size) const;

  /**
   * Retrieves the data buffer of a fixed/var-sized attribute/dimension.
   *
   * @param name The buffer attribute/dimension name. An empty string means
   *     the special default attribute/dimension.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size A pointer to the buffer size to be retrieved.
   * @return Status
   */
  Status get_data_buffer(
      const char* name, void** buffer, uint64_t** buffer_size) const;

  /**
   * Retrieves the data buffer of a fixed or variable-sized dimension label.
   *
   * @param name The name of the the label.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size The size of the buffer.
   */
  void get_label_data_buffer(
      const std::string& name, void** buffer, uint64_t** buffer_size) const;

  /**
   * Retrieves the offset buffer for a variable-sized dimension label.
   *
   * @param name The name of the label.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size The size of the buffer.
   */
  void get_label_offsets_buffer(
      const std::string& name, uint64_t** buffer, uint64_t** buffer_size) const;

  /**
   * Retrieves the offset buffer for a var-sized attribute/dimension.
   *
   * @param name The buffer attribute/dimension name. An empty string means
   * the special default attribute/dimension.
   * @param buffer_off The offsets buffer to be retrieved. This buffer holds
   * the starting offsets of each cell value in the data buffer.
   * @param buffer_off_size A pointer to the buffer size to be retrieved.
   * @return Status
   */
  Status get_offsets_buffer(
      const char* name,
      uint64_t** buffer_off,
      uint64_t** buffer_off_size) const;

  /**
   * Retrieves the validity buffer for a nullable attribute/dimension.
   *
   * @param name The buffer attribute/dimension name. An empty string means
   * the special default attribute/dimension.
   * @param buffer_validity_bytemap The buffer that either have the validity
   * bytemap associated with the input data to be written, or will hold the
   * validity bytemap to be read.
   * @param buffer_validity_bytemap_size In the case of writes, this is the size
   * of `buffer_validity_bytemap` in bytes. In the case of reads, this initially
   * contains the allocated size of `buffer_validity_bytemap`, but after the
   * termination of the query it will contain the size of the useful (read)
   * data in `buffer_validity_bytemap`.
   * @return Status
   */
  Status get_validity_buffer(
      const char* name,
      uint8_t** buffer_validity_bytemap,
      uint64_t** buffer_validity_bytemap_size) const;

  /**
   * This is a deprecated API.
   * Retrieves the buffer and validity bytemap of a fixed-sized, nullable
   * attribute.
   *
   * @param name The buffer attribute name. An empty string means
   *     the special default attribute.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size A pointer to the buffer size to be retrieved.
   * @param buffer The buffer to be retrieved.
   * @param buffer_size A pointer to the buffer size to be retrieved.
   * @return Status
   */
  Status get_buffer_vbytemap(
      const char* name,
      void** buffer,
      uint64_t** buffer_size,
      uint8_t** buffer_validity_bytemap,
      uint64_t** buffer_validity_bytemap_size) const;

  /**
   * This is a deprecated API.
   * Retrieves the offsets, values, and validity bytemap buffers of
   * a var-sized, nullable attribute.
   *
   * @param name The attribute name. An empty string means
   *     the special default attribute.
   * @param buffer_off The offsets buffer to be retrieved.
   * @param buffer_off_size A pointer to the offsets buffer size to be
   * retrieved.
   * @param buffer_val The values buffer to be retrieved.
   * @param buffer_val_size A pointer to the values buffer size to be retrieved.
   * @return Status
   */
  Status get_buffer_vbytemap(
      const char* name,
      uint64_t** buffer_off,
      uint64_t** buffer_off_size,
      void** buffer_val,
      uint64_t** buffer_val_size,
      uint8_t** buffer_validity_bytemap,
      uint64_t** buffer_validity_bytemap_size) const;

  /**
   * Returns the serialization state associated with the given attribute.
   *
   * @param attribute Attribute to get serialization state for
   * @param state Set to point to the serialization state
   * @return Status
   */
  Status get_attr_serialization_state(
      const std::string& attribute, SerializationState::AttrState** state);

  /**
   * Used by serialization to get the map of result sizes
   * @return
   */
  std::unordered_map<std::string, Subarray::ResultSize>
  get_est_result_size_map();

  /**
   * Used by serialization to get the map of max mem sizes
   * @return
   */
  std::unordered_map<std::string, Subarray::MemorySize> get_max_mem_size_map();

  /**
   * Returns `true` if the query has results. Applicable only to read
   * queries (it returns `false` for write queries).
   */
  bool has_results() const;

  /** Initializes the query. */
  Status init();

  /** Returns the first fragment uri. */
  URI first_fragment_uri() const;

  /** Returns the last fragment uri. */
  URI last_fragment_uri() const;

  /** Returns the cell layout. */
  Layout layout() const;

  /**
   * Returns the condition for filtering results in a read query.
   * @return QueryCondition
   */
  const QueryCondition& condition() const;

  /**
   * Returns the update values for an update query.
   * @return UpdateValues
   */
  const std::vector<UpdateValue>& update_values() const;

  /**
   * Returns true if this query requires the use of dimension labels.
   */
  bool uses_dimension_labels() const;

  /** Processes a query. */
  Status process();

  /** Gets the strategy of the query. */
  IQueryStrategy* strategy(bool skip_checks_serialization = false);

  /**
   * Switch the strategy depending on layout. Used by serialization.
   *
   * @param layout New layout
   * @param force_legacy_reader Force use of the legacy reader if the client
   *    requested it.
   * @return Status
   */
  Status reset_strategy_with_layout(Layout layout, bool force_legacy_reader);

  /**
   * Disables checking the global order and coordinate duplicates. Applicable
   * only to writes. This option will supercede the config.
   */
  Status disable_checks_consolidation();

  /**
   * Enables consolidation with timestamps.
   */
  Status set_consolidation_with_timestamps();

  /**
   * Set the processed conditions for writes.
   *
   * @param processed_conditions The processed conditions.
   */
  void set_processed_conditions(std::vector<std::string>& processed_conditions);

  /**
   * Sets the config for the Query
   *
   * This function overrides the config for Query-level parameters only.
   * Semantically, the initial query config is copied from the context
   * config upon initialization. Note that The context config is immutable
   * at the C API level because tiledb_ctx_get_config always returns a copy.
   *
   * Config parameters set here will *only* be applied within the Query.
   */
  Status set_config(const Config& config);

  /**
   * Sets the (zipped) coordinates buffer (set with TILEDB_COORDS as the
   * buffer name).
   *
   * @param buffer The buffer that has the input data to be written.
   * @param buffer_size The size of `buffer` in bytes.
   * @return Status
   */
  Status set_coords_buffer(void* buffer, uint64_t* buffer_size);

  /**
   * Sets the data for a fixed/var-sized attribute/dimension.
   *
   * @param name The attribute/dimension to set the buffer for.
   * @param buffer The buffer that will hold the data to be read.
   * @param buffer_size This initially contains the allocated
   *     size of `buffer`, but after the termination of the function
   *     it will contain the size of the useful (read) data in `buffer`.
   * @param check_null_buffers If true (default), null buffers are not
   * allowed.
   * @return Status
   */
  Status set_data_buffer(
      const std::string& name,
      void* const buffer,
      uint64_t* const buffer_size,
      const bool check_null_buffers = true);

  /**
   * Wrapper to set the internal buffer for a dimension or attribute from a
   * QueryBuffer.
   *
   * This function is intended to be a convenience method for use setting
   * buffers for dimension labels.
   *
   * @WARNING Does not check for or copy validity data.
   *
   * @param name Name of the dimension or attribute to set the buffer for.
   * @param buffer The query buffer to get the data from.
   **/
  void set_dimension_label_buffer(
      const std::string& name, const QueryBuffer& buffer);

  /**
   * Sets the label data buffer for fixed or variable sized dimension labels.
   *
   * @param name The name of the dimension label to set the data buffer for.
   * @param buffer The buffer for the data.
   * @param buffer_size The size of the data.
   * @param check_null_buffers If ``true`` verify buffer and buffersize are not
   *     nullptrs.
   */
  void set_label_data_buffer(
      const std::string& name,
      void* const buffer,
      uint64_t* const buffer_size,
      const bool check_null_buffers = true);

  /**
   * Sets the label offsets buffer for variable sized dimension labels.
   *
   * @param name The name of the dimension label to set the offsets buffer for.
   * @param buffer_offsets The buffer for offsets for the variable size data.
   * @param buffer_offsets_size The size of the offsets buffer.
   * @param check_null_buffers If ``true`` verify buffer and buffersize are not
   *     nullptrs.
   */
  void set_label_offsets_buffer(
      const std::string& name,
      uint64_t* const buffer_offsets,
      uint64_t* const buffer_offsets_size,
      const bool check_null_buffers = true);

  /**
   * Sets the offset buffer for a var-sized attribute/dimension.
   *
   * @param name The attribute/dimension to set the buffer for.
   * @param buffer_offsets The buffer that will hold the data to be read.
   *     This buffer holds the starting offsets of each cell value in
   *     `buffer_val`.
   * @param buffer_offsets_size This initially contains
   *     the allocated size of `buffer_off`, but after the termination of the
   *     function it will contain the size of the useful (read) data in
   *     `buffer_off`.
   * @param check_null_buffers If true (default), null buffers are not
   * allowed.
   * @return Status
   */
  Status set_offsets_buffer(
      const std::string& name,
      uint64_t* const buffer_offsets,
      uint64_t* const buffer_offsets_size,
      const bool check_null_buffers = true);

  /**
   * Sets the validity buffer for nullable attribute/dimension.
   *
   * @param name The attribute/dimension to set the buffer for.
   * @param buffer_validity_bytemap The buffer that either have the validity
   * bytemap associated with the input data to be written, or will hold the
   * validity bytemap to be read.
   * @param buffer_validity_bytemap_size In the case of writes, this is the size
   * of `buffer_validity_bytemap` in bytes. In the case of reads, this initially
   * contains the allocated size of `buffer_validity_bytemap`, but after the
   * termination of the query it will contain the size of the useful (read)
   * data in `buffer_validity_bytemap`.
   * @param check_null_buffers If true (default), null buffers are not
   * allowed.
   * @return Status
   */
  Status set_validity_buffer(
      const std::string& name,
      uint8_t* const buffer_validity_bytemap,
      uint64_t* const buffer_validity_bytemap_size,
      const bool check_null_buffers = true);

  /**
   * Get the config of the query.
   *
   * @return Config from query
   */
  const Config& config() const;

  /**
   * This is a deprecated API.
   * Sets the buffer for a fixed-sized attribute/dimension.
   *
   * @param name The attribute/dimension to set the buffer for.
   * @param buffer The buffer that either have the input data to be written,
   *     or will hold the data to be read.
   * @param buffer_size In the case of writes, this is the size of `buffer`
   *     in bytes. In the case of reads, this initially contains the allocated
   *     size of `buffer`, but after the termination of the query
   *     it will contain the size of the useful (read) data in `buffer`.
   * @param check_null_buffers If true (default), null buffers are not allowed.
   * @return Status
   */
  Status set_buffer(
      const std::string& name,
      void* buffer,
      uint64_t* buffer_size,
      bool check_null_buffers = true);

  /**
   * This is a deprecated API.
   * Sets the buffer for a var-sized attribute/dimension.
   *
   * @param name The attribute/dimension to set the buffer for.
   * @param buffer_off The buffer that either have the input data to be written,
   *     or will hold the data to be read. This buffer holds the starting
   *     offsets of each cell value in `buffer_val`.
   * @param buffer_off_size In the case of writes, it is the size of
   *     `buffer_off` in bytes. In the case of reads, this initially contains
   *     the allocated size of `buffer_off`, but after the termination of the
   *     function it will contain the size of the useful (read) data in
   *     `buffer_off`.
   * @param buffer_val The buffer that either have the input data to be written,
   *     or will hold the data to be read. This buffer holds the actual
   *     var-sized cell values.
   * @param buffer_val_size In the case of writes, it is the size of
   *     `buffer_val` in bytes. In the case of reads, this initially contains
   *     the allocated size of `buffer_val`, but after the termination of the
   *     query it will contain the size of the useful (read) data in
   *     `buffer_val`.
   * @param check_null_buffers If true (default), null buffers are not allowed.
   * @return Status
   */
  Status set_buffer(
      const std::string& name,
      uint64_t* buffer_off,
      uint64_t* buffer_off_size,
      void* buffer_val,
      uint64_t* buffer_val_size,
      bool check_null_buffers = true);

  /**
   * This is a deprecated API.
   * Sets the buffer for a fixed-sized, nullable attribute with a validity
   * bytemap.
   *
   * @param name The attribute to set the buffer for.
   * @param buffer The buffer that either have the input data to be written,
   *     or will hold the data to be read.
   * @param buffer_size In the case of writes, this is the size of `buffer`
   *     in bytes. In the case of reads, this initially contains the allocated
   *     size of `buffer`, but after the termination of the query
   *     it will contain the size of the useful (read) data in `buffer`.
   * @param buffer_validity_bytemap The buffer that either have the validity
   * bytemap associated with the input data to be written, or will hold the
   * validity bytemap to be read.
   * @param buffer_validity_bytemap_size In the case of writes, this is the size
   * of `buffer_validity_bytemap` in bytes. In the case of reads, this initially
   *     contains the allocated size of `buffer_validity_bytemap`, but after the
   *     termination of the query it will contain the size of the useful (read)
   * data in `buffer_validity_bytemap`.
   * @param check_null_buffers If true (default), null buffers are not allowed.
   * @return Status
   */
  Status set_buffer_vbytemap(
      const std::string& name,
      void* buffer,
      uint64_t* buffer_size,
      uint8_t* buffer_validity_bytemap,
      uint64_t* buffer_validity_bytemap_size,
      bool check_null_buffers = true);

  /**
   * This is a deprecated API.
   * Sets the buffer for a var-sized, nullable attribute with a validity
   * bytemap.
   *
   * @param name The attribute to set the buffer for.
   * @param buffer_off The buffer that either have the input data to be written,
   *     or will hold the data to be read. This buffer holds the starting
   *     offsets of each cell value in `buffer_val`.
   * @param buffer_off_size In the case of writes, it is the size of
   *     `buffer_off` in bytes. In the case of reads, this initially contains
   *     the allocated size of `buffer_off`, but after the termination of the
   *     function it will contain the size of the useful (read) data in
   *     `buffer_off`.
   * @param buffer_val The buffer that either have the input data to be written,
   *     or will hold the data to be read. This buffer holds the actual
   *     var-sized cell values.
   * @param buffer_val_size In the case of writes, it is the size of
   *     `buffer_val` in bytes. In the case of reads, this initially contains
   *     the allocated size of `buffer_val`, but after the termination of the
   *     query it will contain the size of the useful (read) data in
   *     `buffer_val`.
   * @param buffer_validity_bytemap The buffer that either have the validity
   * bytemap associated with the input data to be written, or will hold the
   * validity bytemap to be read.
   * @param buffer_validity_bytemap_size In the case of writes, this is the size
   * of `buffer_validity_bytemap` in bytes. In the case of reads, this initially
   *     contains the allocated size of `buffer_validity_bytemap`, but after the
   *     termination of the query it will contain the size of the useful (read)
   * data in `buffer_validity_bytemap`.
   * @param check_null_buffers If true (default), null buffers are not allowed.
   * @return Status
   */
  Status set_buffer_vbytemap(
      const std::string& name,
      uint64_t* buffer_off,
      uint64_t* buffer_off_size,
      void* buffer_val,
      uint64_t* buffer_val_size,
      uint8_t* buffer_validity_bytemap,
      uint64_t* buffer_validity_bytemap_size,
      bool check_null_buffers = true);

  /**
   * Used by serialization to set the estimated result size
   *
   * @param est_result_size map to set
   * @param max_mem_size map to set
   * @return Status
   */
  Status set_est_result_size(
      std::unordered_map<std::string, Subarray::ResultSize>& est_result_size,
      std::unordered_map<std::string, Subarray::MemorySize>& max_mem_size);

  /**
   * Sets the cell layout of the query.
   */
  Status set_layout(Layout layout);

  /**
   * Sets the condition for filtering results in a read query.
   *
   * @param condition The condition object.
   * @return Status
   */
  Status set_condition(const QueryCondition& condition);

  /**
   * Adds an update value for an update query.
   *
   * @param field_name The attribute name.
   * @param update_value The value to set.
   * @param update_value_size The byte size of `update_value`.
   * @return Status
   */
  Status add_update_value(
      const char* field_name,
      const void* update_value,
      uint64_t update_value_size);

  /**
   * Adds ranges to a query initialize with label ranges.
   *
   * @param dim_idx The dimension to update.
   * @param is_point_ranges If ``true`` the data contains point ranges.
   *     Otherwise, it contains standard ranges.
   * @param start Pointer to the start of the range data.
   * @param count Number of total elements in the range data.
   */
  void add_index_ranges_from_label(
      uint32_t dim_idx,
      const bool is_point_range,
      const void* start,
      const uint64_t count);

  /**
   * Set query status, needed for json deserialization
   * @param status
   * @return Status
   */
  void set_status(QueryStatus status);

  /**
   * Sets the query subarray. If it is null, then the subarray will be set to
   * the entire domain.
   *
   * @param subarray The subarray to be set.
   * @return Status
   *
   * @note Setting a subarray for sparse arrays, or for dense arrays
   *     when performing unordered (sparse) writes, has no effect
   *     (will be ingnored).
   */
  Status set_subarray(const void* subarray);

  /** Returns the query subarray. */
  const Subarray* subarray() const;

  /**
   * Sets the query subarray.
   *
   * @param subarray The subarray to be set.
   * @return Status
   *
   * @note Calling set_subarray for sparse arrays, or for dense arrays
   *     when performing unordered (sparse) writes, has no effect.
   */
  Status set_subarray(const tiledb::sm::Subarray& subarray);

  /** Sets the query subarray, without performing any checks. */
  Status set_subarray_unsafe(const Subarray& subarray);

  /** Sets the query subarray, without performing any checks. */
  Status set_subarray_unsafe(const NDRange& subarray);

  /** Submits the query to the storage manager. */
  Status submit();

  /**
   * Submits the query to the storage manager. The query will be
   * processed asynchronously (i.e., in a non-blocking manner).
   * Once the query is completed, the input callback function will
   * be executed using the input callback data.
   */
  Status submit_async(std::function<void(void*)> callback, void* callback_data);

  /** Returns the query status. */
  QueryStatus status() const;

  /** Returns the query status incomplete reason. */
  QueryStatusDetailsReason status_incomplete_reason() const;

  /** Returns the query type. */
  QueryType type() const;

  /** Returns the internal stats object. */
  stats::Stats* stats() const;

  /** Returns the scratch space used for REST requests. */
  shared_ptr<Buffer> rest_scratch() const;

  /** Use the refactored dense reader or not. */
  bool use_refactored_dense_reader(
      const ArraySchema& array_schema, bool all_dense);

  /** Use the refactored sparse global order reader or not. */
  bool use_refactored_sparse_global_order_reader(
      Layout layout, const ArraySchema& array_schema);

  /** Use the refactored sparse unordered with dups reader or not. */
  bool use_refactored_sparse_unordered_with_dups_reader(
      Layout layout, const ArraySchema& array_schema);

  /** Returns if all ranges for this query are non overlapping. */
  tuple<Status, optional<bool>> non_overlapping_ranges();

  /** Returns true if this is a dense query */
  bool is_dense() const;

  /** Returns a reference to the internal WrittenFragmentInfo list */
  std::vector<WrittenFragmentInfo>& get_written_fragment_info();

  /** Called from serialization to mark the query as remote */
  void set_remote_query();

  /**
   * Set a flag to specify we are doing an ordered dimension label read.
   *
   * @param increasing_order Is the query on an array with increasing order? If
   * not assume decreasing order.
   */
  void set_dimension_label_ordered_read(bool increasing_order);

  /**
   * Set the fragment size.
   *
   * @param fragment_size Fragment size.
   */
  void set_fragment_size(uint64_t fragment_size) {
    fragment_size_ = fragment_size;
  }

 private:
  /* ********************************* */
  /*         PRIVATE ATTRIBUTES        */
  /* ********************************* */

  /** A smart pointer to the array the query is associated with.
   * Ensures that the Array object exists as long as the Query object exists. */
  shared_ptr<Array> array_shared_;

  /** The array the query is associated with.
   * Cached copy of array_shared_.get(). */
  Array* array_;

  /** The array schema. */
  shared_ptr<const ArraySchema> array_schema_;

  /** The config for query-level parameters only. */
  Config config_;

  /** A function that will be called upon the completion of an async query. */
  std::function<void(void*)> callback_;

  /** The data input to the callback function. */
  void* callback_data_;

  /** The query type. */
  QueryType type_;

  /** The layout of the cells in the result of the subarray. */
  Layout layout_;

  /** The query subarray (initially the whole domain by default). */
  Subarray subarray_;

  /** The query status. */
  QueryStatus status_;

  /** The storage manager. */
  StorageManager* storage_manager_;

  /** The query strategy. */
  tdb_unique_ptr<IQueryStrategy> strategy_;

  /** The class stats. */
  stats::Stats* stats_;

  /** The class logger. */
  shared_ptr<Logger> logger_;

  /** UID of the logger instance */
  inline static std::atomic<uint64_t> logger_id_ = 0;

  /**
   * Maps attribute/dimension names to their buffers.
   * `TILEDB_COORDS` may be used for the special zipped coordinates
   * buffer.
   * */
  std::unordered_map<std::string, QueryBuffer> buffers_;

  /** Maps label names to their buffers. */
  std::unordered_map<std::string, QueryBuffer> label_buffers_;

  /** Dimension label queries that are part of the main query. */
  tdb_unique_ptr<ArrayDimensionLabelQueries> dim_label_queries_;

  /** Keeps track of the coords data. */
  CoordsInfo coords_info_;

  /** Stores information about the written fragments. */
  std::vector<WrittenFragmentInfo> written_fragment_info_;

  /** The query condition. */
  QueryCondition condition_;

  /** The update values. */
  std::vector<UpdateValue> update_values_;

  /** Set of attributes that have an update value. */
  std::set<std::string> attributes_with_update_value_;

  /** The fragment metadata that this query will focus on. */
  std::vector<shared_ptr<FragmentMetadata>> fragment_metadata_;

  /** The current serialization state. */
  SerializationState serialization_state_;

  /** If the query has coords buffer set or not. */
  bool has_coords_buffer_;

  /** If the query has zipped coords buffer set or not. */
  bool has_zipped_coords_buffer_;

  /** True if at least one separate coordinate buffer is set. */
  bool coord_buffer_is_set_;

  /** True if at least one separate data coordinate buffer is set. */
  bool coord_data_buffer_is_set_;

  /** True if at least one separate offsets coordinate buffer is set. */
  bool coord_offsets_buffer_is_set_;

  /** Keeps track of the name of the data buffer once set. */
  std::string data_buffer_name_;

  /** Keeps track of the name of the offsets buffer once set. */
  std::string offsets_buffer_name_;

  /**
   * If `true`, it will not check if the written coordinates are
   * in the global order or have duplicates. This supercedes the config.
   */
  bool disable_checks_consolidation_;

  /**
   * If `true`, it will enable consolidation with timestamps on the reader.
   */
  bool consolidation_with_timestamps_;

  /* Scratch space used for REST requests. */
  shared_ptr<Buffer> rest_scratch_;

  /* Processed conditions, used for consolidation. */
  std::vector<std::string> processed_conditions_;

  /**
   * Flag to force legacy reader when strategy gets created. This is used by
   * the serialization codebase if a query comes from an older version of the
   * library that doesn't have the refactored readers, we need to run it with
   * the legacy reader.
   */
  bool force_legacy_reader_;

  /**
   * The name of the new fragment to be created for writes.
   *
   * If not set, the fragment name will be created using the latest array
   * timestamp and a generated UUID.
   */
  optional<std::string> fragment_name_;

  /** It tracks if this is a remote query */
  bool remote_query_;

  /** Flag to specify we are doing a dimension label ordered read. */
  bool is_dimension_label_ordered_read_;

  /**
   * Is the dimension label ordered read on an array with increasing order? If
   * not assume decreasing order.
   *
   * NOTE: Only used when is_dimension_label_order_read_ == true.
   */
  bool dimension_label_increasing_;

  /**
   * The desired fragment size. The writer will create a new fragment once
   * this size has been reached.
   *
   * Note: This is only used for global order writes.
   */
  uint64_t fragment_size_;

  /* ********************************* */
  /*           PRIVATE METHODS         */
  /* ********************************* */

  /**
   * Create the strategy.
   *
   * @param skip_checks_serialization Skip checks during serialization.
   */
  Status create_strategy(bool skip_checks_serialization = false);

  Status check_set_fixed_buffer(const std::string& name);

  /** Checks if the buffers names have been appropriately set for the query. */
  Status check_buffer_names();

  /**
   * Internal routine for checking the completeness of all attribute
   * and dimensions buffers. Iteratively searches that all attributes &
   * dimenstions buffers have been set correctly
   * @return Status
   */
  Status check_buffers_correctness();

  /**
   * Returns true if only querying dimension labels.
   *
   * The query will only query dimension labels if all the following are true:
   * 1. At most one dimension buffer is set.
   * 2. No attribute buffers are set.
   * 3. At least one label buffer is set.
   */
  bool only_dim_label_query() const;

  /**
   * This is a deprecated API.
   * Internal routine for setting fixed-sized, nullable attribute buffers with
   * a ValidityVector.
   */
  Status set_buffer(
      const std::string& name,
      void* buffer,
      uint64_t* buffer_size,
      ValidityVector&& validity_vector,
      bool check_null_buffers = true);

  /**
   * This is a deprecated API.
   * Internal routine for setting var-sized, nullable attribute buffers with
   * a ValidityVector.
   */
  Status set_buffer(
      const std::string& name,
      uint64_t* buffer_off,
      uint64_t* buffer_off_size,
      void* buffer_val,
      uint64_t* buffer_val_size,
      ValidityVector&& validity_vector,
      bool check_null_buffers = true);

  /**
   * This is a deprecated API.
   * Internal routine for getting fixed-sized, nullable attribute buffers with
   * a ValidityVector.
   */
  Status get_buffer(
      const char* name,
      void** buffer,
      uint64_t** buffer_size,
      const ValidityVector** validity_vector) const;

  /**
   * This is a deprecated API.
   * Internal routine for getting fixed-sized, nullable attribute buffers with
   * a ValidityVector.
   */
  Status get_buffer(
      const char* name,
      uint64_t** buffer_off,
      uint64_t** buffer_off_size,
      void** buffer_val,
      uint64_t** buffer_val_size,
      const ValidityVector** validity_vector) const;

  /**
   * Check if input buffers are tile aligned. This function should be called
   * only for remote global order writes and it should enforce tile alignment
   * for both dense and sparse arrays.
   */
  Status check_tile_alignment() const;

  /**
   * Check if input buffers are bigger than 5MB. S3 multipart upload
   * requires each part be bigger than 5MB, except the last part.
   * This function should be called only for remote global order writes.
   */
  Status check_buffer_multipart_size() const;

  /**
   * Reset coord buffer markers at end of a global write submit.
   * This will allow for the user to properly set the next write batch.
   */
  void reset_coords_markers();
};

}  // namespace sm
}  // namespace tiledb

#endif  // TILEDB_QUERY_H
