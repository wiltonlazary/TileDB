/**
 * @file   unit-capi-incomplete.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2021 TileDB Inc.
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
 * Tests the C API async queries.
 */

#include <test/support/tdb_catch.h>
#include "test/support/src/helpers.h"
#include "tiledb/sm/c_api/tiledb.h"
#include "tiledb/sm/c_api/tiledb_serialization.h"
#include "tiledb/sm/serialization/query.h"

#include <cstring>
#include <iostream>

using namespace tiledb::test;

/**
 * Tests cases where a read query is incomplete or leads to a buffer
 * overflow.
 */

struct IncompleteFx {
  // Constants
  const char* DENSE_ARRAY_NAME = "test_async_dense";
  const char* SPARSE_ARRAY_NAME = "test_async_sparse";

  // TileDB context
  tiledb_ctx_t* ctx_;

  // Constructors/destructors
  IncompleteFx();
  ~IncompleteFx();

  // Functions
  void create_dense_array();
  void create_sparse_array();
  void write_dense_full(const bool serialized_writes);
  void write_sparse_full(const bool serialized_writes);
  void check_dense_incomplete();
  void check_dense_incomplete_serialized();
  void check_dense_until_complete();
  void check_dense_shrink_buffer_size();
  void check_dense_unsplittable_overflow();
  void check_dense_unsplittable_complete();
  void check_dense_reset_buffers();
  void check_sparse_incomplete();
  void check_sparse_until_complete();
  void check_sparse_unsplittable_overflow();
  void check_sparse_unsplittable_complete();
  void remove_dense_array();
  void remove_sparse_array();
  void remove_array(const std::string& array_name);
  bool is_array(const std::string& array_name);
};

IncompleteFx::IncompleteFx() {
  ctx_ = nullptr;
  REQUIRE(tiledb_ctx_alloc(nullptr, &ctx_) == TILEDB_OK);
}

IncompleteFx::~IncompleteFx() {
  tiledb_ctx_free(&ctx_);
}

void IncompleteFx::create_dense_array() {
  // Create dimensions
  uint64_t dim_domain[] = {1, 4, 1, 4};
  uint64_t tile_extents[] = {2, 2};
  tiledb_dimension_t* d1;
  int rc = tiledb_dimension_alloc(
      ctx_, "d1", TILEDB_UINT64, &dim_domain[0], &tile_extents[0], &d1);
  CHECK(rc == TILEDB_OK);
  tiledb_dimension_t* d2;
  rc = tiledb_dimension_alloc(
      ctx_, "d2", TILEDB_UINT64, &dim_domain[2], &tile_extents[1], &d2);
  CHECK(rc == TILEDB_OK);

  // Create domain
  tiledb_domain_t* domain;
  rc = tiledb_domain_alloc(ctx_, &domain);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_domain_add_dimension(ctx_, domain, d1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_domain_add_dimension(ctx_, domain, d2);
  CHECK(rc == TILEDB_OK);

  // Create attributes
  tiledb_attribute_t* a1;
  rc = tiledb_attribute_alloc(ctx_, "a1", TILEDB_INT32, &a1);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a1, TILEDB_FILTER_LZ4, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a1, 1);
  CHECK(rc == TILEDB_OK);
  tiledb_attribute_t* a2;
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_alloc(ctx_, "a2", TILEDB_CHAR, &a2);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a2, TILEDB_FILTER_GZIP, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a2, TILEDB_VAR_NUM);
  CHECK(rc == TILEDB_OK);
  tiledb_attribute_t* a3;
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_alloc(ctx_, "a3", TILEDB_FLOAT32, &a3);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a3, TILEDB_FILTER_ZSTD, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a3, 2);
  CHECK(rc == TILEDB_OK);

  // Create array schema
  tiledb_array_schema_t* array_schema;
  rc = tiledb_array_schema_alloc(ctx_, TILEDB_DENSE, &array_schema);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_cell_order(ctx_, array_schema, TILEDB_ROW_MAJOR);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_tile_order(ctx_, array_schema, TILEDB_ROW_MAJOR);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_domain(ctx_, array_schema, domain);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a2);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a3);
  CHECK(rc == TILEDB_OK);

  // Check array schema
  rc = tiledb_array_schema_check(ctx_, array_schema);
  CHECK(rc == TILEDB_OK);

  // Create array
  rc = tiledb_array_create(ctx_, DENSE_ARRAY_NAME, array_schema);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_attribute_free(&a1);
  tiledb_attribute_free(&a2);
  tiledb_attribute_free(&a3);
  tiledb_dimension_free(&d1);
  tiledb_dimension_free(&d2);
  tiledb_domain_free(&domain);
  tiledb_array_schema_free(&array_schema);
}

void IncompleteFx::create_sparse_array() {
  // Create dimensions
  uint64_t dim_domain[] = {1, 4, 1, 4};
  uint64_t tile_extents[] = {2, 2};
  tiledb_dimension_t* d1;
  int rc = tiledb_dimension_alloc(
      ctx_, "d1", TILEDB_UINT64, &dim_domain[0], &tile_extents[0], &d1);
  CHECK(rc == TILEDB_OK);
  tiledb_dimension_t* d2;
  rc = tiledb_dimension_alloc(
      ctx_, "d2", TILEDB_UINT64, &dim_domain[2], &tile_extents[1], &d2);
  CHECK(rc == TILEDB_OK);

  // Create domain
  tiledb_domain_t* domain;
  rc = tiledb_domain_alloc(ctx_, &domain);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_domain_add_dimension(ctx_, domain, d1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_domain_add_dimension(ctx_, domain, d2);
  CHECK(rc == TILEDB_OK);

  // Create attributes
  tiledb_attribute_t* a1;
  rc = tiledb_attribute_alloc(ctx_, "a1", TILEDB_INT32, &a1);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a1, TILEDB_FILTER_LZ4, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a1, 1);
  CHECK(rc == TILEDB_OK);
  tiledb_attribute_t* a2;
  rc = tiledb_attribute_alloc(ctx_, "a2", TILEDB_CHAR, &a2);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a2, TILEDB_FILTER_GZIP, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a2, TILEDB_VAR_NUM);
  CHECK(rc == TILEDB_OK);
  tiledb_attribute_t* a3;
  rc = tiledb_attribute_alloc(ctx_, "a3", TILEDB_FLOAT32, &a3);
  CHECK(rc == TILEDB_OK);
  rc = set_attribute_compression_filter(ctx_, a3, TILEDB_FILTER_ZSTD, -1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_attribute_set_cell_val_num(ctx_, a3, 2);
  CHECK(rc == TILEDB_OK);

  // Create array schmea
  tiledb_array_schema_t* array_schema;
  rc = tiledb_array_schema_alloc(ctx_, TILEDB_SPARSE, &array_schema);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_cell_order(ctx_, array_schema, TILEDB_ROW_MAJOR);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_tile_order(ctx_, array_schema, TILEDB_ROW_MAJOR);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_capacity(ctx_, array_schema, 2);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_set_domain(ctx_, array_schema, domain);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a1);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a2);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_schema_add_attribute(ctx_, array_schema, a3);
  CHECK(rc == TILEDB_OK);

  // Check array schema
  rc = tiledb_array_schema_check(ctx_, array_schema);
  CHECK(rc == TILEDB_OK);

  // Create array
  rc = tiledb_array_create(ctx_, SPARSE_ARRAY_NAME, array_schema);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_attribute_free(&a1);
  tiledb_attribute_free(&a2);
  tiledb_attribute_free(&a3);
  tiledb_dimension_free(&d1);
  tiledb_dimension_free(&d2);
  tiledb_domain_free(&domain);
  tiledb_array_schema_free(&array_schema);
}

void IncompleteFx::write_dense_full(const bool serialized_writes) {
  // Set attributes
  const char* attributes[] = {"a1", "a2", "a3"};

  // Prepare cell buffers
  // clang-format off
  int buffer_a1[] = {
      0,  1,  2,  3, 4,  5,  6,  7,
      8,  9,  10, 11, 12, 13, 14, 15
  };
  uint64_t buffer_a2[] = {
      0,  1,  3,  6, 10, 11, 13, 16,
      20, 21, 23, 26, 30, 31, 33, 36
  };
  char buffer_var_a2[] =
      "abbcccdddd"
      "effggghhhh"
      "ijjkkkllll"
      "mnnooopppp";
  float buffer_a3[] = {
      0.1f,  0.2f,  1.1f,  1.2f,  2.1f,  2.2f,  3.1f,  3.2f,
      4.1f,  4.2f,  5.1f,  5.2f,  6.1f,  6.2f,  7.1f,  7.2f,
      8.1f,  8.2f,  9.1f,  9.2f,  10.1f, 10.2f, 11.1f, 11.2f,
      12.1f, 12.2f, 13.1f, 13.2f, 14.1f, 14.2f, 15.1f, 15.2f,
  };
  void* buffers[] = { buffer_a1, buffer_a2, buffer_var_a2, buffer_a3 };
  uint64_t buffer_sizes[] =
  {
      sizeof(buffer_a1),
      sizeof(buffer_a2),
      sizeof(buffer_var_a2)-1,  // No need to store the last '\0' character
      sizeof(buffer_a3)
  };
  // clang-format on

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_WRITE);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_WRITE, &query);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[1], buffers[2], &buffer_sizes[2]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[1], (uint64_t*)buffers[1], &buffer_sizes[1]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[2], buffers[3], &buffer_sizes[3]);
  CHECK(rc == TILEDB_OK);

  if (!serialized_writes) {
    rc = tiledb_query_submit(ctx_, query);
    CHECK(rc == TILEDB_OK);
    rc = tiledb_query_finalize(ctx_, query);
    CHECK(rc == TILEDB_OK);
  } else {
    submit_and_finalize_serialized_query(ctx_, query);
  }

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::write_sparse_full(const bool serialized_writes) {
  // Prepare cell buffers
  int buffer_a1[] = {0, 1, 2, 3, 4, 5, 6, 7};
  uint64_t buffer_a2[] = {0, 1, 3, 6, 10, 11, 13, 16};
  char buffer_var_a2[] = "abbcccddddeffggghhhh";
  float buffer_a3[] = {0.1f,
                       0.2f,
                       1.1f,
                       1.2f,
                       2.1f,
                       2.2f,
                       3.1f,
                       3.2f,
                       4.1f,
                       4.2f,
                       5.1f,
                       5.2f,
                       6.1f,
                       6.2f,
                       7.1f,
                       7.2f};
  uint64_t buffer_coords_dim1[] = {1, 1, 1, 2, 3, 4, 3, 3};
  uint64_t buffer_coords_dim2[] = {1, 2, 4, 3, 1, 2, 3, 4};

  void* buffers[] = {buffer_a1,
                     buffer_a2,
                     buffer_var_a2,
                     buffer_a3,
                     buffer_coords_dim1,
                     buffer_coords_dim2};
  uint64_t buffer_sizes[] = {
      sizeof(buffer_a1),
      sizeof(buffer_a2),
      sizeof(buffer_var_a2) - 1,  // No need to store the last '\0' character
      sizeof(buffer_a3),
      sizeof(buffer_coords_dim1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, SPARSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_WRITE);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  const char* attributes[] = {"a1", "a2", "a3", "d1", "d2"};
  rc = tiledb_query_alloc(ctx_, array, TILEDB_WRITE, &query);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[1], buffers[2], &buffer_sizes[2]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[1], (uint64_t*)buffers[1], &buffer_sizes[1]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[2], buffers[3], &buffer_sizes[3]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[3], buffers[4], &buffer_sizes[4]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[4], buffers[5], &buffer_sizes[4]);
  CHECK(rc == TILEDB_OK);

  if (!serialized_writes) {
    rc = tiledb_query_submit(ctx_, query);
    CHECK(rc == TILEDB_OK);
    rc = tiledb_query_finalize(ctx_, query);
    CHECK(rc == TILEDB_OK);
  } else {
    submit_and_finalize_serialized_query(ctx_, query);
  }

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::remove_array(const std::string& array_name) {
  if (!is_array(array_name))
    return;

  CHECK(tiledb_object_remove(ctx_, array_name.c_str()) == TILEDB_OK);
}

void IncompleteFx::remove_dense_array() {
  remove_array(DENSE_ARRAY_NAME);
}

void IncompleteFx::remove_sparse_array() {
  remove_array(SPARSE_ARRAY_NAME);
}

bool IncompleteFx::is_array(const std::string& array_name) {
  tiledb_object_t type = TILEDB_INVALID;
  REQUIRE(tiledb_object_type(ctx_, array_name.c_str(), &type) == TILEDB_OK);
  return type == TILEDB_ARRAY;
}

void IncompleteFx::check_dense_incomplete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[2];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Free/finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);

  // Check buffer
  int c_buffer_a1[2] = {0, 1};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));
}

void IncompleteFx::check_dense_incomplete_serialized() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[2];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_ROW_MAJOR);
  REQUIRE(rc == TILEDB_OK);

  // Serialize/deserialize into new query using new array (client-side)
  tiledb_buffer_list_t* buff_list;
  REQUIRE(
      tiledb_serialize_query(ctx_, query, TILEDB_CAPNP, 1, &buff_list) ==
      TILEDB_OK);
  tiledb_buffer_t* buff;
  REQUIRE(tiledb_buffer_list_flatten(ctx_, buff_list, &buff) == TILEDB_OK);
  tiledb_array_t* array2;
  rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array2);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array2, TILEDB_READ);
  CHECK(rc == TILEDB_OK);
  tiledb_query_t* query2;
  rc = tiledb_query_alloc(ctx_, array2, TILEDB_READ, &query2);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_deserialize_query(ctx_, buff, TILEDB_CAPNP, 0, query2);
  REQUIRE(rc == TILEDB_OK);

  // Allocate and set buffers on the "server".
  int buffer_a1_server[2];
  void* buffers_server[] = {buffer_a1_server};
  uint64_t buffer_sizes_server[] = {sizeof(buffer_a1_server)};
  rc = tiledb_query_set_data_buffer(
      ctx_, query2, attributes[0], buffers_server[0], &buffer_sizes_server[0]);
  REQUIRE(rc == TILEDB_OK);

  // Submit new query
  rc = tiledb_query_submit(ctx_, query2);
  REQUIRE(rc == TILEDB_OK);

  // Deserialize back to the original query.
  tiledb_buffer_list_t* buff_list2;
  rc = tiledb_serialize_query(ctx_, query2, TILEDB_CAPNP, 0, &buff_list2);
  REQUIRE(rc == TILEDB_OK);
  tiledb_buffer_t* buff2;
  REQUIRE(tiledb_buffer_list_flatten(ctx_, buff_list2, &buff2) == TILEDB_OK);
  // Client-side
  rc = tiledb_deserialize_query(ctx_, buff2, TILEDB_CAPNP, 1, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(status == TILEDB_INCOMPLETE);

  // Close arrays
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_close(ctx_, array2);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_query_free(&query);
  tiledb_query_free(&query2);
  tiledb_array_free(&array);
  tiledb_array_free(&array2);
  tiledb_buffer_free(&buff);
  tiledb_buffer_free(&buff2);
  tiledb_buffer_list_free(&buff_list);
  tiledb_buffer_list_free(&buff_list2);

  // Check buffer
  int c_buffer_a1[2] = {0, 1};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));
}

void IncompleteFx::check_dense_until_complete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[2];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check buffer
  int c_buffer_a1[2] = {0, 1};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Resubmit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_COMPLETED);

  // Check new buffer contents
  int c_buffer_a1_2[2] = {2, 3};
  CHECK(!memcmp(buffer_a1, c_buffer_a1_2, sizeof(c_buffer_a1_2)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));

  // Free/finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_dense_shrink_buffer_size() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[2];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check buffer
  int c_buffer_a1[2] = {0, 1};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Shrink buffer size
  buffer_sizes[0] = sizeof(int);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  CHECK(rc == TILEDB_OK);

  // Resubmit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Check new buffer contents
  CHECK(buffer_sizes[0] == 4);
  CHECK(buffer_a1[0] == 2);

  // Free/finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_dense_unsplittable_overflow() {
  // Initialize a subarray
  const uint64_t subarray[] = {2, 2, 2, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a2"};

  // Prepare the buffers that will store the result
  uint64_t buffer_a2[1];
  char buffer_a2_var[1];
  void* buffers[] = {buffer_a2, buffer_a2_var};
  uint64_t buffer_sizes[] = {sizeof(buffer_a2), sizeof(buffer_a2_var)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[1], &buffer_sizes[1]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[0], (uint64_t*)buffers[0], &buffer_sizes[0]);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  CHECK(rc == TILEDB_OK);

  // Get status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  CHECK(buffer_sizes[0] == 0);
  CHECK(buffer_sizes[1] == 0);

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_dense_unsplittable_complete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 1, 2, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a2"};

  // Prepare the buffers that will store the result
  uint64_t buffer_a2[1];
  char buffer_a2_var[2];
  void* buffers[] = {buffer_a2, buffer_a2_var};
  uint64_t buffer_sizes[] = {sizeof(buffer_a2), sizeof(buffer_a2_var)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[1], &buffer_sizes[1]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[0], (uint64_t*)buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  CHECK(rc == TILEDB_OK);

  // Check buffers
  char c_buffer_a2_var[2] = {'b', 'b'};
  CHECK(!memcmp(buffer_a2_var, c_buffer_a2_var, sizeof(c_buffer_a2_var)));

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_dense_reset_buffers() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[2];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, DENSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check buffer
  int c_buffer_a1[2] = {0, 1};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Resubmit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_COMPLETED);

  // Check new buffer contents
  int c_buffer_a1_2[2] = {2, 3};
  CHECK(!memcmp(buffer_a1, c_buffer_a1_2, sizeof(c_buffer_a1_2)));
  CHECK(buffer_sizes[0] == 2 * sizeof(int));

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_sparse_incomplete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[1];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, SPARSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);

  // Check buffer
  int c_buffer_a1[1] = {0};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == sizeof(int));
}

void IncompleteFx::check_sparse_until_complete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 2, 1, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a1"};

  // Prepare the buffers that will store the result
  int buffer_a1[1];
  void* buffers[] = {buffer_a1};
  uint64_t buffer_sizes[] = {sizeof(buffer_a1)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, SPARSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);

  // Check buffer
  int c_buffer_a1[1] = {0};
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == sizeof(int));

  // Resubmit the query
  rc = tiledb_query_submit(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Check status
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(
      status == (use_refactored_sparse_global_order_reader() ?
                     TILEDB_COMPLETED :
                     TILEDB_INCOMPLETE));

  // Check buffer
  c_buffer_a1[0] = 1;
  CHECK(!memcmp(buffer_a1, c_buffer_a1, sizeof(c_buffer_a1)));
  CHECK(buffer_sizes[0] == sizeof(int));

  /**
   * Old reader needs an extra round here to finish processing all the
   * partitions in the subarray. New reader is done earlier.
   */
  if (!use_refactored_sparse_global_order_reader()) {
    // Submit query
    rc = tiledb_query_submit(ctx_, query);
    REQUIRE(rc == TILEDB_OK);

    // Check status
    rc = tiledb_query_get_status(ctx_, query, &status);
    CHECK(rc == TILEDB_OK);
    CHECK(status == TILEDB_COMPLETED);

    // Check buffer
    CHECK(buffer_sizes[0] == 0);
  }

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_sparse_unsplittable_overflow() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 1, 2, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a2"};

  // Prepare the buffers that will store the result
  uint64_t buffer_a2[1];
  char buffer_a2_var[1];
  void* buffers[] = {buffer_a2, buffer_a2_var};
  uint64_t buffer_sizes[] = {sizeof(buffer_a2), sizeof(buffer_a2_var)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, SPARSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[1], &buffer_sizes[1]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[0], (uint64_t*)buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);

  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  CHECK(rc == TILEDB_OK);
  tiledb_query_status_t status;
  rc = tiledb_query_get_status(ctx_, query, &status);
  CHECK(rc == TILEDB_OK);
  CHECK(status == TILEDB_INCOMPLETE);
  CHECK(buffer_sizes[0] == 0);

  tiledb_query_status_details_t details;
  rc = tiledb_query_get_status_details(ctx_, query, &details);
  CHECK(rc == TILEDB_OK);
  CHECK(details.incomplete_reason == TILEDB_REASON_USER_BUFFER_SIZE);

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

void IncompleteFx::check_sparse_unsplittable_complete() {
  // Initialize a subarray
  const uint64_t subarray[] = {1, 1, 2, 2};

  // Subset over a specific attribute
  const char* attributes[] = {"a2"};

  // Prepare the buffers that will store the result
  uint64_t buffer_a2[1];
  char buffer_a2_var[2];
  void* buffers[] = {buffer_a2, buffer_a2_var};
  uint64_t buffer_sizes[] = {sizeof(buffer_a2), sizeof(buffer_a2_var)};

  // Open array
  tiledb_array_t* array;
  int rc = tiledb_array_alloc(ctx_, SPARSE_ARRAY_NAME, &array);
  CHECK(rc == TILEDB_OK);
  rc = tiledb_array_open(ctx_, array, TILEDB_READ);
  CHECK(rc == TILEDB_OK);

  // Create query
  tiledb_query_t* query;
  rc = tiledb_query_alloc(ctx_, array, TILEDB_READ, &query);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_data_buffer(
      ctx_, query, attributes[0], buffers[1], &buffer_sizes[1]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_offsets_buffer(
      ctx_, query, attributes[0], (uint64_t*)buffers[0], &buffer_sizes[0]);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_subarray(ctx_, query, subarray);
  REQUIRE(rc == TILEDB_OK);
  rc = tiledb_query_set_layout(ctx_, query, TILEDB_GLOBAL_ORDER);
  REQUIRE(rc == TILEDB_OK);

  // Submit query
  rc = tiledb_query_submit(ctx_, query);
  CHECK(rc == TILEDB_OK);

  // Check buffers
  char c_buffer_a2_var[2] = {'b', 'b'};
  CHECK(!memcmp(buffer_a2_var, c_buffer_a2_var, sizeof(c_buffer_a2_var)));

  // Finalize query
  rc = tiledb_query_finalize(ctx_, query);
  REQUIRE(rc == TILEDB_OK);

  // Close array
  rc = tiledb_array_close(ctx_, array);
  CHECK(rc == TILEDB_OK);

  // Clean up
  tiledb_array_free(&array);
  tiledb_query_free(&query);
}

TEST_CASE_METHOD(
    IncompleteFx,
    "C API: Test incomplete read queries, dense",
    "[capi][incomplete][dense-incomplete]") {
  bool serialized_writes = false;
  SECTION("no serialization") {
    serialized_writes = false;
  }
#ifdef TILEDB_SERIALIZATION
  SECTION("serialization enabled global order write") {
    serialized_writes = true;
  }
#endif

  remove_dense_array();
  create_dense_array();
  write_dense_full(serialized_writes);
  check_dense_incomplete();
  check_dense_until_complete();
  check_dense_shrink_buffer_size();
  check_dense_unsplittable_overflow();
  check_dense_unsplittable_complete();
  check_dense_reset_buffers();
  remove_dense_array();
}

TEST_CASE_METHOD(
    IncompleteFx,
    "C API: Test incomplete read queries, sparse",
    "[capi][incomplete][sparse]") {
  bool serialized_writes = false;
  SECTION("no serialization") {
    serialized_writes = false;
  }
#ifdef TILEDB_SERIALIZATION
  SECTION("serialization enabled global order write") {
    serialized_writes = true;
  }
#endif

  remove_sparse_array();
  create_sparse_array();
  write_sparse_full(serialized_writes);
  check_sparse_incomplete();
  check_sparse_until_complete();
  check_sparse_unsplittable_overflow();
  check_sparse_unsplittable_complete();
  remove_sparse_array();
}

#ifdef TILEDB_SERIALIZATION

TEST_CASE_METHOD(
    IncompleteFx,
    "C API: Test incomplete read queries, dense, serialized",
    "[capi][incomplete][dense][serialization]") {
  bool serialized_writes = false;
  SECTION("no serialization") {
    serialized_writes = false;
  }
  SECTION("serialization enabled global order write") {
    serialized_writes = true;
  }

  remove_dense_array();
  create_dense_array();
  write_dense_full(serialized_writes);
  check_dense_incomplete_serialized();
  remove_dense_array();
}

#endif
