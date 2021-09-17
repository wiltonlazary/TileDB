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

#include <cmath>
#include <iostream>
#include <tiledb/tiledb>

using namespace tiledb;

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
    const std::string& array_name,
    const std::vector<test_dim_t<DIM_T>>& test_dims,
    const std::vector<test_attr_t<ATTR_T>>& test_attrs) {
  // Create domain
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
void write(
    const std::string& array_name,
    const std::vector<test_query_buffer_t<ATTR_T>>& test_query_buffers) {
  // Open the array for writing
  Config config;
  config["sm.use_refactored_readers"] = true;
  Context ctx(config);
  Array array(ctx, array_name, TILEDB_WRITE);

  // Create the queries with unordered layouts
  Query query(ctx, array);
  query.set_layout(TILEDB_UNORDERED);

  // Set the query buffers
  for (const auto& test_query_buffer : test_query_buffers) {
    query.set_data_buffer(test_query_buffer.name_, *test_query_buffer.data_);
  }

  // Submit and finalize the queries
  query.submit();
  query.finalize();

  // Close the array
  array.close();
}

template <typename ATTR_T>
void read_array(
    const std::string& array_name,
    const std::vector<test_query_buffer_t<ATTR_T>>& test_query_buffers) {
  Config config;
  config["sm.use_refactored_readers"] = "true";
  Context ctx(config);

  // Open the array for reading
  Array array(ctx, array_name, TILEDB_READ);

  // Print non-empty domain
  auto non_empty_domain = array.non_empty_domain<uint64_t>();
  std::cout << "Non-empty domain: ";
  std::cout << "[" << non_empty_domain[0].second.first << ","
            << non_empty_domain[0].second.second << "], ["
            << non_empty_domain[1].second.first << ","
            << non_empty_domain[1].second.second << "]\n";

  // Create the read query
  Query query(ctx, array, TILEDB_READ);
  query.set_layout(TILEDB_GLOBAL_ORDER);

  // Set the query buffers.
  for (const auto& test_query_buffer : test_query_buffers) {
    query.set_data_buffer(test_query_buffer.name_, *test_query_buffer.data_);
  }

  // Submit the query.
  query.submit();
  query.finalize();

  // Close the array.
  array.close();
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

void array_ordered() {
  Context ctx;

  // Name of array.
  std::string array_ordered("sparse_global_order_reader_ordered_array");

  // Define the dimensions.
  std::vector<test_dim_t<uint64_t>> dims;
  const uint64_t domain_min = 1;
  const uint64_t domain_max = 55;
  const uint64_t tile_extent = 11;
  const std::array<uint64_t, 2> rows_domain = {domain_min, domain_max};
  dims.emplace_back("rows", rows_domain, tile_extent);
  const std::array<uint64_t, 2> cols_domain = {domain_min, domain_max};
  dims.emplace_back("cols", cols_domain, tile_extent);

  // Define the attributes.
  std::vector<test_attr_t<uint64_t>> attrs;
  attrs.emplace_back("a");

  // Define the write query buffers for "a" and
  // dimension query buffers with an unordered write order.
  // Fragment 1
  std::vector<test_query_buffer_t<uint64_t>> write_query_buffers_1;
  std::vector<uint64_t> a_write_buffer_1;
  std::vector<uint64_t> rows_write_buffer_1;
  std::vector<uint64_t> cols_write_buffer_1;
  for (int i = 0; i < 1000; i++) {
    uint64_t coords[2] = {0, 0};
    add_tile_coords(i, dims, coords);
    add_cell_coords(i, dims, coords);

    rows_write_buffer_1.emplace_back(coords[0]);
    cols_write_buffer_1.emplace_back(coords[1]);
    a_write_buffer_1.emplace_back(i);
  }
  write_query_buffers_1.emplace_back("a", &a_write_buffer_1);
  write_query_buffers_1.emplace_back("rows", &rows_write_buffer_1);
  write_query_buffers_1.emplace_back("cols", &cols_write_buffer_1);

  // Fragment 2
  /*std::vector<test_query_buffer_t<uint64_t>> write_query_buffers_2;
  std::vector<uint64_t> a_write_buffer_2;
  std::vector<uint64_t> rows_write_buffer_2;
  std::vector<uint64_t> cols_write_buffer_2;
  for (int i = 1001; i < 2001; i++) {
    uint64_t row = i / domain_max + domain_min;
    uint64_t col = i % domain_max + domain_min;
    uint64_t coords[2] = {row, col};
    auto tile_pos = get_tile_position(coords, dims);
    auto cell_pos = get_cell_position_in_tile(coords, dims);
    auto a = tile_pos * tile_extent * tile_extent + cell_pos;
    rows_write_buffer_2.emplace_back(row);
    cols_write_buffer_2.emplace_back(col);
    a_write_buffer_2.emplace_back(a);
    std::cout << row << " " << col << " " << a << std::endl;
  }
  write_query_buffers_2.emplace_back("a", &a_write_buffer_2);
  write_query_buffers_2.emplace_back("rows", &rows_write_buffer_2);
  write_query_buffers_2.emplace_back("cols", &cols_write_buffer_2);

  // Fragment 3
  std::vector<test_query_buffer_t<uint64_t>> write_query_buffers_3;
  std::vector<uint64_t> a_write_buffer_3;
  std::vector<uint64_t> rows_write_buffer_3;
  std::vector<uint64_t> cols_write_buffer_3;
  for (int i = 2001; i < 3001; i++) {
    uint64_t row = i / domain_max + domain_min;
    uint64_t col = i % domain_max + domain_min;
    uint64_t coords[2] = {row, col};
    auto tile_pos = get_tile_position(coords, dims);
    auto cell_pos = get_cell_position_in_tile(coords, dims);
    auto a = tile_pos * tile_extent * tile_extent + cell_pos;
    rows_write_buffer_2.emplace_back(row);
    cols_write_buffer_2.emplace_back(col);
    a_write_buffer_2.emplace_back(a);
    std::cout << row << " " << col << " " << a << std::endl;
  }
  write_query_buffers_3.emplace_back("a", &a_write_buffer_3);
  write_query_buffers_3.emplace_back("rows", &rows_write_buffer_3);
  write_query_buffers_3.emplace_back("cols", &cols_write_buffer_3);*/

  std::vector<test_query_buffer_t<uint64_t>> read_query_buffers;
  std::vector<uint64_t> data(1000);
  std::vector<uint64_t> coords_rows(1000);
  std::vector<uint64_t> coords_cols(1000);
  read_query_buffers.emplace_back("a", &data);
  read_query_buffers.emplace_back("rows", &coords_rows);
  read_query_buffers.emplace_back("cols", &coords_cols);

  // Create and write the arrays only if they do not exist
  if (Object::object(ctx, array_ordered).type() != Object::Type::Array) {
    create_array(array_ordered, dims, attrs);
    write(array_ordered, write_query_buffers_1);
    // write(array_ordered, write_query_buffers_2);
    // write(array_ordered, write_query_buffers_3);
  }
  read_array(array_ordered, read_query_buffers);
  for (uint64_t i = 0; i < 1000; i++) {
    std::cerr << "data[" << i << "]:" << data[i] << std::endl;
    // if (data[i] != i+1)
    //  std::cerr<<"Data "<<i<<" does not match anticipated value"<<std::endl;
  }
  // validate_read(read_query_buffers, full_data_write_buffer);
}

int main() {
  array_ordered();
  // array_interleaved();
  // array_duplicated();
  return 0;
}