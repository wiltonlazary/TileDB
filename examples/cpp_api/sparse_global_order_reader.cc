/**
 * @file   sparse_global_order_reader.cc
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
 * to it, and read a slice of the data back in the gloabl layout.
 */

#include <iostream>
#include <tiledb/tiledb>

using namespace tiledb;

// Name of arrays.
std::string array_name("sparse_global_order_reader_array");

void create_array() {
  // Create a TileDB config and context with parameters.
  Config config;
  config["sm.use_refactored_readers"] = true;
  Context ctx(config);

  // Create domain
  Domain domain(ctx);
  domain.add_dimension(Dimension::create<int>(ctx, "rows", {{1, 10}}, 1))
      .add_dimension(Dimension::create<int>(ctx, "cols", {{1, 10}}, 1));

  // The array will be sparse.
  ArraySchema schema(ctx, TILEDB_SPARSE);
  schema.set_domain(domain).set_order({{TILEDB_ROW_MAJOR, TILEDB_ROW_MAJOR}});

  // Add a single attribute "a" so each (i,j) cell can store an integer.
  schema.add_attribute(Attribute::create<int>(ctx, "a"));

  // Create the (empty) array on disk.
  Array::create(array_name, schema);
}

void write_array() {
  Config config;
  config["sm.use_refactored_readers"] = true;
  Context ctx(config);

  // Prepare data for writing.
  std::vector<int> coords_rows;
  std::vector<int> coords_cols;
  std::vector<int> data;
  for (int i = 1; i < 101; i += 10) {
    for (int j = 0; j < 10; j++) {
      coords_rows.push_back((i / 10) + 1);
      coords_cols.push_back((j % 10) + 1);
      data.push_back(i + j);
      /*std::cerr<<"row: "<<((i / 10) + 1)<<std::endl;
      std::cerr<<"col: "<<((j % 10) + 1)<<std::endl;
      std::cerr<<"data: "<<i+j<<std::endl;*/
    }
  }

  // Open the array for writing and create the query.
  Array array(ctx, array_name, TILEDB_WRITE);
  Query query(ctx, array);
  query.set_layout(TILEDB_UNORDERED)
      .set_data_buffer("a", data)
      .set_data_buffer("rows", coords_rows)
      .set_data_buffer("cols", coords_cols);

  // Perform the write and close the array.
  query.submit();
  query.finalize();
  array.close();
}

void read_array() {
  Context ctx;

  // Prepare the array for reading
  Array array(ctx, array_name, TILEDB_READ);

  // Print non-empty domain
  auto non_empty_domain = array.non_empty_domain<int>();
  std::cout << "Non-empty domain: ";
  std::cout << "[" << non_empty_domain[0].second.first << ","
            << non_empty_domain[0].second.second << "], ["
            << non_empty_domain[1].second.first << ","
            << non_empty_domain[1].second.second << "]\n";

  // Prepare buffers that will hold the results
  std::vector<int> data(100);
  std::vector<int> coords_rows(10);
  std::vector<int> coords_cols(10);

  // Prepare the query
  Query query(ctx, array, TILEDB_READ);
  query  //.set_subarray(subarray)
      .set_layout(TILEDB_GLOBAL_ORDER)
      .set_data_buffer("a", data)
      .set_data_buffer("rows", coords_rows)
      .set_data_buffer("cols", coords_cols);

  // Submit the query and close the array.
  query.submit();
  array.close();

  // Print out the results.
  auto result_num = (int)query.result_buffer_elements()["a"].second;
  std::cerr << "result_num: " << result_num << std::endl;
  for (int r = 0; r < result_num; r++) {
    int i = coords_rows[r];
    int j = coords_cols[r];
    int a = data[r];
    std::cout << "Cell (" << i << ", " << j << ") has data " << a << "\n";
  }
}

int main() {
  Context ctx;

  // Create and write the array only if it does not exist
  if (Object::object(ctx, array_name).type() != Object::Type::Array) {
    create_array();
    write_array();
  }

  read_array();
  return 0;
}