/**
 * @file tiledb/sm/array_schema/dimension_label_schema.cc
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2022 TileDB, Inc.
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
 */

#include "tiledb/sm/array_schema/dimension_label_schema.h"
#include "tiledb/common/common.h"
#include "tiledb/sm/array_schema/array_schema.h"
#include "tiledb/sm/array_schema/attribute.h"
#include "tiledb/sm/array_schema/dimension.h"
#include "tiledb/sm/array_schema/domain.h"
#include "tiledb/sm/enums/array_type.h"
#include "tiledb/sm/misc/types.h"

using namespace tiledb::common;
using namespace tiledb::type;

namespace tiledb::sm {

DimensionLabelSchema::DimensionLabelSchema(
    LabelOrder label_order,
    Datatype label_type,
    Datatype index_type,
    const void* index_domain,
    const void* index_tile_extent)
    : label_order_(label_order)
    , indexed_array_schema_(make_shared<ArraySchema>(
          HERE(),
          label_order == LabelOrder::UNORDERED_LABELS ? ArrayType::SPARSE :
                                                        ArrayType::DENSE))
    , label_domain_{} {
  // Check the index data type is valid.
  if (!(datatype_is_integer(index_type) || datatype_is_datetime(index_type) ||
        datatype_is_time(index_type))) {
    throw std::invalid_argument(
        "Failed to create dimension label schema; Currently labels are not "
        "support on dimensions with datatype Datatype::" +
        datatype_str(index_type));
  }

  // Check the label data type is valid.
  try {
    ensure_dimension_datatype_is_valid(label_type);
  } catch (...) {
    std::throw_with_nested(std::invalid_argument(
        "Datatype Datatype::" + datatype_str(label_type) +
        " is not a valid dimension datatype."));
  }

  // Create indexed array.
  if (label_order == LabelOrder::UNORDERED_LABELS) {
    throw_if_not_ok(indexed_array_schema_->set_allows_dups(true));
  }
  std::vector<shared_ptr<Dimension>> index_dims{
      make_shared<Dimension>(HERE(), "index", index_type)};
  throw_if_not_ok(index_dims.back()->set_domain(index_domain));
  throw_if_not_ok(index_dims.back()->set_tile_extent(index_tile_extent));
  throw_if_not_ok(indexed_array_schema_->set_domain(make_shared<Domain>(
      HERE(), Layout::ROW_MAJOR, index_dims, Layout::ROW_MAJOR)));
  auto label_attr = make_shared<Attribute>(HERE(), "label", label_type);
  if (label_type == Datatype::STRING_ASCII)
    throw_if_not_ok(label_attr->set_cell_val_num(constants::var_num));
  throw_if_not_ok(indexed_array_schema_->add_attribute(label_attr));
  throw_if_not_ok(indexed_array_schema_->check());
}

DimensionLabelSchema::DimensionLabelSchema(
    LabelOrder label_order, shared_ptr<ArraySchema> indexed_array_schema)
    : label_order_(label_order)
    , indexed_array_schema_(indexed_array_schema)
    , label_domain_{} {
  // Check arrays have one dimension and one attribute
  if (indexed_array_schema->dim_num() != 1)
    throw std::invalid_argument(
        "Invalid dimension label schema; Indexed array must be one "
        "dimensional.");
  if (indexed_array_schema->attribute_num() != 1)
    throw std::invalid_argument(
        "Invalid dimension label schema; Indexed array must have exactly one "
        "attribute.");
  // Check the data type of the index
  auto index_type = indexed_array_schema_->dimension_ptr(0)->type();
  if (!(datatype_is_integer(index_type) || datatype_is_datetime(index_type) ||
        datatype_is_time(index_type)))
    throw std::invalid_argument(
        "Failed to create dimension label schema; Currently labels are not "
        "support on dimensions with datatype Datatype::" +
        datatype_str(index_type));
}

DimensionLabelSchema::DimensionLabelSchema(
    const DimensionLabelSchema* dim_label) {
  label_order_ = dim_label->label_order_;
  indexed_array_schema_ =
      make_shared<ArraySchema>(HERE(), *dim_label->indexed_array_schema_);
  label_domain_ = dim_label->label_domain_;
}

const Dimension* DimensionLabelSchema::index_dimension() const {
  return indexed_array_schema_->dimension_ptr(0);
}

/** Returns the number of values per cell for the index. */
uint32_t DimensionLabelSchema::index_cell_val_num() const {
  return index_dimension()->cell_val_num();
}

/** Returns a reference to the index domain. */
const Range& DimensionLabelSchema::index_domain() const {
  return index_dimension()->domain();
}

/** Returns the index datatype. */
Datatype DimensionLabelSchema::index_type() const {
  return index_dimension()->type();
}

bool DimensionLabelSchema::is_compatible_label(const Dimension* dim) const {
  auto dim0 = index_dimension();
  return dim->type() == dim0->type() &&
         dim->cell_val_num() == dim0->cell_val_num() &&
         dim->domain() == dim0->domain();
}

const Attribute* DimensionLabelSchema::label_attribute() const {
  return indexed_array_schema_->attribute(0);
}

/** Returns the number of values per cell for the index. */
uint32_t DimensionLabelSchema::label_cell_val_num() const {
  return label_attribute()->cell_val_num();
}

/** Returns a reference to the index domain. */
const Range& DimensionLabelSchema::label_domain() const {
  return label_domain_;
}

/** Returns the index datatype. */
Datatype DimensionLabelSchema::label_type() const {
  return label_attribute()->type();
}

}  // namespace tiledb::sm
