#pragma once

#include <iostream>
#include <memory>

#include "storm/models/sparse/Model.h"

namespace storm {
namespace exporter {

/*!
 * Exports a sparse model into the AUT format.
 *
 * @param os           Stream to export to
 * @param sparseModel  Model to export
 */
template<typename ValueType>
void exportModelAsAut(std::ostream& os, std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModel);
}  // namespace exporter
}  // namespace storm
