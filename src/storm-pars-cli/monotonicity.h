#pragma once

#include "storm-cli-utilities/cli.h"
#include "storm-pars/api/region.h"
#include "storm/models/sparse/Model.h"

namespace storm {
// Forward declaration
namespace jani {
class Property;
}

namespace pars {
template<typename ValueType>
void analyzeMonotonicity(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, std::vector<storm::jani::Property> const& properties,
                         std::vector<storm::storage::ParameterRegion<ValueType>> const& regions);

}
}  // namespace storm