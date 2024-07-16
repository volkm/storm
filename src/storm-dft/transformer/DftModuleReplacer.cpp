#include "DftModuleReplacer.h"

#include "storm-dft/builder/DftBuilder.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/macros.h"

namespace storm::dft {
namespace transformer {

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftModuleReplacer<ValueType>::replaceModules(
    storm::dft::storage::DFT<ValueType> const& dft,
    std::vector<std::pair<storm::dft::storage::DftIndependentModule const&, std::shared_ptr<storm::dft::storage::DFT<ValueType> const>>> const&
        moduleReplacements) {
    // Map from module representatives to their new DFT
    std::map<size_t, std::shared_ptr<storm::dft::storage::DFT<ValueType> const>> representativeReplacements;
    // List of elements to replace
    std::set<size_t> replaceElements;
    // Initialize data structures
    for (auto const& replacements : moduleReplacements) {
        auto const& oldModule = replacements.first;
        auto allElements = oldModule.getAllElements();
        auto replaceDft = replacements.second;
        STORM_LOG_THROW(dft.getElement(oldModule.getRepresentative())->name() == replaceDft->getTopLevelElement()->name(),
                        storm::exceptions::InvalidArgumentException, "Replaced module and DFT must have the same toplevel name.");
        representativeReplacements.insert({oldModule.getRepresentative(), replaceDft});
        replaceElements.merge(allElements);
        // After merge, all elements from allElements should be moved to replace Elements
        STORM_LOG_THROW(allElements.empty(), storm::exceptions::InvalidArgumentException, "Modules overlap and cannot be replaced.");
    }

    // Replace each old module by the new DFT
    storm::dft::builder::DftBuilder<ValueType> builder{};
    std::unordered_set<std::string> depInConflict;
    for (auto const id : dft.getAllIds()) {
        auto it = representativeReplacements.find(id);
        if (it != representativeReplacements.end()) {
            // Replace representative of old module by elements from new DFT
            auto const& newModuleDft = it->second;
            for (auto const newDftId : newModuleDft->getAllIds()) {
                auto const newDftElement{newModuleDft->getElement(newDftId)};
                builder.cloneElement(newDftElement);
                // Remember dependency conflict
                if (newDftElement->isDependency() && newModuleDft->isDependencyInConflict(newDftId)) {
                    depInConflict.insert(newDftElement->name());
                }
            }
        } else if (replaceElements.find(id) == replaceElements.end()) {
            // Element is not part of a replaced module -> keep it
            auto const element{dft.getElement(id)};
            builder.cloneElement(element);
            // Remember dependency conflict
            if (element->isDependency() && dft.isDependencyInConflict(id)) {
                depInConflict.insert(element->name());
            }
        }
    }
    builder.setTopLevel(dft.getTopLevelElement()->name());
    // Build new DFT
    auto newDft = std::make_shared<storm::dft::storage::DFT<ValueType>>(builder.build());

    // Update dependency conflicts
    for (size_t id : newDft->getDependencies()) {
        // Set dependencies not in conflict
        if (depInConflict.find(newDft->getElement(id)->name()) == depInConflict.end()) {
            newDft->setDependencyNotInConflict(id);
        }
    }

    return newDft;
}

// Explicitly instantiate the class.
template class DftModuleReplacer<double>;
template class DftModuleReplacer<RationalFunction>;

}  // namespace transformer
}  // namespace storm::dft
