#include "DftModuleTransformer.h"

#include "storm-dft/builder/DftBuilder.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"

namespace storm::dft {
namespace transformer {

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftModuleTransformer<ValueType>::replaceModules(
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

template<typename ValueType>
std::vector<size_t> DftModuleTransformer<ValueType>::findModularisationRewrite(storm::dft::storage::DFT<ValueType> const& dft) {
    for (size_t i = 0; i < dft.nrElements(); ++i) {
        auto const& e = dft.getElement(i);
        if (e->isGate() &&
            (e->type() == storm::dft::storage::elements::DFTElementType::AND || e->type() == storm::dft::storage::elements::DFTElementType::OR)) {
            // Suitable parent gate -> check the independent submodules of the children
            auto const& children = std::static_pointer_cast<storm::dft::storage::elements::DFTGate<ValueType> const>(e)->children();
            for (auto const& child : children) {
                // Compute independent subtree
                auto ISD = std::static_pointer_cast<storm::dft::storage::elements::DFTGate<ValueType>>(child)->independentSubDft(true);
                // Check whether the independent subtree contains some of the other children
                std::vector<size_t> rewrite = {e->id(), child->id()};
                for (size_t isdElemId : ISD) {
                    if (isdElemId == child->id()) {
                        continue;
                    }
                    if (std::find_if(children.begin(), children.end(),
                                     [&isdElemId](std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType>> const& e) {
                                         return e->id() == isdElemId;
                                     }) != children.end()) {
                        // Element in subtree is also a child
                        rewrite.push_back(isdElemId);
                    }
                }
                if (rewrite.size() > 2 && rewrite.size() < children.size() - 1) {
                    // Found children for rewrite
                    return rewrite;
                }
            }
        }
    }
    return {};
}

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftModuleTransformer<ValueType>::prepareModularisation(storm::dft::storage::DFT<ValueType> const& dft) {
    std::vector<size_t> modIdea = DftModuleTransformer<ValueType>::findModularisationRewrite(dft);
    STORM_LOG_DEBUG("Modularisation idea: " << storm::utility::vector::toString(modIdea));

    if (modIdea.empty()) {
        // No rewrite needed
        return std::make_shared<storm::dft::storage::DFT<ValueType>>(dft);
    }

    storm::dft::builder::DftBuilder<ValueType> builder;

    // Copy all other elements which do not change
    for (size_t i = 0; i < dft.nrElements(); ++i) {
        std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> element = dft.getElement(i);
        if (element->id() != modIdea[0]) {
            builder.cloneElement(element);
        }
    }

    // Add rewritten elements
    STORM_LOG_ASSERT(modIdea.size() > 1, "No rewritten elements.");
    std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> child = dft.getElement(modIdea[1]);
    STORM_LOG_ASSERT(child->hasParents(), "Rewritten elements has no parents.");
    STORM_LOG_ASSERT(child->parents().front()->isGate(), "Rewritten element has no parent gate.");
    auto originalParent = dft.getGate(modIdea[0]);
    STORM_LOG_ASSERT(std::find_if(child->parents().begin(), child->parents().end(),
                                  [&originalParent](std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType>> const& p) {
                                      return p->id() == originalParent->id();
                                  }) != child->parents().end(),
                     "Rewritten element has not the same parent");
    std::string newParentName = originalParent->name() + "_module";

    // Accumulate children names
    std::vector<std::string> childrenNames;
    for (size_t i = 1; i < modIdea.size(); ++i) {
        std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> otherChild = dft.getElement(modIdea[i]);
        STORM_LOG_ASSERT(std::find_if(otherChild->parents().begin(), otherChild->parents().end(),
                                      [&originalParent](std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType>> const& p) {
                                          return p->id() == originalParent->id();
                                      }) != otherChild->parents().end(),
                         "Children have not the same father for rewrite " << otherChild->name());

        childrenNames.push_back(otherChild->name());
    }

    // Add element in-between parent and children
    switch (originalParent->type()) {
        case storm::dft::storage::elements::DFTElementType::AND:
            builder.addAndGate(newParentName, childrenNames);
            break;
        case storm::dft::storage::elements::DFTElementType::OR:
            builder.addOrGate(newParentName, childrenNames);
            break;
        default:
            // Other elements are not supported
            STORM_LOG_ASSERT(false, "Dft type " << originalParent->type() << " can not be rewritten.");
            break;
    }

    // Add parent with the new child newParent and all its remaining children
    childrenNames.clear();
    childrenNames.push_back(newParentName);
    for (auto const& child : originalParent->children()) {
        if (std::find(modIdea.begin() + 1, modIdea.end(), child->id()) == modIdea.end()) {
            // Child was not rewritten and must be kept
            childrenNames.push_back(child->name());
        }
    }
    builder.cloneElementWithNewChildren(originalParent, childrenNames);

    builder.setTopLevel(dft.getTopLevelElement()->name());
    // TODO use reference?
    storm::dft::storage::DFT<ValueType> newDft = builder.build();
    STORM_LOG_TRACE(newDft.getElementsString());
    return prepareModularisation(newDft);
}

// Explicitly instantiate the class.
template class DftModuleTransformer<double>;
template class DftModuleTransformer<RationalFunction>;

}  // namespace transformer
}  // namespace storm::dft
