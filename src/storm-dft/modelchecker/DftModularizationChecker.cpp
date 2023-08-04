#include "DftModularizationChecker.h"

#include <sstream>

#include "storm-dft/api/storm-dft.h"
#include "storm-dft/builder/DftBuilder.h"
#include "storm-dft/modelchecker/DFTModelChecker.h"
#include "storm-dft/modelchecker/SftBddChecker.h"
#include "storm-dft/transformations/DftModuleReplacer.h"
#include "storm-dft/transformations/PropertyToBddTransformer.h"
#include "storm-dft/utility/DftModularizer.h"

#include "storm-parsers/api/properties.h"
#include "storm/api/properties.h"

namespace storm::dft {
namespace modelchecker {

template<typename ValueType>
DftModularizationChecker<ValueType>::DftModularizationChecker(std::shared_ptr<storm::dft::storage::DFT<ValueType> const> dft) : dft{dft} {
    // Initialize modules
    storm::dft::utility::DftModularizer<ValueType> modularizer;
    auto topModule = modularizer.computeModules(*dft);
    STORM_LOG_DEBUG("Modularization found the following modules:\n" << topModule.toString(*dft));

    // Gather all dynamic modules
    populateDynamicModules(topModule);
}

template<typename ValueType>
void DftModularizationChecker<ValueType>::populateDynamicModules(storm::dft::storage::DftIndependentModule const& module) {
    if (!module.isStatic()) {
        // Found new dynamic module
        dynamicModules.push_back(module);
    } else if (!module.isFullyStatic()) {
        // Module contains dynamic sub-modules -> recursively visit children
        for (auto const& submodule : module.getSubModules()) {
            populateDynamicModules(submodule);
        }
    }
}

template<typename ValueType>
std::vector<ValueType> DftModularizationChecker<ValueType>::check(FormulaVector const& formulas, size_t chunksize) {
    // Gather all occurring time points
    std::set<ValueType> timepointSet;
    for (auto const& formula : formulas) {
        timepointSet.insert(storm::dft::transformations::PropertyToBddTransformer<ValueType>::getTimebound(*formula));
    }
    std::vector<ValueType> timepoints(timepointSet.begin(), timepointSet.end());

    // Build SFT by replacing dynamic modules
    auto sft = replaceDynamicModules(timepoints);

    auto builder = std::make_shared<storm::dft::builder::BddSftModelBuilder<ValueType>>(sft);
    storm::dft::modelchecker::SftBddChecker checker{builder};
    return checker.check(formulas, chunksize);
}

template<typename ValueType>
std::vector<ValueType> DftModularizationChecker<ValueType>::getProbabilitiesAtTimepoints(std::vector<ValueType> const& timepoints, size_t chunksize) {
    auto sft = replaceDynamicModules(timepoints);
    auto builder = std::make_shared<storm::dft::builder::BddSftModelBuilder<ValueType>>(sft);
    storm::dft::modelchecker::SftBddChecker checker{builder};
    return checker.getProbabilitiesAtTimepoints(timepoints, chunksize);
}

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> DftModularizationChecker<ValueType>::replaceDynamicModules(std::vector<ValueType> const& timepoints) {
    // Store replacements of dynamic module by DFT consisting of single BE
    std::vector<std::pair<storm::dft::storage::DftIndependentModule const&, std::shared_ptr<storm::dft::storage::DFT<ValueType> const>>> replacements;

    // Analyse all dynamic modules and create corresponding DFT consisting of a single BE
    // Instead of building a DFT we could directly replace a module with a BE. However, this approach provides more flexibility.
    for (auto const& mod : dynamicModules) {
        STORM_LOG_DEBUG("Analyse dynamic module " << mod.toString(*dft));
        auto result = analyseDynamicModule(mod, timepoints);
        // Create single BE which has samples corresponding to the previously computed analysis results
        std::map<ValueType, ValueType> probabilities{};
        for (size_t i = 0; i < timepoints.size(); ++i) {
            auto const probability{boost::get<ValueType>(result[i])};
            auto const timebound{timepoints[i]};
            probabilities[timebound] = probability;
        }
        storm::dft::builder::DftBuilder<ValueType> builder;
        std::string const& elementName = dft->getElement(mod.getRepresentative())->name();
        builder.addBasicElementSamples(elementName, probabilities);
        // Create DFT
        builder.setTopLevel(elementName);
        auto beDft = std::make_shared<storm::dft::storage::DFT<ValueType>>(builder.build());
        STORM_LOG_ASSERT(beDft->nrElements() == 1 && beDft->getTopLevelElement()->isBasicElement(), "DFT should be single BE");
        STORM_LOG_DEBUG("Replace module " << dft->getElement(mod.getRepresentative())->toString() << " by DFT with single BE: " << beDft->getElementsString());
        replacements.push_back({mod, beDft});
    }

    // Perform replacements
    auto sft = storm::dft::transformations::DftModuleReplacer<ValueType>::replaceModules(*dft, replacements);
    STORM_LOG_ASSERT(sft->getDependencies().empty(), "SFT should have no dependencies.");
    STORM_LOG_DEBUG("Remaining static FT: " << sft->getElementsString());
    return sft;
}

template<typename ValueType>
typename storm::dft::modelchecker::DFTModelChecker<ValueType>::dft_results DftModularizationChecker<ValueType>::analyseDynamicModule(
    storm::dft::storage::DftIndependentModule const& module, std::vector<ValueType> const& timepoints) {
    STORM_LOG_ASSERT(!module.isStatic() && !module.isFullyStatic(), "Module should be dynamic.");
    STORM_LOG_ASSERT(!dft->getElement(module.getRepresentative())->isBasicElement(), "Dynamic module should not be a single BE.");

    auto subDft = module.getSubtree(*dft);

    // Create properties
    std::stringstream propertyStream{};
    for (auto const timebound : timepoints) {
        propertyStream << "Pmin=? [F<=" << timebound << "\"failed\"];";
    }
    auto const props{storm::api::extractFormulasFromProperties(storm::api::parseProperties(propertyStream.str()))};

    storm::dft::modelchecker::DFTModelChecker<ValueType> modelchecker(true);
    // Check DFT
    // TODO activate symred and other optimizations
    return std::move(modelchecker.check(subDft, props, false, false, {}));
}

// Explicitly instantiate the class.
template class DftModularizationChecker<double>;

}  // namespace modelchecker
}  // namespace storm::dft
