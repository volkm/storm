#include "monotonicity.h"

#include "storm-pars-cli/print.h"
#include "storm-pars/analysis/MonotonicityHelper.h"
#include "storm-pars/api/region.h"
#include "storm-pars/settings/modules/MonotonicitySettings.h"
#include "storm-pars/settings/modules/ParametricSettings.h"
#include "storm-pars/settings/modules/RegionSettings.h"
#include "storm-pars/utility/parametric.h"
#include "storm/api/verification.h"
#include "storm/exceptions/BaseException.h"
#include "storm/io/file.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/models/ModelBase.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/storage/jani/Property.h"
#include "storm/utility/Engine.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/initialize.h"
#include "storm/utility/macros.h"

namespace storm::pars {
template<typename ValueType>
void analyzeMonotonicity(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, std::vector<storm::jani::Property> const& properties,
                         std::vector<storm::storage::ParameterRegion<ValueType>> const& regions) {
    std::ofstream outfile;
    auto monSettings = storm::settings::getModule<storm::settings::modules::MonotonicitySettings>();

    if (monSettings.isExportMonotonicitySet()) {
        storm::utility::openFile(monSettings.getExportMonotonicityFilename(), outfile);
    }
    std::vector<std::shared_ptr<storm::logic::Formula const>> formulas = storm::api::extractFormulasFromProperties(properties);
    storm::utility::Stopwatch monotonicityWatch(true);
    STORM_LOG_THROW(regions.size() <= 1, storm::exceptions::InvalidArgumentException, "Monotonicity analysis only allowed on single region");
    if (!monSettings.isMonSolutionSet()) {
        auto monotonicityHelper = storm::analysis::MonotonicityHelper<ValueType, double>(
            model, formulas, regions, monSettings.getNumberOfSamples(), storm::settings::getModule<storm::settings::modules::GeneralSettings>().getPrecision(),
            monSettings.isDotOutputSet());
        if (monSettings.isExportMonotonicitySet()) {
            monotonicityHelper.checkMonotonicityInBuild(outfile, monSettings.isUsePLABoundsSet(), monSettings.getDotOutputFilename());
        } else {
            monotonicityHelper.checkMonotonicityInBuild(std::cout, monSettings.isUsePLABoundsSet(), monSettings.getDotOutputFilename());
        }
    } else {
        // Checking monotonicity based on solution function

        auto parametricSettings = storm::settings::getModule<storm::settings::modules::ParametricSettings>();
        auto regionSettings = storm::settings::getModule<storm::settings::modules::RegionSettings>();

        std::function<std::unique_ptr<storm::modelchecker::CheckResult>(std::shared_ptr<storm::logic::Formula const> const& formula)> verificationCallback;
        std::function<void(std::unique_ptr<storm::modelchecker::CheckResult> const&)> postprocessingCallback;

        // Check the given set of regions with or without refinement
        verificationCallback = [&](std::shared_ptr<storm::logic::Formula const> const& formula) {
            std::unique_ptr<storm::modelchecker::CheckResult> result =
                storm::api::verifyWithSparseEngine<ValueType>(model, storm::api::createTask<ValueType>(formula, true));
            return result;
        };

        for (auto& property : properties) {
            auto result = verificationCallback(property.getRawFormula())->asExplicitQuantitativeCheckResult<ValueType>().getValueVector();
            ValueType valuation;

            auto states = model->getInitialStates();
            for (auto state : states) {
                valuation += result[state];
            }

            storm::analysis::MonotonicityResult<storm::RationalFunctionVariable> monRes;
            for (auto& var : storm::models::sparse::getProbabilityParameters(*model)) {
                auto res = storm::analysis::MonotonicityChecker<ValueType>::checkDerivative(valuation.derivative(var), regions[0]);

                if (res.first && res.second) {
                    monRes.addMonotonicityResult(var, analysis::MonotonicityResult<storm::RationalFunctionVariable>::Monotonicity::Constant);
                } else if (res.first) {
                    monRes.addMonotonicityResult(var, analysis::MonotonicityResult<storm::RationalFunctionVariable>::Monotonicity::Incr);
                } else if (res.second) {
                    monRes.addMonotonicityResult(var, analysis::MonotonicityResult<storm::RationalFunctionVariable>::Monotonicity::Decr);
                } else {
                    monRes.addMonotonicityResult(var, analysis::MonotonicityResult<storm::RationalFunctionVariable>::Monotonicity::Not);
                }
            }
            if (monSettings.isExportMonotonicitySet()) {
                outfile << monRes.toString();
            } else {
                STORM_PRINT(monRes.toString());
            }
        }
    }

    if (monSettings.isExportMonotonicitySet()) {
        storm::utility::closeFile(outfile);
    }

    monotonicityWatch.stop();
    STORM_PRINT("\nTotal time for monotonicity checking: " << monotonicityWatch << ".\n\n");
    return;
}

template void analyzeMonotonicity(std::shared_ptr<storm::models::sparse::Model<storm::RationalFunction>> const& model,
                                  std::vector<storm::jani::Property> const& properties,
                                  std::vector<storm::storage::ParameterRegion<storm::RationalFunction>> const& regions);
}  // namespace storm::pars