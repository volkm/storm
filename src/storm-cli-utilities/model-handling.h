#pragma once

#include <type_traits>

#include "storm-cli-utilities/print.h"
#include "storm-counterexamples/settings/modules/CounterexampleGeneratorSettings.h"
#include "storm/api/bisimulation.h"
#include "storm/builder/BuilderOptions.h"
#include "storm/environment/Environment.h"
#include "storm/exceptions/InvalidSettingsException.h"
#include "storm/modelchecker/results/CheckResult.h"
#include "storm/models/symbolic/StandardRewardModel.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/AbstractionSettings.h"
#include "storm/settings/modules/BuildSettings.h"
#include "storm/settings/modules/DebugSettings.h"
#include "storm/settings/modules/IOSettings.h"
#include "storm/settings/modules/TransformationSettings.h"
#include "storm/storage/jani/Property.h"
#include "storm/utility/Engine.h"
#include "storm/utility/Stopwatch.h"

namespace storm::cli {

struct SymbolicInput {
    // The symbolic model description.
    boost::optional<storm::storage::SymbolicModelDescription> model;

    // The original properties to check.
    std::vector<storm::jani::Property> properties;

    // The preprocessed properties to check (in case they needed amendment).
    boost::optional<std::vector<storm::jani::Property>> preprocessedProperties;
};

void parseSymbolicModelDescription(storm::settings::modules::IOSettings const& ioSettings, SymbolicInput& input);

void parseProperties(storm::settings::modules::IOSettings const& ioSettings, SymbolicInput& input,
                     boost::optional<std::set<std::string>> const& propertyFilter);

SymbolicInput parseSymbolicInputQvbs(storm::settings::modules::IOSettings const& ioSettings);

SymbolicInput parseSymbolicInput();

struct ModelProcessingInformation {
    // The engine to use
    storm::utility::Engine engine;

    // If set, bisimulation will be applied.
    bool applyBisimulation;

    // If set, a transformation to Jani will be enforced
    bool transformToJani;

    // Which data type is to be used for numbers ...
    enum class ValueType { FinitePrecision, Exact, Parametric };
    ValueType buildValueType;         // ... during model building
    ValueType verificationValueType;  // ... during model verification

    // The Dd library to be used
    storm::dd::DdType ddType;

    // The environment used during model checking
    storm::Environment env;

    // A flag which is set to true, if the settings were detected to be compatible.
    // If this is false, it could be that the query can not be handled.
    bool isCompatible;
};

void getModelProcessingInformationAutomatic(SymbolicInput const& input, ModelProcessingInformation& mpi);

/*!
 * Sets the model processing information based on the given input.
 * Finding the right model processing information might require a conversion to jani.
 * In this case, the jani conversion is stored in the transformedJaniInput pointer (unless it is null)
 */
ModelProcessingInformation getModelProcessingInformation(SymbolicInput const& input, std::shared_ptr<SymbolicInput> const& transformedJaniInput = nullptr);

void ensureNoUndefinedPropertyConstants(std::vector<storm::jani::Property> const& properties);

std::pair<SymbolicInput, ModelProcessingInformation> preprocessSymbolicInput(SymbolicInput const& input);

void exportSymbolicInput(SymbolicInput const& input);

std::vector<std::shared_ptr<storm::logic::Formula const>> createFormulasToRespect(std::vector<storm::jani::Property> const& properties);

template<storm::dd::DdType DdType, typename ValueType>
std::shared_ptr<storm::models::ModelBase> buildModelDd(SymbolicInput const& input);

storm::builder::BuilderOptions createBuildOptionsSparseFromSettings(SymbolicInput const& input);

template<typename ValueType>
std::shared_ptr<storm::models::ModelBase> buildModelSparse(SymbolicInput const& input, storm::builder::BuilderOptions const& options);

template<typename ValueType>
std::shared_ptr<storm::models::ModelBase> buildModelExplicit(storm::settings::modules::IOSettings const& ioSettings,
                                                             storm::settings::modules::BuildSettings const& buildSettings);

template<storm::dd::DdType DdType, typename ValueType>
std::shared_ptr<storm::models::ModelBase> buildModel(SymbolicInput const& input, storm::settings::modules::IOSettings const& ioSettings,
                                                     ModelProcessingInformation const& mpi);

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> preprocessSparseMarkovAutomaton(
    std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> const& model);

template<typename ValueType>
std::shared_ptr<storm::models::sparse::Model<ValueType>> preprocessSparseModelBisimulation(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, SymbolicInput const& input,
    storm::settings::modules::BisimulationSettings const& bisimulationSettings);

template<typename ValueType>
std::pair<std::shared_ptr<storm::models::sparse::Model<ValueType>>, bool> preprocessSparseModel(
    std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<typename ValueType>
void exportSparseModel(std::shared_ptr<storm::models::sparse::Model<ValueType>> const& model, SymbolicInput const& input);

template<storm::dd::DdType DdType, typename ValueType>
void exportDdModel(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model, SymbolicInput const&);

template<storm::dd::DdType DdType, typename ValueType>
void exportModel(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input);

template<storm::dd::DdType DdType, typename ValueType>
typename std::enable_if<DdType != storm::dd::DdType::Sylvan && !std::is_same<ValueType, double>::value, std::shared_ptr<storm::models::Model<ValueType>>>::type
preprocessDdMarkovAutomaton(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model) {
    return model;
}

template<storm::dd::DdType DdType, typename ValueType>
typename std::enable_if<DdType == storm::dd::DdType::Sylvan || std::is_same<ValueType, double>::value, std::shared_ptr<storm::models::Model<ValueType>>>::type
preprocessDdMarkovAutomaton(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model) {
    auto ma = model->template as<storm::models::symbolic::MarkovAutomaton<DdType, ValueType>>();
    if (!ma->isClosed()) {
        return std::make_shared<storm::models::symbolic::MarkovAutomaton<DdType, ValueType>>(ma->close());
    } else {
        return model;
    }
}

template<storm::dd::DdType DdType, typename ValueType, typename ExportValueType = ValueType>
std::shared_ptr<storm::models::Model<ExportValueType>> preprocessDdModelBisimulation(
    std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model, SymbolicInput const& input,
    storm::settings::modules::BisimulationSettings const& bisimulationSettings, ModelProcessingInformation const& mpi) {
    STORM_LOG_WARN_COND(!bisimulationSettings.isWeakBisimulationSet(),
                        "Weak bisimulation is currently not supported on DDs. Falling back to strong bisimulation.");

    auto quotientFormat = bisimulationSettings.getQuotientFormat();
    if (mpi.engine == storm::utility::Engine::DdSparse && quotientFormat != storm::dd::bisimulation::QuotientFormat::Sparse &&
        bisimulationSettings.isQuotientFormatSetFromDefaultValue()) {
        STORM_LOG_INFO("Setting bisimulation quotient format to 'sparse'.");
        quotientFormat = storm::dd::bisimulation::QuotientFormat::Sparse;
    }
    STORM_LOG_INFO("Performing bisimulation minimization...");
    return storm::api::performBisimulationMinimization<DdType, ValueType, ExportValueType>(
        model, createFormulasToRespect(input.properties), storm::storage::BisimulationType::Strong, bisimulationSettings.getSignatureMode(), quotientFormat);
}

template<storm::dd::DdType DdType, typename ValueType, typename ExportValueType = ValueType>
std::pair<std::shared_ptr<storm::models::ModelBase>, bool> preprocessDdModel(std::shared_ptr<storm::models::symbolic::Model<DdType, ValueType>> const& model,
                                                                             SymbolicInput const& input, ModelProcessingInformation const& mpi) {
    auto bisimulationSettings = storm::settings::getModule<storm::settings::modules::BisimulationSettings>();
    std::pair<std::shared_ptr<storm::models::Model<ValueType>>, bool> intermediateResult = std::make_pair(model, false);

    if (model->isOfType(storm::models::ModelType::MarkovAutomaton)) {
        intermediateResult.first = preprocessDdMarkovAutomaton(intermediateResult.first->template as<storm::models::symbolic::Model<DdType, ValueType>>());
        intermediateResult.second = true;
    }

    std::unique_ptr<std::pair<std::shared_ptr<storm::models::Model<ExportValueType>>, bool>> result;
    auto symbolicModel = intermediateResult.first->template as<storm::models::symbolic::Model<DdType, ValueType>>();
    if (mpi.applyBisimulation) {
        std::shared_ptr<storm::models::Model<ExportValueType>> newModel =
            preprocessDdModelBisimulation<DdType, ValueType, ExportValueType>(symbolicModel, input, bisimulationSettings, mpi);
        result = std::make_unique<std::pair<std::shared_ptr<storm::models::Model<ExportValueType>>, bool>>(newModel, true);
    } else {
        result = std::make_unique<std::pair<std::shared_ptr<storm::models::Model<ExportValueType>>, bool>>(
            symbolicModel->template toValueType<ExportValueType>(), !std::is_same<ValueType, ExportValueType>::value);
    }

    if (result && result->first->isSymbolicModel() && mpi.engine == storm::utility::Engine::DdSparse) {
        // Mark as changed.
        result->second = true;

        std::shared_ptr<storm::models::symbolic::Model<DdType, ExportValueType>> symbolicModel =
            result->first->template as<storm::models::symbolic::Model<DdType, ExportValueType>>();
        std::vector<std::shared_ptr<storm::logic::Formula const>> formulas;
        for (auto const& property : input.properties) {
            formulas.emplace_back(property.getRawFormula());
        }
        result->first = storm::api::transformSymbolicToSparseModel(symbolicModel, formulas);
        STORM_LOG_THROW(result, storm::exceptions::NotSupportedException, "The translation to a sparse model is not supported for the given model type.");
    }

    return *result;
}

template<storm::dd::DdType DdType, typename BuildValueType, typename ExportValueType = BuildValueType>
std::pair<std::shared_ptr<storm::models::ModelBase>, bool> preprocessModel(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input,
                                                                           ModelProcessingInformation const& mpi) {
    storm::utility::Stopwatch preprocessingWatch(true);

    std::pair<std::shared_ptr<storm::models::ModelBase>, bool> result = std::make_pair(model, false);
    if (model->isSparseModel()) {
        result = preprocessSparseModel<BuildValueType>(result.first->as<storm::models::sparse::Model<BuildValueType>>(), input, mpi);
    } else {
        STORM_LOG_ASSERT(model->isSymbolicModel(), "Unexpected model type.");
        result =
            preprocessDdModel<DdType, BuildValueType, ExportValueType>(result.first->as<storm::models::symbolic::Model<DdType, BuildValueType>>(), input, mpi);
    }

    preprocessingWatch.stop();

    if (result.second) {
        STORM_PRINT("\nTime for model preprocessing: " << preprocessingWatch << ".\n\n");
    }
    return result;
}

template<typename ValueType>
void generateCounterexamples(std::shared_ptr<storm::models::ModelBase> const&, SymbolicInput const&);

template<typename ValueType>
void printResult(std::unique_ptr<storm::modelchecker::CheckResult> const& result, storm::jani::Property const& property,
                 storm::utility::Stopwatch* watch = nullptr);

using VerificationCallbackType = std::function<std::unique_ptr<storm::modelchecker::CheckResult>(std::shared_ptr<storm::logic::Formula const> const& formula,
                                                                                                 std::shared_ptr<storm::logic::Formula const> const& states)>;
using PostprocessingCallbackType = std::function<void(std::unique_ptr<storm::modelchecker::CheckResult> const&)>;

struct PostprocessingIdentity {
    void operator()(std::unique_ptr<storm::modelchecker::CheckResult> const&) {
        // Intentionally left empty.
    }
};

/*!
 * Verifies the given formula plus a filter formula to identify relevant states and warns the user in case of issues
 * @param formula the formula to check
 * @param statesFilter a second formula that identifies the relevant states (needs to be qualitative)
 * @param verificationCallback Function to perform the actual verification task for a given formula plus a filter formula to identify relevant states
 */
template<typename ValueType>
std::unique_ptr<storm::modelchecker::CheckResult> verifyProperty(std::shared_ptr<storm::logic::Formula const> const& formula,
                                                                 std::shared_ptr<storm::logic::Formula const> const& statesFilter,
                                                                 VerificationCallbackType const& verificationCallback);

/*!
 * Verifies all (potentially preprocessed) properties given in `input`.
 * @param input Where the properties are read from
 * @param verificationCallback Function to perform the actual verification task for a given formula plus a filter formula to identify relevant states
 * @param postprocessingCallback Function that processes the verification result, such as e.g. output to a file
 */
template<typename ValueType>
void verifyProperties(SymbolicInput const& input, VerificationCallbackType const& verificationCallback,
                      PostprocessingCallbackType const& postprocessingCallback = PostprocessingIdentity());

/*!
 * Computes values for each state (such as the steady-state probability distribution).
 * If one or more formulas are given, they serve as filter to identify which states are relevant.
 * @param description A string description of what is to be computed used for sensible command line output.
 * @param computationCallback A function that performs the actual computation of the state values
 * @param input Where the properties are read from
 * @param verificationCallback Function to perform the actual verification task for a given formula plus a filter formula to identify relevant states
 * @param postprocessingCallback Function that processes the verification result, such as e.g. output to a file
 */
template<typename ValueType>
void computeStateValues(std::string const& description, std::function<std::unique_ptr<storm::modelchecker::CheckResult>()> const& computationCallback,
                        SymbolicInput const& input, VerificationCallbackType const& verificationCallback,
                        PostprocessingCallbackType const& postprocessingCallback = PostprocessingIdentity());

std::vector<storm::expressions::Expression> parseConstraints(storm::expressions::ExpressionManager const& expressionManager,
                                                             std::string const& constraintsString);

std::vector<std::vector<storm::expressions::Expression>> parseInjectedRefinementPredicates(storm::expressions::ExpressionManager const& expressionManager,
                                                                                           std::string const& refinementPredicatesString);

template<storm::dd::DdType DdType, typename ValueType>
void verifyWithAbstractionRefinementEngine(SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<typename ValueType>
void verifyWithExplorationEngine(SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<typename ValueType>
void verifyWithSparseEngine(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
void verifyWithHybridEngine(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
void verifyWithDdEngine(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
void verifyWithAbstractionRefinementEngine(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input,
                                           ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
typename std::enable_if<DdType != storm::dd::DdType::CUDD || std::is_same<ValueType, double>::value, void>::type verifySymbolicModel(
    std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
typename std::enable_if<DdType == storm::dd::DdType::CUDD && !std::is_same<ValueType, double>::value, void>::type verifySymbolicModel(
    std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename ValueType>
void verifyModel(std::shared_ptr<storm::models::ModelBase> const& model, SymbolicInput const& input, ModelProcessingInformation const& mpi);

template<storm::dd::DdType DdType, typename BuildValueType, typename VerificationValueType = BuildValueType>
std::shared_ptr<storm::models::ModelBase> buildPreprocessModelWithValueTypeAndDdlib(SymbolicInput const& input, ModelProcessingInformation const& mpi) {
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    auto buildSettings = storm::settings::getModule<storm::settings::modules::BuildSettings>();
    std::shared_ptr<storm::models::ModelBase> model;
    if (!buildSettings.isNoBuildModelSet()) {
        model = buildModel<DdType, BuildValueType>(input, ioSettings, mpi);
    }

    if (model) {
        model->printModelInformationToStream(std::cout);
    }

    STORM_LOG_THROW(model || input.properties.empty(), storm::exceptions::InvalidSettingsException, "No input model.");

    if (model) {
        auto preprocessingResult = preprocessModel<DdType, BuildValueType, VerificationValueType>(model, input, mpi);
        if (preprocessingResult.second) {
            model = preprocessingResult.first;
            model->printModelInformationToStream(std::cout);
        }
    }
    return model;
}

template<storm::dd::DdType DdType, typename BuildValueType, typename VerificationValueType = BuildValueType>
std::shared_ptr<storm::models::ModelBase> buildPreprocessExportModelWithValueTypeAndDdlib(SymbolicInput const& input, ModelProcessingInformation const& mpi) {
    auto model = buildPreprocessModelWithValueTypeAndDdlib<DdType, BuildValueType, VerificationValueType>(input, mpi);
    if (model) {
        exportModel<DdType, BuildValueType>(model, input);
    }
    return model;
}

template<storm::dd::DdType DdType, typename BuildValueType, typename VerificationValueType = BuildValueType>
void processInputWithValueTypeAndDdlib(SymbolicInput const& input, ModelProcessingInformation const& mpi) {
    auto abstractionSettings = storm::settings::getModule<storm::settings::modules::AbstractionSettings>();
    auto counterexampleSettings = storm::settings::getModule<storm::settings::modules::CounterexampleGeneratorSettings>();

    // For several engines, no model building step is performed, but the verification is started right away.
    if (mpi.engine == storm::utility::Engine::AbstractionRefinement &&
        abstractionSettings.getAbstractionRefinementMethod() == storm::settings::modules::AbstractionSettings::Method::Games) {
        verifyWithAbstractionRefinementEngine<DdType, VerificationValueType>(input, mpi);
    } else if (mpi.engine == storm::utility::Engine::Exploration) {
        verifyWithExplorationEngine<VerificationValueType>(input, mpi);
    } else {
        std::shared_ptr<storm::models::ModelBase> model =
            buildPreprocessExportModelWithValueTypeAndDdlib<DdType, BuildValueType, VerificationValueType>(input, mpi);
        if (model) {
            if (counterexampleSettings.isCounterexampleSet()) {
                generateCounterexamples<VerificationValueType>(model, input);
            } else {
                verifyModel<DdType, VerificationValueType>(model, input, mpi);
            }
        }
    }
}

template<typename ValueType>
void processInputWithValueType(SymbolicInput const& input, ModelProcessingInformation const& mpi) {
    if (mpi.ddType == storm::dd::DdType::CUDD) {
        STORM_LOG_ASSERT(mpi.verificationValueType == ModelProcessingInformation::ValueType::FinitePrecision &&
                             mpi.buildValueType == ModelProcessingInformation::ValueType::FinitePrecision && (std::is_same<ValueType, double>::value),
                         "Unexpected value type for Dd library cudd.");
        processInputWithValueTypeAndDdlib<storm::dd::DdType::CUDD, double>(input, mpi);
    } else {
        STORM_LOG_ASSERT(mpi.ddType == storm::dd::DdType::Sylvan, "Unknown DD library.");
        if (mpi.buildValueType == mpi.verificationValueType) {
            processInputWithValueTypeAndDdlib<storm::dd::DdType::Sylvan, ValueType>(input, mpi);
        } else {
            // Right now, we only require (buildType == Exact and verificationType == FinitePrecision).
            // We exclude all other combinations to safe a few template instantiations.
            STORM_LOG_THROW((std::is_same<ValueType, double>::value) && mpi.buildValueType == ModelProcessingInformation::ValueType::Exact,
                            storm::exceptions::InvalidArgumentException, "Unexpected combination of buildValueType and verificationValueType");
            processInputWithValueTypeAndDdlib<storm::dd::DdType::Sylvan, storm::RationalNumber, double>(input, mpi);
        }
    }
}
}  // namespace storm::cli
