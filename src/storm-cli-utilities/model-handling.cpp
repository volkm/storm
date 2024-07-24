#include "model-handling.h"

#include "storm-parsers/parser/ExpressionParser.h"
#include "storm/settings/modules/HintSettings.h"
#include "storm/storage/Qvbs.h"
#include "storm/storage/jani/localeliminator/JaniLocalEliminator.h"
#include "storm/utility/AutomaticSettings.h"

namespace storm {
namespace cli {

void parseSymbolicModelDescription(storm::settings::modules::IOSettings const& ioSettings, SymbolicInput& input) {
    auto buildSettings = storm::settings::getModule<storm::settings::modules::BuildSettings>();
    if (ioSettings.isPrismOrJaniInputSet()) {
        storm::utility::Stopwatch modelParsingWatch(true);
        if (ioSettings.isPrismInputSet()) {
            input.model =
                storm::api::parseProgram(ioSettings.getPrismInputFilename(), buildSettings.isPrismCompatibilityEnabled(), !buildSettings.isNoSimplifySet());
        } else {
            boost::optional<std::vector<std::string>> propertyFilter;
            if (ioSettings.isJaniPropertiesSet()) {
                if (ioSettings.areJaniPropertiesSelected()) {
                    propertyFilter = ioSettings.getSelectedJaniProperties();
                } else {
                    propertyFilter = boost::none;
                }
            } else {
                propertyFilter = std::vector<std::string>();
            }
            auto janiInput = storm::api::parseJaniModel(ioSettings.getJaniInputFilename(), propertyFilter);
            input.model = std::move(janiInput.first);
            if (ioSettings.isJaniPropertiesSet()) {
                input.properties = std::move(janiInput.second);
            }
        }
        modelParsingWatch.stop();
        STORM_PRINT("Time for model input parsing: " << modelParsingWatch << ".\n\n");
    }
}

void parseProperties(storm::settings::modules::IOSettings const& ioSettings, SymbolicInput& input,
                     boost::optional<std::set<std::string>> const& propertyFilter) {
    if (ioSettings.isPropertySet()) {
        std::vector<storm::jani::Property> newProperties;
        if (input.model) {
            newProperties = storm::api::parsePropertiesForSymbolicModelDescription(ioSettings.getProperty(), input.model.get(), propertyFilter);
        } else {
            newProperties = storm::api::parseProperties(ioSettings.getProperty(), propertyFilter);
        }

        input.properties.insert(input.properties.end(), newProperties.begin(), newProperties.end());
    }
}

SymbolicInput parseSymbolicInputQvbs(storm::settings::modules::IOSettings const& ioSettings) {
    // Parse the model input
    SymbolicInput input;
    storm::storage::QvbsBenchmark benchmark(ioSettings.getQvbsModelName());
    STORM_PRINT_AND_LOG(benchmark.getInfo(ioSettings.getQvbsInstanceIndex(), ioSettings.getQvbsPropertyFilter()));
    storm::utility::Stopwatch modelParsingWatch(true);
    auto janiInput = storm::api::parseJaniModel(benchmark.getJaniFile(ioSettings.getQvbsInstanceIndex()), ioSettings.getQvbsPropertyFilter());
    input.model = std::move(janiInput.first);
    input.properties = std::move(janiInput.second);
    modelParsingWatch.stop();
    STORM_PRINT("Time for model input parsing: " << modelParsingWatch << ".\n\n");

    // Parse additional properties
    boost::optional<std::set<std::string>> propertyFilter = storm::api::parsePropertyFilter(ioSettings.getPropertyFilter());
    parseProperties(ioSettings, input, propertyFilter);

    // Substitute constant definitions
    auto constantDefinitions = input.model.get().parseConstantDefinitions(benchmark.getConstantDefinition(ioSettings.getQvbsInstanceIndex()));
    input.model = input.model.get().preprocess(constantDefinitions);
    if (!input.properties.empty()) {
        input.properties = storm::api::substituteConstantsInProperties(input.properties, constantDefinitions);
    }

    return input;
}

SymbolicInput parseSymbolicInput() {
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    if (ioSettings.isQvbsInputSet()) {
        return parseSymbolicInputQvbs(ioSettings);
    } else {
        // Parse the property filter, if any is given.
        boost::optional<std::set<std::string>> propertyFilter = storm::api::parsePropertyFilter(ioSettings.getPropertyFilter());

        SymbolicInput input;
        parseSymbolicModelDescription(ioSettings, input);
        parseProperties(ioSettings, input, propertyFilter);
        return input;
    }
}

void getModelProcessingInformationAutomatic(SymbolicInput const& input, ModelProcessingInformation& mpi) {
    auto hints = storm::settings::getModule<storm::settings::modules::HintSettings>();

    STORM_LOG_THROW(input.model.is_initialized(), storm::exceptions::InvalidArgumentException, "Automatic engine requires a JANI input model.");
    STORM_LOG_THROW(input.model->isJaniModel(), storm::exceptions::InvalidArgumentException, "Automatic engine requires a JANI input model.");
    std::vector<storm::jani::Property> const& properties =
        input.preprocessedProperties.is_initialized() ? input.preprocessedProperties.get() : input.properties;
    STORM_LOG_THROW(!properties.empty(), storm::exceptions::InvalidArgumentException, "Automatic engine requires a property.");
    STORM_LOG_WARN_COND(properties.size() == 1,
                        "Automatic engine does not support decisions based on multiple properties. Only the first property will be considered.");

    storm::utility::AutomaticSettings as;
    if (hints.isNumberStatesSet()) {
        as.predict(input.model->asJaniModel(), properties.front(), hints.getNumberStates());
    } else {
        as.predict(input.model->asJaniModel(), properties.front());
    }

    mpi.engine = as.getEngine();
    if (as.enableBisimulation()) {
        mpi.applyBisimulation = true;
    }
    if (as.enableExact() && mpi.verificationValueType == ModelProcessingInformation::ValueType::FinitePrecision) {
        mpi.verificationValueType = ModelProcessingInformation::ValueType::Exact;
    }
    STORM_PRINT_AND_LOG("Automatic engine picked the following settings: \n"
                        << "\tengine=" << mpi.engine << std::boolalpha << "\t bisimulation=" << mpi.applyBisimulation
                        << "\t exact=" << (mpi.verificationValueType != ModelProcessingInformation::ValueType::FinitePrecision) << std::noboolalpha << '\n');
}

/*!
 * Sets the model processing information based on the given input.
 * Finding the right model processing information might require a conversion to jani.
 * In this case, the jani conversion is stored in the transformedJaniInput pointer (unless it is null)
 */
ModelProcessingInformation getModelProcessingInformation(SymbolicInput const& input, std::shared_ptr<SymbolicInput> const& transformedJaniInput) {
    ModelProcessingInformation mpi;
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    auto coreSettings = storm::settings::getModule<storm::settings::modules::CoreSettings>();
    auto generalSettings = storm::settings::getModule<storm::settings::modules::GeneralSettings>();
    auto bisimulationSettings = storm::settings::getModule<storm::settings::modules::BisimulationSettings>();

    // Set the engine.
    mpi.engine = coreSettings.getEngine();

    // Set whether bisimulation is to be used.
    mpi.applyBisimulation = generalSettings.isBisimulationSet();

    // Set the value type used for numeric values
    if (generalSettings.isParametricSet()) {
        mpi.verificationValueType = ModelProcessingInformation::ValueType::Parametric;
    } else if (generalSettings.isExactSet()) {
        mpi.verificationValueType = ModelProcessingInformation::ValueType::Exact;
    } else {
        mpi.verificationValueType = ModelProcessingInformation::ValueType::FinitePrecision;
    }
    auto originalVerificationValueType = mpi.verificationValueType;

    // Since the remaining settings could depend on the ones above, we need apply the automatic engine now.
    bool useAutomatic = input.model.is_initialized() && mpi.engine == storm::utility::Engine::Automatic;
    if (useAutomatic) {
        if (input.model->isJaniModel()) {
            // This can potentially overwrite the settings above, but will not overwrite settings that were explicitly set by the user (e.g. we will not disable
            // bisimulation or disable exact arithmetic)
            getModelProcessingInformationAutomatic(input, mpi);
        } else {
            // Transform Prism to jani first
            STORM_LOG_ASSERT(input.model->isPrismProgram(), "Unexpected type of input.");
            SymbolicInput janiInput;
            janiInput.properties = input.properties;
            storm::prism::Program const& prog = input.model.get().asPrismProgram();
            auto modelAndProperties = prog.toJani(input.preprocessedProperties.is_initialized() ? input.preprocessedProperties.get() : input.properties);
            janiInput.model = modelAndProperties.first;
            if (!modelAndProperties.second.empty()) {
                janiInput.preprocessedProperties = std::move(modelAndProperties.second);
            }
            // This can potentially overwrite the settings above, but will not overwrite settings that were explicitly set by the user (e.g. we will not disable
            // bisimulation or disable exact arithmetic)
            getModelProcessingInformationAutomatic(janiInput, mpi);
            if (transformedJaniInput) {
                // We cache the transformation result.
                *transformedJaniInput = std::move(janiInput);
            }
        }
    }

    // Check whether these settings are compatible with the provided input.
    if (input.model) {
        auto checkCompatibleSettings = [&mpi, &input] {
            switch (mpi.verificationValueType) {
                case ModelProcessingInformation::ValueType::Parametric:
                    return storm::utility::canHandle<storm::RationalFunction>(
                        mpi.engine, input.preprocessedProperties.is_initialized() ? input.preprocessedProperties.get() : input.properties, input.model.get());
                case ModelProcessingInformation::ValueType::Exact:
                    return storm::utility::canHandle<storm::RationalNumber>(
                        mpi.engine, input.preprocessedProperties.is_initialized() ? input.preprocessedProperties.get() : input.properties, input.model.get());
                    break;
                case ModelProcessingInformation::ValueType::FinitePrecision:
                    return storm::utility::canHandle<double>(
                        mpi.engine, input.preprocessedProperties.is_initialized() ? input.preprocessedProperties.get() : input.properties, input.model.get());
            }
            return false;
        };
        mpi.isCompatible = checkCompatibleSettings();
        if (!mpi.isCompatible) {
            if (useAutomatic) {
                bool useExact = mpi.verificationValueType != ModelProcessingInformation::ValueType::FinitePrecision;
                STORM_LOG_WARN("The settings picked by the automatic engine (engine="
                               << mpi.engine << ", bisim=" << mpi.applyBisimulation << ", exact=" << useExact
                               << ") are incompatible with this model. Falling back to default settings.");
                mpi.engine = storm::utility::Engine::Sparse;
                mpi.applyBisimulation = false;
                mpi.verificationValueType = originalVerificationValueType;
                // Retry check with new settings
                mpi.isCompatible = checkCompatibleSettings();
            }
        }
    } else {
        // If there is no input model, nothing has to be done, actually
        mpi.isCompatible = true;
    }

    // Set whether a transformation to jani is required or necessary
    mpi.transformToJani = ioSettings.isPrismToJaniSet();
    if (input.model) {
        auto builderType = storm::utility::getBuilderType(mpi.engine);
        bool transformToJaniForDdMA = (builderType == storm::builder::BuilderType::Dd) &&
                                      (input.model->getModelType() == storm::storage::SymbolicModelDescription::ModelType::MA) && (!input.model->isJaniModel());
        STORM_LOG_WARN_COND(mpi.transformToJani || !transformToJaniForDdMA,
                            "Dd-based model builder for Markov Automata is only available for JANI models, automatically converting the input model.");
        mpi.transformToJani |= transformToJaniForDdMA;
    }

    // Set the Valuetype used during model building
    mpi.buildValueType = mpi.verificationValueType;
    if (bisimulationSettings.useExactArithmeticInDdBisimulation()) {
        if (storm::utility::getBuilderType(mpi.engine) == storm::builder::BuilderType::Dd && mpi.applyBisimulation) {
            if (mpi.buildValueType == ModelProcessingInformation::ValueType::FinitePrecision) {
                mpi.buildValueType = ModelProcessingInformation::ValueType::Exact;
            }
        } else {
            STORM_LOG_WARN("Requested using exact arithmetic in Dd bisimulation but no dd bisimulation is applied.");
        }
    }

    // Set the Dd library
    mpi.ddType = coreSettings.getDdLibraryType();
    if (mpi.ddType == storm::dd::DdType::CUDD && coreSettings.isDdLibraryTypeSetFromDefaultValue()) {
        if (!(mpi.buildValueType == ModelProcessingInformation::ValueType::FinitePrecision &&
              mpi.verificationValueType == ModelProcessingInformation::ValueType::FinitePrecision)) {
            STORM_LOG_INFO("Switching to DD library sylvan to allow for rational arithmetic.");
            mpi.ddType = storm::dd::DdType::Sylvan;
        }
    }
    return mpi;
}

void ensureNoUndefinedPropertyConstants(std::vector<storm::jani::Property> const& properties) {
    // Make sure there are no undefined constants remaining in any property.
    for (auto const& property : properties) {
        std::set<storm::expressions::Variable> usedUndefinedConstants = property.getUndefinedConstants();
        if (!usedUndefinedConstants.empty()) {
            std::vector<std::string> undefinedConstantsNames;
            for (auto const& constant : usedUndefinedConstants) {
                undefinedConstantsNames.emplace_back(constant.getName());
            }
            STORM_LOG_THROW(
                false, storm::exceptions::InvalidArgumentException,
                "The property '" << property << " still refers to the undefined constants " << boost::algorithm::join(undefinedConstantsNames, ",") << ".");
        }
    }
}

std::pair<SymbolicInput, ModelProcessingInformation> preprocessSymbolicInput(SymbolicInput const& input) {
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();

    SymbolicInput output = input;

    // Preprocess properties (if requested)
    if (ioSettings.isPropertiesAsMultiSet()) {
        STORM_LOG_THROW(!input.properties.empty(), storm::exceptions::InvalidArgumentException,
                        "Can not translate properties to multi-objective formula because no properties were specified.");
        output.properties = {storm::api::createMultiObjectiveProperty(output.properties)};
    }

    // Substitute constant definitions in symbolic input.
    std::string constantDefinitionString = ioSettings.getConstantDefinitionString();
    std::map<storm::expressions::Variable, storm::expressions::Expression> constantDefinitions;
    if (output.model) {
        constantDefinitions = output.model.get().parseConstantDefinitions(constantDefinitionString);
        output.model = output.model.get().preprocess(constantDefinitions);
    }
    if (!output.properties.empty()) {
        output.properties = storm::api::substituteConstantsInProperties(output.properties, constantDefinitions);
    }
    ensureNoUndefinedPropertyConstants(output.properties);
    auto transformedJani = std::make_shared<SymbolicInput>();
    ModelProcessingInformation mpi = getModelProcessingInformation(output, transformedJani);

    // Check whether conversion for PRISM to JANI is requested or necessary.
    if (output.model && output.model.get().isPrismProgram()) {
        if (mpi.transformToJani) {
            if (transformedJani->model) {
                // Use the cached transformation if possible
                output = std::move(*transformedJani);
            } else {
                storm::prism::Program const& model = output.model.get().asPrismProgram();
                auto modelAndProperties = model.toJani(output.properties);

                output.model = modelAndProperties.first;

                if (!modelAndProperties.second.empty()) {
                    output.preprocessedProperties = std::move(modelAndProperties.second);
                }
            }
        }
    }

    if (output.model && output.model.get().isJaniModel()) {
        storm::jani::ModelFeatures supportedFeatures = storm::api::getSupportedJaniFeatures(storm::utility::getBuilderType(mpi.engine));
        storm::api::simplifyJaniModel(output.model.get().asJaniModel(), output.properties, supportedFeatures);

        const auto& buildSettings = storm::settings::getModule<storm::settings::modules::BuildSettings>();
        if (buildSettings.isLocationEliminationSet()) {
            auto locationHeuristic = buildSettings.getLocationEliminationLocationHeuristic();
            auto edgesHeuristic = buildSettings.getLocationEliminationEdgesHeuristic();
            output.model->setModel(storm::jani::JaniLocalEliminator::eliminateAutomatically(output.model.get().asJaniModel(), output.properties,
                                                                                            locationHeuristic, edgesHeuristic));
        }
    }

    return {output, mpi};
}

void exportSymbolicInput(SymbolicInput const& input) {
    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    if (input.model && input.model.get().isJaniModel()) {
        storm::storage::SymbolicModelDescription const& model = input.model.get();
        if (ioSettings.isExportJaniDotSet()) {
            storm::api::exportJaniModelAsDot(model.asJaniModel(), ioSettings.getExportJaniDotFilename());
        }
    }
}

std::vector<std::shared_ptr<storm::logic::Formula const>> createFormulasToRespect(std::vector<storm::jani::Property> const& properties) {
    std::vector<std::shared_ptr<storm::logic::Formula const>> result = storm::api::extractFormulasFromProperties(properties);

    for (auto const& property : properties) {
        if (!property.getFilter().getStatesFormula()->isInitialFormula()) {
            result.push_back(property.getFilter().getStatesFormula());
        }
    }

    return result;
}

storm::builder::BuilderOptions createBuildOptionsSparseFromSettings(SymbolicInput const& input) {
    auto buildSettings = storm::settings::getModule<storm::settings::modules::BuildSettings>();
    storm::builder::BuilderOptions options(createFormulasToRespect(input.properties), input.model.get());
    options.setBuildChoiceLabels(options.isBuildChoiceLabelsSet() || buildSettings.isBuildChoiceLabelsSet());
    options.setBuildStateValuations(options.isBuildStateValuationsSet() || buildSettings.isBuildStateValuationsSet());
    options.setBuildAllLabels(options.isBuildAllLabelsSet() || buildSettings.isBuildAllLabelsSet());
    options.setBuildObservationValuations(options.isBuildObservationValuationsSet() || buildSettings.isBuildObservationValuationsSet());
    bool buildChoiceOrigins = options.isBuildChoiceOriginsSet() || buildSettings.isBuildChoiceOriginsSet();
    if (storm::settings::manager().hasModule(storm::settings::modules::CounterexampleGeneratorSettings::moduleName)) {
        auto counterexampleGeneratorSettings = storm::settings::getModule<storm::settings::modules::CounterexampleGeneratorSettings>();
        if (counterexampleGeneratorSettings.isCounterexampleSet()) {
            buildChoiceOrigins |= counterexampleGeneratorSettings.isMinimalCommandSetGenerationSet();
        }
    }
    options.setBuildChoiceOrigins(buildChoiceOrigins);

    if (buildSettings.isApplyNoMaximumProgressAssumptionSet()) {
        options.setApplyMaximalProgressAssumption(false);
    }

    if (buildSettings.isExplorationChecksSet()) {
        options.setExplorationChecks();
    }
    options.setReservedBitsForUnboundedVariables(buildSettings.getBitsForUnboundedVariables());

    options.setAddOutOfBoundsState(buildSettings.isBuildOutOfBoundsStateSet());
    if (buildSettings.isBuildFullModelSet()) {
        options.clearTerminalStates();
        options.setApplyMaximalProgressAssumption(false);
        options.setBuildAllLabels(true);
        options.setBuildAllRewardModels(true);
    }

    if (buildSettings.isAddOverlappingGuardsLabelSet()) {
        options.setAddOverlappingGuardsLabel(true);
    }

    auto ioSettings = storm::settings::getModule<storm::settings::modules::IOSettings>();
    if (ioSettings.isComputeExpectedVisitingTimesSet() || ioSettings.isComputeSteadyStateDistributionSet()) {
        options.clearTerminalStates();
    }
    return options;
}

void printComputingCounterexample(storm::jani::Property const& property) {
    STORM_PRINT("Computing counterexample for property " << *property.getRawFormula() << " ...\n");
}

void printCounterexample(std::shared_ptr<storm::counterexamples::Counterexample> const& counterexample, storm::utility::Stopwatch* watch) {
    if (counterexample) {
        STORM_PRINT(*counterexample << '\n');
        if (watch) {
            STORM_PRINT("Time for computation: " << *watch << ".\n");
        }
    } else {
        STORM_PRINT(" failed.\n");
    }
}

void printModelCheckingProperty(storm::jani::Property const& property) {
    STORM_PRINT("\nModel checking property \"" << property.getName() << "\": " << *property.getRawFormula() << " ...\n");
}

std::vector<storm::expressions::Expression> parseConstraints(storm::expressions::ExpressionManager const& expressionManager,
                                                             std::string const& constraintsString) {
    std::vector<storm::expressions::Expression> constraints;

    std::vector<std::string> constraintsAsStrings;
    boost::split(constraintsAsStrings, constraintsString, boost::is_any_of(","));

    storm::parser::ExpressionParser expressionParser(expressionManager);
    std::unordered_map<std::string, storm::expressions::Expression> variableMapping;
    for (auto const& variableTypePair : expressionManager) {
        variableMapping[variableTypePair.first.getName()] = variableTypePair.first;
    }
    expressionParser.setIdentifierMapping(variableMapping);

    for (auto const& constraintString : constraintsAsStrings) {
        if (constraintString.empty()) {
            continue;
        }

        storm::expressions::Expression constraint = expressionParser.parseFromString(constraintString);
        STORM_LOG_TRACE("Adding special (user-provided) constraint " << constraint << ".");
        constraints.emplace_back(constraint);
    }

    return constraints;
}

std::vector<std::vector<storm::expressions::Expression>> parseInjectedRefinementPredicates(storm::expressions::ExpressionManager const& expressionManager,
                                                                                           std::string const& refinementPredicatesString) {
    std::vector<std::vector<storm::expressions::Expression>> injectedRefinementPredicates;

    storm::parser::ExpressionParser expressionParser(expressionManager);
    std::unordered_map<std::string, storm::expressions::Expression> variableMapping;
    for (auto const& variableTypePair : expressionManager) {
        variableMapping[variableTypePair.first.getName()] = variableTypePair.first;
    }
    expressionParser.setIdentifierMapping(variableMapping);

    std::vector<std::string> predicateGroupsAsStrings;
    boost::split(predicateGroupsAsStrings, refinementPredicatesString, boost::is_any_of(";"));

    if (!predicateGroupsAsStrings.empty()) {
        for (auto const& predicateGroupString : predicateGroupsAsStrings) {
            if (predicateGroupString.empty()) {
                continue;
            }

            std::vector<std::string> predicatesAsStrings;
            boost::split(predicatesAsStrings, predicateGroupString, boost::is_any_of(":"));

            if (!predicatesAsStrings.empty()) {
                injectedRefinementPredicates.emplace_back();
                for (auto const& predicateString : predicatesAsStrings) {
                    storm::expressions::Expression predicate = expressionParser.parseFromString(predicateString);
                    STORM_LOG_TRACE("Adding special (user-provided) refinement predicate " << predicateString << ".");
                    injectedRefinementPredicates.back().emplace_back(predicate);
                }

                STORM_LOG_THROW(!injectedRefinementPredicates.back().empty(), storm::exceptions::InvalidArgumentException,
                                "Expecting non-empty list of predicates to inject for each (mentioned) refinement step.");

                // Finally reverse the list, because we take the predicates from the back.
                std::reverse(injectedRefinementPredicates.back().begin(), injectedRefinementPredicates.back().end());
            }
        }

        // Finally reverse the list, because we take the predicates from the back.
        std::reverse(injectedRefinementPredicates.begin(), injectedRefinementPredicates.end());
    }

    return injectedRefinementPredicates;
}

}  // namespace cli
}  // namespace storm
