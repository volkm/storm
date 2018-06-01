#include "storm/modelchecker/multiobjective/multiObjectiveModelChecking.h"

#include "storm/utility/macros.h"
#include "storm/environment/modelchecker/MultiObjectiveModelCheckerEnvironment.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/modelchecker/multiobjective/preprocessing/SparseMultiObjectivePreprocessor.h"
#include "storm/modelchecker/multiobjective/pcaa/SparsePcaaAchievabilityQuery.h"
#include "storm/modelchecker/multiobjective/pcaa/SparsePcaaQuantitativeQuery.h"
#include "storm/modelchecker/multiobjective/pcaa/SparsePcaaParetoQuery.h"
#include "storm/modelchecker/multiobjective/constraintbased/SparseCbAchievabilityQuery.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/CoreSettings.h"
#include "storm/utility/Stopwatch.h"

#include "storm/exceptions/InvalidArgumentException.h"


namespace storm {
    namespace modelchecker {
        namespace multiobjective {
            
            template<typename SparseModelType>
            std::unique_ptr<CheckResult> performMultiObjectiveModelChecking(Environment const& env, SparseModelType const& model, storm::logic::MultiObjectiveFormula const& formula) {
                storm::utility::Stopwatch swTotal(true);
                storm::utility::Stopwatch swPreprocessing(true);
                STORM_LOG_ASSERT(model.getInitialStates().getNumberOfSetBits() == 1, "Multi-objective Model checking on model with multiple initial states is not supported.");
                
                // If we consider an MA, ensure that it is closed
                if(model.isOfType(storm::models::ModelType::MarkovAutomaton)) {
                    STORM_LOG_THROW(dynamic_cast<storm::models::sparse::MarkovAutomaton<typename SparseModelType::ValueType> const *>(&model)->isClosed(), storm::exceptions::InvalidArgumentException, "Unable to check multi-objective formula on non-closed Markov automaton.");
                }
                
                // Preprocess the model
                auto preprocessorResult = preprocessing::SparseMultiObjectivePreprocessor<SparseModelType>::preprocess(model, formula);
                swPreprocessing.stop();
                if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isShowStatisticsSet()) {
                    STORM_PRINT_AND_LOG("Preprocessing done in " << swPreprocessing << " seconds." << std::endl << " Result: " << preprocessorResult << std::endl);
                } else {
                    STORM_LOG_INFO("Preprocessing done in " << swPreprocessing << " seconds." << std::endl << " Result: " << preprocessorResult << std::endl);
                }
                
                
                // Invoke the analysis
                storm::utility::Stopwatch swAnalysis(true);
                std::unique_ptr<CheckResult> result;
                MultiObjectiveMethod method = env.modelchecker().multi().getMethod();
                switch (method) {
                    case MultiObjectiveMethod::Pcaa:
                    {
                        std::unique_ptr<SparsePcaaQuery<SparseModelType, storm::RationalNumber>> query;
                        switch (preprocessorResult.queryType) {
                            case preprocessing::SparseMultiObjectivePreprocessorResult<SparseModelType>::QueryType::Achievability:
                                query = std::unique_ptr<SparsePcaaQuery<SparseModelType, storm::RationalNumber>> (new SparsePcaaAchievabilityQuery<SparseModelType, storm::RationalNumber>(preprocessorResult));
                                break;
                            case preprocessing::SparseMultiObjectivePreprocessorResult<SparseModelType>::QueryType::Quantitative:
                                query = std::unique_ptr<SparsePcaaQuery<SparseModelType, storm::RationalNumber>> (new SparsePcaaQuantitativeQuery<SparseModelType, storm::RationalNumber>(preprocessorResult));
                                break;
                            case preprocessing::SparseMultiObjectivePreprocessorResult<SparseModelType>::QueryType::Pareto:
                                query = std::unique_ptr<SparsePcaaQuery<SparseModelType, storm::RationalNumber>> (new SparsePcaaParetoQuery<SparseModelType, storm::RationalNumber>(preprocessorResult));
                                break;
                            default:
                                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The multi-objective query type is not supported for the selected solution method '" << toString(method) << "'.");
                                break;
                        }
        
                        result = query->check(env);
                                        
                        if (env.modelchecker().multi().isExportPlotSet()) {
                            query->exportPlotOfCurrentApproximation(env);
                        }
                        break;
                    }
                    case MultiObjectiveMethod::ConstraintBased:
                    {
                        std::unique_ptr<SparseCbQuery<SparseModelType>> query;
                        switch (preprocessorResult.queryType) {
                            case preprocessing::SparseMultiObjectivePreprocessorResult<SparseModelType>::QueryType::Achievability:
                                query = std::unique_ptr<SparseCbQuery<SparseModelType>> (new SparseCbAchievabilityQuery<SparseModelType>(preprocessorResult));
                                break;
                            default:
                                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The multi-objective query type is not supported for the selected solution method '" << toString(method) << "'.");
                                break;
                        }
                        
                        result = query->check(env);
                     
                        if (env.modelchecker().multi().isExportPlotSet()) {
                            STORM_LOG_ERROR("Can not export plot for the constrained based solver.");
                        }
                        break;
                    }
                    default:
                        STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "The multi-objective solution method '" << toString(method) << "' is not supported.");
                }
                swAnalysis.stop();
                
                swTotal.stop();
                if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isShowStatisticsSet()) {
                    STORM_PRINT_AND_LOG("Solving multi-objective query took " << swTotal << " seconds (consisting of " << swPreprocessing << " seconds for preprocessing and " << swAnalysis << " seconds for analyzing the preprocessed model)." << std::endl);
                }

                return result;
            }
     
            
            template std::unique_ptr<CheckResult> performMultiObjectiveModelChecking<storm::models::sparse::Mdp<double>>(Environment const& env, storm::models::sparse::Mdp<double> const& model, storm::logic::MultiObjectiveFormula const& formula);
            template std::unique_ptr<CheckResult> performMultiObjectiveModelChecking<storm::models::sparse::MarkovAutomaton<double>>(Environment const& env, storm::models::sparse::MarkovAutomaton<double> const& model, storm::logic::MultiObjectiveFormula const& formula);
            template std::unique_ptr<CheckResult> performMultiObjectiveModelChecking<storm::models::sparse::Mdp<storm::RationalNumber>>(Environment const& env, storm::models::sparse::Mdp<storm::RationalNumber> const& model, storm::logic::MultiObjectiveFormula const& formula);
            template std::unique_ptr<CheckResult> performMultiObjectiveModelChecking<storm::models::sparse::MarkovAutomaton<storm::RationalNumber>>(Environment const& env, storm::models::sparse::MarkovAutomaton<storm::RationalNumber> const& model, storm::logic::MultiObjectiveFormula const& formula);

        }
    }
}
