#include "analysis.h"

#include "storm-dft/builder/BddSftModelBuilder.h"
#include "storm-dft/modelchecker/DFTASFChecker.h"
#include "storm-dft/modelchecker/DftModularizationChecker.h"
#include "storm-dft/modelchecker/SftBddChecker.h"
#include "storm-dft/utility/FDEPConflictFinder.h"
#include "storm-dft/utility/FailureBoundFinder.h"

namespace storm::dft {
namespace api {

template<typename ValueType>
std::pair<uint64_t, uint64_t> computeBEFailureBounds(storm::dft::storage::DFT<ValueType> const& dft, bool useSMT, double solverTimeout) {
    uint64_t lowerBEBound = storm::dft::utility::FailureBoundFinder::getLeastFailureBound(dft, useSMT, solverTimeout);
    uint64_t upperBEBound = storm::dft::utility::FailureBoundFinder::getAlwaysFailedBound(dft, useSMT, solverTimeout);
    return std::make_pair(lowerBEBound, upperBEBound);
}

template<typename ValueType>
bool computeDependencyConflicts(storm::dft::storage::DFT<ValueType>& dft, bool useSMT, double solverTimeout) {
    // Initialize which DFT elements have dynamic behavior
    dft.setDynamicBehaviorInfo();

    std::vector<std::pair<uint64_t, uint64_t>> fdepConflicts =
        storm::dft::utility::FDEPConflictFinder<ValueType>::getDependencyConflicts(dft, useSMT, solverTimeout);

    for (auto const& pair : fdepConflicts) {
        STORM_LOG_DEBUG("Conflict between " << dft.getElement(pair.first)->name() << " and " << dft.getElement(pair.second)->name());
    }

    // Set the conflict map of the dft
    std::set<size_t> conflict_set;
    for (auto const& conflict : fdepConflicts) {
        conflict_set.insert(size_t(conflict.first));
        conflict_set.insert(size_t(conflict.second));
    }
    for (size_t depId : dft.getDependencies()) {
        if (!conflict_set.count(depId)) {
            dft.setDependencyNotInConflict(depId);
        }
    }
    return !fdepConflicts.empty();
}

template<>
void analyzeDFTBdd(std::shared_ptr<storm::dft::storage::DFT<double>> const& dft, bool exportToDot, std::string const& filename, bool calculateMttf,
                   double mttfPrecision, double mttfStepsize, storm::dft::utility::MTTFApproximationAlgorithm mttfAlgorithm, bool calculateMCS,
                   bool calculateProbability, bool useModularisation, std::string const& importanceMeasureName, std::vector<double> const& timepoints,
                   std::vector<std::shared_ptr<storm::logic::Formula const>> const& properties, std::vector<std::string> const& additionalRelevantEventNames,
                   size_t chunksize) {
    if (calculateMttf) {
        switch (mttfAlgorithm) {
            case storm::dft::utility::MTTFApproximationAlgorithm::Proceeding:
                std::cout << "The numerically approximated MTTF is " << storm::dft::utility::MTTFHelperProceeding(dft, mttfStepsize, mttfPrecision) << '\n';
                break;
            case storm::dft::utility::MTTFApproximationAlgorithm::VariableChange:
                std::cout << "The numerically approximated MTTF is " << storm::dft::utility::MTTFHelperVariableChange(dft, mttfStepsize) << '\n';
                break;
            default:
                STORM_LOG_THROW(false, storm::exceptions::InvalidArgumentException, "MTTF approximation algorithm not known.");
        }
    }

    if (useModularisation && calculateProbability) {
        storm::dft::modelchecker::DftModularizationChecker<double> checker{dft};
        if (chunksize == 1) {
            for (auto const& timebound : timepoints) {
                auto const probability{checker.getProbabilityAtTimebound(timebound)};
                std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
            }
        } else {
            auto const probabilities{checker.getProbabilitiesAtTimepoints(timepoints, chunksize)};
            for (size_t i{0}; i < timepoints.size(); ++i) {
                auto const timebound{timepoints[i]};
                auto const probability{probabilities[i]};
                std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
            }
        }
        if (!properties.empty()) {
            auto const probabilities{checker.check(properties, chunksize)};
            for (size_t i{0}; i < probabilities.size(); ++i) {
                std::cout << "Property \"" << properties.at(i)->toString() << "\" has result " << probabilities.at(i) << '\n';
            }
        }
        return;
    } else {
        STORM_LOG_THROW(dft->nrDynamicElements() == 0, storm::exceptions::NotSupportedException,
                        "DFT is dynamic. "
                        "Bdds can only be used on static fault trees. "
                        "Try modularisation.");
    }

    auto builder = std::make_shared<storm::dft::builder::BddSftModelBuilder<double>>(dft);
    builder->getBddManager().execute([&]() {
        storm::dft::utility::RelevantEvents relevantEvents{additionalRelevantEventNames.begin(), additionalRelevantEventNames.end()};
        storm::dft::modelchecker::SftBddChecker checker{builder};

        if (exportToDot) {
            builder->buildBdds(relevantEvents);
            builder->exportToDot(builder->getBddForTopLevelElement(), filename);
        }

        if (calculateMCS) {
            auto const minimalCutSets{checker.getMinimalCutSetsAsIndices()};

            std::cout << "{\n";
            for (auto const& minimalCutSet : minimalCutSets) {
                std::cout << '{';
                for (auto const& beIndex : minimalCutSet) {
                    std::cout << builder->getName(beIndex) << ' ';
                }
                std::cout << "},\n";
            }
            std::cout << "}\n";
        }

        if (calculateProbability) {
            if (chunksize == 1) {
                for (auto const& timebound : timepoints) {
                    auto const probability{checker.getProbabilityAtTimebound(timebound)};
                    std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
                }
            } else {
                auto const probabilities{checker.getProbabilitiesAtTimepoints(timepoints, chunksize)};
                for (size_t i{0}; i < timepoints.size(); ++i) {
                    auto const timebound{timepoints[i]};
                    auto const probability{probabilities[i]};
                    std::cout << "System failure probability at timebound " << timebound << " is " << probability << '\n';
                }
            }

            if (!properties.empty()) {
                auto const probabilities{checker.check(properties, chunksize, relevantEvents)};
                for (size_t i{0}; i < probabilities.size(); ++i) {
                    std::cout << "Property \"" << properties.at(i)->toString() << "\" has result " << probabilities.at(i) << '\n';
                }
            }
        }

        if (importanceMeasureName != "" && timepoints.size() == 1) {
            auto const bes{dft->getBasicElements()};
            std::vector<double> values{};
            if (importanceMeasureName == "MIF") {
                values = checker.getAllBirnbaumFactorsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "CIF") {
                values = checker.getAllCIFsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "DIF") {
                values = checker.getAllDIFsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "RAW") {
                values = checker.getAllRAWsAtTimebound(timepoints[0]);
            }
            if (importanceMeasureName == "RRW") {
                values = checker.getAllRRWsAtTimebound(timepoints[0]);
            }

            for (size_t i{0}; i < bes.size(); ++i) {
                std::cout << importanceMeasureName << " for the basic event " << bes[i]->name() << " at timebound " << timepoints[0] << " is " << values[i]
                          << '\n';
            }
        } else if (importanceMeasureName != "") {
            auto const bes{dft->getBasicElements()};
            std::vector<std::vector<double>> values{};
            if (importanceMeasureName == "MIF") {
                values = checker.getAllBirnbaumFactorsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "CIF") {
                values = checker.getAllCIFsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "DIF") {
                values = checker.getAllDIFsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "RAW") {
                values = checker.getAllRAWsAtTimepoints(timepoints, chunksize);
            }
            if (importanceMeasureName == "RRW") {
                values = checker.getAllRRWsAtTimepoints(timepoints, chunksize);
            }
            for (size_t i{0}; i < bes.size(); ++i) {
                for (size_t j{0}; j < timepoints.size(); ++j) {
                    std::cout << importanceMeasureName << " for the basic event " << bes[i]->name() << " at timebound " << timepoints[j] << " is "
                              << values[i][j] << '\n';
                }
            }
        }
    });
}

template<>
void analyzeDFTBdd(std::shared_ptr<storm::dft::storage::DFT<storm::RationalFunction>> const& dft, bool exportToDot, std::string const& filename,
                   bool calculateMttf, double mttfPrecision, double mttfStepsize, storm::dft::utility::MTTFApproximationAlgorithm mttfAlgorithm,
                   bool calculateMCS, bool calculateProbability, bool useModularisation, std::string const& importanceMeasureName,
                   std::vector<double> const& timepoints, std::vector<std::shared_ptr<storm::logic::Formula const>> const& properties,
                   std::vector<std::string> const& additionalRelevantEventNames, size_t chunksize) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "BDD analysis is not supported for this data type.");
}

template<>
void analyzeDFTSMT(storm::dft::storage::DFT<double> const& dft, bool printOutput) {
    uint64_t solverTimeout = 10;

    storm::dft::modelchecker::DFTASFChecker smtChecker(dft);
    smtChecker.toSolver();
    // Removed bound computation etc. here
    smtChecker.setSolverTimeout(solverTimeout);
    smtChecker.checkTleNeverFailed();
    smtChecker.unsetSolverTimeout();
}

template<>
void analyzeDFTSMT(storm::dft::storage::DFT<storm::RationalFunction> const& dft, bool printOutput) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Analysis by SMT not supported for this data type.");
}

// Explicitly instantiate methods
// Double
template std::pair<uint64_t, uint64_t> computeBEFailureBounds(storm::dft::storage::DFT<double> const&, bool, double);
template bool computeDependencyConflicts(storm::dft::storage::DFT<double>&, bool, double);

// Rational function
template std::pair<uint64_t, uint64_t> computeBEFailureBounds(storm::dft::storage::DFT<storm::RationalFunction> const&, bool, double);
template bool computeDependencyConflicts(storm::dft::storage::DFT<storm::RationalFunction>&, bool, double);

}  // namespace api
}  // namespace storm::dft
