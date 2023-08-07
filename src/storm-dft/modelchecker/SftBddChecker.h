#pragma once

#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "storm-dft/builder/BddSftModelBuilder.h"
#include "storm-dft/storage/DFT.h"
#include "storm-dft/storage/SylvanBddManager.h"
#include "storm/adapters/EigenAdapter.h"
#include "storm/storage/PairHash.h"

namespace storm::dft {
namespace modelchecker {

/*!
 * SFT checker based on BDDs.
 */
template<typename ValueType>
class SftBddChecker {
   public:
    using Bdd = sylvan::Bdd;
    using BEPointer = std::shared_ptr<storm::dft::storage::elements::DFTBE<ValueType> const>;

    /*!
     * Constructor.
     *
     * @param builder BDD builder for SFT.
     */
    SftBddChecker(std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType>> builder);

    /*!
     * Get the set of minimal cut sets (MCS) by the element names.
     *
     * @return List of MCS given by the names of the SFT elements.
     */
    std::vector<std::vector<std::string>> getMinimalCutSets();

    /*!
     * Get the set of minimal cut sets (MCS) by the element indices in the BDD manager.
     *
     * @return List of MCS given by indices of the SFT elements in the BDD manager.
     */
    std::vector<std::vector<uint32_t>> getMinimalCutSetsAsIndices();

    /*!
     * Compute probability that the top level event fails within the given time bound.
     *
     * @param timebound Time bound.
     * @return Probability at time bound.
     */
    ValueType getProbabilityAtTimebound(ValueType timebound) {
        return getProbabilityAtTimebound(builder->getOrCreateBddForTopLevelElement(), timebound);
    }

    /*!
     * Compute probability that the given BDD fails within the given time bound.
     *
     * @param bdd The BDD that represents an event in the DFT.
     * @param timebound Time bound.
     * @return Probability at time bound.
     */
    ValueType getProbabilityAtTimebound(Bdd bdd, ValueType timebound) const;

    /*!
     * Compute probabilities that the top level event fails within the given time points.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of probabilities corresponding to given time points.
     */
    std::vector<ValueType> getProbabilitiesAtTimepoints(std::vector<ValueType> const &timepoints, size_t const chunksize = 0) {
        return getProbabilitiesAtTimepoints(builder->getOrCreateBddForTopLevelElement(), timepoints, chunksize);
    }

    /*!
     * Compute probabilities that the given BDD fails within the given time points.
     *
     * @param bdd The BDD that represents an event in the DFT.
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of probabilities corresponding to given time points.
     */
    std::vector<ValueType> getProbabilitiesAtTimepoints(Bdd bdd, std::vector<ValueType> const &timepoints, size_t chunksize = 0) const;

    /*!
     * Check the given formulas.
     *
     * @param formulas Formulas to check.
     * @param chunksize Splits to time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @param relevantEvents Additional relevant events which should be considered while building the BDDs.
     * @return List of values corresponding to the given formulas.
     */
    std::vector<ValueType> check(std::vector<std::shared_ptr<storm::logic::Formula const>> const &formulas, size_t const chunksize = 0,
                                 storm::dft::utility::RelevantEvents relevantEvents = {});

    /*!
     * Compute Birnbaum importance index for given BE at given time bound.
     *
     * @param be BE.
     * @param timebound Time bound.
     * @return Birnbaum factor.
     */
    ValueType getBirnbaumFactorAtTimebound(BEPointer be, ValueType timebound);

    /*!
     * Compute Birnbaum importance indices for all BEs at given time bound.
     * The method is faster than looping over getBirnbaumFactorAtTimebound().
     *
     * @param timebound Time bound.
     * @return List of Birnbaum factors corresponding to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getAllBirnbaumFactorsAtTimebound(ValueType timebound);

    /*!
     * Compute Birnbaum imortance index for given BE at given time points.
     *
     * @param be BE.
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of Birnbaum factors corresponding to given time points.
     */
    std::vector<ValueType> getBirnbaumFactorsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute Birnbaum importance indices for all BEs at given time points.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Birnbaum factors for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<std::vector<ValueType>> getAllBirnbaumFactorsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute the Critical Importance Factor (CIF) for given BE at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param be BE.
     * @param timebound Time bound.
     * @return Critical Importance Factor.
     */
    ValueType getCIFAtTimebound(BEPointer be, ValueType timebound);

    /*!
     * Compute Critical Importance Factors (CIF) for all BEs at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     * The method is faster than looping over getCIFAtTimebound().
     *
     * @param timebound Time bound.
     * @return List of Critical Importance Factors corresponding to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getAllCIFsAtTimebound(ValueType timebound);

    /*!
     * Compute Critical Importance Factors (CIF) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Critical Importance Factors for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getCIFsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute Critical Importance Factors (CIF) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Critical Importance Factors for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<std::vector<ValueType>> getAllCIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute the Diagnostic Importance Factor (DIF) for given BE at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param be BE.
     * @param timebound Time bound.
     * @return Diagnostic Importance Factor.
     */
    ValueType getDIFAtTimebound(BEPointer be, ValueType timebound);

    /*!
     * Compute Diagnostic Importance Factors (DIF) for all BEs at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     * The method is faster than looping over getDIFAtTimebound().
     *
     * @param timebound Time bound.
     * @return List of Diagnostic Importance Factors corresponding to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getAllDIFsAtTimebound(ValueType timebound);

    /*!
     * Compute Diagnostic Importance Factors (DIF) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Diagnostic Importance Factors for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getDIFsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute Diagnostic Importance Factors (DIF) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Diagnostic Importance Factors for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<std::vector<ValueType>> getAllDIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute the Risk Achievement Worth (RAW) for given BE at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param be BE.
     * @param timebound Time bound.
     * @return Risk Achievement Worth.
     */
    ValueType getRAWAtTimebound(BEPointer be, ValueType timebound);

    /*!
     * Compute Risk Achievement Worth (RAW) for all BEs at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     * The method is faster than looping over getRAWAtTimebound().
     *
     * @param timebound Time bound.
     * @return List of Risk Achievement Worth corresponding to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getAllRAWsAtTimebound(ValueType timebound);

    /*!
     * Compute Risk Achievement Worth (RAW) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Risk Achievement Worth for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getRAWsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute Risk Achievement Worth (RAW) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Risk Achievement Worth for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<std::vector<ValueType>> getAllRAWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute the Risk Reduction Worth (RRW) for given BE at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param be BE.
     * @param timebound Time bound.
     * @return Risk Reduction Worth.
     */
    ValueType getRRWAtTimebound(BEPointer be, ValueType timebound);

    /*!
     * Compute Risk Reduction Worth (RRW) for all BEs at given time bound.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     * The method is faster than looping over getRRWAtTimebound().
     *
     * @param timebound Time bound.
     * @return List of Risk Reduction Worth corresponding to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getAllRRWsAtTimebound(ValueType timebound);

    /*!
     * Compute Risk Reduction Worth (RRW) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Risk Reduction Worth for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<ValueType> getRRWsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /*!
     * Compute Risk Reduction Worth (RRW) for all BEs at given time points.
     * Details on the definition are given in https://doi.org/10.1016/S0951-8320(01)00004-7.
     *
     * @param timepoints List of time points.
     * @param chunksize Splits the time-points array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of list of Risk Reduction Worth for all time points per BE.
     * The list corresponds to the order given by dft->getBasicElements().
     */
    std::vector<std::vector<ValueType>> getAllRRWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

   private:
    /*!
     * Recursively compute probability that the BDD is true given the probabilities that the variables are true.
     *
     * @param bdd BDD.
     * @param indexToProbability Mapping from each variable in the BDD to a probability.
     * @param bddToProbability Cache for common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @return Probability for given BDD.
     */
    ValueType recursiveProbability(Bdd const bdd, std::map<uint32_t, ValueType> const &indexToProbability,
                                   std::map<uint64_t, ValueType> &bddToProbability) const;

    /*!
     * Recursively compute Birnbaum importance factor for the given variable.
     *
     * @param variableIndex Index of variable.
     * @param bdd BDD.
     * @param indexToProbability Mapping from each variable in the BDD to a probability.
     * @param bddToProbability Cache for probabilities of common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @param bddToBirnbaumFactor Cache for Birnbaum factor of common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @return Birnbaum importance factor for given variable.
     */
    ValueType recursiveBirnbaumFactor(uint32_t const variableIndex, Bdd const bdd, std::map<uint32_t, ValueType> const &indexToProbability,
                                      std::map<uint64_t, ValueType> &bddToProbability, std::map<uint64_t, ValueType> &bddToBirnbaumFactor) const;

    /*!
     * Recursively compute the (multiple) probabilities that the BDD is true given the probabilities that the variables are true.
     *
     * @param chunksize The width of the Eigen arrays.
     * @param bdd BDD.
     * @param indexToProbabilities Mapping from each variable in the BDD to probabilities.
     * @param bddToProbabilities Cache for probabilities of common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @return Probabilities for the given BDD.
     * @note Great care was made to ensure that all pointers are valid elements in bddToProbabilities.
     */
    Eigen::ArrayXd const *recursiveProbabilities(size_t const chunksize, Bdd const bdd, std::map<uint32_t, Eigen::ArrayXd> const &indexToProbabilities,
                                                 std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToProbabilities) const;

    /*!
     * Recursively compute (multiple) Birnbaum importance factors for the given variable.
     *
     * @param chunksize The width of the Eigen arrays.
     * @param variableIndex Index of variable.
     * @param bdd BDD.
     * @param indexToProbabilities Mapping from each variable in the BDD to probabilities.
     * @param bddToProbabilities Cache for probabilities of common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @param bddToBirnbaumFactors  Cache for Birnbaum factor of common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @return Birnbaum importance factors for given variable.
     */
    Eigen::ArrayXd const *recursiveBirnbaumFactors(size_t const chunksize, uint32_t const variableIndex, Bdd const bdd,
                                                   std::map<uint32_t, Eigen::ArrayXd> const &indexToProbabilities,
                                                   std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToProbabilities,
                                                   std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToBirnbaumFactors) const;

    /*!
     * Recursively traverse the given BDD and return the minimal cut sets.
     *
     * @param bdd Current BDD.
     * @param buffer Reference to a vector that is used as a stack. The stack temporarily stores the positive variables encountered.
     * @param minimalCutSets Reference to a set of minimal cut sets. It will be populated by the function and contains the MCS in the end.
     */
    void recursiveMCS(Bdd const bdd, std::vector<uint32_t> &buffer, std::vector<std::vector<uint32_t>> &minimalCutSets) const;

    template<typename FuncType>
    void chunkCalculationTemplate(std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func) const;

    template<typename FuncType>
    ValueType getImportanceMeasureAtTimebound(BEPointer be, ValueType timebound, FuncType func);

    template<typename FuncType>
    std::vector<ValueType> getAllImportanceMeasuresAtTimebound(ValueType timebound, FuncType func);

    template<typename FuncType>
    std::vector<ValueType> getImportanceMeasuresAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func);

    template<typename FuncType>
    std::vector<std::vector<ValueType>> getAllImportanceMeasuresAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func);

    std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType>> builder;
};

}  // namespace modelchecker
}  // namespace storm::dft
