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

    /*!
     * Constructor.
     * @param builder BDD builder for SFT.
     */
    SftBddChecker(std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType>> builder);

    /*!
     * Get the set of minimal cut sets (MCS) by the element names.
     * @return List of MCS given by the names of the SFT elements.
     */
    std::vector<std::vector<std::string>> getMinimalCutSets();

    /*!
     * Get the set of minimal cut sets (MCS) by the element indices in the BDD manager.
     * @return List of MCS given by indices of the SFT elements in the BDD manager.
     */
    std::vector<std::vector<uint32_t>> getMinimalCutSetsAsIndices();

    /*!
     * Compute probability that the top level event fails within the given timebound.
     * @param timebound Time bound.
     * @return Probability at time bound.
     */
    ValueType getProbabilityAtTimebound(ValueType timebound) {
        return getProbabilityAtTimebound(builder->getBddForTopLevelElement(), timebound);
    }

    /*!
     * Compute probability that the given BDD fails within the given timebound.
     * @param bdd The BDD that represents an event in the DFT.
     * @param timebound Time bound.
     * @return Probability at time bound.
     */
    ValueType getProbabilityAtTimebound(Bdd bdd, ValueType timebound) const;

    /*!
     * Compute probabilities that the top level event fails within the given time points.
     * @param timepoints List of time points.
     * @param chunksize Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of probabilities corresponding to given timepoints.
     */
    std::vector<ValueType> getProbabilitiesAtTimepoints(std::vector<ValueType> const &timepoints, size_t const chunksize = 0) {
        return getProbabilitiesAtTimepoints(builder->getBddForTopLevelElement(), timepoints, chunksize);
    }

    /*!
     * Compute probabilities that the given BDD fails within the given time points.
     * @param bdd The BDD that represents an event in the DFT.
     * @param timepoints List of time points.
     * @param chunksize Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     * @return List of probabilities corresponding to given timepoints.
     */
    std::vector<ValueType> getProbabilitiesAtTimepoints(Bdd bdd, std::vector<ValueType> const &timepoints, size_t chunksize = 0) const;

    /**
     * \return
     * The birnbaum importance factor of the given basic event
     * at the given timebound
     */
    ValueType getBirnbaumFactorAtTimebound(std::string const &beName, ValueType timebound);

    /**
     * \return
     * The birnbaum importance factor of all basic events
     * at the given timebound
     *
     * \note
     * Sorted after the order of dft->getBasicElements.
     * Faster than looping over getBirnbaumFactorAtTimebound.
     */
    std::vector<ValueType> getAllBirnbaumFactorsAtTimebound(ValueType timebound);

    /**
     * \return
     * The birnbaum importance factors of the given basic event
     *
     * \param timepoints
     * Array of timebounds to calculate the factors for.
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<ValueType> getBirnbaumFactorsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The birnbaum importance factors of all basic events
     *
     * \param timepoints
     * Array of timebounds to calculate the factors for.
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<std::vector<ValueType>> getAllBirnbaumFactorsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The Critical importance factor of the given basic event
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     */
    ValueType getCIFAtTimebound(std::string const &beName, ValueType timebound);

    /**
     * \return
     * The Critical importance factor of all basic events
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \note
     * Sorted after the order of dft->getBasicElements.
     * Faster than looping over getCIFAtTimebound.
     */
    std::vector<ValueType> getAllCIFsAtTimebound(ValueType timebound);

    /**
     * \return
     * The Critical importance factor of the given basic event
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<ValueType> getCIFsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The Critical importance factor of all basic events
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<std::vector<ValueType>> getAllCIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The Diagnostic importance factor of the given basic event
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     */
    ValueType getDIFAtTimebound(std::string const &beName, ValueType timebound);

    /**
     * \return
     * The Diagnostic importance factor of all basic events
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \note
     * Sorted after the order of dft->getBasicElements.
     * Faster than looping over getDIFAtTimebound.
     */
    std::vector<ValueType> getAllDIFsAtTimebound(ValueType timebound);

    /**
     * \return
     * The Diagnostic importance factor of the given basic event
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<ValueType> getDIFsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The Diagnostic importance factor of all basic events
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<std::vector<ValueType>> getAllDIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The risk achievement worth of the given basic event
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     */
    ValueType getRAWAtTimebound(std::string const &beName, ValueType timebound);

    /**
     * \return
     * The risk achievement worth of all basic events
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \note
     * Sorted after the order of dft->getBasicElements.
     * Faster than looping over getRAWAtTimebound.
     */
    std::vector<ValueType> getAllRAWsAtTimebound(ValueType timebound);

    /**
     * \return
     * The risk achievement worth of the given basic event
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<ValueType> getRAWsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The risk achievement worth of all basic events
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<std::vector<ValueType>> getAllRAWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The risk reduction worth of the given basic event
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     */
    ValueType getRRWAtTimebound(std::string const &beName, ValueType timebound);

    /**
     * \return
     * The risk reduction worth of all basic events
     * at the given timebound as defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \note
     * Sorted after the order of dft->getBasicElements.
     * Faster than looping over getRRWAtTimebound.
     */
    std::vector<ValueType> getAllRRWsAtTimebound(ValueType timebound);

    /**
     * \return
     * The risk reduction worth of the given basic event
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<ValueType> getRRWsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize = 0);

    /**
     * \return
     * The risk reduction worth of all basic events
     * defined in
     * 10.1016/S0951-8320(01)00004-7
     *
     * \param chunksize
     * Splits the timepoints array into chunksize chunks.
     * A value of 0 represents to calculate the whole array at once.
     */
    std::vector<std::vector<ValueType>> getAllRRWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize = 0);

   private:
    /*!
     * Recursively compute probability that the BDD is true given the probabilities that the variables are true.
     * @param bdd BDD.
     * @param indexToProbability Mapping from each variable in the BDD to a probability.
     * @param bddToProbability Cache for common sub-BDDs. It is either empty or filled from an earlier call with an ancestor BDD.
     * @return Probability for given BDD.
     */
    ValueType recursiveProbability(Bdd const bdd, std::map<uint32_t, ValueType> const &indexToProbability,
                                   std::map<uint64_t, ValueType> &bddToProbability) const;

    /*!
     * Recursively compute Birnbaum importance factor for the given variable.
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
     * @param bdd Current BDD.
     * @param buffer Reference to a vector that is used as a stack. The stack temporarily stores the positive variables encountered.
     * @param minimalCutSets Reference to a set of minimal cut sets. It will be populated by the function and contains the MCS in the end.
     */
    void recursiveMCS(Bdd const bdd, std::vector<uint32_t> &buffer, std::vector<std::vector<uint32_t>> &minimalCutSets) const;

    template<typename FuncType>
    void chunkCalculationTemplate(std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func) const;

    template<typename FuncType>
    ValueType getImportanceMeasureAtTimebound(std::string const &beName, ValueType timebound, FuncType func);

    template<typename FuncType>
    std::vector<ValueType> getAllImportanceMeasuresAtTimebound(ValueType timebound, FuncType func);

    template<typename FuncType>
    std::vector<ValueType> getImportanceMeasuresAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize,
                                                             FuncType func);

    template<typename FuncType>
    std::vector<std::vector<ValueType>> getAllImportanceMeasuresAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func);

    std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType>> builder;
};

}  // namespace modelchecker
}  // namespace storm::dft
