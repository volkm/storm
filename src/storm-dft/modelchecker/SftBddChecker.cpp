#include <gmm/gmm_std.h>

#include <memory>
#include <vector>

#include "storm-dft/builder/BddSftModelBuilder.h"
#include "storm-dft/modelchecker/SftBddChecker.h"
#include "storm-dft/transformer/PropertyToBddTransformer.h"

namespace storm::dft {
namespace modelchecker {

template<typename ValueType>
SftBddChecker<ValueType>::SftBddChecker(std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType>> builder) : builder{builder} {}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::check(std::vector<std::shared_ptr<storm::logic::Formula const>> const &formulas, size_t const chunksize,
                                                       storm::dft::utility::RelevantEvents relevantEvents) {
    // Mark all events from formulas as relevant (even from formulas which cannot be handled)
    relevantEvents.insertNamesFromProperties(formulas.begin(), formulas.end());
    // Build relevant BDDs
    builder->buildBdds(relevantEvents);

    // Create BDDs for formulas
    std::vector<Bdd> bdds{};
    bdds.reserve(formulas.size());
    std::map<uint64_t, Bdd> bddMap{};
    for (auto const &formula : formulas) {
        auto const &bdd = storm::dft::transformer::PropertyToBddTransformer<ValueType>::translate(*formula, builder);
        bdds.push_back(bdd);
        bddMap[bdd.GetBDD()] = bdd;
    }

    std::map<uint64_t, std::vector<double>> bddToReversedTimepoints{};
    // A vector of timepoints is necessary as formula-BDDs can occur multiple times.
    // The vector is reversed as it later allows to pop the results from the back which is more efficient.
    for (size_t i = 0; i < bdds.size(); ++i) {
        auto const reversedIndex{bdds.size() - 1 - i};
        auto const &bdd{bdds[reversedIndex]};
        auto const &formula{formulas[reversedIndex]};
        auto const timebound{storm::dft::transformer::PropertyToBddTransformer<ValueType>::getTimebound(*formula)};

        bddToReversedTimepoints[bdd.GetBDD()].push_back(timebound);
    }

    std::map<uint64_t, std::vector<double>> bddToReversedProbabilities{};
    for (auto const &pair : bddToReversedTimepoints) {
        auto const bdd{bddMap.at(pair.first)};
        bddToReversedProbabilities[pair.first] = getProbabilitiesAtTimepoints(bdd, pair.second, chunksize);
    }

    std::vector<ValueType> rval{};
    rval.reserve(bdds.size());
    for (size_t i = 0; i < bdds.size(); ++i) {
        auto const &bdd{bdds[i]};
        auto &tmpVec{bddToReversedProbabilities.at(bdd.GetBDD())};
        rval.push_back(tmpVec.back());
        tmpVec.pop_back();
    }

    return rval;
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::recursiveProbability(Bdd const bdd, std::map<uint32_t, ValueType> const &indexToProbability,
                                                         std::map<uint64_t, ValueType> &bddToProbability) const {
    if (bdd.isOne()) {
        return storm::utility::one<ValueType>();
    } else if (bdd.isZero()) {
        return storm::utility::zero<ValueType>();
    }

    // Check if result was already computed before.
    auto const it{bddToProbability.find(bdd.GetBDD())};
    if (it != bddToProbability.end()) {
        return it->second;
    }

    auto const currentVar{bdd.TopVar()};
    auto const currentProbability{indexToProbability.at(currentVar)};

    auto const thenProbability{recursiveProbability(bdd.Then(), indexToProbability, bddToProbability)};
    auto const elseProbability{recursiveProbability(bdd.Else(), indexToProbability, bddToProbability)};

    // P(Ite(x, f1, f2)) = P(x) * P(f1) + P(!x) * P(f2)
    auto const probability{currentProbability * thenProbability + (storm::utility::one<ValueType>() - currentProbability) * elseProbability};
    bddToProbability[bdd.GetBDD()] = probability;
    return probability;
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::recursiveBirnbaumFactor(uint32_t const variableIndex, Bdd const bdd,
                                                            std::map<uint32_t, ValueType> const &indexToProbability,
                                                            std::map<uint64_t, ValueType> &bddToProbability,
                                                            std::map<uint64_t, ValueType> &bddToBirnbaumFactor) const {
    if (bdd.isTerminal()) {
        return storm::utility::zero<ValueType>();
    }

    // Check if result was already computed before.
    auto const it{bddToBirnbaumFactor.find(bdd.GetBDD())};
    if (it != bddToBirnbaumFactor.end()) {
        return it->second;
    }

    auto const currentVar{bdd.TopVar()};
    auto const currentProbability{indexToProbability.at(currentVar)};

    ValueType birnbaumFactor{storm::utility::zero<ValueType>()};

    if (currentVar > variableIndex) {
        return storm::utility::zero<ValueType>();
    } else if (currentVar == variableIndex) {
        auto const thenProbability{recursiveProbability(bdd.Then(), indexToProbability, bddToProbability)};
        auto const elseProbability{recursiveProbability(bdd.Else(), indexToProbability, bddToProbability)};
        birnbaumFactor = thenProbability - elseProbability;
    } else if (currentVar < variableIndex) {
        auto const thenBirnbaumFactor{recursiveBirnbaumFactor(variableIndex, bdd.Then(), indexToProbability, bddToProbability, bddToBirnbaumFactor)};
        auto const elseBirnbaumFactor{recursiveBirnbaumFactor(variableIndex, bdd.Else(), indexToProbability, bddToProbability, bddToBirnbaumFactor)};

        birnbaumFactor = currentProbability * thenBirnbaumFactor + (storm::utility::one<ValueType>() - currentProbability) * elseBirnbaumFactor;
    }

    bddToBirnbaumFactor[bdd.GetBDD()] = birnbaumFactor;
    return birnbaumFactor;
}

template<typename ValueType>
Eigen::ArrayXd const *SftBddChecker<ValueType>::recursiveProbabilities(
    size_t const chunksize, Bdd const bdd, std::map<uint32_t, Eigen::ArrayXd> const &indexToProbabilities,
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToProbabilities) const {
    auto const bddId{bdd.GetBDD()};

    // Check if result was already computed before.
    auto const it{bddToProbabilities.find(bddId)};
    if (it != bddToProbabilities.end() && it->second.first) {
        return &it->second.second;
    }

    auto &bddToProbabilitiesElement{bddToProbabilities[bddId]};
    if (bdd.isOne() || bdd.isZero()) {
        bddToProbabilitiesElement.first = true;
        bddToProbabilitiesElement.second = Eigen::ArrayXd::Constant(chunksize, bdd.isOne() ? 1 : 0);
        return &bddToProbabilitiesElement.second;
    }

    auto const &thenProbabilities{*recursiveProbabilities(chunksize, bdd.Then(), indexToProbabilities, bddToProbabilities)};
    auto const &elseProbabilities{*recursiveProbabilities(chunksize, bdd.Else(), indexToProbabilities, bddToProbabilities)};

    auto const currentVar{bdd.TopVar()};
    auto const &currentProbabilities{indexToProbabilities.at(currentVar)};

    // P(Ite(x, f1, f2)) = P(x) * P(f1) + P(!x) * P(f2)
    bddToProbabilitiesElement.first = true;
    bddToProbabilitiesElement.second = currentProbabilities * thenProbabilities + (1 - currentProbabilities) * elseProbabilities;
    return &bddToProbabilitiesElement.second;
}

template<typename ValueType>
Eigen::ArrayXd const *SftBddChecker<ValueType>::recursiveBirnbaumFactors(
    size_t const chunksize, uint32_t const variableIndex, Bdd const bdd, std::map<uint32_t, Eigen::ArrayXd> const &indexToProbabilities,
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToProbabilities,
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> &bddToBirnbaumFactors) const {
    auto const bddId{bdd.GetBDD()};

    // Check if result was already computed before.
    auto const it{bddToBirnbaumFactors.find(bddId)};
    if (it != bddToBirnbaumFactors.end() && it->second.first) {
        return &it->second.second;
    }

    auto &bddToBirnbaumFactorsElement{bddToBirnbaumFactors[bddId]};
    if (bdd.isTerminal() || bdd.TopVar() > variableIndex) {
        // Return zero vector
        bddToBirnbaumFactorsElement.second = Eigen::ArrayXd::Constant(chunksize, 0);
        return &bddToBirnbaumFactorsElement.second;
    }

    auto const currentVar{bdd.TopVar()};
    auto const &currentProbabilities{indexToProbabilities.at(currentVar)};

    if (currentVar == variableIndex) {
        auto const &thenProbabilities{*recursiveProbabilities(chunksize, bdd.Then(), indexToProbabilities, bddToProbabilities)};
        auto const &elseProbabilities{*recursiveProbabilities(chunksize, bdd.Else(), indexToProbabilities, bddToProbabilities)};

        bddToBirnbaumFactorsElement.first = true;
        bddToBirnbaumFactorsElement.second = thenProbabilities - elseProbabilities;
        return &bddToBirnbaumFactorsElement.second;
    }

    // currentVar < variableIndex
    auto const &thenBirnbaumFactors{
        *recursiveBirnbaumFactors(chunksize, variableIndex, bdd.Then(), indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};
    auto const &elseBirnbaumFactors{
        *recursiveBirnbaumFactors(chunksize, variableIndex, bdd.Else(), indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};

    bddToBirnbaumFactorsElement.first = true;
    bddToBirnbaumFactorsElement.second = currentProbabilities * thenBirnbaumFactors + (1 - currentProbabilities) * elseBirnbaumFactors;
    return &bddToBirnbaumFactorsElement.second;
}

template<typename ValueType>
void SftBddChecker<ValueType>::recursiveMCS(Bdd const bdd, std::vector<uint32_t> &buffer, std::vector<std::vector<uint32_t>> &minimalCutSets) const {
    if (bdd.isOne()) {
        minimalCutSets.push_back(buffer);
    } else if (!bdd.isZero()) {
        auto const currentVar{bdd.TopVar()};

        buffer.push_back(currentVar);
        recursiveMCS(bdd.Then(), buffer, minimalCutSets);
        buffer.pop_back();
        recursiveMCS(bdd.Else(), buffer, minimalCutSets);
    }
}

template<typename ValueType>
std::vector<std::vector<std::string>> SftBddChecker<ValueType>::getMinimalCutSets() {
    std::vector<std::vector<uint32_t>> mcsIndices{getMinimalCutSetsAsIndices()};

    std::vector<std::vector<std::string>> mcs{};
    mcs.reserve(mcsIndices.size());
    while (!mcsIndices.empty()) {
        std::vector<std::string> tmp{};
        tmp.reserve(mcsIndices.back().size());
        for (auto const &beIndex : mcsIndices.back()) {
            tmp.push_back(builder->getName(beIndex));
        }
        mcs.push_back(std::move(tmp));
        mcsIndices.pop_back();
    }

    return mcs;
}

template<typename ValueType>
std::vector<std::vector<uint32_t>> SftBddChecker<ValueType>::getMinimalCutSetsAsIndices() {
    auto const bdd{builder->getOrCreateBddForTopLevelElement().Minsol()};

    std::vector<std::vector<uint32_t>> mcs{};
    std::vector<uint32_t> buffer{};
    recursiveMCS(bdd, buffer, mcs);

    return mcs;
}

template<>
template<typename FuncType>
void SftBddChecker<double>::chunkCalculationTemplate(std::vector<double> const &timepoints, size_t chunksize, FuncType func) const {
    if (chunksize == 0) {
        chunksize = timepoints.size();
    }

    // Caches
    std::map<uint32_t, Eigen::ArrayXd> indexToProbabilities{};

    // The current timepoints we calculate with
    Eigen::ArrayXd timepointsArray{chunksize};

    for (size_t currentIndex = 0; currentIndex < timepoints.size(); currentIndex += chunksize) {
        auto const sizeLeft{timepoints.size() - currentIndex};
        if (sizeLeft < chunksize) {
            chunksize = sizeLeft;
            timepointsArray = Eigen::ArrayXd{chunksize};
        }

        // Update current timepoints
        for (size_t i = 0; i < chunksize; ++i) {
            timepointsArray(i) = timepoints[i + currentIndex];
        }

        // Update the probabilities of the basic elements
        for (auto const &[be, beIndex] : builder->getBeVariables()) {
            // Vectorize known BETypes
            if (be->beType() == storm::dft::storage::elements::BEType::EXPONENTIAL) {
                // Vectorize exponential distribution p(T <= t) = 1 - exp(-lambda*t)
                auto const failureRate{std::static_pointer_cast<storm::dft::storage::elements::BEExponential<double> const>(be)->activeFailureRate()};
                indexToProbabilities[beIndex] = 1 - (-failureRate * timepointsArray).exp();
            } else {
                // Fallback to getUnreliability()
                auto probabilities{timepointsArray};
                for (size_t i = 0; i < chunksize; ++i) {
                    probabilities(i) = be->getUnreliability(timepointsArray(i));
                }
                indexToProbabilities[beIndex] = probabilities;
            }
        }

        func(chunksize, timepointsArray, indexToProbabilities);
    }
}

template<typename ValueType>
template<typename FuncType>
void SftBddChecker<ValueType>::chunkCalculationTemplate(std::vector<ValueType> const &timepoints, size_t chunksize, FuncType func) const {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Chunk calculation not supported for ValueTypes other than double.");
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getProbabilityAtTimebound(Bdd bdd, ValueType timebound) const {
    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &[be, beIndex] : builder->getBeVariables()) {
        indexToProbability[beIndex] = be->getUnreliability(timebound);
    }

    std::map<uint64_t, ValueType> bddToProbability{};
    auto const probability{recursiveProbability(bdd, indexToProbability, bddToProbability)};
    return probability;
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getProbabilitiesAtTimepoints(Bdd bdd, std::vector<ValueType> const &timepoints, size_t chunksize) const {
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToProbabilities{};
    std::vector<ValueType> resultProbabilities{};
    resultProbabilities.reserve(timepoints.size());

    chunkCalculationTemplate(timepoints, chunksize, [&](auto const currentChunksize, auto const &timepointsArray, auto const &indexToProbabilities) {
        // Invalidate BDD cache
        for (auto &i : bddToProbabilities) {
            i.second.first = false;
        }

        // Great care was made so that the pointer returned is always valid
        // and points to an element in bddToProbabilities
        auto const &probabilitiesArray{*recursiveProbabilities(currentChunksize, bdd, indexToProbabilities, bddToProbabilities)};

        // Update result probabilities
        for (size_t i = 0; i < currentChunksize; ++i) {
            resultProbabilities.push_back(probabilitiesArray(i));
        }
    });

    return resultProbabilities;
}

template<typename ValueType>
template<typename FuncType>
ValueType SftBddChecker<ValueType>::getImportanceMeasureAtTimebound(BEPointer be, ValueType timebound, FuncType func) {
    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &[beIt, beItIndex] : builder->getBeVariables()) {
        indexToProbability[beItIndex] = beIt->getUnreliability(timebound);
    }

    auto const bdd{builder->getOrCreateBddForTopLevelElement()};
    auto const beIndex{builder->getIndex(be)};
    std::map<uint64_t, ValueType> bddToProbability{};
    std::map<uint64_t, ValueType> bddToBirnbaumFactor{};
    auto const probability{recursiveProbability(bdd, indexToProbability, bddToProbability)};
    auto const birnbaumFactor{recursiveBirnbaumFactor(beIndex, bdd, indexToProbability, bddToProbability, bddToBirnbaumFactor)};
    auto const &beProbability{indexToProbability[beIndex]};

    return func(beProbability, probability, birnbaumFactor);
}

template<typename ValueType>
template<typename FuncType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllImportanceMeasuresAtTimebound(ValueType timebound, FuncType func) {
    auto const bdd{builder->getOrCreateBddForTopLevelElement()};

    std::vector<ValueType> resultVector{};
    resultVector.reserve(builder->getBeVariables().size());

    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &[be, beIndex] : builder->getBeVariables()) {
        indexToProbability[beIndex] = be->getUnreliability(timebound);
    }
    std::map<uint64_t, ValueType> bddToProbability{};

    auto const probability{recursiveProbability(bdd, indexToProbability, bddToProbability)};

    for (auto const &[be, beIndex] : builder->getBeVariables()) {
        std::map<uint64_t, ValueType> bddToBirnbaumFactor{};
        auto const birnbaumFactor{recursiveBirnbaumFactor(beIndex, bdd, indexToProbability, bddToProbability, bddToBirnbaumFactor)};
        auto const &beProbability{indexToProbability[beIndex]};
        resultVector.push_back(func(beProbability, probability, birnbaumFactor));
    }
    return resultVector;
}

template<typename ValueType>
template<typename FuncType>
std::vector<ValueType> SftBddChecker<ValueType>::getImportanceMeasuresAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize,
                                                                                   FuncType func) {
    auto const bdd{builder->getOrCreateBddForTopLevelElement()};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToProbabilities{};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToBirnbaumFactors{};
    std::vector<ValueType> resultVector{};
    resultVector.reserve(timepoints.size());

    chunkCalculationTemplate(timepoints, chunksize, [&](auto const currentChunksize, auto const &timepointsArray, auto const &indexToProbabilities) {
        // Invalidate BDD caches
        for (auto &i : bddToProbabilities) {
            i.second.first = false;
        }
        for (auto &i : bddToBirnbaumFactors) {
            i.second.first = false;
        }

        // Great care was made so that the pointer returned is always valid
        auto const index{builder->getIndex(be)};
        auto const &probabilitiesArray{*recursiveProbabilities(currentChunksize, bdd, indexToProbabilities, bddToProbabilities)};
        auto const &birnbaumFactorsArray{
            *recursiveBirnbaumFactors(currentChunksize, index, bdd, indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};

        auto const &beProbabilitiesArray{indexToProbabilities.at(index)};
        auto const ImportanceMeasureArray{func(beProbabilitiesArray, probabilitiesArray, birnbaumFactorsArray)};

        // Update result probabilities
        for (size_t i = 0; i < currentChunksize; ++i) {
            resultVector.push_back(ImportanceMeasureArray(i));
        }
    });

    return resultVector;
}

template<typename ValueType>
template<typename FuncType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllImportanceMeasuresAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize,
                                                                                                   FuncType func) {
    auto const bdd{builder->getOrCreateBddForTopLevelElement()};
    auto const basicElements{builder->getBeVariables()};

    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToProbabilities{};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToBirnbaumFactors{};
    std::vector<std::vector<ValueType>> resultVector{};
    resultVector.resize(basicElements.size());
    for (auto &i : resultVector) {
        i.reserve(timepoints.size());
    }

    chunkCalculationTemplate(timepoints, chunksize, [&](auto const currentChunksize, auto const &timepointsArray, auto const &indexToProbabilities) {
        // Invalidate BDD cache
        for (auto &i : bddToProbabilities) {
            i.second.first = false;
        }

        auto const &probabilitiesArray{*recursiveProbabilities(currentChunksize, bdd, indexToProbabilities, bddToProbabilities)};

        size_t resultIndex = 0;
        for (auto const &[be, beIndex] : basicElements) {
            STORM_LOG_ASSERT(beIndex == resultIndex, "Indices should be equal.");
            // Invalidate bdd cache
            for (auto &i : bddToBirnbaumFactors) {
                i.second.first = false;
            }

            // Great care was made so that the pointer returned is always
            // valid and points to an element in bddToProbabilities
            auto const &birnbaumFactorsArray{
                *recursiveBirnbaumFactors(currentChunksize, beIndex, bdd, indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};

            auto const &beProbabilitiesArray{indexToProbabilities.at(beIndex)};

            auto const ImportanceMeasureArray{func(beProbabilitiesArray, probabilitiesArray, birnbaumFactorsArray)};

            // Update result probabilities
            for (size_t i = 0; i < currentChunksize; ++i) {
                resultVector[resultIndex].push_back(ImportanceMeasureArray(i));
            }
            ++resultIndex;
        }
    });

    return resultVector;
}

namespace {

struct BirnbaumFunctor {
    template<typename ValueType>
    constexpr auto operator()(ValueType const &beProbability, ValueType const &probability, ValueType const &birnbaumFactor) const {
        return birnbaumFactor;
    }
};

struct CIFFunctor {
    template<typename ValueType>
    constexpr ValueType operator()(ValueType const &beProbability, ValueType const &probability, ValueType const &birnbaumFactor) const {
        return (beProbability / probability) * birnbaumFactor;
    }
};

struct DIFFunctor {
    template<typename ValueType>
    constexpr ValueType operator()(ValueType const &beProbability, ValueType const &probability, ValueType const &birnbaumFactor) const {
        return beProbability + (beProbability * (1 - beProbability) * birnbaumFactor) / probability;
    }
};

struct RAWFunctor {
    template<typename ValueType>
    constexpr ValueType operator()(ValueType const &beProbability, ValueType const &probability, ValueType const &birnbaumFactor) const {
        return 1 + ((1 - beProbability) * birnbaumFactor) / probability;
    }
};

struct RRWFunctor {
    template<typename ValueType>
    constexpr ValueType operator()(ValueType const &beProbability, ValueType const &probability, ValueType const &birnbaumFactor) const {
        return probability / (probability - beProbability * birnbaumFactor);
    }
};

}  // namespace

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getBirnbaumFactorAtTimebound(BEPointer be, ValueType timebound) {
    return getImportanceMeasureAtTimebound(be, timebound, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllBirnbaumFactorsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getBirnbaumFactorsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(be, timepoints, chunksize, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllBirnbaumFactorsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, BirnbaumFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getCIFAtTimebound(BEPointer be, ValueType timebound) {
    return getImportanceMeasureAtTimebound(be, timebound, CIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllCIFsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, CIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getCIFsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(be, timepoints, chunksize, CIFFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllCIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, CIFFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getDIFAtTimebound(BEPointer be, ValueType timebound) {
    return getImportanceMeasureAtTimebound(be, timebound, DIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllDIFsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, DIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getDIFsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(be, timepoints, chunksize, DIFFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllDIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, DIFFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getRAWAtTimebound(BEPointer be, ValueType timebound) {
    return getImportanceMeasureAtTimebound(be, timebound, RAWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllRAWsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, RAWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getRAWsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(be, timepoints, chunksize, RAWFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllRAWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, RAWFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getRRWAtTimebound(BEPointer be, ValueType timebound) {
    return getImportanceMeasureAtTimebound(be, timebound, RRWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllRRWsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, RRWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getRRWsAtTimepoints(BEPointer be, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(be, timepoints, chunksize, RRWFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllRRWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, RRWFunctor{});
}

// Explicitly instantiate the class
template class SftBddChecker<double>;

}  // namespace modelchecker
}  // namespace storm::dft
