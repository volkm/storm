#include <gmm/gmm_std.h>

#include <memory>
#include <vector>

#include "storm-dft/builder/BddSftModelBuilder.h"
#include "storm-dft/modelchecker/SftBddChecker.h"
#include "storm-dft/transformations/PropertyToBddTransformer.h"

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
    // Init BDD-based checker
    auto checker = std::make_shared<storm::dft::modelchecker::SftBddChecker<ValueType>>(builder);

    // Create BDDs for formulas
    std::vector<Bdd> bdds{};
    bdds.reserve(formulas.size());
    std::map<uint64_t, Bdd> bddMap{};
    for (auto const &formula : formulas) {
        auto const &bdd = storm::dft::transformations::PropertyToBddTransformer<ValueType>::translate(*formula, builder);
        bdds.push_back(bdd);
        bddMap[bdd.GetBDD()] = bdd;
    }

    std::map<uint64_t, std::vector<double>> bddToReversedTimepoints{};
    // A vector of timepoints is necessary as formula-BDDs can occur multiple times.
    // The vector is reversed as it later allows to pop the results from the back which is more efficient.
    for (size_t i{0}; i < bdds.size(); ++i) {
        auto const reversedIndex{bdds.size() - 1 - i};
        auto const &bdd{bdds[reversedIndex]};
        auto const &formula{formulas[reversedIndex]};
        auto const timebound{storm::dft::transformations::PropertyToBddTransformer<ValueType>::getTimebound(*formula)};

        bddToReversedTimepoints[bdd.GetBDD()].push_back(timebound);
    }

    std::map<uint64_t, std::vector<double>> bddToReversedProbabilities{};
    for (auto const &pair : bddToReversedTimepoints) {
        auto const bdd{bddMap.at(pair.first)};
        bddToReversedProbabilities[pair.first] = checker->getProbabilitiesAtTimepoints(bdd, pair.second, chunksize);
    }

    std::vector<ValueType> rval{};
    rval.reserve(bdds.size());
    for (size_t i{0}; i < bdds.size(); ++i) {
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
    auto const it{bddToProbabilities.find(bddId)};
    if (it != bddToProbabilities.end() && it->second.first) {
        return &it->second.second;
    }

    auto &bddToProbabilitiesElement{bddToProbabilities[bddId]};
    if (bdd.isOne()) {
        bddToProbabilitiesElement.first = true;
        bddToProbabilitiesElement.second = Eigen::ArrayXd::Constant(chunksize, 1);
        return &bddToProbabilitiesElement.second;
    } else if (bdd.isZero()) {
        bddToProbabilitiesElement.first = true;
        bddToProbabilitiesElement.second = Eigen::ArrayXd::Constant(chunksize, 0);
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
    auto const it{bddToBirnbaumFactors.find(bddId)};
    if (it != bddToBirnbaumFactors.end() && it->second.first) {
        return &it->second.second;
    }

    auto &bddToBirnbaumFactorsElement{bddToBirnbaumFactors[bddId]};
    if (bdd.isTerminal() || bdd.TopVar() > variableIndex) {
        // return vector 0;
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
std::vector<std::vector<std::string>> SftBddChecker<ValueType>::getMinimalCutSets() {
    std::vector<std::vector<uint32_t>> mcsIndices{getMinimalCutSetsAsIndices()};

    std::vector<std::vector<std::string>> mcs{};
    mcs.reserve(mcsIndices.size());
    while (!mcsIndices.empty()) {
        std::vector<std::string> tmp{};
        tmp.reserve(mcsIndices.back().size());
        for (auto const &be : mcsIndices.back()) {
            tmp.push_back(builder->getSylvanBddManager().getName(be));
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

    // caches
    auto const basicElements{builder->getSft()->getBasicElements()};
    std::map<uint32_t, Eigen::ArrayXd> indexToProbabilities{};

    // The current timepoints we calculate with
    Eigen::ArrayXd timepointsArray{chunksize};

    for (size_t currentIndex{0}; currentIndex < timepoints.size(); currentIndex += chunksize) {
        auto const sizeLeft{timepoints.size() - currentIndex};
        if (sizeLeft < chunksize) {
            chunksize = sizeLeft;
            timepointsArray = Eigen::ArrayXd{chunksize};
        }

        // Update current timepoints
        for (size_t i{currentIndex}; i < currentIndex + chunksize; ++i) {
            timepointsArray(i - currentIndex) = timepoints[i];
        }

        // Update the probabilities of the basic elements
        for (auto const &be : basicElements) {
            auto const beIndex{builder->getSylvanBddManager().getIndex(be->name())};
            // Vectorize known BETypes
            // fallback to getUnreliability() otherwise
            if (be->beType() == storm::dft::storage::elements::BEType::EXPONENTIAL) {
                auto const failureRate{std::static_pointer_cast<storm::dft::storage::elements::BEExponential<double>>(be)->activeFailureRate()};

                // exponential distribution
                // p(T <= t) = 1 - exp(-lambda*t)
                indexToProbabilities[beIndex] = 1 - (-failureRate * timepointsArray).exp();
            } else {
                auto probabilities{timepointsArray};
                for (size_t i{0}; i < chunksize; ++i) {
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
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Chunk calculation not supported for datatypes other than double.");
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getProbabilityAtTimebound(Bdd bdd, ValueType timebound) const {
    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &be : builder->getSft()->getBasicElements()) {
        auto const currentIndex{builder->getSylvanBddManager().getIndex(be->name())};
        indexToProbability[currentIndex] = be->getUnreliability(timebound);
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

        // Update result Probabilities
        for (size_t i{0}; i < currentChunksize; ++i) {
            resultProbabilities.push_back(probabilitiesArray(i));
        }
    });

    return resultProbabilities;
}

template<typename ValueType>
template<typename FuncType>
ValueType SftBddChecker<ValueType>::getImportanceMeasureAtTimebound(std::string const &beName, ValueType timebound, FuncType func) {
    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &be : builder->getSft()->getBasicElements()) {
        auto const currentIndex{builder->getSylvanBddManager().getIndex(be->name())};
        indexToProbability[currentIndex] = be->getUnreliability(timebound);
    }

    auto const bdd{builder->getOrCreateBddForTopLevelElement()};
    auto const index{builder->getSylvanBddManager().getIndex(beName)};
    std::map<uint64_t, ValueType> bddToProbability{};
    std::map<uint64_t, ValueType> bddToBirnbaumFactor{};
    auto const probability{recursiveProbability(bdd, indexToProbability, bddToProbability)};
    auto const birnbaumFactor{recursiveBirnbaumFactor(index, bdd, indexToProbability, bddToProbability, bddToBirnbaumFactor)};
    auto const &beProbability{indexToProbability[index]};

    return func(beProbability, probability, birnbaumFactor);
}

template<typename ValueType>
template<typename FuncType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllImportanceMeasuresAtTimebound(ValueType timebound, FuncType func) {
    auto const bdd{builder->getOrCreateBddForTopLevelElement()};

    std::vector<ValueType> resultVector{};
    resultVector.reserve(builder->getSft()->getBasicElements().size());

    std::map<uint32_t, ValueType> indexToProbability{};
    for (auto const &be : builder->getSft()->getBasicElements()) {
        auto const currentIndex{builder->getSylvanBddManager().getIndex(be->name())};
        indexToProbability[currentIndex] = be->getUnreliability(timebound);
    }
    std::map<uint64_t, ValueType> bddToProbability{};

    auto const probability{recursiveProbability(bdd, indexToProbability, bddToProbability)};

    for (auto const &be : builder->getSft()->getBasicElements()) {
        auto const index{builder->getSylvanBddManager().getIndex(be->name())};
        std::map<uint64_t, ValueType> bddToBirnbaumFactor{};
        auto const birnbaumFactor{recursiveBirnbaumFactor(index, bdd, indexToProbability, bddToProbability, bddToBirnbaumFactor)};
        auto const &beProbability{indexToProbability[index]};
        resultVector.push_back(func(beProbability, probability, birnbaumFactor));
    }
    return resultVector;
}

template<typename ValueType>
template<typename FuncType>
std::vector<ValueType> SftBddChecker<ValueType>::getImportanceMeasuresAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints,
                                                                                   size_t chunksize, FuncType func) {
    auto const bdd{builder->getOrCreateBddForTopLevelElement()};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToProbabilities{};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToBirnbaumFactors{};
    std::vector<ValueType> resultVector{};
    resultVector.reserve(timepoints.size());

    chunkCalculationTemplate(timepoints, chunksize, [&](auto const currentChunksize, auto const &timepointsArray, auto const &indexToProbabilities) {
        // Invalidate bdd caches
        for (auto &i : bddToProbabilities) {
            i.second.first = false;
        }
        for (auto &i : bddToBirnbaumFactors) {
            i.second.first = false;
        }

        // Great care was made so that the pointer returned is always valid
        auto const index{builder->getSylvanBddManager().getIndex(beName)};
        auto const &probabilitiesArray{*recursiveProbabilities(currentChunksize, bdd, indexToProbabilities, bddToProbabilities)};
        auto const &birnbaumFactorsArray{
            *recursiveBirnbaumFactors(currentChunksize, index, bdd, indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};

        auto const &beProbabilitiesArray{indexToProbabilities.at(index)};
        auto const ImportanceMeasureArray{func(beProbabilitiesArray, probabilitiesArray, birnbaumFactorsArray)};

        // Update result Probabilities
        for (size_t i{0}; i < currentChunksize; ++i) {
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
    auto const basicElements{builder->getSft()->getBasicElements()};

    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToProbabilities{};
    std::unordered_map<uint64_t, std::pair<bool, Eigen::ArrayXd>> bddToBirnbaumFactors{};
    std::vector<std::vector<ValueType>> resultVector{};
    resultVector.resize(builder->getSft()->getBasicElements().size());
    for (auto &i : resultVector) {
        i.reserve(timepoints.size());
    }

    chunkCalculationTemplate(timepoints, chunksize, [&](auto const currentChunksize, auto const &timepointsArray, auto const &indexToProbabilities) {
        // Invalidate bdd cache
        for (auto &i : bddToProbabilities) {
            i.second.first = false;
        }

        auto const &probabilitiesArray{*recursiveProbabilities(currentChunksize, bdd, indexToProbabilities, bddToProbabilities)};

        for (size_t basicElementIndex{0}; basicElementIndex < basicElements.size(); ++basicElementIndex) {
            auto const &be{basicElements[basicElementIndex]};
            // Invalidate bdd cache
            for (auto &i : bddToBirnbaumFactors) {
                i.second.first = false;
            }

            // Great care was made so that the pointer returned is always
            // valid and points to an element in bddToProbabilities
            auto const index{builder->getSylvanBddManager().getIndex(be->name())};
            auto const &birnbaumFactorsArray{
                *recursiveBirnbaumFactors(currentChunksize, index, bdd, indexToProbabilities, bddToProbabilities, bddToBirnbaumFactors)};

            auto const &beProbabilitiesArray{indexToProbabilities.at(index)};

            auto const ImportanceMeasureArray{func(beProbabilitiesArray, probabilitiesArray, birnbaumFactorsArray)};

            // Update result Probabilities
            for (size_t i{0}; i < currentChunksize; ++i) {
                resultVector[basicElementIndex].push_back(ImportanceMeasureArray(i));
            }
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
ValueType SftBddChecker<ValueType>::getBirnbaumFactorAtTimebound(std::string const &beName, ValueType timebound) {
    return getImportanceMeasureAtTimebound(beName, timebound, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllBirnbaumFactorsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getBirnbaumFactorsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints,
                                                                                size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(beName, timepoints, chunksize, BirnbaumFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllBirnbaumFactorsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, BirnbaumFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getCIFAtTimebound(std::string const &beName, ValueType timebound) {
    return getImportanceMeasureAtTimebound(beName, timebound, CIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllCIFsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, CIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getCIFsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(beName, timepoints, chunksize, CIFFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllCIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, CIFFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getDIFAtTimebound(std::string const &beName, ValueType timebound) {
    return getImportanceMeasureAtTimebound(beName, timebound, DIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllDIFsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, DIFFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getDIFsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(beName, timepoints, chunksize, DIFFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllDIFsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, DIFFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getRAWAtTimebound(std::string const &beName, ValueType timebound) {
    return getImportanceMeasureAtTimebound(beName, timebound, RAWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllRAWsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, RAWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getRAWsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(beName, timepoints, chunksize, RAWFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllRAWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, RAWFunctor{});
}

template<typename ValueType>
ValueType SftBddChecker<ValueType>::getRRWAtTimebound(std::string const &beName, ValueType timebound) {
    return getImportanceMeasureAtTimebound(beName, timebound, RRWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getAllRRWsAtTimebound(ValueType timebound) {
    return getAllImportanceMeasuresAtTimebound(timebound, RRWFunctor{});
}

template<typename ValueType>
std::vector<ValueType> SftBddChecker<ValueType>::getRRWsAtTimepoints(std::string const &beName, std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getImportanceMeasuresAtTimepoints(beName, timepoints, chunksize, RRWFunctor{});
}

template<typename ValueType>
std::vector<std::vector<ValueType>> SftBddChecker<ValueType>::getAllRRWsAtTimepoints(std::vector<ValueType> const &timepoints, size_t chunksize) {
    return getAllImportanceMeasuresAtTimepoints(timepoints, chunksize, RRWFunctor{});
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

// Explicitly instantiate the class.
template class SftBddChecker<double>;
// template class SftBddChecker<storm::RationalFunction>;

}  // namespace modelchecker
}  // namespace storm::dft
