#include "BddSftModelBuilder.h"

#include "storm/exceptions/NotSupportedException.h"

namespace storm::dft {
namespace builder {

template<typename ValueType>
BddSftModelBuilder<ValueType>::BddSftModelBuilder(std::shared_ptr<storm::dft::storage::DFT<ValueType>> sft) : sft{std::move(sft)} {
    STORM_LOG_THROW(this->sft->isStatic(), storm::exceptions::NotSupportedException, "Fault tree is not static and cannot be translated into BDDs.");

    // Create variables for the BEs
    for (auto const& be : this->sft->getBasicElements()) {
        // Filter constantBeTrigger
        if (be->name() != "constantBeTrigger") {
            variables.push_back(this->sylvanBddManager.createVariable(be->name()));
        }
    }
}
template<typename ValueType>
void BddSftModelBuilder<ValueType>::buildBdds(storm::dft::utility::RelevantEvents relevantEvents) {
    this->relevantEvents = relevantEvents;
    relevantEventBdds.clear();
    relevantEventBdds[getSft()->getTopLevelElement()->name()] = translate(sft->getTopLevelElement());
}
template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd const& BddSftModelBuilder<ValueType>::getBddForElement(std::string const& element) const {
    STORM_LOG_THROW(relevantEventBdds.count(element) > 0, storm::exceptions::InvalidArgumentException, "BDD for element '" << element << "' was not built.");
    return relevantEventBdds.at(element);
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd const& BddSftModelBuilder<ValueType>::getBddForTopLevelElement() const {
    return this->getBddForElement(getSft()->getTopLevelElement()->name());
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd const& BddSftModelBuilder<ValueType>::getOrCreateBddForTopLevelElement() {
    std::string topName = getSft()->getTopLevelElement()->name();
    if (relevantEventBdds.count(topName) == 0) {
        // Build BDDs
        storm::dft::utility::RelevantEvents relevantEvents{};
        this->buildBdds(relevantEvents);
    }
    return this->getBddForElement(topName);
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd BddSftModelBuilder<ValueType>::translate(
    std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> element) {
    auto isRelevant{relevantEvents.isRelevant(element->name())};

    // Check if BDD already exists
    if (isRelevant) {
        auto const it{relevantEventBdds.find(element->name())};
        if (it != relevantEventBdds.end()) {
            return it->second;
        }
    }

    Bdd elemBdd;
    if (element->isGate()) {
        elemBdd = translate(std::static_pointer_cast<storm::dft::storage::elements::DFTGate<ValueType> const>(element));
    } else if (element->isBasicElement()) {
        elemBdd = translate(std::static_pointer_cast<storm::dft::storage::elements::DFTBE<ValueType> const>(element));
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Element of type '" << element->typestring() << "' cannot be translated into BDD.");
        elemBdd = sylvanBddManager.getZero();
    }

    if (isRelevant) {
        // Store BDD
        relevantEventBdds[element->name()] = elemBdd;
    }
    return elemBdd;
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd BddSftModelBuilder<ValueType>::translate(
    std::shared_ptr<storm::dft::storage::elements::DFTBE<ValueType> const> be) {
    return sylvanBddManager.getPositiveLiteral(be->name());
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd BddSftModelBuilder<ValueType>::translate(
    std::shared_ptr<storm::dft::storage::elements::DFTGate<ValueType> const> gate) {
    STORM_LOG_THROW(gate->isStaticElement(), storm::exceptions::NotSupportedException, "Dynamic gate '" << gate->name() << "' cannot be translated to BDD.");
    switch (gate->type()) {
        case storm::dft::storage::elements::DFTElementType::AND: {
            // used only in conjunctions therefore neutral element -> 1
            auto tmpBdd{sylvanBddManager.getOne()};
            for (auto const& child : gate->children()) {
                tmpBdd &= translate(child);
            }
            return tmpBdd;
        }
        case storm::dft::storage::elements::DFTElementType::OR: {
            // used only in disjunctions therefore neutral element -> 0
            auto tmpBdd{sylvanBddManager.getZero()};
            for (auto const& child : gate->children()) {
                tmpBdd |= translate(child);
            }
            return tmpBdd;
        }
        case storm::dft::storage::elements::DFTElementType::VOT: {
            auto vot = std::static_pointer_cast<storm::dft::storage::elements::DFTVot<ValueType> const>(gate);
            std::vector<Bdd> bdds;
            bdds.reserve(vot->children().size());
            for (auto const& child : vot->children()) {
                bdds.push_back(translate(child));
            }
            auto const tmpBdd{translateVot(0, vot->threshold(), bdds)};
            return tmpBdd;
        }
        default:
            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Gate of type '" << gate->typestring() << "' cannot be translated to BDD.");
            return sylvanBddManager.getZero();
    }
}

template<typename ValueType>
typename BddSftModelBuilder<ValueType>::Bdd BddSftModelBuilder<ValueType>::translateVot(size_t const currentIndex, size_t const threshold,
                                                                                        std::vector<Bdd> const& bdds) const {
    if (threshold == 0) {
        return sylvanBddManager.getOne();
    } else if (currentIndex >= bdds.size()) {
        return sylvanBddManager.getZero();
    }

    auto const notChosenBdd{translateVot(currentIndex + 1, threshold, bdds)};
    auto const chosenBdd{translateVot(currentIndex + 1, threshold - 1, bdds)};

    return bdds[currentIndex].Ite(chosenBdd, notChosenBdd);
}

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType> const> BddSftModelBuilder<ValueType>::getSft() const {
    return sft;
}

template<typename ValueType>
storm::dft::storage::SylvanBddManager const& BddSftModelBuilder<ValueType>::getSylvanBddManager() const {
    return sylvanBddManager;
}

// Explicitly instantiate the class.
template class BddSftModelBuilder<double>;
template class BddSftModelBuilder<storm::RationalFunction>;

}  // namespace builder
}  // namespace storm::dft
