#pragma once

#include <memory>
#include <vector>

#include "storm-dft/storage/DFT.h"
#include "storm-dft/storage/SylvanBddManager.h"
#include "storm-dft/utility/RelevantEvents.h"

namespace storm::dft {
namespace builder {

/*!
 * Build a BDD for an SFT.
 */
template<typename ValueType>
class BddSftModelBuilder {
   public:
    using Bdd = sylvan::Bdd;

    /*!
     * Constructor.
     * @param sft Static fault tree.
     */
    BddSftModelBuilder(std::shared_ptr<storm::dft::storage::DFT<ValueType>> sft);

    /*!
     * Build BDDs for all relevant events in the SFT.
     *
     * @param relevantEvents Relevant events for which the BDDs should be created.
     */
    void buildBdds(storm::dft::utility::RelevantEvents relevantEvents);

    /*!
     * Get the BDD representing the function of the given SFT element.
     * Assumes that buildBdds() was called before.
     *
     * @return BDD for element.
     */
    Bdd const& getBddForElement(std::string const& element) const;

    /*!
     * Get the BDD representing the function of the top level element.
     * Calls buildBdds() if the BDD for the top level element is not already present.
     *
     * @return BDD for top level element.
     */
    Bdd const& getBddForTopLevelElement();

    /*!
     * Return SFT.
     * @return SFT.
     */
    std::shared_ptr<storm::dft::storage::DFT<ValueType> const> getSft() const;

    /*!
     * Return manager for Sylvan BDDs.
     * @return Manager.
     */
    storm::dft::storage::SylvanBddManager const& getSylvanBddManager() const;

   private:
    /*!
     * Translate an SFT element into a BDD.
     * @param element DFT element.
     * @return BDD.
     */
    Bdd translate(std::shared_ptr<storm::dft::storage::elements::DFTElement<ValueType> const> element);

    /**
     * Translate a DFT basic element gate into a BDD.
     *
     * \note This is the recursion anchor.
     *
     */
    Bdd translate(std::shared_ptr<storm::dft::storage::elements::DFTBE<ValueType> const> be);

    /*!
     * Translate an SFT gate into a BDD.
     * @param gate SFT gate.
     * @return BDD.
     */
    Bdd translate(std::shared_ptr<storm::dft::storage::elements::DFTGate<ValueType> const> gate);

    /*!
     * Helper function to translate a VOTing gate into a BDD using the Shannon decomposition.
     *
     * @param currentIndex Index of the child currently considered.
     * @param threshold How many children need to be chosen.
     * @param bdds Reference to the BDDs corresponding to all children of the VOTing gate.
     * @return BDD representation of VOTing gate.
     */
    Bdd translateVot(size_t currentIndex, size_t threshold, std::vector<Bdd> const& bdds) const;

    std::shared_ptr<storm::dft::storage::DFT<ValueType>> sft;
    storm::dft::utility::RelevantEvents relevantEvents;
    storm::dft::storage::SylvanBddManager sylvanBddManager;
    std::vector<uint32_t> variables{};
    std::map<std::string, Bdd> relevantEventBdds{};
};

}  // namespace builder
}  // namespace storm::dft
