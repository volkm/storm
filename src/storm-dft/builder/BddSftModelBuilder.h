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
    using BEPointer = std::shared_ptr<storm::dft::storage::elements::DFTBE<ValueType> const>;

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
     * Assumes that buildBdds() was called before.
     *
     * @return BDD for top level element.
     */
    Bdd const& getBddForTopLevelElement() const;

    /*!
     * Get the BDD representing the function of the top level element.
     * Calls buildBdds() to create the BDD for the top level element if it was not already present.
     *
     * @return BDD for top level element.
     */
    Bdd const& getOrCreateBddForTopLevelElement();

    /*!
     * Get index of variable for given BE.
     * @param be Basic event.
     * @return Index.
     */
    uint32_t getIndex(std::shared_ptr<storm::dft::storage::elements::DFTBE<ValueType> const> be) const;

    /*!
     * Get name of BE corresponding to given variable index.
     * @param index Variable index.
     * @return Name of BE.
     */
    std::string getName(uint32_t const index) const;

    /*!
     * Get list of BEs and the index of their corresponding BDD variables.
     * @return List of tuples (BE, index of corresponding BDD variable).
     */
    std::vector<std::pair<BEPointer, uint32_t>> const& getBeVariables() const;

    /*!
     * Get BDD manager associated with this builder.
     * @return BDD manager.
     */
    storm::dft::storage::SylvanBddManager& getBddManager();

    /*!
     * Export the given BDD to a file in the dot format.
     * @param bdd BDD.
     * @param filename File to export to.
     */
    void exportToDot(sylvan::Bdd const& bdd, std::string const& filename);

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

    storm::dft::storage::SylvanBddManager sylvanBddManager;
    std::shared_ptr<storm::dft::storage::DFT<ValueType>> sft;
    storm::dft::utility::RelevantEvents relevantEvents;
    std::vector<std::pair<BEPointer, uint32_t>> beVariables{};
    std::map<std::string, Bdd> relevantEventBdds{};
};

}  // namespace builder
}  // namespace storm::dft
