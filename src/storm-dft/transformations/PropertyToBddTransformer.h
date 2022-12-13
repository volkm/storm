#pragma once

#include "storm-dft/modelchecker/SftBddChecker.h"
#include "storm-dft/storage/DFT.h"
#include "storm/logic/AtomicLabelFormula.h"
#include "storm/logic/BinaryBooleanStateFormula.h"
#include "storm/logic/BoundedUntilFormula.h"
#include "storm/logic/ProbabilityOperatorFormula.h"
#include "storm/logic/StateFormula.h"
#include "storm/logic/UnaryBooleanStateFormula.h"

namespace storm::dft {
namespace transformations {

/*!
 * Transformer from logic properties to BDDs.
 * Used for SFT analysis via BDDs.
 * The transformer requires that the BDDs for the relevant labels have already been constructed.
 */
template<typename ValueType>
class PropertyToBddTransformer {
    using Bdd = typename storm::dft::builder::BddSftModelBuilder<ValueType>::Bdd;
    using BuilderPointer = typename std::shared_ptr<storm::dft::builder::BddSftModelBuilder<ValueType> const>;

   public:
    /*!
     * Check whether a given formula is supported by the BDD analysis.
     * Supported formulas are of type 'P=? [F op t phi]' where op is in {<=, <, =} and phi is a state formula.
     * @param formula Formula.
     * @return True iff the formula is supported.
     */
    static bool canHandle(storm::logic::Formula const& formula);

    /*!
     * Translate formula to BDD.
     * Note the that BDDs for the relevant labels have to be already present.
     * @param formula Formula.
     * @return BDD representing the given formula.
     */
    static Bdd translate(storm::logic::Formula const& formula, BuilderPointer builder);

    /*!
     * Get time bound from a (bounded until) formula.
     * @param formula Formula.
     * @return Time bound.
     */
    static double getTimebound(storm::logic::Formula const& formula);

   private:
    /*!
     * Translate state formula to BDD.
     * @param stateFormula State formula.
     * @param builder BDD builder for SFTs.
     * @param enableNegation Whether negation is enabled. Note that negation only works with timepoints but not time bounds.
     * @return BDD.
     */
    static Bdd translate(storm::logic::StateFormula const& stateFormula, BuilderPointer builder, bool enableNegation);

    /*!
     * Translate binary Boolean state formula to BDD.
     *
     * @param stateFormula Binary Boolean state formula.
     * @param builder BDD builder for SFTs.
     * * @param enableNegation Whether negation is enabled. Note that negation only works with timepoints but not time bounds.
     * @return BDD.
     */
    static Bdd translate(storm::logic::BinaryBooleanStateFormula const& stateFormula, BuilderPointer builder, bool enableNegation);

    /*!
     * Translate unary Boolean state formula to BDD.
     * @param stateFormula Unary Boolean state formula.
     * @param builder BDD builder for SFTs.
     * @param enableNot Whether negation is enabled. Note that negation only works with timepoints but not time bounds.
     * @return BDD.
     */
    static Bdd translate(storm::logic::UnaryBooleanStateFormula const& stateFormula, BuilderPointer builder, bool enableNegation);

    /*!
     * Translate atomic label formula to BDD.
     * @param atomicFormula Atomic label formula.
     * @param builder BDD builder for SFTs.
     * @return BDD.
     */
    static Bdd translate(storm::logic::AtomicLabelFormula const& atomicFormula, BuilderPointer builder);
};

}  // namespace transformations
}  // namespace storm::dft
