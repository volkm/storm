#include "PropertyToBddTransformer.h"

namespace storm::dft {
namespace transformations {

template<typename ValueType>
bool PropertyToBddTransformer<ValueType>::canHandle(storm::logic::Formula const& formula) {
    // TODO replace by declaring supported fragment
    STORM_LOG_THROW(formula.isProbabilityOperatorFormula(), storm::exceptions::NotSupportedException,
                    "Formula '" << formula << "' must be a ProbabilityOperatorFormulas.");
    auto const& subFormula = formula.asProbabilityOperatorFormula().getSubformula();
    STORM_LOG_THROW(subFormula.isBoundedUntilFormula(), storm::exceptions::NotSupportedException,
                    "Subformula of '" << formula << "' must be a bounded until formula.");
    auto const& boundedUntil = subFormula.asBoundedUntilFormula();
    STORM_LOG_THROW(boundedUntil.getLeftSubformula().isTrueFormula(), storm::exceptions::NotSupportedException,
                    "Left side of until formula '" << formula << "' must be true.");
    STORM_LOG_THROW(boundedUntil.getRightSubformula().isStateFormula(), storm::exceptions::NotSupportedException,
                    "Right side of until formula '" << formula << "' must be state formula.");

    STORM_LOG_THROW(boundedUntil.hasUpperBound(), storm::exceptions::NotSupportedException, "Until formula '" << formula << "' must be bounded.");
    if (boundedUntil.hasUpperBound() && boundedUntil.hasLowerBound()) {
        // Check if '[F=t phi]' was used.
        STORM_LOG_THROW(boundedUntil.getUpperBound().evaluateAsDouble() == boundedUntil.getLowerBound().evaluateAsDouble(),
                        storm::exceptions::NotSupportedException, "Only upper bounds are supported for formula '" << formula << "'.");
    }
    return true;
}

template<typename ValueType>
typename PropertyToBddTransformer<ValueType>::Bdd PropertyToBddTransformer<ValueType>::translate(storm::logic::Formula const& formula, BuilderPointer builder) {
    if (!canHandle(formula)) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Formula '" << formula << "' cannot be translated to a BDD.");
    }
    // Formula can be handled -> no checks necessary
    auto const& boundedUntil = formula.asProbabilityOperatorFormula().getSubformula().asBoundedUntilFormula();
    auto const& stateFormula = boundedUntil.getRightSubformula().asStateFormula();

    bool hasSameBounds = false;
    if (boundedUntil.hasLowerBound() && boundedUntil.hasUpperBound()) {
        if (boundedUntil.getUpperBound().evaluateAsDouble() == boundedUntil.getLowerBound().evaluateAsDouble()) {
            hasSameBounds = true;
        }
    }

    return translate(stateFormula, builder, hasSameBounds);
}

template<typename ValueType>
double PropertyToBddTransformer<ValueType>::getTimebound(storm::logic::Formula const& formula) {
    auto const& boundedUntil = formula.asProbabilityOperatorFormula().getSubformula().asBoundedUntilFormula();
    return boundedUntil.getUpperBound().evaluateAsDouble();
}

template<typename ValueType>
typename PropertyToBddTransformer<ValueType>::Bdd PropertyToBddTransformer<ValueType>::translate(storm::logic::StateFormula const& stateFormula,
                                                                                                 BuilderPointer builder, bool enableNegation) {
    if (stateFormula.isBinaryBooleanStateFormula()) {
        return translate(stateFormula.asBinaryBooleanStateFormula(), builder, enableNegation);
    } else if (stateFormula.isAtomicLabelFormula()) {
        return translate(stateFormula.asAtomicLabelFormula(), builder);
    } else if (stateFormula.isUnaryBooleanStateFormula()) {
        return translate(stateFormula.asUnaryBooleanStateFormula(), builder, enableNegation);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Illegal StateFormula '" << stateFormula << "'.");
    }
}

template<typename ValueType>
typename PropertyToBddTransformer<ValueType>::Bdd PropertyToBddTransformer<ValueType>::translate(storm::logic::BinaryBooleanStateFormula const& stateFormula,
                                                                                                 BuilderPointer builder, bool enableNegation) {
    auto const leftBdd{translate(stateFormula.getLeftSubformula().asStateFormula(), builder, enableNegation)};
    auto const rightBdd{translate(stateFormula.getRightSubformula().asStateFormula(), builder, enableNegation)};

    if (stateFormula.isAnd()) {
        return leftBdd & rightBdd;
    } else if (stateFormula.isOr()) {
        return leftBdd | rightBdd;
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Illegal BinaryStateFormula '" << stateFormula << "'.");
    }
}
template<typename ValueType>
typename PropertyToBddTransformer<ValueType>::Bdd PropertyToBddTransformer<ValueType>::translate(storm::logic::UnaryBooleanStateFormula const& stateFormula,
                                                                                                 BuilderPointer builder, bool enableNegation) {
    if (!enableNegation) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException,
                        "Illegal UnaryStateFormula '" << stateFormula << "'. Can only use negation with a formula of the form 'P=? [F=t phi]'");
    }
    auto const subBdd{translate(stateFormula.getSubformula().asStateFormula(), builder, enableNegation)};
    if (stateFormula.isNot()) {
        return !subBdd;
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Illegal UnaryStateFormula '" << stateFormula << "'.");
    }
}
template<typename ValueType>
typename PropertyToBddTransformer<ValueType>::Bdd PropertyToBddTransformer<ValueType>::translate(storm::logic::AtomicLabelFormula const& atomicFormula,
                                                                                                 BuilderPointer builder) {
    auto const label{atomicFormula.getLabel()};
    if (label == "failed") {
        return builder->getBddForTopLevelElement();
    } else if (boost::ends_with(label, "_failed")) {
        auto const name{label.substr(0, label.size() - 7)};
        return builder->getBddForElement(name);
    } else {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Label '" << label << "' is not a supported label.");
    }
}

// Explicitly instantiate the class.
template class PropertyToBddTransformer<double>;
template class PropertyToBddTransformer<RationalFunction>;

}  // namespace transformations
}  // namespace storm::dft
