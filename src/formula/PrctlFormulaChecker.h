#ifndef STORM_FORMULA_PRCTLFORMULACHECKER_H_
#define STORM_FORMULA_PRCTLFORMULACHECKER_H_

#include "src/formula/AbstractFormulaChecker.h"
#include "src/formula/Formulas.h"

#include <iostream>

namespace storm {
namespace formula {

/*!
 *	@brief Checks formulas if they are within PRCTL.
 *
 *	This class implements AbstractFormulaChecker to check if a given formula
 *	is part of PRCTL logic.
 */
template <class T>
class PrctlFormulaChecker : public AbstractFormulaChecker<T> {
	public:
		/*!
		 *	Implementation of AbstractFormulaChecker::conforms() using code
		 *	looking exactly like the sample code given there.
		 *
		 *	We currently allow And, Ap, Eventually, Not, Or,
		 *	ProbabilisticNoBoundOperator.
		 */
		virtual bool conforms(const AbstractFormula<T>* formula) const {
			if (
					dynamic_cast<const And<T>*>(formula) ||
					dynamic_cast<const Ap<T>*>(formula) ||
					dynamic_cast<const Eventually<T>*>(formula) ||
					dynamic_cast<const Not<T>*>(formula) ||
					dynamic_cast<const Or<T>*>(formula) ||
					dynamic_cast<const ProbabilisticNoBoundOperator<T>*>(formula)
				) {
				return formula->conforms(*this);
			}
			return false;
		}
};

} // namespace formula
} // namespace storm

#endif