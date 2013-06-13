#ifndef STORM_MODELS_ABSTRACTDETERMINISTICMODEL_H_
#define STORM_MODELS_ABSTRACTDETERMINISTICMODEL_H_

#include "AbstractModel.h"

#include <memory>
#include <sstream>

namespace storm {

namespace models {

/*!
 *	@brief	Base class for all deterministic model classes.
 *
 *	This is base class defines a common interface for all deterministic models.
 */
template<class T>
class AbstractDeterministicModel: public AbstractModel<T> {

	public:
		/*! Constructs an abstract determinstic model from the given parameters.
		 * @param transitionMatrix The matrix representing the transitions in the model.
		 * @param stateLabeling The labeling that assigns a set of atomic
		 * propositions to each state.
		 * @param stateRewardVector The reward values associated with the states.
		 * @param transitionRewardMatrix The reward values associated with the transitions of the model.
		 */
		AbstractDeterministicModel(std::shared_ptr<storm::storage::SparseMatrix<T>> transitionMatrix,
			std::shared_ptr<storm::models::AtomicPropositionsLabeling> stateLabeling,
			std::shared_ptr<std::vector<T>> stateRewardVector, std::shared_ptr<storm::storage::SparseMatrix<T>> transitionRewardMatrix)
			: AbstractModel<T>(transitionMatrix, stateLabeling, stateRewardVector, transitionRewardMatrix) {
		}

		/*!
		 * Destructor.
		 */
		virtual ~AbstractDeterministicModel() {
			// Intentionally left empty.
		}

		/*!
		 * Copy Constructor.
		 */
		AbstractDeterministicModel(AbstractDeterministicModel const& other) : AbstractModel<T>(other) {
			// Intentionally left empty.
		}
    
        /*!
         * Returns an iterator to the successors of the given state.
         *
         * @param state The state for which to return the iterator.
         * @return An iterator to the successors of the given state.
         */
        virtual typename storm::storage::SparseMatrix<T>::ConstIndexIterator constStateSuccessorIteratorBegin(uint_fast64_t state) const {
            return this->transitionMatrix->constColumnIteratorBegin(state);
        }
    
        /*!
         * Returns an iterator pointing to the element past the successors of the given state.
         *
         * @param state The state for which to return the iterator.
         * @return An iterator pointing to the element past the successors of the given state.
         */
        virtual typename storm::storage::SparseMatrix<T>::ConstIndexIterator constStateSuccessorIteratorEnd(uint_fast64_t state) const {
            return this->transitionMatrix->constColumnIteratorEnd(state);
        }
    
        virtual void writeDotToStream(std::ostream& outStream, bool includeLabeling = true, storm::storage::BitVector const* subsystem = nullptr, std::vector<T> const* firstValue = nullptr, std::vector<T> const* secondValue = nullptr, std::vector<uint_fast64_t> const* stateColoring = nullptr, std::vector<std::string> const* colors = nullptr, std::vector<uint_fast64_t>* scheduler = nullptr, bool finalizeOutput = true) const override {
            AbstractModel<T>::writeDotToStream(outStream, includeLabeling, subsystem, firstValue, secondValue, stateColoring, colors, scheduler, false);
            
            for (auto const& transition : *this->transitionMatrix) {
                if (transition.value() != storm::utility::constGetZero<T>()) {
                    if (subsystem == nullptr || subsystem->get(transition.column())) {
                        outStream << "\t" << transition.row() << " -> " << transition.column() << " [ label= \"" << transition.value() << "\" ];" << std::endl;
                    }
                }
            }
                        
            if (finalizeOutput) {
                outStream << "}" << std::endl;
            }
        }
};

} // namespace models
} // namespace storm

#endif /* STORM_MODELS_ABSTRACTDETERMINISTICMODEL_H_ */
