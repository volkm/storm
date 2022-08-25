#pragma once

#include "DFTGate.h"

namespace storm::dft {
namespace storage {
namespace elements {

/*!
 * OR gate.
 * Fails if at least one child has failed.
 */
template<typename ValueType>
class DFTOr : public DFTGate<ValueType> {
   public:
    /*!
     * Constructor.
     * @param id Id.
     * @param name Name.
     * @param children Children.
     */
    DFTOr(size_t id, std::string const& name, std::vector<std::shared_ptr<DFTElement<ValueType>>> const& children = {})
        : DFTGate<ValueType>(id, name, children) {
        // Intentionally empty
    }

    std::shared_ptr<DFTElement<ValueType>> clone() const override {
        return std::shared_ptr<DFTElement<ValueType>>(new DFTOr<ValueType>(this->id(), this->name(), {}));
    }

    storm::dft::storage::elements::DFTElementType type() const override {
        return storm::dft::storage::elements::DFTElementType::OR;
    }

    void checkFails(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        STORM_LOG_ASSERT(this->hasFailedChild(state), "No failed child.");
        if (state.isOperational(this->mId)) {
            this->fail(state, queues);
        }
    }

    void checkFailsafe(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        for (auto const& child : this->children()) {
            if (!state.isFailsafe(child->id())) {
                return;
            }
        }
        // All children are failsafe
        this->failsafe(state, queues);
    }

    void checkRepairs(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        if (state.hasFailed(this->mId)) {
            for (auto const& child : this->children()) {
                if (!state.hasFailed(child->id())) {
                    this->repair(state, queues);
                    return;
                }
            }
        }
    }

   private:
    /*!
     * Check whether it has a failed child.
     * @param state Current state of DFT.
     * @return True iff failed child exists.
     */
    bool hasFailedChild(storm::dft::storage::DFTState<ValueType> const& state) const {
        for (auto const& child : this->children()) {
            if (state.hasFailed(child->id())) {
                return true;
            }
        }
        return false;
    }
};

}  // namespace elements
}  // namespace storage
}  // namespace storm::dft
