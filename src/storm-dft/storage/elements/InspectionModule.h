#pragma once

#include "DFTChildren.h"

namespace storm::dft {
namespace storage {
namespace elements {

/*!
 * Inspection modules for repairs of DFT elements.
 */
template<typename ValueType>
class InspectionModule : public DFTChildren<ValueType> {
    using DFTElementPointer = std::shared_ptr<DFTElement<ValueType>>;
    using DFTElementVector = std::vector<DFTElementPointer>;

   public:
    /*!
     * Constructor.
     * @param id Id.
     * @param name Name.
     * @param children Children.
     */
    InspectionModule(size_t id, std::string const& name, ValueType rate, unsigned phases, DFTElementVector const& children = {})
        : DFTChildren<ValueType>(id, name, children), mRate(rate), mPhases(phases) {
        // Intentionally left empty.
    }

    /*!
     * Destructor
     */
    virtual ~InspectionModule() = default;

    std::shared_ptr<DFTElement<ValueType>> clone() const override {
        return std::shared_ptr<DFTElement<ValueType>>(new InspectionModule<ValueType>(this->id(), this->name(), this->rate(), this->phases(), {}));
    }

    storm::dft::storage::elements::DFTElementType type() const override {
        return storm::dft::storage::elements::DFTElementType::INSPECTION;
    }

    std::string typestring() const override {
        std::stringstream stream;
        stream << "INSPECTION (" << this->phases() << " phases with rate " << this->rate() << ")";
        return stream.str();
    }

    /*!
     * Return rate between states.
     * @return Rate.
     */
    ValueType const& rate() const {
        return mRate;
    }

    /*!
     * Return number of phases (also called the shape).
     * @return Number of phases.
     */
    unsigned phases() const {
        return mPhases;
    }

    bool isInspectionModule() const override {
        return true;
    }

    void extendSpareModule(std::set<size_t>& elementsInSpareModule) const override {
        // Do nothing
    }

    void checkFails(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        STORM_LOG_ASSERT(queues.failurePropagationDone(), "Failure propagation not finished.");
        // TODO check for repairs
        for (auto const& child : this->children()) {
            if (state.hasFailed(child->id())) {
                // child->repair(state, queues);
            }
        }
    }

    void checkFailsafe(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        // Do nothing
    }

    bool checkDontCareAnymore(storm::dft::storage::DFTState<ValueType>& state,
                              storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        return false;
    }

   protected:
    void fail(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        // Do nothing
    }

    void failsafe(storm::dft::storage::DFTState<ValueType>& state, storm::dft::storage::DFTStateSpaceGenerationQueues<ValueType>& queues) const override {
        // Do nothing
    }

   private:
    ValueType mRate;
    unsigned mPhases;
};

}  // namespace elements
}  // namespace storage
}  // namespace storm::dft
