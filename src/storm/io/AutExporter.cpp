#include "AutExporter.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm {
namespace exporter {

template<typename ValueType>
void exportModelAsAut(std::ostream& os, std::shared_ptr<storm::models::sparse::Model<ValueType>> sparseModel) {
    STORM_LOG_THROW(sparseModel->getType() == storm::models::ModelType::Dtmc || sparseModel->getType() == storm::models::ModelType::Mdp,
                    storm::exceptions::NotSupportedException, "Model type " << sparseModel->getType() << " is not supported for AUT export.");

    // Write header
    STORM_LOG_THROW(sparseModel->getInitialStates().getNumberOfSetBits() == 1, storm::exceptions::NotSupportedException,
                    "Cannot export model with multiple initial states.");
    os << "des(" << *sparseModel->getInitialStates().begin() << ", " << sparseModel->getNumberOfChoices() << ", " << sparseModel->getNumberOfStates() << ")\n";

    // Iterate over states and export outgoing transitions
    storm::storage::SparseMatrix<ValueType> const& matrix = sparseModel->getTransitionMatrix();
    for (typename storm::storage::SparseMatrix<ValueType>::index_type group = 0; group < matrix.getRowGroupCount(); ++group) {
        os << "(" << group << ",";

        // Write transitions
        typename storm::storage::SparseMatrix<ValueType>::index_type start = matrix.hasTrivialRowGrouping() ? group : matrix.getRowGroupIndices()[group];
        typename storm::storage::SparseMatrix<ValueType>::index_type end = matrix.hasTrivialRowGrouping() ? group + 1 : matrix.getRowGroupIndices()[group + 1];
        // Iterate over all actions
        for (typename storm::storage::SparseMatrix<ValueType>::index_type row = start; row < end; ++row) {
            // Write choice
            os << "\"";
            if (sparseModel->hasChoiceLabeling()) {
                if (sparseModel->getChoiceLabeling().getLabelsOfChoice(row).empty()) {
                    os << "tau";
                }
                bool first = true;
                for (auto const& label : sparseModel->getChoiceLabeling().getLabelsOfChoice(row)) {
                    if (!first) {
                        os << ",";
                        first = false;
                    }
                    os << label;
                }
            } else {
                os << "tau";
            }
            os << "\",";

            // Write distribution
            for (auto it = matrix.begin(row), itNext = matrix.begin(row); it != matrix.end(row); ++it) {
                ++itNext;
                if (itNext == matrix.end(row)) {
                    // Last element
                    os << it->getColumn();
                } else {
                    os << it->getColumn() << " " << it->getValue() << " ";
                }
            }
            os << ")\n";
        }
    }
}

// Template instantiations
template void exportModelAsAut<double>(std::ostream& os, std::shared_ptr<storm::models::sparse::Model<double>> sparseModel);
template void exportModelAsAut<storm::RationalNumber>(std::ostream& os, std::shared_ptr<storm::models::sparse::Model<storm::RationalNumber>> sparseModel);

}  // namespace exporter
}  // namespace storm
