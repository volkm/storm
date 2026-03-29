#include "storm-config.h"
#include "test/storm_gtest.h"

#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/exceptions/InvalidStateException.h"
#include "storm/exceptions/OutOfRangeException.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/FlexibleSparseMatrix.h"
#include "storm/storage/SparseMatrix.h"

TEST(FlexibleSparseMatrix, Insert) {
    // Build sparse matrix first
    storm::storage::SparseMatrixBuilder<double> matrixBuilder(5, 3, 11, true, true, 3);
    matrixBuilder.newRowGroup(0);
    matrixBuilder.addNextValue(0, 1, 0.7);
    matrixBuilder.addNextValue(0, 2, 0.3);
    matrixBuilder.addNextValue(1, 0, 0.2);
    matrixBuilder.addNextValue(1, 1, 0.4);
    matrixBuilder.addNextValue(1, 2, 0.4);
    matrixBuilder.newRowGroup(2);
    matrixBuilder.addNextValue(2, 0, 0.4);
    matrixBuilder.addNextValue(2, 1, 0.4);
    matrixBuilder.addNextValue(2, 2, 0.2);
    matrixBuilder.addNextValue(3, 0, 1.0);
    matrixBuilder.newRowGroup(4);
    matrixBuilder.addNextValue(4, 0, 0.5);
    matrixBuilder.addNextValue(4, 2, 0.5);
    auto matrix = matrixBuilder.build();
    ASSERT_EQ(5ul, matrix.getRowCount());
    ASSERT_EQ(3ul, matrix.getColumnCount());
    ASSERT_EQ(11ul, matrix.getEntryCount());
    ASSERT_EQ(3ul, matrix.getRowGroupCount());

    auto flexibleMatrix = storm::storage::FlexibleSparseMatrix<double>(matrix);
    EXPECT_EQ(5ul, flexibleMatrix.getRowCount());
    EXPECT_EQ(3ul, flexibleMatrix.getColumnCount());
    EXPECT_EQ(11ul, flexibleMatrix.getNonzeroEntryCount());
    EXPECT_EQ(3ul, flexibleMatrix.getRowGroupCount());

    auto newRow = flexibleMatrix.addNewRowGroup();
    EXPECT_EQ(5ul, newRow);
    storm::storage::FlexibleSparseMatrix<double>::row_type& row = flexibleMatrix.getRow(newRow);
    row.push_back(storm::storage::MatrixEntry<storm::storage::FlexibleSparseMatrix<double>::index_type, double>(3, 1));
    flexibleMatrix.updateDimensions();
    EXPECT_EQ(6ul, flexibleMatrix.getRowCount());
    EXPECT_EQ(4ul, flexibleMatrix.getColumnCount());
    EXPECT_EQ(12ul, flexibleMatrix.getNonzeroEntryCount());
    EXPECT_EQ(4ul, flexibleMatrix.getRowGroupCount());

    newRow = flexibleMatrix.addNewRowsToRowGroup(1, 2);
    EXPECT_EQ(4ul, newRow);
    EXPECT_EQ(8ul, flexibleMatrix.getRowCount());
    EXPECT_EQ(4ul, flexibleMatrix.getColumnCount());
    EXPECT_EQ(12ul, flexibleMatrix.getNonzeroEntryCount());
    EXPECT_EQ(4ul, flexibleMatrix.getRowGroupCount());
    auto rowGroups = flexibleMatrix.getRowGroupIndices();
    EXPECT_EQ(0ul, rowGroups[0]);
    EXPECT_EQ(2ul, rowGroups[1]);
    EXPECT_EQ(6ul, rowGroups[2]);
    EXPECT_EQ(7ul, rowGroups[3]);
    EXPECT_EQ(8ul, rowGroups[4]);

    auto matrix2 = flexibleMatrix.createSparseMatrix();
    EXPECT_EQ(8ul, matrix2.getRowCount());
    EXPECT_EQ(4ul, matrix2.getColumnCount());
    EXPECT_EQ(12ul, matrix2.getNonzeroEntryCount());
    EXPECT_EQ(4ul, matrix2.getRowGroupCount());
    for (size_t i = 0; i < flexibleMatrix.getRowCount(); ++i) {
        auto rowF = flexibleMatrix.getRow(i);
        auto rowM = matrix2.getRow(i);
        for (auto [itF, itM] = std::tuple{rowF.begin(), rowM.begin()}; itF != rowF.end() && itM != rowM.end(); ++itF, ++itM) {
            EXPECT_EQ(itF->getColumn(), itM->getColumn());
            EXPECT_EQ(itF->getValue(), itM->getValue());
        }
    }
}
