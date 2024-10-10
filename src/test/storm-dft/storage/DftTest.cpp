#include "storm-config.h"
#include "test/storm_gtest.h"

#include "storm-dft/api/storm-dft.h"

namespace {

TEST(DftTest, Subtree) {
    std::string file = STORM_TEST_RESOURCES_DIR "/dft/modules2.dft";
    std::shared_ptr<storm::dft::storage::DFT<double>> dft = storm::dft::api::loadDFTGalileoFile<double>(file);

    size_t root = dft->getIndex("F3");
    auto subtree = dft->getSubtree(root);
    EXPECT_EQ(subtree.getTopLevelElement()->name(), "F3");
    EXPECT_EQ(subtree.nrElements(), 3ul);
    EXPECT_TRUE(subtree.existsName("x2"));
    EXPECT_TRUE(subtree.existsName("x3"));
    EXPECT_FALSE(subtree.existsName("x1"));

    root = dft->getIndex("F4");
    subtree = dft->getSubtree(root);
    EXPECT_EQ(subtree.getTopLevelElement()->name(), "F4");
    EXPECT_EQ(subtree.nrElements(), 7ul);
    EXPECT_TRUE(subtree.existsName("F5"));
    EXPECT_TRUE(subtree.existsName("F6"));
    EXPECT_TRUE(subtree.existsName("x4"));
    EXPECT_TRUE(subtree.existsName("x5"));
    EXPECT_TRUE(subtree.existsName("x6"));
    EXPECT_TRUE(subtree.existsName("x7"));
    EXPECT_FALSE(subtree.existsName("S1"));

    file = STORM_TEST_RESOURCES_DIR "/dft/spare_two_modules.dft";
    dft = storm::dft::api::loadDFTGalileoFile<double>(file);
    root = dft->getIndex("B");
    subtree = dft->getSubtree(root);
    EXPECT_EQ(subtree.getTopLevelElement()->name(), "B");
    EXPECT_EQ(subtree.nrElements(), 4ul);
    EXPECT_TRUE(subtree.existsName("B"));
    EXPECT_TRUE(subtree.existsName("K"));
    EXPECT_TRUE(subtree.existsName("I"));
    EXPECT_TRUE(subtree.existsName("J"));
}

}  // namespace
