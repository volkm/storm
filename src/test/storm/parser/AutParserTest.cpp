#include "storm-config.h"
#include "test/storm_gtest.h"

#include "storm-parsers/parser/AutParser.h"
#include "storm/models/sparse/Mdp.h"

TEST(AutParserTest, TwoDiceParsing) {
    std::shared_ptr<storm::models::sparse::Model<double>> modelPtr =
        storm::parser::AutParser<double>::parseModel(STORM_TEST_RESOURCES_DIR "/dtmc/two_dice.aut");

    // Test if parsed correctly.
    ASSERT_EQ(storm::models::ModelType::Dtmc, modelPtr->getType());
    ASSERT_EQ(27ul, modelPtr->getNumberOfStates());
    ASSERT_EQ(54ul, modelPtr->getNumberOfTransitions());
    ASSERT_TRUE(modelPtr->hasLabel("init"));
    ASSERT_EQ(1ul, modelPtr->getInitialStates().getNumberOfSetBits());
    ASSERT_EQ(26ul, *modelPtr->getInitialStates().begin());
    ASSERT_TRUE(modelPtr->hasLabel("s21"));
    ASSERT_EQ(1ul, modelPtr->getStates("s21").getNumberOfSetBits());
    ASSERT_EQ(21ul, *modelPtr->getStates("s21").begin());
}

TEST(AutParserTest, DiningPhilosophersParsing) {
    std::shared_ptr<storm::models::sparse::Model<double>> modelPtr =
        storm::parser::AutParser<double>::parseModel(STORM_TEST_RESOURCES_DIR "/mdp/dining_philosophers.aut");

    // Test if parsed correctly.
    ASSERT_EQ(storm::models::ModelType::Mdp, modelPtr->getType());
    ASSERT_EQ(10ul, modelPtr->getNumberOfStates());
    ASSERT_EQ(12ul, modelPtr->getNumberOfTransitions());
    ASSERT_TRUE(modelPtr->hasLabel("init"));
    ASSERT_EQ(1ul, modelPtr->getInitialStates().getNumberOfSetBits());
    ASSERT_EQ(0ul, *modelPtr->getInitialStates().begin());
    ASSERT_TRUE(modelPtr->hasLabel("s2"));
    ASSERT_EQ(1ul, modelPtr->getStates("s2").getNumberOfSetBits());
    ASSERT_EQ(2ul, *modelPtr->getStates("s2").begin());
}