#include "storm-config.h"
#include "test/storm_gtest.h"

#include "storm-dft/api/storm-dft.h"
#include "storm-dft/transformations/PropertyToBddTransformer.h"
#include "storm-parsers/api/properties.h"
#include "storm/api/properties.h"

namespace {

// Helper functions
std::shared_ptr<storm::dft::builder::BddSftModelBuilder<double>> createBdds(std::string const& file) {
    auto dft = storm::dft::api::loadDFTGalileoFile<double>(file);
    auto relevantEvents = storm::dft::api::computeRelevantEvents<double>(*dft, {}, {"all"});
    auto builder = std::make_shared<storm::dft::builder::BddSftModelBuilder<double>>(dft);
    builder->buildBdds(relevantEvents);
    return builder;
}

typename storm::dft::builder::BddSftModelBuilder<double>::Bdd translate(std::shared_ptr<storm::dft::builder::BddSftModelBuilder<double>> builder,
                                                                        std::string const& formula) {
    auto const formulas{storm::api::extractFormulasFromProperties(storm::api::parseProperties(formula))};
    return storm::dft::transformations::PropertyToBddTransformer<double>::translate(*formulas[0], builder);
}

TEST(PropertyToBddTransformerTest, AndOrFormula) {
    auto builder = createBdds(STORM_TEST_RESOURCES_DIR "/dft/bdd/AndOrTest.dft");
    auto bdd = translate(builder, "P=? [F <= 1 \"failed\"]");
    auto bdd2 = translate(builder, "P=? [F <= 1 \"F_failed\"]");
    auto bdd3 = translate(builder, "P=? [F <= 1 \"F1_failed\" & \"F2_failed\"]");
    auto bdd4 = translate(builder, "P=? [F  = 1 !\"failed\"]");
    auto bdd5 = translate(builder, "P=? [F  = 1 !\"F1_failed\"]");
    auto bdd6 = translate(builder, "P=? [F  = 1 !\"F2_failed\"]");
    auto bdd7 = translate(builder, "P=? [F <= 1 \"F1_failed\"]");
    auto bdd8 = translate(builder, "P=? [F <= 1 \"F2_failed\"]");

    EXPECT_EQ(bdd, bdd2);
    EXPECT_EQ(bdd.GetBDD(), bdd2.GetBDD());
    EXPECT_EQ(bdd2, bdd3);
    EXPECT_EQ(bdd2.GetBDD(), bdd3.GetBDD());
    EXPECT_EQ(bdd3, bdd);
    EXPECT_EQ(bdd3.GetBDD(), bdd.GetBDD());

    EXPECT_EQ(bdd.GetBDD(), (!bdd4).GetBDD());
    EXPECT_NE(bdd.GetBDD(), bdd4.GetBDD());

    EXPECT_EQ(bdd7.GetBDD(), (!bdd5).GetBDD());
    EXPECT_NE(bdd7.GetBDD(), bdd5.GetBDD());
    EXPECT_EQ(bdd8.GetBDD(), (!bdd6).GetBDD());
    EXPECT_NE(bdd8.GetBDD(), bdd6.GetBDD());

    EXPECT_NE(bdd4.GetBDD(), bdd5.GetBDD());

    EXPECT_EQ(bdd.GetShaHash(), "fc1e9a418e3c207e81ffa7fde7768f027b6996732c4216c1ed5de6861dbc86ae");
    EXPECT_EQ(bdd2.GetShaHash(), "fc1e9a418e3c207e81ffa7fde7768f027b6996732c4216c1ed5de6861dbc86ae");
    EXPECT_EQ(bdd3.GetShaHash(), "fc1e9a418e3c207e81ffa7fde7768f027b6996732c4216c1ed5de6861dbc86ae");
    EXPECT_EQ(bdd4.GetShaHash(), "fc1e9a418e3c207e81ffa7fde7768f027b6996732c4216c1ed5de6861dbc86ae");
    EXPECT_EQ(bdd5.GetShaHash(), "c5cf2304417926961c3e1ce1d876fc2886ece1365fd946bfd3e1abd71401696d");
    EXPECT_EQ(bdd6.GetShaHash(), "a4f129fa27c6cd32625b088811d4b12f8059ae0547ee035c083deed9ef9d2c59");
    EXPECT_EQ(bdd7.GetShaHash(), "c5cf2304417926961c3e1ce1d876fc2886ece1365fd946bfd3e1abd71401696d");
    EXPECT_EQ(bdd8.GetShaHash(), "a4f129fa27c6cd32625b088811d4b12f8059ae0547ee035c083deed9ef9d2c59");
}

TEST(PropertyToBddTransformerTest, Formulas) {
    auto builder = createBdds(STORM_TEST_RESOURCES_DIR "/dft/bdd/ImportanceTest.dft");
    auto bdd = translate(builder, "P=? [F=1 \"failed\"]");
    auto bdd2 = translate(builder, "P=? [F=1 \"F1_failed\" & \"F2_failed\" & \"x1_failed\"]");
    auto bdd3 = translate(builder, "P=? [F=1 \"F_failed\" & !\"x1_failed\"]");
    auto bdd4 = translate(builder, "P=? [F=1 \"F_failed\" & !\"x2_failed\" & (\"x5_failed\" | \"x6_failed\")]");
    auto bdd5 = translate(builder, "P=? [F=1 \"F_failed\" & !(\"x2_failed\" | (!\"x5_failed\" & !\"x6_failed\"))]");

    EXPECT_EQ(bdd, bdd2);
    EXPECT_EQ(bdd.GetBDD(), bdd2.GetBDD());

    EXPECT_EQ(bdd3, sylvan::Bdd::bddZero());
    EXPECT_TRUE(bdd3.isZero());

    EXPECT_EQ(0ul, bdd3.SatCount(6));
    EXPECT_EQ(1ul, bdd3.NodeCount());

    EXPECT_EQ(6ul, bdd4.SatCount(6));
    EXPECT_EQ(6ul, bdd4.NodeCount());

    EXPECT_EQ(bdd4, bdd5);
    EXPECT_EQ(bdd4.GetBDD(), bdd5.GetBDD());
}

}  // namespace