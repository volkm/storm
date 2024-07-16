#include "transformation.h"

#include "storm-dft/transformer/DftToGspnTransformer.h"
#include "storm-dft/transformer/DftTransformer.h"
#include "storm-dft/utility/DftValidator.h"
#include "storm-gspn/builder/JaniGSPNBuilder.h"
#include "storm/storage/jani/Model.h"

namespace storm::dft {
namespace api {

template<typename ValueType>
std::pair<bool, std::string> isWellFormed(storm::dft::storage::DFT<ValueType> const& dft, bool validForMarkovianAnalysis) {
    std::stringstream stream;
    bool wellFormed = false;
    if (validForMarkovianAnalysis) {
        wellFormed = storm::dft::utility::DftValidator<ValueType>::isDftValidForMarkovianAnalysis(dft, stream);
    } else {
        wellFormed = storm::dft::utility::DftValidator<ValueType>::isDftWellFormed(dft, stream);
    }
    return std::pair<bool, std::string>(wellFormed, stream.str());
}

template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> applyTransformations(storm::dft::storage::DFT<ValueType> const& dft, bool uniqueBE, bool binaryFDEP,
                                                                          bool exponentialDistributions) {
    std::shared_ptr<storm::dft::storage::DFT<ValueType>> transformedDft = std::make_shared<storm::dft::storage::DFT<ValueType>>(dft);
    if (exponentialDistributions && !storm::dft::transformer::DftTransformer<ValueType>::hasOnlyExponentialDistributions(*transformedDft)) {
        transformedDft = storm::dft::transformer::DftTransformer<ValueType>::transformExponentialDistributions(*transformedDft);
    }
    if (uniqueBE && !storm::dft::transformer::DftTransformer<ValueType>::hasUniqueFailedBE(*transformedDft)) {
        transformedDft = storm::dft::transformer::DftTransformer<ValueType>::transformUniqueFailedBE(*transformedDft);
    }
    if (binaryFDEP && storm::dft::transformer::DftTransformer<ValueType>::hasNonBinaryDependency(*transformedDft)) {
        transformedDft = storm::dft::transformer::DftTransformer<ValueType>::transformBinaryDependencies(*transformedDft);
    }
    STORM_LOG_ASSERT(!exponentialDistributions || storm::dft::transformer::DftTransformer<ValueType>::hasOnlyExponentialDistributions(*transformedDft),
                     "DFT still has non-exponential distributions.");
    STORM_LOG_ASSERT(!uniqueBE || storm::dft::transformer::DftTransformer<ValueType>::hasUniqueFailedBE(*transformedDft), "DFT still has multiple failed BEs.");
    STORM_LOG_ASSERT(!binaryFDEP || !storm::dft::transformer::DftTransformer<ValueType>::hasNonBinaryDependency(*transformedDft),
                     "DFT still has non-binary dependencies.");
    return transformedDft;
}

template<>
std::pair<std::shared_ptr<storm::gspn::GSPN>, uint64_t> transformToGSPN(storm::dft::storage::DFT<double> const& dft, bool disableDC, bool mergeDCFailed,
                                                                        bool extendPriorities, bool smartTransformation) {
    // Set Don't Care elements
    std::set<uint64_t> dontCareElements;
    if (!disableDC) {
        // Insert all elements as Don't Care elements
        for (std::size_t i = 0; i < dft.nrElements(); i++) {
            dontCareElements.insert(dft.getElement(i)->id());
        }
    }

    // Transform to GSPN
    storm::dft::transformer::DftToGspnTransformer<double> gspnTransformator(dft);
    auto priorities = gspnTransformator.computePriorities(extendPriorities);
    gspnTransformator.transform(priorities, dontCareElements, smartTransformation, mergeDCFailed, extendPriorities);
    std::shared_ptr<storm::gspn::GSPN> gspn(gspnTransformator.obtainGSPN());
    return std::make_pair(gspn, gspnTransformator.toplevelFailedPlaceId());
}

std::pair<std::shared_ptr<storm::jani::Model>, std::vector<storm::jani::Property>> transformToJani(storm::gspn::GSPN const& gspn,
                                                                                                   uint64_t toplevelFailedPlace) {
    // Build Jani model
    storm::builder::JaniGSPNBuilder builder(gspn);
    std::shared_ptr<storm::jani::Model> model(builder.build("dft_gspn"));

    // Build properties
    std::shared_ptr<storm::expressions::ExpressionManager> const& exprManager = gspn.getExpressionManager();
    storm::jani::Variable const& topfailedVar = builder.getPlaceVariable(toplevelFailedPlace);
    storm::expressions::Expression targetExpression = exprManager->integer(1) == topfailedVar.getExpressionVariable().getExpression();
    // Add variable for easier access to 'failed' state
    builder.addTransientVariable(model.get(), "failed", targetExpression);
    auto failedFormula = std::make_shared<storm::logic::AtomicExpressionFormula>(targetExpression);
    std::vector<storm::jani::Property> properties = builder.getStandardProperties(model.get(), failedFormula, "Failed", "a failed state", true);
    return std::make_pair(model, properties);
}

template<>
std::pair<std::shared_ptr<storm::gspn::GSPN>, uint64_t> transformToGSPN(storm::dft::storage::DFT<storm::RationalFunction> const&, bool, bool, bool, bool) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Transformation to GSPN not supported for this data type.");
}

// Explicitly instantiate methods
// Double
template std::pair<bool, std::string> isWellFormed(storm::dft::storage::DFT<double> const&, bool);
template std::shared_ptr<storm::dft::storage::DFT<double>> applyTransformations(storm::dft::storage::DFT<double> const&, bool, bool binaryFDEP, bool);

// Rational function
template std::pair<bool, std::string> isWellFormed(storm::dft::storage::DFT<storm::RationalFunction> const&, bool);
template std::shared_ptr<storm::dft::storage::DFT<storm::RationalFunction>> applyTransformations(storm::dft::storage::DFT<storm::RationalFunction> const&, bool,
                                                                                                 bool binaryFDEP, bool);
}  // namespace api
}  // namespace storm::dft
