#pragma once

#include "storm-dft/storage/DFT.h"

namespace storm {

namespace jani {
class Property;
class Model;
}  // namespace jani
namespace gspn {
class GSPN;
}  // namespace gspn

namespace dft {
namespace api {

/*!
 * Check whether the DFT is well-formed.
 *
 * @param dft DFT.
 * @param validForMarkovianAnalysis If true, additional checks are performed to check whether the DFT is valid for analysis via Markov models.
 * @return Pair where the first entry is true iff the DFT is well-formed. The second entry contains the error messages for ill-formed parts.
 */
template<typename ValueType>
std::pair<bool, std::string> isWellFormed(storm::dft::storage::DFT<ValueType> const& dft, bool validForMarkovianAnalysis = true);

/*!
 * Apply transformations for DFT.
 *
 * @param dft DFT.
 * @param uniqueBE Flag whether a unique constant failed BE is created.
 * @param binaryFDEP Flag whether all dependencies should be binary (only one dependent child).
 * @param exponentialDistributions Flag whether distributions should be transformed to exponential distributions (if possible).
 * @return Transformed DFT.
 */
template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> applyTransformations(storm::dft::storage::DFT<ValueType> const& dft, bool uniqueBE, bool binaryFDEP,
                                                                          bool exponentialDistributions);

/*!
 * Apply transformations to make DFT feasible for analysis via Markov models.
 *
 * @param dft DFT.
 * @return Transformed DFT.
 */
template<typename ValueType>
std::shared_ptr<storm::dft::storage::DFT<ValueType>> prepareForMarkovAnalysis(storm::dft::storage::DFT<ValueType> const& dft) {
    return storm::dft::api::applyTransformations(dft, true, true, true);
}

/*!
 * Transform DFT to GSPN.
 *
 * @param dft DFT.
 * @param disableDC Flag indicating if Don't Care propagation for DFT elements which should be disabled.
 * @param mergeDCFailed Flag indicating if Don't Care places and Failed places should be merged.
 * @param extendPriorities Flag indicating if the extended priority calculation is used.
 * @param smartTransformation Flag indicating if smart semantics should be used.
 *                            Smart semantics will only generate necessary parts of the GSPNs.
 * @return Pair of GSPN and id of failed place corresponding to the top level element.
 */
template<typename ValueType>
std::pair<std::shared_ptr<storm::gspn::GSPN>, uint64_t> transformToGSPN(storm::dft::storage::DFT<ValueType> const& dft, bool disableDC, bool mergeDCFailed,
                                                                        bool extendPriorities, bool smartTransformation);

/*!
 * Transform GSPN to Jani model.
 *
 * @param gspn GSPN.
 * @param toplevelFailedPlace Id of the failed place in the GSPN for the top level element in the DFT.
 * @return Pair of JANI model and Jani properties.
 */
std::pair<std::shared_ptr<storm::jani::Model>, std::vector<storm::jani::Property>> transformToJani(storm::gspn::GSPN const& gspn, uint64_t toplevelFailedPlace);

}  // namespace api
}  // namespace dft
}  // namespace storm