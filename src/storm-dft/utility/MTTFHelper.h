#pragma once

#include <memory>

namespace storm::dft {

namespace storage {
// Forward declaration
template<typename ValueType>
class DFT;
}  // namespace storage

namespace utility {

/*!
 * Algorithms to approximate the MTTF.
 */
enum class MTTFApproximationAlgorithm { Proceeding, VariableChange };

/*!
 * Numerically approximate the MTTF of the given DFT
 * by integrating 1 - cdf(dft) with Simpson's rule.
 * @param dft DFT.
 * @param stepsize Step size.
 * @param precision Precision for convergence.
 * @return Approximate MTTF.
 */
double MTTFHelperProceeding(std::shared_ptr<storm::dft::storage::DFT<double> const> const dft, double const stepsize = 1e-10, double const precision = 1e-12);

/*!
 * Numerically approximate the MTTF of the given DFT
 * by integrating 1 - cdf(dft) by changing the variable
 * such that the interval (0,oo) becomes (0,1).
 * @param dft DFT.
 * @param stepsize Step size.
 * @return Approximate MTTF.
 */
double MTTFHelperVariableChange(std::shared_ptr<storm::dft::storage::DFT<double> const> const dft, double const stepsize = 1e-6);

}  // namespace utility
}  // namespace storm::dft
