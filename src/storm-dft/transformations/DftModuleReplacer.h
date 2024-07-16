#pragma once

#include "storm-dft/storage/DFT.h"
#include "storm-dft/storage/DftModule.h"

namespace storm::dft {
namespace transformations {

/*!
 * Replaces DFT modules.
 */
template<typename ValueType>
class DftModuleReplacer {
   public:
    /*!
     * Replace independent modules in DFT by new (sub-)DFTs.
     * @param dft DFT.
     * @param moduleReplacements List of the module replacements given as pairs <oldModule, new DFT>.
     *                           The representative of oldModule and the TLE of the new DFT must have the same name.
     * @return DFT where each old module has been replaced with its new DFT.
     */
    static std::shared_ptr<storm::dft::storage::DFT<ValueType>> replaceModules(
        storm::dft::storage::DFT<ValueType> const& dft,
        std::vector<std::pair<storm::dft::storage::DftIndependentModule const&, std::shared_ptr<storm::dft::storage::DFT<ValueType> const>>> const&
            moduleReplacements);
};

}  // namespace transformations
}  // namespace storm::dft
