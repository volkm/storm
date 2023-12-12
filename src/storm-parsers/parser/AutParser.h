#pragma once

#include "storm/models/sparse/Model.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/Distribution.h"
#include "storm/storage/sparse/ModelComponents.h"

namespace storm {
namespace parser {

template<typename T>
class ValueParser;

/*!
 *	Parser for models in the AUT format (also called Aldebaran format).
 *	@see https://cadp.inria.fr/man/aut.html for the original description supported by CADP.
 *	@see https://www.mcrl2.org/web/user_manual/tools/lts.html for the probabilistic extension supported by mCRL2.
 */
template<typename ValueType, typename RewardModelType = models::sparse::StandardRewardModel<ValueType>>
class AutParser {
   public:
    /*!
     * Load a model in AUT format from a file and create the model.
     *
     * @param filename The AUT file to be parsed.
     *
     * @return A sparse model
     */
    static std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> parseModel(std::string const& filename);

   private:
    /*!
     * Parse states and return transition matrix.
     *
     * @param file Input file stream.
     * @param stateSize No. of states
     * @param nrChoices No. of choices (corresponds to no. of transtions).
     * @param initDistribution Distribution of initial state.
     * @param valueParser Value parser.
     *
     * @return Transition matrix.
     */
    static std::shared_ptr<storm::storage::sparse::ModelComponents<ValueType, RewardModelType>> parseStates(
        std::istream& file, size_t stateSize, size_t nrChoices, storm::storage::Distribution<ValueType> const& initDistribution,
        ValueParser<ValueType> const& valueParser);

    /*!
     * Parse a distribution of the form (s_0 p_0 s_1 ... p_{n-1} s_n).
     *
     * @param distributionString String representation of distribution.
     * @param valueParser Value parser.
     * @return Distribution.
     */
    static storm::storage::Distribution<ValueType> parseDistribution(std::string const& distributionString, ValueParser<ValueType> const& valueParser);
};

}  // namespace parser
}  // namespace storm
