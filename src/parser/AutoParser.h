#ifndef STORM_PARSER_AUTOPARSER_H_
#define STORM_PARSER_AUTOPARSER_H_

#include "src/models/AbstractModel.h"

#include <string>

namespace storm {

	/*!
	 * Contains all file parsers and helper classes.
	 *
	 * This namespace contains everything needed to load data files (like
	 * atomic propositions, transition systems, formulas, etc.) including
	 * methods for efficient file access (see MappedFile).
	 */
	namespace parser {

		/*!
		 * This class automatically chooses the correct parser for the given files and returns the corresponding model.
		 * The choice of the parser is made using the model hint at the beginning of the given transition file.
		 */
		class AutoParser {
		public:

			/*!
			 * Checks the given files and parses the model within these files.
			 *
			 * This parser analyzes the format hint in the first line of the transition file.
			 * If this is a valid format, it will use the parser for this format, otherwise it will throw an exception.
			 *
			 * When the files are parsed successfully, a shared pointer owning the resulting model is returned.
			 * The concrete model can be obtained using the as<Type>() member of the AbstractModel class.
			 *
			 * @note The number of states of the model is determined by the transitions file.
			 *       The labeling file may therefore not contain labels of states that are not contained in the transitions file.
			 *
			 * @param transitionsFilename The path and name of the file containing the transitions of the model.
			 * @param labelingFilename The path and name of the file containing the labels for the states of the model.
			 * @param stateRewardFilename The path and name of the file that contains the state reward of the model. This file is optional.
			 * @param transitionRewardFilename The path and name of the file that contains the transition rewards of the model. This file is optional.
			 * @return A shared_ptr containing the resulting model.
			 */
			static std::shared_ptr<storm::models::AbstractModel<double>> parseModel(std::string const & transitionsFilename,
                                                                                    std::string const & labelingFilename,
                                                                                    std::string const & stateRewardFilename = "",
                                                                                    std::string const & transitionRewardFilename = "");

		private:
            // Define the maximal length of a hint in the file.
            static constexpr uint_fast64_t hintLength = 10;

			/*!
			 *	Opens the given file and parses the file format hint.
			 *
			 *	@param filename The path and name of the file that is to be analysed.
			 *	@return The type of the model as an enum value.
			 */
			static storm::models::ModelType analyzeHint(std::string const & filename);
		};
		
	} // namespace parser
} // namespace storm

#endif /* STORM_PARSER_AUTOPARSER_H_ */
