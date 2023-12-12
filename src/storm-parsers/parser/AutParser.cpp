#include "storm-parsers/parser/AutParser.h"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <string>

#include "storm-parsers/parser/ValueParser.h"
#include "storm/exceptions/AbortException.h"
#include "storm/exceptions/WrongFormatException.h"
#include "storm/io/file.h"
#include "storm/storage/Distribution.h"
#include "storm/utility/SignalHandler.h"
#include "storm/utility/builder.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm {
namespace parser {

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::models::sparse::Model<ValueType, RewardModelType>> AutParser<ValueType, RewardModelType>::parseModel(std::string const& filename) {
    // Load file
    STORM_LOG_INFO("Reading from file " << filename);
    std::ifstream file;
    storm::utility::openFile(filename, file);
    std::string line;

    // Initialize
    std::shared_ptr<storm::storage::sparse::ModelComponents<ValueType, RewardModelType>> modelComponents;

    // Parse header
    while (storm::utility::getline(file, line)) {
        if (line.empty() || boost::starts_with(line, "//")) {
            continue;
        }

        if (boost::starts_with(line, "des (")) {
            // Parse header of form "des (<first_state>, <nr_of_transitions>, <nr_of_states>)"
            std::string substr = line.substr(5, line.length() - 5 - 1);
            std::vector<std::string> header;
            boost::split(header, substr, boost::is_any_of(","));
            STORM_LOG_THROW(header.size() == 3, storm::exceptions::WrongFormatException, "Invalid header definition");

            ValueParser<ValueType> valueParser;
            storm::storage::Distribution<ValueType> initDist = parseDistribution(header[0], valueParser);
            size_t nrChoices = parseNumber<size_t>(header[1]);
            size_t nrStates = parseNumber<size_t>(header[2]);

            // Construct model components
            modelComponents = parseStates(file, nrStates, nrChoices, initDist, valueParser);
            break;
        } else {
            STORM_LOG_THROW(false, storm::exceptions::WrongFormatException, "Could not parse line '" << line << "'.");
        }
    }
    // Done parsing
    storm::utility::closeFile(file);

    // Build model
    if (modelComponents->transitionMatrix.hasTrivialRowGrouping()) {
        return storm::utility::builder::buildModelFromComponents(storm::models::ModelType::Dtmc, std::move(*modelComponents));
    } else {
        return storm::utility::builder::buildModelFromComponents(storm::models::ModelType::Mdp, std::move(*modelComponents));
    }
}

template<typename ValueType, typename RewardModelType>
std::shared_ptr<storm::storage::sparse::ModelComponents<ValueType, RewardModelType>> AutParser<ValueType, RewardModelType>::parseStates(
    std::istream& file, size_t stateSize, size_t nrChoices, storm::storage::Distribution<ValueType> const& initDistribution,
    ValueParser<ValueType> const& valueParser) {
    // Offset if dedicated initial state is needed
    size_t initOffset = initDistribution.size() > 1 ? 1 : 0;
    // Initialize
    auto modelComponents = std::make_shared<storm::storage::sparse::ModelComponents<ValueType, RewardModelType>>();
    modelComponents->stateLabeling = storm::models::sparse::StateLabeling(stateSize + initOffset);
    modelComponents->choiceLabeling = storm::models::sparse::ChoiceLabeling(nrChoices + initOffset);
    // Temporary data structure to store transitions as states can have an arbitrary order in the file
    std::vector<std::vector<std::pair<std::string, storm::storage::Distribution<ValueType>>>> transitions(stateSize);

    // Iterate over all lines
    std::string line;
    uint64_t lineNumber = 0;
    bool deterministic = true;

    while (storm::utility::getline(file, line)) {
        lineNumber++;
        if (boost::starts_with(line, "//")) {
            continue;
        }
        if (line.empty()) {
            continue;
        }

        // Parse line
        STORM_LOG_TRACE("Parsing line no " << lineNumber << " : " << line);
        STORM_LOG_ASSERT(line.front() == '(' && line.back() == ')', "Invalid transition definition: " << line << ".");
        std::string substr = line.substr(1, line.length() - 2);

        // Parse state
        auto posEndState = substr.find(',');
        STORM_LOG_THROW(posEndState != std::string::npos, storm::exceptions::WrongFormatException, "Invalid transition definition: " << line << ".");
        size_t state = parseNumber<size_t>(substr.substr(0, posEndState));
        STORM_LOG_THROW(state <= stateSize, storm::exceptions::WrongFormatException, "More states detected than declared.");

        // Parse choice
        auto posBeginDistribution = substr.find_last_of(',');
        STORM_LOG_THROW(posBeginDistribution != std::string::npos, storm::exceptions::WrongFormatException, "Invalid transition definition: " << line << ".");
        STORM_LOG_THROW(substr[posEndState + 1] == '\"' && substr[posBeginDistribution - 1] == '\"', storm::exceptions::WrongFormatException,
                        "Invalid transition definition: " << line << ".");
        std::string choice = substr.substr(posEndState + 2, posBeginDistribution - posEndState - 3);

        // Parse distribution
        storm::storage::Distribution distribution = parseDistribution(substr.substr(posBeginDistribution + 1), valueParser);
        transitions[state].push_back(std::make_pair(choice, distribution));
        if (deterministic && transitions[state].size() > 1) {
            deterministic = false;
        }

        if (storm::utility::resources::isTerminate()) {
            std::cout << "Parsed " << lineNumber << " lines before abort.\n";
            STORM_LOG_THROW(false, storm::exceptions::AbortException, "Aborted in state space exploration.");
            break;
        }

    }  // end state iteration
    STORM_LOG_TRACE("Finished parsing");

    // Create actual data structures now that all transitions are known
    // Build sparse matrix
    storm::storage::SparseMatrixBuilder<ValueType> builder =
        storm::storage::SparseMatrixBuilder<ValueType>(nrChoices, stateSize, 0, false, !deterministic, stateSize);
    size_t row = 0;
    for (size_t state = 0; state < stateSize; ++state) {
        // Add states
        if (!deterministic) {
            builder.newRowGroup(row);
        }

        for (auto const& [choice, distribution] : transitions[state]) {
            // Add choices
            for (auto const& entry : distribution) {
                // Add transitions
                STORM_LOG_THROW(entry.first < stateSize, storm::exceptions::WrongFormatException,
                                "For state " << state << ", target state " << entry.first << " is greater than state size " << stateSize);
                STORM_LOG_TRACE("Transition " << state << " -> " << entry.first << ": " << entry.second);
                builder.addNextValue(row, entry.first, entry.second);
            }
            // Add choice label
            if (!modelComponents->choiceLabeling.value().containsLabel(choice)) {
                modelComponents->choiceLabeling.value().addLabel(choice);
            }
            modelComponents->choiceLabeling.value().addLabelToChoice(choice, row);
            ++row;
        }

        // Add custom state label identifying the state id
        std::string label = "s" + std::to_string(state);
        modelComponents->stateLabeling.addLabel(label);
        modelComponents->stateLabeling.addLabelToState(label, state);
    }

    STORM_LOG_THROW(
        nrChoices == 0 || builder.getLastRow() + 1 == nrChoices, storm::exceptions::WrongFormatException,
        "Number of transitions detected (at least " << builder.getLastRow() + 1 << ") does not match number of transitions declared (" << nrChoices << ").");

    // Add a dedicated initial state if needed
    modelComponents->stateLabeling.addLabel("init");
    if (initOffset > 0) {
        size_t initState = stateSize;
        if (!deterministic) {
            builder.newRowGroup(row);
        }
        for (auto const& entry : initDistribution) {
            // Add transitions
            STORM_LOG_THROW(entry.first < stateSize, storm::exceptions::WrongFormatException,
                            "For initial distribution, target state " << entry.first << " is greater than state size " << stateSize);
            STORM_LOG_TRACE("Initial transition " << initState << " -> " << entry.first << ": " << entry.second);
            builder.addNextValue(row, entry.first, entry.second);
        }
        modelComponents->stateLabeling.addLabelToState("init", initState);
    } else {
        size_t initState = initDistribution.begin()->first;
        modelComponents->stateLabeling.addLabelToState("init", initState);
    }

    // Build transition matrix
    modelComponents->transitionMatrix = builder.build(row + initOffset, stateSize + initOffset, stateSize + initOffset);
    STORM_LOG_TRACE("Built matrix");
    return modelComponents;
}

template<typename ValueType, typename RewardModelType>
storm::storage::Distribution<ValueType> AutParser<ValueType, RewardModelType>::parseDistribution(std::string const& distributionString,
                                                                                                 ValueParser<ValueType> const& valueParser) {
    std::vector<std::string> distValues;
    boost::split(distValues, distributionString, boost::is_any_of(" "));
    STORM_LOG_THROW(distValues.size() % 2 == 1, storm::exceptions::WrongFormatException, "Invalid distribution definition " << distributionString << ".");

    ValueType sum = storm::utility::zero<ValueType>();
    ValueType probability = storm::utility::zero<ValueType>();
    storm::storage::Distribution<ValueType> distribution;
    distribution.reserve(distValues.size() / 2 + 1);
    for (size_t i = 0; i + 1 < distValues.size(); i += 2) {
        probability = valueParser.parseValue(distValues[i + 1]);
        sum += probability;
        distribution.addProbability(parseNumber<size_t>(distValues[i]), probability);
    }
    ValueType remaining = storm::utility::one<ValueType>() - sum;
    STORM_LOG_ASSERT(!storm::utility::isZero(remaining), "Last state in " << distributionString << " should be reached with non-zero probability.");
    distribution.addProbability(parseNumber<size_t>(distValues.back()), remaining);
    return distribution;
}

// Template instantiations.
template class AutParser<double>;
template class AutParser<storm::RationalNumber>;

}  // namespace parser
}  // namespace storm
