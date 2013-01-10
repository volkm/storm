/*
 * Program.h
 *
 *  Created on: 04.01.2013
 *      Author: chris
 */

#ifndef PROGRAM_H_
#define PROGRAM_H_

#include "src/ir/expressions/ConstantExpression.h"
#include "Module.h"
#include "RewardModel.h"

namespace storm {

namespace ir {

class Program {
public:

	enum ModelType {UNDEFINED, DTMC, CTMC, MDP, CTMDP};

	Program() : modelType(UNDEFINED), booleanUndefinedConstantExpressions(), integerUndefinedConstantExpressions(), doubleUndefinedConstantExpressions(), modules(), rewards() {

	}

	Program(ModelType modelType, std::map<std::string, std::shared_ptr<storm::ir::expressions::BooleanConstantExpression>> booleanUndefinedConstantExpressions, std::map<std::string, std::shared_ptr<storm::ir::expressions::IntegerConstantExpression>> integerUndefinedConstantExpressions, std::map<std::string, std::shared_ptr<storm::ir::expressions::DoubleConstantExpression>> doubleUndefinedConstantExpressions, std::vector<storm::ir::Module> modules, std::map<std::string, storm::ir::RewardModel> rewards)
		: modelType(modelType), booleanUndefinedConstantExpressions(booleanUndefinedConstantExpressions), integerUndefinedConstantExpressions(integerUndefinedConstantExpressions), doubleUndefinedConstantExpressions(doubleUndefinedConstantExpressions), modules(modules), rewards(rewards) {

	}

	Program(ModelType modelType, std::map<std::string, std::shared_ptr<storm::ir::expressions::BooleanConstantExpression>> booleanUndefinedConstantExpressions, std::map<std::string, std::shared_ptr<storm::ir::expressions::IntegerConstantExpression>> integerUndefinedConstantExpressions, std::map<std::string, std::shared_ptr<storm::ir::expressions::DoubleConstantExpression>> doubleUndefinedConstantExpressions, std::vector<storm::ir::Module> modules)
		: modelType(modelType), booleanUndefinedConstantExpressions(booleanUndefinedConstantExpressions), integerUndefinedConstantExpressions(integerUndefinedConstantExpressions), doubleUndefinedConstantExpressions(doubleUndefinedConstantExpressions), modules(modules), rewards() {

	}


	Program(ModelType modelType, std::vector<storm::ir::Module> modules) : modelType(modelType), booleanUndefinedConstantExpressions(), integerUndefinedConstantExpressions(), doubleUndefinedConstantExpressions(), modules(modules), rewards() {

	}

	std::string toString() {
		std::string result = "";
		switch (modelType) {
		case UNDEFINED: result += "undefined\n\n"; break;
		case DTMC: result += "dtmc\n\n"; break;
		case CTMC: result += "ctmc\n\n"; break;
		case MDP: result += "mdp\n\n"; break;
		case CTMDP: result += "ctmdp\n\n"; break;
		}
		for (auto element : booleanUndefinedConstantExpressions) {
			result += "const bool " + element.first + ";\n";
		}
		for (auto element : integerUndefinedConstantExpressions) {
			result += "const int " + element.first + ";\n";
		}
		for (auto element : doubleUndefinedConstantExpressions) {
			result += "const double " + element.first + ";\n";
		}
		result += "\n";
		for (auto mod : modules) {
			result += mod.toString();
			result += "\n";
		}

		for (auto rewardModel : rewards) {
			result += rewardModel.second.toString();
			result +="\n";
		}

		return result;
	}

	void addBooleanUndefinedConstantExpression(std::string constantName, std::shared_ptr<storm::ir::expressions::BooleanConstantExpression> constantExpression) {
		booleanUndefinedConstantExpressions[constantName] = constantExpression;
	}

	void addIntegerUndefinedConstantExpression(std::string constantName, std::shared_ptr<storm::ir::expressions::IntegerConstantExpression> constantExpression) {
		integerUndefinedConstantExpressions[constantName] = constantExpression;
	}

	void addDoubleUndefinedConstantExpression(std::string constantName, std::shared_ptr<storm::ir::expressions::DoubleConstantExpression> constantExpression) {
		doubleUndefinedConstantExpressions[constantName] = constantExpression;
	}

private:
	ModelType modelType;
	std::map<std::string, std::shared_ptr<storm::ir::expressions::BooleanConstantExpression>> booleanUndefinedConstantExpressions;
	std::map<std::string, std::shared_ptr<storm::ir::expressions::IntegerConstantExpression>> integerUndefinedConstantExpressions;
	std::map<std::string, std::shared_ptr<storm::ir::expressions::DoubleConstantExpression>> doubleUndefinedConstantExpressions;
	std::vector<storm::ir::Module> modules;
	std::map<std::string, storm::ir::RewardModel> rewards;
};

}

}

#endif /* PROGRAM_H_ */
