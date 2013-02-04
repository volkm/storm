#ifndef STORM_PARSER_NONDETTRAPARSER_H_
#define STORM_PARSER_NONDETTRAPARSER_H_

#include "src/storage/SparseMatrix.h"

#include "src/parser/Parser.h"
#include "src/utility/OsDetection.h"

#include <utility>
#include <memory>
#include <vector>

namespace storm {
namespace parser {
	
typedef boost::bimap<uint_fast64_t, std::pair<uint_fast64_t,std::string>> RowStateMapping;

/*!
 *	@brief	Load a nondeterministic transition system from file and create a
 *	sparse adjacency matrix whose entries represent the weights of the edges
 */
class NonDeterministicSparseTransitionParser : public Parser {
	public:
		NonDeterministicSparseTransitionParser(std::string const &filename);
		
		inline std::shared_ptr<storm::storage::SparseMatrix<double>> getMatrix() const {
			return this->matrix;
		}
		
		inline std::shared_ptr<std::vector<uint_fast64_t>> getRowMapping() const {
			return this->rowMapping;
		}
	
	private:
		std::shared_ptr<storm::storage::SparseMatrix<double>> matrix;
		std::shared_ptr<std::vector<uint_fast64_t>> rowMapping;
		
		uint_fast64_t firstPass(char* buf, uint_fast64_t& choices, int_fast64_t& maxnode);
	
};
		
} // namespace parser
} // namespace storm

#endif /* STORM_PARSER_NONDETTRAPARSER_H_ */
