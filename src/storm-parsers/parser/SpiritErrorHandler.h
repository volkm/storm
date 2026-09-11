#pragma once

#include <memory>
#include <string>
#include <utility>

#include "storm-parsers/parser/SpiritParserDefinitions.h"
#include "storm/exceptions/WrongFormatException.h"
#include "storm/utility/macros.h"

namespace storm {
namespace parser {

// Collects a deferred parsing error message.
//
// Semantic actions and error handlers invoked by Boost.Spirit must not throw: on some toolchains
// (notably macOS x86_64) a C++ exception cannot be unwound out through Spirit's heavily-inlined
// Qi/Phoenix template frames and hits std::terminate instead of the caller's catch. Instead they
// record the message here and signal failure through Spirit's own channels (qi::_pass = false, or
// an on_error handler returning qi::fail). The driver raises the WrongFormatException from normal
// control flow once phrase_parse has returned.
//
// Sink to capture error messages. Recording is last-wins.
struct SpiritErrorSink {
    bool hasError = false;
    std::string message;

    void record(std::string msg) {
        hasError = true;
        message = std::move(msg);
    }

    void reset() {
        hasError = false;
        message.clear();
    }
};

// Functor used for displaying error information.
struct SpiritErrorHandler {
    typedef qi::error_handler_result result_type;

    // Sink for recording errors. If null, it falls back to logging.
    std::shared_ptr<SpiritErrorSink> sink;

    template<typename T1, typename T2, typename T3, typename T4>
    qi::error_handler_result operator()(T1 b, T2 e, T3 where, T4 const& what) const {
        auto lineStart = boost::spirit::get_line_start(b, where);
        auto lineEnd = std::find(where, e, '\n');
        std::string line(lineStart, lineEnd);

        std::stringstream stream;
        stream << "Parsing error at " << get_line(where) << ":" << boost::spirit::get_column(lineStart, where) << ": "
               << " expecting " << what << ", here:\n";
        stream << "\t" << line << '\n';
        auto caretColumn = boost::spirit::get_column(lineStart, where);
        stream << "\t" << std::string(caretColumn - 1, ' ') << "^\n";

        if (sink) {
            sink->record(stream.str());
        } else {
            STORM_LOG_ERROR(stream.str());
        }
        return qi::fail;
    }
};
}  // namespace parser
}  // namespace storm
