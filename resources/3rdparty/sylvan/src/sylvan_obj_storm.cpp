#include "storm_wrapper.h"
#include "sylvan_storm_rational_number.h"
#include "sylvan_storm_rational_function.h"

/*********************************************
 Functions added to sylvan's Sylvan class.
 *********************************************/
void
Sylvan::initCustomMtbdd()
{
    sylvan_init_mt();
    sylvan_storm_rational_number_init();
    sylvan_storm_rational_function_init();
}

/*********************************************
 Functions added to sylvan's Bdd class.
 *********************************************/

Mtbdd
Bdd::toDoubleMtbdd() const {
    LACE_ME;
    return mtbdd_bool_to_double(bdd);
}

Mtbdd
Bdd::toInt64Mtbdd() const {
    LACE_ME;
    return mtbdd_bool_to_int64(bdd);
}

Mtbdd
Bdd::toStormRationalNumberMtbdd() const {
    LACE_ME;
    return mtbdd_bool_to_storm_rational_number(bdd);
}

#if defined(SYLVAN_HAVE_CARL) || defined(STORM_HAVE_CARL)
Mtbdd
Bdd::toStormRationalFunctionMtbdd() const {
    LACE_ME;
    return mtbdd_bool_to_storm_rational_function(bdd);
}
#endif

Mtbdd
Bdd::Ite(Mtbdd const& thenDd, Mtbdd const& elseDd) const {
    LACE_ME;
    return mtbdd_ite(bdd, thenDd.GetMTBDD(), elseDd.GetMTBDD());
}

Bdd
Bdd::ExistAbstractRepresentative(const BddSet& cube) const {
    LACE_ME;
    return sylvan_existsRepresentative(bdd, cube.set.bdd);
}

/*********************************************
 Functions added to sylvan's Mtbdd class.
 *********************************************/

Bdd
Mtbdd::NotZero() const
{
    LACE_ME;
    return mtbdd_not_zero(mtbdd);
}

size_t
Mtbdd::CountLeaves() const {
    LACE_ME;
    return mtbdd_leafcount(mtbdd);
}

double
Mtbdd::NonZeroCount(size_t variableCount) const {
    LACE_ME;
    return mtbdd_non_zero_count(mtbdd, variableCount);
}

bool
Mtbdd::isValid() const {
    LACE_ME;
    return mtbdd_test_isvalid(mtbdd) == 1;
}

void
Mtbdd::PrintDot(FILE *out) const {
    mtbdd_fprintdot(out, mtbdd);
}

std::string
Mtbdd::GetShaHash() const {
    char buf[65];
    mtbdd_getsha(mtbdd, buf);
    return std::string(buf);
}

Mtbdd
Mtbdd::Minus(const Mtbdd &other) const
{
    LACE_ME;
    return mtbdd_minus(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::Divide(const Mtbdd &other) const
{
    LACE_ME;
    return mtbdd_divide(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::Equals(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_equals(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::Less(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_less_as_bdd(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::LessOrEqual(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_less_or_equal_as_bdd(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::AbstractMinRepresentative(const BddSet &variables) const
{
    LACE_ME;
    return mtbdd_min_abstract_representative(mtbdd, variables.set.bdd);
}

Bdd
Mtbdd::AbstractMaxRepresentative(const BddSet &variables) const
{
    LACE_ME;
    return mtbdd_max_abstract_representative(mtbdd, variables.set.bdd);
}

Mtbdd
Mtbdd::Pow(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_pow(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::Mod(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_mod(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::Logxy(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_logxy(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::Floor() const {
    LACE_ME;
    return mtbdd_floor(mtbdd);
}

Mtbdd
Mtbdd::Ceil() const {
    LACE_ME;
    return mtbdd_ceil(mtbdd);
}

Mtbdd
Mtbdd::Minimum() const {
    LACE_ME;
    return mtbdd_minimum(mtbdd);
}

Mtbdd
Mtbdd::Maximum() const {
    LACE_ME;
    return mtbdd_maximum(mtbdd);
}

bool
Mtbdd::EqualNorm(const Mtbdd& other, double epsilon) const {
    LACE_ME;
    return mtbdd_equal_norm_d(mtbdd, other.mtbdd, epsilon);
}

bool
Mtbdd::EqualNormRel(const Mtbdd& other, double epsilon) const {
    LACE_ME;
    return mtbdd_equal_norm_rel_d(mtbdd, other.mtbdd, epsilon);
}

// Functions for Mtbdds over rational numbers.
Mtbdd
Mtbdd::stormRationalNumberTerminal(storm::RationalNumber const& value)
{
    storm_rational_number_ptr ptr = (storm_rational_number_ptr)(&value);
    return mtbdd_storm_rational_number(ptr);
}

Bdd
Mtbdd::EqualsRN(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_equals_rational_number(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::LessRN(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_number_less(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::LessOrEqualRN(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_number_less_or_equal(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinRN(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_number_min(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MaxRN(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_number_max(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::PlusRN(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_number_plus(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinusRN(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_number_minus(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::TimesRN(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_number_times(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::DivideRN(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_number_divide(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::FloorRN() const {
    LACE_ME;
    return sylvan_storm_rational_number_floor(mtbdd);
}

Mtbdd
Mtbdd::CeilRN() const {
    LACE_ME;
    return sylvan_storm_rational_number_ceil(mtbdd);
}

Mtbdd
Mtbdd::PowRN(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_number_pow(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinimumRN() const {
    LACE_ME;
    return sylvan_storm_rational_number_minimum(mtbdd);
}

Mtbdd
Mtbdd::MaximumRN() const {
    LACE_ME;
    return sylvan_storm_rational_number_maximum(mtbdd);
}

Mtbdd
Mtbdd::AndExistsRN(const Mtbdd &other, const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_number_and_exists(mtbdd, other.mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractPlusRN(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_number_abstract_plus(mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractMinRN(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_number_abstract_min(mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractMaxRN(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_number_abstract_max(mtbdd, variables.set.bdd);
}

Bdd
Mtbdd::BddThresholdRN(storm::RationalNumber const& rn) const {
    LACE_ME;
    return sylvan_storm_rational_number_threshold(mtbdd, (storm_rational_number_ptr)&rn);
}

Bdd
Mtbdd::BddStrictThresholdRN(storm::RationalNumber const& rn) const {
    LACE_ME;
    return sylvan_storm_rational_number_strict_threshold(mtbdd, (storm_rational_number_ptr)&rn);
}

bool
Mtbdd::EqualNormRN(const Mtbdd& other, storm::RationalNumber const& epsilon) const {
    LACE_ME;
    return sylvan_storm_rational_number_equal_norm_d(mtbdd, other.mtbdd, (storm_rational_number_ptr)&epsilon);
}

bool
Mtbdd::EqualNormRelRN(const Mtbdd& other, storm::RationalNumber const& epsilon) const {
    LACE_ME;
    return sylvan_storm_rational_number_equal_norm_rel_d(mtbdd, other.mtbdd, (storm_rational_number_ptr)&epsilon);
}

Mtbdd
Mtbdd::ToDoubleRN() const {
    LACE_ME;
    return sylvan_storm_rational_number_to_double(mtbdd);
}

// Functions for Mtbdds over rational functions.
#if defined(SYLVAN_HAVE_CARL) || defined(STORM_HAVE_CARL)
Mtbdd
Mtbdd::stormRationalFunctionTerminal(storm::RationalFunction const& value)
{
    storm_rational_function_ptr ptr = (storm_rational_function_ptr)(&value);
    return mtbdd_storm_rational_function(ptr);
}

Bdd
Mtbdd::EqualsRF(const Mtbdd& other) const {
    LACE_ME;
    return mtbdd_equals_rational_function(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::LessRF(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_function_less(mtbdd, other.mtbdd);
}

Bdd
Mtbdd::LessOrEqualRF(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_function_less_or_equal(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinRF(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_function_min(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MaxRF(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_function_max(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::PlusRF(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_function_plus(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinusRF(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_function_minus(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::TimesRF(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_function_times(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::DivideRF(const Mtbdd &other) const
{
    LACE_ME;
    return sylvan_storm_rational_function_divide(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::FloorRF() const {
    LACE_ME;
    return sylvan_storm_rational_function_floor(mtbdd);
}

Mtbdd
Mtbdd::CeilRF() const {
    LACE_ME;
    return sylvan_storm_rational_function_ceil(mtbdd);
}

Mtbdd
Mtbdd::PowRF(const Mtbdd& other) const {
    LACE_ME;
    return sylvan_storm_rational_function_pow(mtbdd, other.mtbdd);
}

Mtbdd
Mtbdd::MinimumRF() const {
    LACE_ME;
    return sylvan_storm_rational_function_minimum(mtbdd);
}

Mtbdd
Mtbdd::MaximumRF() const {
    LACE_ME;
    return sylvan_storm_rational_function_maximum(mtbdd);
}

Mtbdd
Mtbdd::AndExistsRF(const Mtbdd &other, const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_function_and_exists(mtbdd, other.mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractPlusRF(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_function_abstract_plus(mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractMinRF(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_function_abstract_min(mtbdd, variables.set.bdd);
}

Mtbdd Mtbdd::AbstractMaxRF(const BddSet &variables) const {
    LACE_ME;
    return sylvan_storm_rational_function_abstract_max(mtbdd, variables.set.bdd);
}

Bdd
Mtbdd::BddThresholdRF(storm::RationalFunction const& rf) const {
    LACE_ME;
    return sylvan_storm_rational_function_threshold(mtbdd, (storm_rational_function_ptr)&rf);
}

Bdd
Mtbdd::BddStrictThresholdRF(storm::RationalFunction const& rf) const {
    LACE_ME;
    return sylvan_storm_rational_function_strict_threshold(mtbdd, (storm_rational_function_ptr)&rf);
}

bool
Mtbdd::EqualNormRF(const Mtbdd& other, storm::RationalFunction const& epsilon) const {
    LACE_ME;
    return sylvan_storm_rational_function_equal_norm_d(mtbdd, other.mtbdd, (storm_rational_number_ptr)&epsilon);
}

bool
Mtbdd::EqualNormRelRF(const Mtbdd& other, storm::RationalFunction const& epsilon) const {
    LACE_ME;
    return sylvan_storm_rational_function_equal_norm_rel_d(mtbdd, other.mtbdd, (storm_rational_number_ptr)&epsilon);
}

Mtbdd
Mtbdd::ToDoubleRF() const {
    LACE_ME;
    return sylvan_storm_rational_function_to_double(mtbdd);
}
#endif











