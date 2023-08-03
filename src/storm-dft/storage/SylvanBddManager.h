#pragma once

#include <map>

#include "storm/adapters/sylvan.h"
#include "storm/storage/dd/sylvan/InternalSylvanDdManager.h"
#include "storm/utility/macros.h"

namespace storm::dft {
namespace storage {

/*!
 * Simple Manager for Sylvan.
 * Note that the manager is not thread safe.
 * Initializes Sylvan which can only be initialized once per program invocation.
 */
class SylvanBddManager {
   public:
    /*!
     * Initialize Sylvan.
     * Internally Sylvan is initialized by a InternalSylvanDdManager. This ensures compatibility.
     */
    SylvanBddManager() = default;

    // We can only initialize Sylvan once therefore no copy semantics
    SylvanBddManager(SylvanBddManager const&) = delete;

    SylvanBddManager(SylvanBddManager&&) = default;
    SylvanBddManager& operator=(SylvanBddManager&&) = default;

    /*!
     * Destructor.
     */
    ~SylvanBddManager() = default;

    /*!
     * All code that manipulates DDs shall be called through this function.
     * This is generally needed to set-up the correct context.
     * Specifically for sylvan, this is required to make sure that DD-manipulating code is executed as a LACE task.
     * Example usage: `manager->execute([&]() { bar = foo(arg1,arg2); }`
     *
     * @param f the function that is executed
     */
    void execute(std::function<void()> const& f) const {
        internalManager.execute(f);
    }

    /**
     * Create a variable with the given name.
     *
     * @param name Variable name.
     * @return Index of the variable.
     */
    uint32_t createVariable(std::string const& name) {
        nameToIndex[name] = nextFreeVariableIndex;
        indexToName[nextFreeVariableIndex] = name;
        sylvan::Bdd::bddVar(nextFreeVariableIndex);
        return nextFreeVariableIndex++;
    }

    /*!
     * Get positive literal for variable.
     * @param name Variable name.
     * @return BDD representation of the variable: ITE(x, 1, 0).
     */
    sylvan::Bdd getPositiveLiteral(std::string const& name) const {
        return getPositiveLiteral(getIndex(name));
    }

    /*!
     * Get negative literal for variable.
     * @param name Variable name.
     * @return Negative BDD representation of the variable: ITE(x, 0, 1).
     */
    sylvan::Bdd getNegativeLiteral(std::string const& name) const {
        return getNegativeLiteral(getIndex(name));
    }

    /*!
     * Get positive literal for variable.
     * @param index Variable index.
     * @return BDD representation of the variable: ITE(x, 1, 0).
     */
    sylvan::Bdd getPositiveLiteral(uint32_t const index) const {
        return sylvan::Bdd::bddVar(index);
    }

    /*!
     * Get negative literal for variable.
     * @param index Variable index.
     * @return Negative BDD representation of the variable: ITE(x, 0, 1).
     */
    sylvan::Bdd getNegativeLiteral(uint32_t const index) const {
        return !sylvan::Bdd::bddVar(index);
    }

    /*!
     * Get BDD representation of the constant function 1.
     * @return BDD representing 1.
     */
    sylvan::Bdd getOne() const {
        return sylvan::Bdd::bddOne();
    }

    /*!
     * Get BDD representation of the constant function 0.
     * @return BDD representing 0.
     */
    sylvan::Bdd getZero() const {
        return sylvan::Bdd::bddZero();
    }

    /*!
     * Get index of variable with given name.
     * @param name Variable name.
     * @return Index.
     */
    uint32_t getIndex(std::string const& name) const {
        return nameToIndex.at(name);
    }

    /*!
     * Get name of variable at given index.
     * @param index Variable index.
     * @return Name.
     */
    std::string getName(uint32_t const index) const {
        return indexToName.at(index);
    }

   private:
    storm::dd::InternalDdManager<storm::dd::DdType::Sylvan> internalManager{};
    uint32_t nextFreeVariableIndex{0};

    std::map<std::string, uint32_t> nameToIndex{};
    std::map<uint32_t, std::string> indexToName{};
};

}  // namespace storage
}  // namespace storm::dft
