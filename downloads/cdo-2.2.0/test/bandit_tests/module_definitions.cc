#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <iostream>
#include "../../src/modules.h"

bool is_alias(const std::string &);

const std::string alias_name = "alias1-1";
void *
testFunction(void *test)
{
  return test;
}
static const module_t module_with_alias = { "aliasModule",
                                            testFunction,
                                            {},
                                            { "oper1-1", "oper1-2" },
                                            1,
                                            0,
                                            1,
                                            1,
                                            NoRestriction,
                                            { Alias(alias_name, "oper1-1") } };

using namespace snowhouse;

go_bandit([]() {
  bandit::describe("Testing for registered module with alias", [&]() {
    bandit::it("has registered the module",
               [&]() { AssertThat(get_modules(), Is().OfLength(1)); });
    bandit::it(
        "has added both operators to module map while not adding the alias",
        [&]() { AssertThat(get_module_map(), Is().OfLength(2)); });
    bandit::it("has registered the alias",
               [&]() { AssertThat(get_aliases(), Is().OfLength(1)); });
    bandit::it("has registered the correct original name", [&]() {
      AssertThat(get_aliases()[alias_name], Is().EqualTo("oper1-1"));
    });
  });
  bandit::describe("Testing the interface for alias and name retreval", [&]() {
    bandit::it("gets the right name for alias", [&]() {
      AssertThat(get_original(alias_name.c_str()), Is().EqualTo("oper1-1"));
    });
    bandit::it("extracts the operator name", [&]() {
      AssertThat(get_original("oper1-1"), Is().EqualTo("oper1-1"));
    });
    bandit::it("recognizes an alias",
               [&]() { AssertThat(is_alias(alias_name), IsTrue()); });
  });
});

int
main(int argc, char **argv)
{

  int result = bandit::run(argc, argv);

  return result;
}
