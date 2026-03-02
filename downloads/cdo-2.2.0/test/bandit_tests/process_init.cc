
#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <iostream>
#include <memory>

#include "../../src/modules.h"
#include "../../src/process.h"
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
  cdo::progname = "process_init_test";
  std::vector<std::string> arguements;
  auto shared_process = std::make_shared<Process>(0, alias_name, arguements);

  bandit::it(
      "assignes the origial of the alias as operatorName in init of process",
      [&]() {
        AssertThat(shared_process->operatorName, Is().EqualTo("oper1-1"));
      });

  bandit::it("creates the right prompt, discarding the alias in favor of the original name", [&]() {
    AssertThat(shared_process->prompt,
               Is().EqualTo(std::string(cdo::progname) + "    " + shared_process->operatorName));
  });
});

int
main(int argc, char **argv)
{

  int result = bandit::run(argc, argv);

  return result;
}
