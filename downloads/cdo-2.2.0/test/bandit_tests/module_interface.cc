
#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <iostream>
#include "../../src/modules.h"

using namespace snowhouse;

go_bandit([]() {
  bandit::describe(
      "Testing operator name extraction from command string", [&]() {
        bandit::it("returns the commandString if it and has no '-' and has no "
                   "arguments",
                   [&]() {
                     AssertThat(extract_operator_name("test"),
                                Is().EqualTo("test"));
                   });

        bandit::it("returns the commandString without the '-' while it has no "
                   "arguments",
                   [&]() {
                     AssertThat(extract_operator_name("-test"),
                                Is().EqualTo("test"));
                   });

        bandit::it(
            "returns the commandString without the '-' while it has arguments",
            [&]() {
              AssertThat(extract_operator_name("-test,arg"),
                         Is().EqualTo("test"));
            });

        bandit::it(
            "returns the commandString while it has arguments and has no '-'",
            [&]() {
              AssertThat(extract_operator_name("test,arg"),
                         Is().EqualTo("test"));
            });
      });

  bandit::describe(
      "Testing operator name and argument extraction from command string",
      [&]() {
        bandit::it("returns the right name and argument with single "
                   "argument and usage of '-'",
                   [&]() {
                     std::string operName;
                     std::string operArgument;
                     extract_name_and_argument("-test,arg", operName,
                                               operArgument);
                     AssertThat(operName, Is().EqualTo("test"));
                     AssertThat(operArgument, Is().EqualTo("arg"));
                   });
        bandit::it("does not cut off multiple arguments with no '-'", [&]() {
          std::string operName;
          std::string operArgument;
          extract_name_and_argument("test,arg,arg2,arg3", operName,
                                    operArgument);
          AssertThat(operName, Is().EqualTo("test"));
          AssertThat(operArgument, Is().EqualTo("arg,arg2,arg3"));
        });
        bandit::it(
            "works with operators that have no arguments and are written "
            "without '-'",
            [&]() {
              std::string operName;
              std::string operArgument;
              extract_name_and_argument("test", operName, operArgument);
              AssertThat(operName, Is().EqualTo("test"));
              AssertThat(operArgument, Is().EqualTo(""));
            });
        bandit::it("works with operators that have no arguments "
                   "and are written with '-'",
                   [&]() {
                     std::string operName;
                     std::string operArgument;
                     extract_name_and_argument("-test", operName, operArgument);
                     AssertThat(operName, Is().EqualTo("test"));
                     AssertThat(operArgument, Is().EqualTo(""));
                   });
      });
});

int
main(int argc, char **argv)
{

  int result = bandit::run(argc, argv);

  return result;
}
