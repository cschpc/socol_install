#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <vector>
#include <string.h>

#include "../../src/param_conversion.h"

void
cdoTestExit()
{
  exit(EXIT_FAILURE);
}
std::string testContext = "SeltimeTest";
const char *
getTestContext()
{
  return testContext.c_str();
}

go_bandit([]() {
  bandit::describe("testing split_intstring", []() {
    //----------------------------------------------------------------
    bandit::it("converts simple start to end string", [&]() {
      std::vector<std::string> test = { "0/4" };
      std::vector<int> expected_result = { 0, 1, 2, 3, 4 };
      std::vector<int> result = cdo_argv_to_int(test);
      AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
    });
    //----------------------------------------------------------------
    bandit::it("converts input with negative start", [&]() {
      std::vector<std::string> test = { "-2/2" };
      std::vector<int> expected_result = { -2, -1, 0, 1, 2 };
      std::vector<int> result = cdo_argv_to_int(test);
      AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
    });
    //----------------------------------------------------------------
    bandit::it("converts input with positive start and negative end", [&]() {
      std::vector<std::string> test = { "2/-2" };
      std::vector<int> expected_result = { 2, 1, 0, -1, -2 };
      std::vector<int> result = cdo_argv_to_int(test);
      AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
    });
    //----------------------------------------------------------------
    bandit::it("converts string with increment of 2", [&]() {
      std::vector<std::string> test = { "0/6/2" };
      std::vector<int> expected_result = { 0, 2, 4, 6 };
      std::vector<int> result = cdo_argv_to_int(test);
      AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
    });
    //----------------------------------------------------------------
    bandit::it("converts string with negative start and increment of 3", [&]() {
      std::vector<std::string> test = { "-4/4/3" };
      std::vector<int> expected_result = { -4, -1, 2 };
      std::vector<int> result = cdo_argv_to_int(test);
      AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
    });
    //----------------------------------------------------------------
    bandit::it(
        "converts string positve start and negative end and negative increment",
        [&]() {
          std::vector<std::string> test = { "2/-2/-1" };
          std::vector<int> expected_result = { 2, 1, 0, -1, -2 };
          std::vector<int> result = cdo_argv_to_int(test);
          AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
        });
    //----------------------------------------------------------------
    bandit::it("returns empty array if start is positive end end is negative "
               "and positive increment",
               [&]() {
                 std::vector<std::string> test = { "2/-2/1" };
                 std::vector<int> expected_result = {};
                 std::vector<int> result = cdo_argv_to_int(test);
                 AssertThat(result,
                            snowhouse::Is().EqualToContainer(expected_result));
               });
    //----------------------------------------------------------------
    bandit::it("returns empty array if start is negative end end is positive "
               "and positive increment",
               [&]() {
                 std::vector<std::string> test = { "-2/2/-1" };
                 std::vector<int> expected_result = {};
                 std::vector<int> result = cdo_argv_to_int(test);
                 AssertThat(result,
                            snowhouse::Is().EqualToContainer(expected_result));
               });
    //----------------------------------------------------------------
    bandit::it(
        "returns empty array if start is larger than end and inc is positive "
        "and positive increment",
        [&]() {
          std::vector<std::string> test = { "2/-2/1" };
          std::vector<int> expected_result = {};
          std::vector<int> result = cdo_argv_to_int(test);
          AssertThat(result, snowhouse::Is().EqualToContainer(expected_result));
        });
    //----------------------------------------------------------------
  });
});

int
main(int argc, char **argv)
{
  int result = bandit::run(argc, argv);

  return result;
}

