#include "bandit/bandit/bandit.h"
#include "../../src/cdo_options.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include "../../src/util_string.h"
#include <string>
#include <tuple>

using namespace snowhouse;
void
cdoExit()
{
}
go_bandit([]() {
  //==============================================================================

  bandit::describe(
      "Testing function 'tokenize_comma_seperated_int_list'", []() {
        bandit::describe("Testing tokenization with proper input", []() {
          std::string to_be_tokenized{ "1,-2,3,4" };
          std::vector<std::string> expected_tokens = { "1", "-2", "3", "4" };
          std::tuple<bool, std::vector<std::string>> result
              = tokenize_comma_seperated_int_list(to_be_tokenized);
          bandit::it("returns true on success",
                     [&]() { AssertThat(std::get<0>(result), Equals(true)); });
          bandit::it("returns the right tokens", [&]() {
            AssertThat(std::get<1>(result), EqualsContainer(expected_tokens));
          });
        });

        bandit::describe(
            "Testing tokenization with input containing non int", []() {
              std::vector<std::string> expected_tokens = {};
              std::string to_be_tokenized{ "1,a,3,4" };
              std::tuple<bool, std::vector<std::string>> result
                  = tokenize_comma_seperated_int_list(to_be_tokenized);

              bandit::it("returns false on failure", [&]() {
                AssertThat(std::get<0>(result), Equals(false));
              });
              bandit::it("returns a empty list", [&]() {
                AssertThat(std::get<1>(result),
                           EqualsContainer(expected_tokens));
              });
            });
        bandit::describe(
            "Testing tokenization with input containing a float value", []() {
              std::vector<std::string> expected_tokens = {};
              std::string to_be_tokenized{ "1,1.2,3" };
              std::tuple<bool, std::vector<std::string>> result
                  = tokenize_comma_seperated_int_list(to_be_tokenized);

              bandit::it("returns false", [&]() {
                AssertThat(std::get<0>(result), Equals(false));
              });
              bandit::it("returns empty list", [&]() {
                AssertThat(std::get<1>(result),
                           EqualsContainer(expected_tokens));
              });
            });
        bandit::describe(
            "Testing tokenization with input containing a single value", []() {
              std::vector<std::string> expected_tokens = { "1" };
              std::string to_be_tokenized{ "1" };
              std::tuple<bool, std::vector<std::string>> result
                  = tokenize_comma_seperated_int_list(to_be_tokenized);

              bandit::it("returns true", [&]() {
                AssertThat(std::get<0>(result), Equals(true));
              });
              bandit::it("returns the token  '1'", [&]() {
                AssertThat(std::get<1>(result),
                           EqualsContainer(expected_tokens));
              });
            });
      });

  bandit::describe("Testing function 'tokenize_comma_seperated_int_list'",
                   []() {

                   });
});
int
main(int argc, char **argv)
{
  int result = bandit::run(argc, argv);

  return result;
}
