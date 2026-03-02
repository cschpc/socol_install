
#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <vector>
#include <string.h>

#include "../../src/param_conversion.h"
#include "../../src/cdo_output.h"

/*WARNING inclusion of cc files !!! */
#include "../../src/Select.cc"
#include "../../src/Seltime.cc"
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

/* pre declarations for used functions */
//std::vector<int> get_season_list(std::vector<std::string> &seasonString);
//std::vector<double> get_date_list(std::vector<std::string> &arguments);

go_bandit([]() {
  cdo::set_exit_function(cdoTestExit);
  cdo::set_context_function(getTestContext);

  bandit::describe("generating season lists", [&]() {
    std::vector<std::vector<std::string>> test_input = { { "1", "2", "3", "4" },
                                                         { "2", "4" },
                                                         { "JFM" },
                                                         { "MAM" },
                                                         { "NDJF" },
                                                         { "ANN" } };
    std::vector<std::vector<int>> test_ref
        = { { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
            { 3, 4, 5, 9, 10, 11 },
            { 1, 2, 3 },
            { 3, 4, 5 },
            { 1, 2, 11, 12 },
            { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } };
    bandit::it("has the right size and values", [&]() {
      for (size_t t = 0; t < test_input.size(); t++)
        {
          // std::cout << "test t: " << t << std::endl;
          std::vector<int> res = get_season_list(test_input[t]);
          bandit::it("accecpts input ", [&]() {
            AssertThat(res.size(), snowhouse::Equals(test_ref[t].size()));
            for (size_t r = 0; r < test_ref[t].size(); r++)
              {
                AssertThat(res[r], snowhouse::Equals(test_ref[t][r]));
              }
          });
        }
    });
  });

  bandit::describe("generating date lists", [&]() {
    std::vector<std::vector<std::string>> test_input
        = { { "2005-12-24", "2006-06-24" },
            { "2006-12-24" },
            { "2006-12-24T13:00:37" },
            { "2006-12-24T10:37:00", "-02006-12-24T00:02:00" },
            { "-", "2006-12-24T00:37:22" },
            { "2006-12-24T00:37:22", "-" } };

    std::vector<std::vector<double>> test_ref
        = { { 20051224.0, 20060624.999 },
            { 20061224.0, 20061224.999 },
            { 20061224.130037 },
            { 20061224.1037, -20061224.000200 },
            { -99999999999.0, 20061224.003722 },
            { 20061224.003722, 99999999999.0 } };
    for (size_t t = 0; t < test_input.size(); t++)
      {
        std::vector<double> res = get_date_list(test_input[t]);

        std::string message
            = "test: " + std::to_string(t + 1) + " has the right size";
        bandit::it(message, [&]() {
          AssertThat(res.size(), snowhouse::Equals(test_ref[t].size()));
          message = "test " + std::to_string(t + 1) + " has the right values";
          bandit::it(message, [&]() {
            std::string a;
            for (size_t r = 0; r < res.size() || r < test_ref[t].size(); r++)
              {
                AssertThat(res[r],
                           snowhouse::EqualsWithDelta(test_ref[t][r], 0.001));
              }
          });
        });
      }
  });
});

int
main(int argc, char **argv)
{
  int result = bandit::run(argc, argv);

  return result;
}
