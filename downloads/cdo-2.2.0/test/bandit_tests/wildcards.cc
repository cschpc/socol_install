#include "bandit/bandit/bandit.h"

#ifdef HAVE_CONFIG_H
#include "../../src/config.h"
#endif

#ifdef HAVE_WORDEXP_H
#include "../../src/util_wildcards.h"
#include <fstream>

//==========================================================================
go_bandit([]() {

  // File name rules
  const std::string prefix = "testFile";
  const std::string suffix = ".banditTestFile";
  std::vector<std::string> wildcards = {"*", "?"};
  std::vector<std::string> testArgv;
  std::string toBeExpanded;

  const int fileCount = 5;

  // +1 because there has to be another string at the binning of testArgv since
  // the first entry will be ignored by expand wildcards

  // storage for temp files;
  std::ofstream files[fileCount];
  std::vector<std::string> fileNames;

  //--------------------------------------------------------------------------
  std::string fileName;
  // Creating test files
  for (int fileID = 0; fileID < fileCount; fileID++) {
    fileName = prefix + std::to_string(fileID) + suffix;
    files[fileID].open(fileName);
    fileNames.push_back(fileName);
    files[fileID].close();
  }
  // see comment at definition of expectedFileCnt
  testArgv.push_back("wildcards");
  // setting the actual filenames
  for (unsigned int wildCardIDX = 0; wildCardIDX < wildcards.size(); wildCardIDX++) {
    testArgv.push_back(prefix + wildcards[wildCardIDX] + suffix);
  }

  //for ignoring non wildcard elements as '[' and ']'
  int k = 2;
  //--------------------------------------------------------------------------
  bandit::describe("Expanding wildcard", [&]() {
    std::vector<std::string> expandedWildCards = expand_wild_cards(testArgv);
    for (unsigned int j = 0; j < wildcards.size(); j++) {
      for (int i = 0; i < fileCount; i++) {
        int argvIdx = k + i + (fileCount * j);
        bandit::it("has expanded the filename", [&]() {
          AssertThat(expandedWildCards[argvIdx],
                     snowhouse::Equals(fileNames[i]));
        });
      }
      k += 2;
    }
    //--------------------------------------------------------------------------
  });
});

int main(int argc, char **argv) {
  int result = bandit::run(argc, argv);
  return result;
}
#else
#define SKIPPED_TEST 77
int main(int argc, char **argv) {
  int result = SKIPPED_TEST;
  return result;
}
#endif
