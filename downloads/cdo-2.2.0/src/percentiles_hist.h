#ifndef PERCENTILES_HIST_H_
#define PERCENTILES_HIST_H_

#include <cassert>

#include "field.h"

struct Histogram
{
  void *ptr = nullptr;
  float min = 0.0f;
  float max = 0.0f;
  float step = 0.0f;
  int nsamp = 0;
  int capacity = 0;
  short nbins = 0;
  bool isUint32 = false;
};

// clang-format off
class  // HistogramSet
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
HistogramSet
// clang-format on
{
private:
  int nvars = 0;
  int nsteps = 0;
  std::vector<int> var_nlevels;
  std::vector<size_t> var_nhists;
  std::vector<std::vector<std::vector<Histogram>>> histograms;

  void
  init()
  {
    var_nlevels.resize(nvars, 0);
    var_nhists.resize(nvars, 0);
    histograms.resize(nvars);
  }

public:
  HistogramSet() {}

  explicit HistogramSet(int _nvars) : nvars(_nvars)
  {
    assert(nvars > 0);
    init();
  }

  HistogramSet(int _nvars, int _nsteps) : nvars(_nvars), nsteps(_nsteps)
  {
    assert(nvars > 0);
    init();
  }

  void
  create(int _nvars, int _nsteps = 0)
  {
    nvars = _nvars;
    nsteps = _nsteps;
    assert(nvars > 0);
    init();
  }

  ~HistogramSet()
  {
    for (auto varID = nvars; varID-- > 0;)
      {
        const auto nhists = this->var_nhists[varID];
        for (auto levelID = this->var_nlevels[varID]; levelID-- > 0;)
          {
            for (auto histID = nhists; histID-- > 0;) free(this->histograms[varID][levelID][histID].ptr);
          }
      }
  }

  void createVarLevels(int varID, int nlevels, size_t nhists);
  void defVarLevelBounds(int varID, int levelID, const Field &field1, const Field &field2);
  int addVarLevelValues(int varID, int levelID, const Field &field);
  int subVarLevelValues(int varID, int levelID, const Field &field);
  void getVarLevelPercentiles(Field &field, int varID, int levelID, double p);
  // void reset(int varID, int levelID); // unused
};

#endif /* PERCENTILES_HIST_H_ */
