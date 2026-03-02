#include <cstdio>
#include <algorithm>

#include "progress.h"

namespace progress
{

bool stdoutIsTerminal = false;
bool silentMode = false;

static bool ps_lhead = false;
static int ps_nch = 0;
static int ps_cval = -1;

static const char *context = "";
const char *(*getContext)(void) = nullptr;

/**
 * parameter p_context:
 *          will be displayed in status message and indicates sub process for
 *          which the progress is displayed.
 */
void
set_context_function(const char *(*func)(void) )
{
  getContext = func;
}

void
init()
{
  ps_lhead = false;
  ps_nch = 0;
  ps_cval = -1;

  if (getContext != nullptr) context = getContext();
}

void
init(const char *p_context)
{
  ps_lhead = false;
  ps_nch = 0;
  ps_cval = -1;

  context = p_context;
}

void
update(double offset, double refval, double curval)
{
  if (silentMode) return;
  if (!stdoutIsTerminal) return;

  offset = std::min(std::max(offset, 0.), 1.);
  refval = std::min(std::max(refval, 0.), 1.);
  curval = std::min(std::max(curval, 0.), 1.);

  const int ival = (offset + refval * curval) * 100;

  if (ps_cval == -1)
    {
      ps_nch = fprintf(stdout, "%s: %3d%%", context, 0);
      fflush(stdout);
      ps_lhead = true;
    }

  if (ival != ps_cval)
    {
      ps_cval = ival;
      fprintf(stdout, "\b\b\b\b%3d%%", ps_cval);
      fflush(stdout);
    }

  if (ps_cval == 100 && ps_lhead)
    {
      ps_lhead = false;
      while (ps_nch--) fprintf(stdout, "\b \b");
      fflush(stdout);
    }
}

}  // namespace progress
