/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

/*
   This module contains the following operators:
*/

#include <cdi.h>

#include <climits>
#include <string>
#include <vector>
#include <iterator>

#include "cdo_options.h"
#include "process_int.h"
#include "counter.h"
#include "varray.h"

static bool mode_auto = false;
static bool mode_demo = false;

static VarList varList;

static CdiStreamID gl_streamID = 0;
static int gl_vlistID = 0;
static int gl_varID = 0;
static int gl_nvars = 0;
static int gl_levelID = 0;
static int gl_tsID1 = 0;
static int gl_tsID2 = 0;
static int gl_ntsteps = -1;
static Varray<double> gl_data;

static int Done = 0;

static int com_help(const std::string &);
static int com_list(const std::string &);
static int com_quit(const std::string &);
static int com_stat(const std::string &);
static int com_set(const std::string &);
static int com_vars(const std::string &);

struct command_t
{
  int (*func)(const std::string &);  // Function to call to do the job.
  const std::string name;            // User printable name of the function.
  const std::string doc;             // Documentation for this function.
};

static const std::vector<command_t> commands = { { com_help, "help", "Display this text" },
                                                 { com_help, "?", "Synonym for 'help'" },
                                                 { com_list, "list", "List files in DIR" },
                                                 { com_quit, "quit", "Quit using CDO" },
                                                 { com_stat, "stat", "Statistic for selected field" },
                                                 { com_set, "set", "set variables" },
                                                 { com_vars, "vars", "list variables" } };
static const int ncommands = commands.size();

// Return non-zero if ARG is a valid argument for CALLER, else print an error message and return zero.
/*
static int
valid_argument(char *caller, char *arg)
{
  if (!arg || !*arg)
    {
      fprintf(stderr, "%s: Argument required.\n", caller);
      return 0;
    }
  return 1;
}
*/
// Print out help for ARG, or for all of the commands if ARG is not present.
static int
com_help(const std::string &arg)
{
  int printed = 0;

  for (int i = 0; i < ncommands; ++i)
    {
      if (arg.empty() || (arg == commands[i].name))
        {
          printf("%s\t\t%s.\n", commands[i].name.c_str(), commands[i].doc.c_str());
          printed++;
        }
    }

  if (!printed)
    {
      printf("No commands match '%s'. Possibilties are:\n", arg.c_str());
      for (int i = 0; i < ncommands; ++i)
        {
          if (printed == 6)  // Print in six columns
            {
              printed = 0;
              printf("\n");
            }
          printf("%s\t", commands[i].name.c_str());
          printed++;
        }

      if (printed) printf("\n");
    }

  return 0;
}

// List the file(s) named in arg
static int
com_list(const std::string &arg)
{
  (void) (arg);

  return 0;
}

// The user wishes to quit using this program. Just set DONE non-zero.
static int
com_quit(const std::string &arg)
{
  (void) (arg);

  Done = 1;

  return 0;
}

static int
com_stat(const std::string &arg)
{
  (void) (arg);  // unused

  auto name = varList[gl_varID].name.c_str();

  auto tsID2 = gl_tsID2;
  if (tsID2 == -1) tsID2 = (gl_ntsteps != -1) ? gl_ntsteps - 1 : INT_MAX - 1;

  for (int tsID = gl_tsID1; tsID <= tsID2; ++tsID)
    {
      const auto nrecs = streamInqTimestep(gl_streamID, tsID);
      if (nrecs == 0)
        {
          if (gl_ntsteps == -1)
            gl_ntsteps = tsID;
          else
            fprintf(stderr, "Timestep %d out of range!\n", tsID + 1);
          break;
        }
      else
        {
          cdo::Counter counter;
          counter.start();

          const auto nlevels = varList[gl_varID].nlevels;
          const auto gridsize = varList[gl_varID].gridsize;
          const auto missval = varList[gl_varID].missval;

          const auto levelID = (nlevels > 1) ? gl_levelID : 0;

          size_t nmiss;
          streamReadVarSlice(gl_streamID, gl_varID, levelID, gl_data.data(), &nmiss);

          const auto mmm = varray_min_max_mean_mv(gridsize, gl_data, missval);

          counter.stop();
          fprintf(stdout, "%s:  z=%d  t=%d  size=%zu nmiss=%zu  min=%.5g mean=%.5g max=%.5g [%.2fs]\n", name, levelID + 1, tsID + 1,
                  gridsize, nmiss, mmm.min, mmm.mean, mmm.max, counter.cputime());
        }
    }

  return 0;
}

static inline void
paramWarning(const char *name, const char *cstring, const char *endptr)
{
  fprintf(stdout, "%s parameter >%s< contains invalid character at position %d!\n", name, cstring, (int) (endptr - cstring + 1));
}

int
param2int(const std::string &str)
{
  char *endptr = nullptr;
  const int ival = (int) strtol(str.c_str(), &endptr, 10);
  if (*endptr) paramWarning("Integer", str.c_str(), endptr);
  return ival;
}

static int
check_tsID(int tsID)
{
  if (gl_ntsteps >= 0 && (tsID < -1 || tsID >= gl_ntsteps))
    {
      fprintf(stdout, "t=%d out of range (max=%d)!\n", (tsID >= 0) ? tsID + 1 : tsID, gl_ntsteps);
      return 1;
    }

  return 0;
}

static void
set_tsID(int tsID1, int tsID2)
{
  if (check_tsID(tsID1)) return;
  if (check_tsID(tsID2)) return;
  gl_tsID1 = tsID1;
  gl_tsID2 = tsID2;

  if (mode_auto) com_stat("");
}

static int
check_levelID(int levelID)
{
  const auto nlevels = varList[gl_varID].nlevels;
  if (levelID >= nlevels)
    {
      fprintf(stdout, "z=%d out of range (max=%d)!\n", levelID + 1, nlevels);
      return 1;
    }

  return 0;
}

static void
set_levelID(int levelID)
{
  if (check_levelID(levelID)) return;
  gl_levelID = levelID;

  if (mode_auto) com_stat("");
}

static int
check_varID(int varID)
{
  if (varID < 0 || varID >= (int) varList.size())
    {
      fprintf(stdout, "varID out of range (max=%d)!\n", (int) varList.size());
      return 1;
    }

  return 0;
}

static void
set_varID(int varID)
{
  if (check_varID(varID)) return;
  gl_varID = varID;

  if (mode_auto) com_stat("");
}

static void
set_t(const std::vector<std::string> &argv)
{
  auto argc = argv.size();
  if (argc == 1)
    {
      fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
      return;
    }
  else if (argc > 3)
    {
      fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
      return;
    }

  auto t1 = param2int(argv[1]);
  auto t2 = (argc == 3) ? param2int(argv[2]) : t1;
  set_tsID((t1 > 0) ? t1 - 1 : t1, (t2 > 0) ? t2 - 1 : t2);
}

static void
set_z(const std::vector<std::string> &argv)
{
  auto argc = argv.size();
  if (argc == 1)
    {
      fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
      return;
    }
  else if (argc > 2)
    {
      fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
      return;
    }

  auto z = param2int(argv[1]);
  set_levelID(z - 1);
}

static void
set_var(const std::vector<std::string> &argv)
{
  auto argc = argv.size();
  if (argc == 1)
    {
      fprintf(stdout, "  set %s: Too few arguments\n", argv[0].c_str());
      return;
    }
  else if (argc > 2)
    {
      fprintf(stdout, "  set %s: Too many arguments\n", argv[0].c_str());
      return;
    }

  auto &name = argv[1];
  auto lfound = false;
  for (int varID = 0; varID < gl_nvars; ++varID)
    {
      if (name == varList[varID].name)
        {
          lfound = true;
          set_varID(varID);
          break;
        }
    }

  if (!lfound)
    {
      fprintf(stdout, "  set %s: Variable name <%s> not found!\n", argv[0].c_str(), name.c_str());
      return;
    }
}

static int
com_set(const std::string &arg)
{
  printf("com_set: >%s<\n", arg.c_str());
  if (arg.empty())
    {
      fprintf(stdout, "  command set: argument missing!\n");
      return -1;
    }

  std::istringstream iss(arg);
  std::vector<std::string> argv(std::istream_iterator<std::string>{ iss }, std::istream_iterator<std::string>());

  // for (int i = 0; i < argv.size(); ++i) printf(">%s<\n", argv[i].c_str());

  if (argv[0] == "t") set_t(argv);
  if (argv[0] == "z")
    set_z(argv);
  else if (argv[0] == "var")
    set_var(argv);

  return 0;
}

static int
com_vars(const std::string &arg)
{
  char paramstr[32];

  printf("com_vars: %s %d\n", arg.c_str(), gl_nvars);

  for (int varID = 0; varID < gl_nvars; ++varID)
    {
      cdiParamToString(varList[varID].param, paramstr, sizeof(paramstr));

      fprintf(stdout, "varID=%3d, param=%s, name=%s, longname=\"%s\", units=\"%s\"\n", varID + 1, paramstr,
              varList[varID].name.c_str(), varList[varID].longname.c_str(), varList[varID].units.c_str());
    }

  return 0;
}

/* Look up NAME as the name of a command, and return a pointer to that command. Return a nullptr pointer if NAME isn't a command
 * name.
 */
static const command_t *
find_command(const std::string &name)
{
  for (int i = 0; i < ncommands; ++i)
    if (name == commands[i].name) return &commands[i];

  return (command_t *) nullptr;
}

// Execute a command line.
static int
execute_line(const std::string &line)
{
  if (line.empty()) return 0;

  // Isolate the command word.
  int i = 0;
  while (line[i] && isspace(line[i])) i++;
  int pos = i;
  while (line[i] && !isspace(line[i])) i++;
  int count = i;

  const auto word = line.substr(pos, count);

  const command_t *command = find_command(word);
  if (!command)
    {
      fprintf(stderr, "%s: No such command!\n", word.c_str());
      return -1;
    }
  // Get argument to command, if any.
  while (isspace(line[i])) i++;

  pos = i;
  auto args = line.substr(pos);

  // Call the function.
  return (*(command->func))(args);
}

std::string
trim(const std::string &str, const std::string &chars = "\t\n\v\f\r ")
{
  const auto first = str.find_first_not_of(chars);
  if (std::string::npos == first) return str;
  const auto last = str.find_last_not_of(chars);

  return str.substr(first, (last - first + 1));
}

extern "C" size_t getPeakRSS();

static std::string
peakRSS_string()
{
  char memstring[32] = { "" };
  size_t memmax = getPeakRSS();
  if (memmax)
    {
      size_t muindex = 0;
      static const char *mu[] = { "B", "KB", "MB", "GB", "TB", "PB" };
      static const size_t nmu = sizeof(mu) / sizeof(char *);
      while (memmax > 9999 && muindex < nmu - 1)
        {
          memmax /= 1024;
          muindex++;
        }
      std::snprintf(memstring, sizeof(memstring), "%zu%s", memmax, mu[muindex]);
    }

  return memstring;
}

static void
read_line(const std::string &prompt, std::string &line)
{
  fputs(prompt.c_str(), stdout);
  if (Options::cdoVerbose)
    {
      fputs(" [", stdout);
      fputs(peakRSS_string().c_str(), stdout);
      fputs("]", stdout);
    }
  fputs("> ", stdout);
  fflush(stdout);

  std::getline(std::cin, line);
}

static void
command_init()
{
  gl_vlistID = streamInqVlist(gl_streamID);

  const auto taxisID = vlistInqTaxis(gl_vlistID);
  (void) (taxisID);  // unused

  auto ntsteps = vlistNtsteps(gl_vlistID);
  if (ntsteps == 0) ntsteps = 1;
  gl_ntsteps = ntsteps;

  const auto gridsizemax = vlistGridsizeMax(gl_vlistID);
  gl_data.resize(gridsizemax);

  gl_nvars = vlistNvars(gl_vlistID);
  varListInit(varList, gl_vlistID);

  set_varID(0);
}

static void
run_demo()
{
  gl_tsID1 = 0;
  gl_tsID2 = -1;
  for (int varID = 0; varID < gl_nvars; ++varID)
    {
      const auto nlevels = varList[varID].nlevels;
      for (int levelID = 0; levelID < nlevels; ++levelID)
        {
          set_varID(varID);
          set_levelID(levelID);
          com_stat(varList[varID].name);
        }
    }
}

void *
Command(void *process)
{
  cdo_initialize(process);

  if (cdo_operator_argc() == 1)
    {
      if (cdo_operator_argv(0) == "auto")
        mode_auto = true;
      else if (cdo_operator_argv(0) == "demo")
        mode_demo = true;
      else
        cdo_abort("Unsupported parameter: %s", cdo_operator_argv(0));
    }

  gl_streamID = streamOpenRead(cdo_get_stream_name(0));

  command_init();

  if (mode_demo) { run_demo(); }
  else
    {
      // Loop reading and executing lines until the user quits.
      const std::string prompt = "cdo cmd";
      while (!Done)
        {
          std::string line;
          read_line(prompt, line);
          execute_line(trim(line));
        }
    }

  streamClose(gl_streamID);

  cdo_finish();

  return nullptr;
}
