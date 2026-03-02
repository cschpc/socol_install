#include "cdo_query.h"
#include "cdo_options.h"
#include "cdo_output.h"
#include "param_conversion.h"
#include "util_string.h"

static std::vector<int>
decode_int_parameter(const std::string &param)
{
  const std::string separator("/to/");

  std::vector<int> values;
  auto n = param.find(separator);
  if (n == std::string::npos) { values.push_back(parameter_to_int(param)); }
  else
    {
      auto v1 = parameter_to_int(param.substr(0, n));
      auto v2 = parameter_to_int(param.substr(n + separator.size()));
      if (v2 < v1) cdo_abort("Second parameter %d muss be greater than first parameter %d!", v2, v1);
      auto numVals = (v2 - v1 + 1);
      values.resize(numVals);
      for (int i = 0; i < numVals; ++i) values[i] = v1 + i;
    }

  return values;
}

static std::vector<size_t>
decode_cell_parameter(const std::string &param)
{
  const std::string separator("/to/");

  std::vector<size_t> cells;
  auto n = param.find(separator);
  if (n == std::string::npos) { cells.push_back(parameter_to_size_t(param)); }
  else
    {
      auto v1 = parameter_to_size_t(param.substr(0, n));
      auto v2 = parameter_to_size_t(param.substr(n + separator.size()));
      if (v2 < v1) cdo_abort("Second parameter %zu muss be greater than first parameter %zu!", v2, v1);
      cells.push_back(v1);
      auto numVals = (v2 - v1 + 1);
      if (numVals > 1) cells.push_back(numVals);
    }

  return cells;
}

std::string
set_query_parameter(const KVList &kvlist, CdiQuery *query)
{
  std::string path;

  for (const auto &kv : kvlist)
    {
      const auto &key = kv.key;
      int numValues = kv.nvalues;
      if (numValues < 1) cdo_abort("Missing value for parameter key >%s<!", key);

      // clang-format off
      if (key == "name")
        {
          std::vector<char*> queryNames(numValues);
          for (int i = 0; i < numValues; ++i)
            {
              queryNames[i] = (char*)kv.values[i].c_str();
            }
          cdiQuerySetNames(query, queryNames.size(), queryNames.data());
        }
      else if (key == "step")
        {
          std::vector<int> querySteps;
          for (int i = 0; i < numValues; ++i)
            {
              auto steps = decode_int_parameter(kv.values[i]);
              auto numSteps = querySteps.size();
              querySteps.resize(numSteps + steps.size());
              for (size_t k = 0; k < steps.size(); ++k) querySteps[numSteps + k] = steps[k];
            }
          cdiQuerySetStepidx(query, querySteps.size(), querySteps.data());
        }
      else if (key == "cell")
        {
          if (numValues > 1) cdo_abort("Too many values for key=cell (maxvalues=1 or range=start/to/end)");
          auto queryCells = decode_cell_parameter(kv.values[0]);
          cdiQuerySetCellidx(query, queryCells.size(), queryCells.data());
        }
      else if (key == "path")
        {
          path = kv.values[0];
          for (int i = 1; i < numValues; ++i) path += "," + kv.values[i];
        }
      else cdo_abort("Invalid parameter key >%s<!", key);
      // clang-format on
    }

  return path;
}

std::string
set_query_parameter(const std::string &params, CdiQuery *query)
{
  auto paramsArgv = split_with_seperator(params, ',');

  KVList kvlist;
  kvlist.name = "QUERY";
  if (kvlist.parse_arguments(paramsArgv.size(), paramsArgv) != 0) cdo_abort("Parse error!");
  if (Options::cdoVerbose) kvlist.print();

  return set_query_parameter(kvlist, query);
}
