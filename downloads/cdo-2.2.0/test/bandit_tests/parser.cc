#include "bandit/bandit/bandit.h"
// BANDIT NEEDS TO BE INCLUDED FIRST!!!

#include <iostream>
#include <vector>
#include <memory>
#include <string>
#include <iterator>

#include "../../src/parser.h"
#include "../../src/node.h"
#include "../../src/modules.h"
#include "../../src/process_int.h"
#include "../../src/cdo_options.h"
#include "test_module_list.h"
#include "../../src/cdo_node_attach_exception.h"
#include "../../src/cdo_exception.h"

// Util stuff tests start at marker: TESTS
using namespace snowhouse;

// Required functions
void
cdoExit()
{
}
const char *
process_inq_prompt(void)
{
  static const char *context = "cdo_test";
  return context;
}

const char *
process_inq_err_prompt(void)
{
  return "cdo (Abort): ";
}

unsigned int
getRequiredElements(std::vector<std::string> input)
{
  unsigned int numNodes = 0;
  for (auto x : input)
    {
      if (x != "[" && x != "]" && x != ":")
        {
          numNodes += 1;
        }
    }
  return numNodes;
}

unsigned int
getNumChildren(std::shared_ptr<Node> root)
{
  if (root->children.size() == 0) return 1;
  unsigned int sum = 1;
  for (auto &c : root->children)
    {
      sum += getNumChildren(c);
    }
  return sum;
}

void
check(std::string description, std::vector<std::string> in, std::string out)
{
  bandit::it(description, [&]() {
    std::string err_location = " ";
    std::string node_structure = "";
    unsigned numChildrenExpected = -1337;
    try
      {
        auto res = Parser::parse(in, process_inq_prompt);
        node_structure = res[0]->to_string();
        numChildrenExpected = getNumChildren(res[0]);
        AssertThat(node_structure, Equals(out));
        AssertThat(res, Is().OfLength(1));
        AssertThat(numChildrenExpected, Equals(getRequiredElements(in)));
      }
    catch (CdoException &exp)
      {
        err_location = "thrown from: " + exp.file + ":" + exp.line;
        AssertionException(exp.what() + std::string("\n") + err_location);
      }
  });
}

void
checkNegative(std::string description, std::vector<std::string> in,
              std::string expected_err_msg)
{
  bandit::it(description, [&, in]() {
    AssertThrows(CdoSyntaxError, Parser::_parse(in, process_inq_prompt));
    AssertThat(LastException<CdoSyntaxError>().what(), Contains(expected_err_msg));
  });
}

void
checkApply(std::string description, std::vector<std::string> in,
           std::string out, unsigned int numChildren)
{
  bandit::it(description, [&]() {
    std::string node_structure;
    std::string err_location = "";
    unsigned numChildrenExpected;
    try
      {
        auto res = Parser::parse(in, process_inq_prompt);
        if (res.size() > 0)
          {
            node_structure = res[0]->to_string();
            numChildrenExpected = getNumChildren(res[0]);
          }
        else
          {
            node_structure = "";
            numChildrenExpected = -1337;
          }
        AssertThat(numChildrenExpected, Equals(numChildren));
        AssertThat(node_structure, Equals(out));
        AssertThat(res, Is().OfLength(1));
      }
    catch (CdoException &exp)
      {
        err_location = "thrown from: " + exp.file + ":" + exp.line;
        AssertionException(exp.what() + std::string("\n") + err_location);
      }
  });
}

// TESTS
go_bandit([]() {
  //==============================================================================
  cdo::progname = "cdo_bandit_test";
  cdo::set_exit_function(cdoExit);
  cdo::set_context_function(process_inq_prompt);

  // cdo::set_debug(1024);
  //-----------------------------Test_01------------------------------------------
  //------------------------------------------------------------------------------
  //
  bandit::describe("Parser", [&]() {
    bandit::describe("The Core Functionality", [&]() {
      check("handles a single operator with in and output files",
            { "-in1_out1", "in", "out" }, "out [ in1_out1 [ in ] ]");
      check("handles a single operator with 2 in files and a single output "
            "file",
            { "-in2_out1", "infile1", "infile2", "out" },
            "out [ in2_out1 [ infile1 infile2 ] ]");
      check("handles a single operator with 0 input and a single output",
            { "-in0_out1", "out" }, "out [ in0_out1 ]");

      check("handles a single operator with 1 input and no output",
            { "-in1_out0", "in" }, "in1_out0 [ in ]");

      check("handles a single operator with variable input and no output",
            { "-inVariable_out0", "infile1", "infile2", "infile3" },
            "inVariable_out0 [ infile1 infile2 infile3 ]");

      check("handles a single operator with 2 input and no output",
            { "-in2_out0", "infile1", "infile2" },
            "in2_out0 [ infile1 infile2 ]");

      check("handles a operator with obase feature "
            "subgroups",
            { "-in2_outObase", "f1", "f2", "obase" },
            "obase [ in2_outObase [ f1 f2 ] ]");

      check("handles a operator that takes only files ",
            { "files_only", "f1", "out" }, "out [ files_only [ f1 ] ]");
    });
    bandit::describe("Variable input Functionality", [&]() {
      check("handles a single operator with variable inputs (1) and a single "
            "outfile",
            { "-inVariable_out1", "infile1", "out" },
            "out [ inVariable_out1 [ infile1 ] ]");

      check("handles a single operator with variable inputs (2) and a single "
            "outfile",
            { "-inVariable_out1", "infile1", "infile2", "out" },
            "out [ inVariable_out1 [ infile1 infile2 ] ]");

      check("handles a single operator with variable inputs (3) and a single "
            "outfile",
            { "-inVariable_out1", "infile1", "infile2", "infile3", "out" },
            "out [ inVariable_out1 [ infile1 infile2 infile3 ] ]");

      check("handles nested subgroups in subgrups",
            { "-inVariable_out1", "[", "[", "infile1", "infile2", "]", "[",
              "infile3", "infile4", "]", "]", "out" },
            "out [ inVariable_out1 [ infile1 infile2 infile3 infile4 ] ]");
    });

    bandit::describe("Subgroup Functionality", [&]() {
      check("handles a multiple nested variable input operators while using "
            "subgroups",
            { "-inVariable_out1", "[", "infile1", "-inVariable_out1", "infile2",
              "infile3", "]", "out" },
            "out [ inVariable_out1 [ infile1 inVariable_out1 [ infile2 infile3 "
            "] ] "
            "]");

      check("handles mixed input (files then operators)",
            { "-inVariable_out1", "[", "infile1", "-in0_out1", "]", "out" },
            "out [ inVariable_out1 [ infile1 in0_out1 ] ]");

      check("handles mixed input (operators then files)",
            { "-inVariable_out1", "[", "-in0_out1", "infile1", "]", "out" },
            "out [ inVariable_out1 [ in0_out1 infile1 ] ]");

      check("handles only files",
            { "-inVariable_out1", "[", "infile1", "infile2", "]", "out" },
            "out [ inVariable_out1 [ infile1 infile2 ] ]");

      check("handles operators with no input",
            { "-inVariable_out1", "[", "-in0_out1", "-in0_out1", "]", "out" },
            "out [ inVariable_out1 [ in0_out1 in0_out1 ] ]");

      check("handles multiple variable input operators",
            { "-inVariable_out1", "[", "-inVariable_out1", "[", "-in0_out1",
              "-in0_out1", "]", "-inVariable_out1", "[", "infile1", "infile2",
              "]", "]", "outfile" },
            "outfile [ inVariable_out1 [ inVariable_out1 [ in0_out1 in0_out1 ] "
            "inVariable_out1 [ infile1 infile2 ] ] ]");
    });

    bandit::describe("Apply feature", [&]() {
      checkApply(
          "handles merge with multiple groups as input",
          { "-inVariable_out1", "[", "[", "-in1_out1",
            "E5ml00_1H_2000-01-01_129-1", "]", "[", "-in1_out1",
            "E5ml00_1H_2000-01-01_152-1", "]", "]", "tmp1" },
          "tmp1 [ inVariable_out1 [ in1_out1 [ E5ml00_1H_2000-01-01_129-1 ] "
          "in1_out1 [ E5ml00_1H_2000-01-01_152-1 ] ] ]",
          6);

      checkApply("handles apply with chains as argument",
                 { "-inVariable_out0", "[", "-in1_out1", "-in1_out1", ":", "f1",
                   "f2", "]" },
                 "inVariable_out0 [ in1_out1 [ in1_out1 [ f1 ] ] in1_out1 [ "
                 "in1_out1 [ "
                 "f2 ] ] ]",
                 7);

      // DO NOT SOURROUND THE APPLY WITH \" does not work
      checkApply(
          "handles the old way apply worked",
          { "-inVariable_out1", "-apply,-in1_out1", "[", "infile1", "infile2",
            "infile3", "]", "out" },
          "out [ inVariable_out1 [ in1_out1 [ infile1 ] in1_out1 [ infile2 ] "
          "in1_out1 [ infile3 ] ] ]",
          8);

      checkApply(
          "Apply symbol '[ : ]' works",
          { "-inVariable_out1", "[", "-in1_out1", ":", "infile1", "infile2",
            "infile3", "]", "out" },
          "out [ inVariable_out1 [ in1_out1 [ infile1 ] in1_out1 [ infile2 ] "
          "in1_out1 [ infile3 ] ] ]",
          8);
    });
  });

  bandit::describe("Error handling:", [&]() {
    bandit::describe("Syntax Errors:", [&]() {
      std::vector<std::string> argv = { "in", "infile2" };
      checkNegative("fails on file at pos 1", argv,
                                    Parser::err_msg_oper_not_found(argv[0]));

      argv = { "in", "in1_out1", "infile2", "out" };
      checkNegative(
          "fails on file at pos 1 with other operators follwing", argv,
          Parser::err_msg_oper_not_found(argv[0]));

      argv = { "in2_out1",      "[",  "infile1", "infile2", "]",
               "file_too_much", "out" };
      checkNegative(
          "aborts with an unattached file after subgroup", argv,
          Parser::errmsg_unprocessed_inputs);

      checkNegative(
          "aborts detects too much inputs",
          { "-in2_out1", "-in0_out1", "-in0_out1", "-in0_out1", "out" },
          Parser::errmsg_unprocessed_inputs);

      argv = { "only_a_file" };
      checkNegative(
          "aborts when no in- and output are present (file)", argv,
          Parser::err_msg_oper_not_found(argv[0]));

      argv = { "-in1_out1" };
      checkNegative(
          "aborts when no in- and output are present (in1 out1)", argv,
          Parser::errmsg_missing_outputs);
      argv = { "-in1_out0" };
      checkNegative(
          "aborts when no in- and output are present (in1 out0)", argv,
          Parser::errmsg_missing_inputs);

      argv = { "-in0_out1" };
      checkNegative(
          "aborts when no in- and output are present (in0 out1)", argv,
          Parser::errmsg_missing_outputs);

      checkNegative(
          "error on a single operator with variable inputs (0) and a single "
          "outfile",
          { "-inVariable_out1", "out" }, Parser::errmsg_missing_inputs);

      checkNegative(
          "multiple variable input operators are not allowed without "
          "grouping "
          "feature ",
          { "-inVariable_out1", "-inVariable_out1", "infile1", "infile2",
            "out" },
          Parser::errmsg_multiple_variable);
    });
    bandit::describe("Apply Errors:", [&]() {
      checkNegative("detects arguments with multiple inputs ",
                                    { "inVariable_out1", "[", "in2_out1", ":",
                                      "infile1", "infile2", "infile3", "]",
                                      "out" },
                                    Parser::errmsg_only_1_to_1_operators);

      checkNegative("apply detects if a file is in front",
                                    { "inVariable_out1", "WRONG", "[",
                                      "-in1_out1", ":", "infile1", "infile2",
                                      "infile3", "]", "out" },
                                    Parser::errmsg_mixed_input);

      checkNegative("apply detects if a operator is in front",
                                    { "inVariable_out1", "-in0_out1", "[",
                                      "-in1_out1", ":", "infile1", "infile2",
                                      "]", "out" },
                                    Parser::errmsg_mixed_input);

      checkNegative(
          "apply detects if a subgroup returns to many for the target",
          { "inVariable_out1", "-in1_out1", "[", "-in1_out1", ":", "infile1",
            "infile2", "]", "out" },
          errmsg_node_to_many_inputs);

      checkNegative(
          "aborts when apply is used as first operator",
          { "-apply,-in1_out1", "[", "-in0out1", "]", "out" },
          Parser::errmsg_apply_in_first_pos);

      checkNegative(
          "aborts when old apply argument contains unkown operator",
          { "-inVariable_out1", "-apply,NOOPER", "[", "f1", "f2", "]", "out" },
          Parser::err_msg_oper_not_found("NOOPER"));
      checkNegative(
          "aborts when old apply has no inputs",
          { "-inVariable_out1", "-apply,-in1_out1", "out" },
          Parser::errmsg_apply_requires_bracket);

      checkNegative(
          "aborts when old apply has no inputs but '[ ]'",
          { "-inVariable_out1", "-apply,-in1_out1", "[", "]", "out" },
          Parser::errmsg_apply_missing_argument);

      checkNegative(
          "aborts when new apply has no inputs",
          { "-inVariable_out1", "[", "-in1_out1", ":", "]", "out" },
          Parser::errmsg_apply_missing_argument);

      checkNegative("aborts when old apply has no inputs and "
                                    "variable input has no output",
                                    { "-inVariable_out0", "-apply,-in1_out1" },
                                    Parser::errmsg_apply_requires_bracket);
      checkNegative(
          "aborts when old apply has the '[' but nothing else",
          { "-inVariable_out0", "-apply,-in1_out1", "[" },
          Parser::errmsg_bracket_not_closed);
      checkNegative("detects a missing duplicate bracket",
                                    { "inVariable_out1", "[", "in2_out1", ":",
                                      "infile1", "infile2", "infile3", "out" },
                                    Parser::errmsg_bracket_not_closed);
      checkNegative(
          "apply only allows chains with single in and output",
          { "inVariable_out1", "[", "-in2_out0", "-in0_out1", "-in1_out1", ":",
            "infile1", "infile2", "]", "out" },
          Parser::errmsg_only_1_to_1_operators);
    });
    bandit::describe("SubGroups Errors:", [&]() {
      checkNegative(
          "handles nested subgroups in subgrups with too many brackets",
          { "-inVariable_out1", "[", "[", "infile1", "infile2", "]", "[", "[",
            "infile3", "infile4", "]", "]", "out" },
          Parser::errmsg_bracket_not_closed);
      checkNegative(
          "empty [ ] are not ignored",
          { "-in2_out0", "infile1", "[", "]", "infile2", "[", "]" },
          Parser::errmsg_empty_subgroup);

      checkNegative(
          "missing ']' bracket detected",
          { "inVariable_out1", "[", "infile1", "infile2", "infile3", "out" },
          Parser::errmsg_bracket_not_closed);

      checkNegative(
          "missing '[' bracket detected",
          { "inVariable_out1", "infile1", "infile2", "infile3", "]", "out" },
          Parser::errmsg_missing_sub_group);

      checkNegative(
          "handles nested subgroups in subgrups with too many '[' brackets",
          { "-inVariable_out1", "[", "[", "infile1", "infile2", "]", "[", "[",
            "infile3", "infile4", "]", "]", "out" },
          Parser::errmsg_bracket_not_closed);

      checkNegative(
          "handles nested subgroups in subgrups with too many ']' brackets",
          { "-inVariable_out1", "[", "[", "infile1", "infile2", "]", "]", "[",
            "infile3", "infile4", "]", "]", "out" },
          Parser::errmsg_missing_sub_group);
    });

    bandit::describe("Errors handled by Node instead of Parser", [&]() {
      checkNegative("abort when only_file operator has pipe",
                                    { "files_only", "-in0_out1", "out" },
                                    errmsg_node_only_accepts_files);
      checkNegative(
          "using 0 output operators as input for other node is caught",
          { "-in1_out1", "-in1_out0", "out" }, errmsg_node_no_output);

      checkNegative(
          "sub groups cannot overflow",
          { "-in2_out1", "[", "infile1", "infile2", "infile3", "]", "out" },
          errmsg_node_to_many_inputs);

      checkNegative("detects malformed e.g. subgroup operators without inptu",
                                    { "-inVariable_out1", "[", "[", "-in2_out1",
                                      "]", "[", "infile2", "infile3", "]","]",
                                      "out" },
                                    Parser::errmsg_malformed_subgroup);
    });
  });
});

//==============================================================================
#define EXCEPTION_EXTRA_INFO = 1;
int
main(int argc, char **argv)
{
  std::vector<char *> argv_v(argv, argv + argc);
  std::string reporter = "--reporter=spec";
  argv_v.push_back(&reporter[0]);
  int result = bandit::run(argc + 1, argv_v.data());

  return result;
}
