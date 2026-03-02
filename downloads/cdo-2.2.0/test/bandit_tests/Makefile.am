noinst_PROGRAMS = wildcards Seltime_test param_conversion_test module_definitions process_init module_interface parser string_util

tests_CPPFLAGS = -I$(top_srcdir)/test/bandit_tests -I$(top_srcdir)/test/bandit_tests/bandit -I$(top_srcdir)/libcdi/src -std=c++20
tests_LDADD = $(top_builddir)/src/libcdo.la $(top_builddir)/libcdi/src/libcdi.la
tests_LDADD   += $(top_builddir)/src/lib/healpix/libhealpix.la

parser_SOURCES = parser.cc
parser_CPPFLAGS = $(tests_CPPFLAGS)
parser_LDADD = $(tests_LDADD)

param_conversion_test_SOURCES = param_conversion_test.cc
param_conversion_test_CPPFLAGS = $(tests_CPPFLAGS)
param_conversion_test_LDADD = $(tests_LDADD)

Seltime_test_SOURCES = Seltime_test.cc
Seltime_test_CPPFLAGS = $(tests_CPPFLAGS)
Seltime_test_LDADD = $(tests_LDADD)

wildcards_SOURCES = wildcards.cc
wildcards_CPPFLAGS = $(tests_CPPFLAGS)
wildcards_LDADD = $(tests_LDADD)

module_definitions_SOURCES = module_definitions.cc
module_definitions_CPPFLAGS = $(tests_CPPFLAGS)
module_definitions_LDADD = $(tests_LDADD)

process_init_SOURCES = process_init.cc
process_init_CPPFLAGS = $(tests_CPPFLAGS)
process_init_LDADD = $(tests_LDADD)

module_interface_SOURCES = module_interface.cc
module_interface_CPPFLAGS = $(tests_CPPFLAGS)
module_interface_LDADD = $(tests_LDADD)

string_util_SOURCES = util_string.cc
string_util_CPPFLAGS = $(tests_CPPFLAGS)
string_util_LDADD = $(tests_LDADD)

LD_LIBRARY_PATH = $(tests_LDADD)
TESTS = ./parser ./wildcards ./Seltime_test ./param_conversion_test ./module_definitions ./process_init ./module_interface ./string_util
