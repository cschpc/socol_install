/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#include <cstdio>
#include <vector>

extern int CDO_parser_errorno;

enum CoordIndex
{
  TIMESTEP = 0,
  DATE,
  TIME,
  DELTAT,
  DAY,
  MONTH,
  YEAR,
  SECOND,
  MINUTE,
  HOUR,
  LEN
};

enum class NodeEnum
{
  typeCon,
  typeVar,
  typeFun,
  typeOpr,
  typeCmd
};

// commands
struct cmdNodeType
{
  char *cmdName;  // command name
  char *varName;  // variable name
};

// constants
struct conNodeType
{
  double value;  // value of constant
};

// variables
struct varNodeType
{
  char *name;  // variable name
};

// functions
struct funNodeType
{
  char *name;                 // function name
  int nops;                   // number of operands
  struct nodeType *op[1];     // operands (expandable)
};

// operators
struct oprNodeType
{
  int oper;                   // operator
  int nops;                   // number of operands
  struct nodeType *op[1];     // operands (expandable)
};

enum class ParamType
{
  UNDEFINED,
  VAR,
  CONST
};

// parameter
struct ParamEntry
{
  ParamType type = ParamType::UNDEFINED;
  bool isValid = false;
  bool select = false;
  bool remove = false;
  bool hasMV = false;
  int coord = 0;
  int gridID = -1;
  int zaxisID = -1;
  int datatype = -1;
  int steptype = -1;
  size_t ngp = 0;
  size_t nlat = 0;
  size_t nlev = 0;
  size_t nmiss = 0;
  char *name = nullptr;
  char *longname = nullptr;
  char *stdname = nullptr;
  char *units = nullptr;
  double *data = nullptr;
  double *weight = nullptr;
  double missval = 0.0;
};

struct nodeType
{
  ParamEntry param;
  NodeEnum type;  // type of node
  bool isTmpObj;

  // union must be last entry in nodeType because oprNodeType adn funNodeType may dynamically increase
  union
  {
    conNodeType con;  // constants
    varNodeType var;  // variables
    cmdNodeType cmd;  // commands
    funNodeType fun;  // functions
    oprNodeType opr;  // operators
  } u;
};

struct CoordType
{
  std::vector<double> data;
  std::vector<char> units;
  std::vector<char> longname;
  size_t size;
  int coord;
  int cdiID;
  bool needed;
};

struct ParseParamType
{
  std::vector<bool> needed;
  std::vector<CoordType> coords;
  std::vector<ParamEntry> params;
  int maxparams;
  int nparams;
  int cnparams;  // current number of valid params
  int nvars1;
  int ncoords;
  int maxCoords;
  int tsID;
  int pointID;
  int zonalID;
  int surfaceID;
  bool init;
  bool debug;
};

typedef union
{
  double cvalue;   // constant value
  char *varnm;     // variable name
  char *fname;     // function name
  nodeType *nPtr;  // node pointer
} yysType;

#define YYSTYPE yysType
#define YY_EXTRA_TYPE ParseParamType *

#define YY_DECL int yylex(YYSTYPE *yylval_param, void *yyscanner)
YY_DECL;

int yyparse(ParseParamType &parseArg, void *);
void yyerror(const ParseParamType &parseArg, void *scanner, const char *errstr);

int yylex_init(void **);
int yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);

nodeType *expr_run(nodeType *p, ParseParamType &parseArg);
int params_get_coord_ID(const ParseParamType &parseArg, int coord, int cdiID);
