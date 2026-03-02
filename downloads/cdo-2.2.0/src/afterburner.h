#ifndef AFTERBURNER_H
#define AFTERBURNER_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "dmemory.h"
#include "process_int.h"
#include "transform.h"
#include "cdo_task.h"

#define MaxLevel 1024

#define MaxCodes 277

#define S_ECHAM5 1

struct Date
{
  int yr;
  int mo;
  int dy;
  int hr;
  int mn;
};

struct AfterControl
{
  int Verbose;

  int Mean;
  int MeanCount0;
  int MeanCount;
  int Multi;
  int Nfiles;
  int TermCount;

  int OutputInterval;
  int EndOfInterval;

  int AnalysisData;  // 0 = ECHAM Data, 1 = ECMWF Spectral Analyses
  int DayIn;         // day increment of infiles if Multi = true
  int Debug;
  bool extrapolate;
  bool szip;

  int istreamID;
  CdoStreamID ostreamID;
  CdoStreamID ostreamID2;
  int ivlistID;
  int ovlistID;
  int ovlistID2;
  int taxisID;
  int taxisID2;

  struct Date NextDate;
  struct Date NewDate;
  struct Date OldDate;
  struct Date StartDate;

  int nvct;
  double *vct;

  int *vertIndex;
  size_t *pnmiss;
  double *Orography;
  double *p_of_height;

  int Type;
  int unitsel;

  int Fouriers;
  int Latitudes;
  int Longitudes;
  int HalfLevels;
  int Gaussian;
  int Spectral;

  int Truncation;
  int Waves;

  long Dim3FC, Dim3SP, Dim3GP;
  long DimFC, DimGP, DimSP;
  long DimSP_half;

  double *poli;
  double *pold;
  double *pdev;
  double *pol2;
  double *pol3;

  double *dv2uv_f1;
  double *dv2uv_f2;

  int NumCodesRequest;

  int NumLevel;
  int NumLevelFound;
  int NumLevelRequest;
  double LevelRequest[MaxLevel];

  double *rcoslat;
  double *coslat;
  double *DerivationFactor;
  double *Field;
};

struct Variable
{
  int needed0;   // var needed for process
  int needed;    // var needed for process
  int selected;  // var selected for output
  int detected;  // var detected in input
  int comp;      // compute var if selected and not detected
  int sfit;
  int hlev;
  int plev;
  int ivarID;
  int ovarID;   // 1st variable ID
  int ovarID2;  // 2nd variable ID used for variance
  int tableID;
  int igridID;
  int ogridID;
  int izaxisID;
  int ozaxisID;
  size_t nmiss0;
  size_t nmiss;
  double missval;
  double *spectral;
  double *spectral0;
  double *fourier;
  double *hybrid;
  double *hybrid0;
  double *height;
  double *grid;
  double *grid0;
  double *mean;
  double *variance;
  int *samp;
};

// clang-format off
#define    LOW_CLOUD   34
#define    MID_CLOUD   35
#define    HIH_CLOUD   36
#define    LOW_WATER   37  /* not used ?   */
#define    MID_WATER   38  /* not used ?   */
#define    HIH_WATER   39  /* not used ?   */
#define    ALL_WATER   40  /* not used ?   */

#define GEOPOTENTIAL  129
#define  TEMPERATURE  130
#define       U_WIND  131
#define       V_WIND  132
#define     HUMIDITY  133
#define           PS  134
#define        OMEGA  135
#define    VORTICITY  138
#define           TS  139
#define       STREAM  148
#define      VELOPOT  149
#define          SLP  151
#define         LNPS  152
#define   DIVERGENCE  155
#define GEOPOTHEIGHT  156
#define    RHUMIDITY  157

#define   SW_BOT_CLF  189  /* not used ?   */
#define   LW_BOT_CLF  190  /* not used ?   */
#define   SW_TOP_CLF  191  /* not used ?   */
#define   LW_TOP_CLF  192  /* not used ?   */
#define  NET_TOP_CLF  193  /* not computed */

#define    WINDSPEED  259
#define       PRECIP  260
#define      NET_TOP  261
#define      NET_BOT  262
#define     NET_HEAT  263
#define    NET_WATER  264
#define       SW_CLF  265
#define       LW_CLF  266
#define      NET_CLF  267
#define       SW_ATM  268
#define       LW_ATM  269
#define      NET_ATM  270
#define  SURF_RUNOFF  271
#define        DPSDX  273
#define        DPSDY  274
#define  FRESH_WATER  275
#define      PS_PROG  276  // PS for prognostic timestep
#define   HALF_PRESS  277
#define   FULL_PRESS  278
#define       THETAH  279
#define       THETAF  280
// clang-format on

void after_gp2sp(const AfterControl &globs, struct Variable *vars, int ccode);
void after_GP2FC(double *gp, double *fc, long nlat, long nlon, long nlev, long nfc);
void after_FC2GP(double *fc, double *gp, long nlat, long nlon, long nlev, long nfc);
void after_FCrh2FCsh(const AfterControl &globs, struct Variable *vars);
void after_SPuv2SPdv(const AfterControl &globs, struct Variable *vars);
void after_FCsh2FCrh(const AfterControl &globs, struct Variable *vars);

void after_EchamCompGP(const AfterControl &globs, struct Variable *vars);
void after_processPL(AfterControl &globs, struct Variable *vars);
void after_processML(AfterControl &globs, struct Variable *vars);

void after_AnalysisAddRecord(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                             size_t nmiss);
double *after_get_dataptr(struct Variable *vars, int code, int gridID, int zaxisID, int levelID);
void after_EchamAddRecord(const AfterControl *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID,
                          size_t nmiss);

void after_AnalysisDependencies(struct Variable *vars, int ncodes);
void after_EchamDependencies(struct Variable *vars, int ncodes, int type, int source);

void after_legini_setup(AfterControl &globs, struct Variable *vars);

template <typename... Args>
void
afterAbort(const std::string &format, Args const &...args)
{
  extern cdo::Task *afterReadTask;
  extern bool afterReadAsync;
  if (afterReadAsync && afterReadTask) afterReadTask->wait();
  cdo_abort(format, args...);
}

#endif /* AFTERBURNER_H */
