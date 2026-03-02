/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef MODULE_DEFINITIONS_H
#define MODULE_DEFINITIONS_H /* \cond */

// clang-format off

void   *Adisit(void *argument);
#define AdisitOperators                   {"adisit", "adipot"}

void   *Afterburner(void *argument);
#define AfterburnerOperators              {"after"}

void   *Arith(void *argument);
#define ArithOperators                    {"add", "sub", "mul", "div", "min", "max", "atan2", "setmiss"}

void   *Arithc(void *argument);
#define ArithcOperators                   {"addc", "subc", "mulc", "divc", "minc", "maxc", "mod"}

void   *Arithdays(void *argument);
#define ArithdaysOperators                {"muldpm", "divdpm", "muldpy", "divdpy", "muldoy"}

void   *Arithlat(void *argument);
#define ArithlatOperators                 {"mulcoslat", "divcoslat"}

void   *Bitrounding(void *argument);
#define BitroundingOperators              {"bitrounding"}

void   *Cat(void *argument);
#define CatOperators                      {"cat"}

void   *CDItest(void *argument);
#define CDItestOperators                  {"ncopy"}

void   *CDIread(void *argument);
#define CDIreadOperators                  {"cdiread"}

void   *CDIwrite(void *argument);
#define CDIwriteOperators                 {"cdiwrite"}

void   *Change(void *argument);
#define ChangeOperators                   {"chcode", "chtabnum", "chparam", "chname", "chunit", "chlevel", "chlevelc", "chlevelv", "chltype"}

void   *Change_e5slm(void *argument);
#define Change_e5slmOperators             {"change_e5slm", "change_e5lsm", "change_e5mask"}

void   *Cloudlayer(void *argument);
#define CloudlayerOperators               {"cloudlayer"}

void   *CMOR(void *argument);
#define CMOROperators                     {"cmor"}

void   *CMOR_lite(void *argument);
#define CMORliteOperators                 {"cmorlite"}

void   *CMOR_table(void *argument);
#define CMORtableOperators                {"dump_cmor_table", "conv_cmor_table"}

void   *Collgrid(void *argument);
#define CollgridOperators                 {"collgrid"}

void   *Command(void *argument);
#define CommandOperators                  {"command", "com", "cmd"}

void   *Comp(void *argument);
#define CompOperators                     {"eq", "ne", "le", "lt", "ge", "gt"}

void   *Compc(void *argument);
#define CompcOperators                    {"eqc", "nec", "lec", "ltc", "gec", "gtc"}

void   *Complextorect(void *argument);
#define ComplextorectOperators            {"complextorect", "complextopol"}

void   *Cond(void *argument);
#define CondOperators                     {"ifthen", "ifnotthen"}

void   *Cond2(void *argument);
#define Cond2Operators                    {"ifthenelse"}

void   *Condc(void *argument);
#define CondcOperators                    {"ifthenc", "ifnotthenc"}

void   *Consecstat(void *argument);
#define ConsecstatOperators               {"consects", "consecsum"}

void   *Copy(void *argument);
#define CopyOperators                     {"copy", "clone", "szip"}

void   *DCW_util(void *argument);
#define DCW_utilOperators                 {"dcw"}

void   *Dayarith(void *argument);
#define DayarithOperators                 {"dayadd", "daysub", "daymul", "daydiv"}

void   *Deltat(void *argument);
#define DeltatOperators                   {"deltat", "timederivative"}

void   *Deltime(void *argument);
#define DeltimeOperators                  {"delday", "del29feb"}

void   *Depth(void *argument);
#define DepthOperators                    {"zsdepth"}

void   *Derivepar(void *argument);
#define DeriveparOperators                {"gheight", "gheighthalf", "sealevelpressure"}

void   *Detrend(void *argument);
#define DetrendOperators                  {"detrend"}

void   *Diff(void *argument);
#define DiffOperators                     {"diff", "diffp", "diffn", "diffc"}

void   *Distgrid(void *argument);
#define DistgridOperators                 {"distgrid"}

void   *Duplicate(void *argument);
#define DuplicateOperators                {"duplicate"}


void   *Echam5ini(void *argument);
#define Echam5iniOperators                {"import_e5ml", "export_e5ml"}

void   *Enlarge(void *argument);
#define EnlargeOperators                  {"enlarge"}

void   *Enlargegrid(void *argument);
#define EnlargegridOperators              {"enlargegrid"}

void   *Ensstat(void *argument);
#define EnsstatOperators                  {"ensrange", "ensmin", "ensmax", "enssum", "ensmean", "ensavg", "ensvar", "ensvar1", "ensstd", "ensstd1", "ensskew", "enskurt", "ensmedian", "enspctl"}

void   *Ensstat3(void *argument);
#define Ensstat3Operators                 {"ensrkhistspace", "ensrkhisttime", "ensroc"}

void   *Ensval(void *argument);
#define EnsvalOperators                   {"enscrps", "ensbrs"}

void   *Eofcoeff(void *argument);
#define EofcoeffOperators                 {"eofcoeff"}

void   *Eofcoeff3d(void *argument);
#define Eofcoeff3dOperators               {"eofcoeff3d"}

void   *EOFs(void *argument);
#define EOFsOperators                     {"eof", "eofspatial", "eoftime"}

void   *EOF3d(void *argument);
#define EOF3dOperators                    {"eof3d","eof3dspatial","eof3dtime"}

void   *EstFreq(void *argument);
#define EstFreqOperators                  {"estfreq"}

void   *Expr(void *argument);
#define ExprOperators                     {"expr", "exprf", "aexpr", "aexprf"}

void   *FC(void *argument);
#define FCOperators                       {"fc2sp", "sp2fc", "fc2gp", "gp2fc", "fourier2grid", "grid2fourier"}

void   *Filedes(void *argument);
#define FiledesOperators                  {"filedes", "griddes", "griddes2", "zaxisdes", "vct", "vct2", "codetab", "vlist", "partab", "partab2", "spartab"}

void   *Fillmiss(void *argument);
#define FillmissOperators                 {"fillmiss", "fillmiss2"}

void   *Filter(void *argument);
#define FilterOperators                   {"bandpass", "highpass", "lowpass"}

void   *Fldrms(void *argument);
#define FldrmsOperators                   {"fldrms"}

void   *Fldstat(void *argument);
#define FldstatOperators                  {"fldrange", "fldmin", "fldmax", "fldsum", "fldint", "fldmean", "fldavg", "fldstd", "fldstd1", "fldvar", "fldvar1", "fldskew", "fldkurt", "fldmedian", "fldcount", "fldpctl"}

void   *Fldstat2(void *argument);
#define FldcorOperators                   {"fldcor"}
#define FldcovarOperators                 {"fldcovar"}

void   *Fourier(void *argument);
#define FourierOperators                  {"fourier"}

void   *Gengrid(void *argument);
#define GengridOperators                  {"gengrid"}

void   *Gradsdes(void *argument);
#define GradsdesOperators                 {"gradsdes", "dumpmap"}

void   *Gridboxstat(void *argument);
#define GridboxstatOperators              {"gridboxrange", "gridboxmin", "gridboxmax", "gridboxsum", "gridboxmean", "gridboxavg", "gridboxstd", "gridboxstd1", "gridboxvar", "gridboxvar1", "gridboxskew", "gridboxkurt", "gridboxmedian"}

void   *Gridcell(void *argument);
#define GridcellOperators                 {"gridarea", "gridweights", "gridmask", "griddx", "griddy", "gridcellidx"}

void   *Gridsearch(void *argument);
#define GridsearchOperators               {"testpointsearch", "testcellsearch"}

void   *Harmonic(void *argument);
#define HarmonicOperators                 {"harmonic"}

void   *Healpix(void *argument);
#define HealpixOperators                  {"hpupgrade", "hpdegrade"}

void   *Histogram(void *argument);
#define HistogramOperators                {"histcount", "histsum", "histmean", "histfreq"}

void   *Importamsr(void *argument);
#define ImportamsrOperators               {"import_amsr"}

void   *Importbinary(void *argument);
#define ImportbinaryOperators             {"import_binary"}

void   *Importcmsaf(void *argument);
#define ImportcmsafOperators              {"import_cmsaf"}

void   *Importobs(void *argument);
#define ImportobsOperators                {"import_obs"}

void   *Importfv3grid(void *argument);
#define Importfv3gridOperators            {"import_fv3grid"}

void   *Info(void *argument);
#define InfoOperators                     {"info", "infop", "infon", "infoc", "vinfon", "xinfon", "map"}

void   *Input(void *argument);
#define InputOperators                    {"input", "inputsrv", "inputext"}

void   *Intgrid(void *argument);
#define IntgridOperators                  {"intgridbil", "intgriddis", "intgridnn", "interpolate", "boxavg", "thinout"}

void   *Intgridtraj(void *argument);
#define IntgridtrajOperators              {"intgridtraj"}

void   *Intlevel(void *argument);
#define IntlevelOperators                 {"intlevel", "intlevelx"}

void   *Intlevel3d(void *argument);
#define Intlevel3dOperators               {"intlevel3d", "intlevelx3d"}

void   *Inttime(void *argument);
#define InttimeOperators                  {"inttime"}

void   *Intntime(void *argument);
#define IntntimeOperators                 {"intntime"}

void   *Intyear(void *argument);
#define IntyearOperators                  {"intyear"}

void   *Invert(void *argument);
#define InvertOperators                   {"invertlat", "invertlon", "invertlatdes", "invertlondes", "invertlatdata", "invertlondata"}

void   *Invertlev(void *argument);
#define InvertlevOperators                {"invertlev"}

void   *Lic(void *argument);
#define LicOperators                      {"lic"}

void   *Longinfo(void *argument);
#define LonginfoOperators                 {"linfo"}

void   *MapReduce(void *argument);
#define MapReduceOperators                {"reducegrid"}

void   *Maskbox(void *argument);
#define MaskboxOperators                  {"masklonlatbox", "maskindexbox"}
#define MaskregionOperators               {"maskregion"}

void   *Mastrfu(void *argument);
#define MastrfuOperators                  {"mastrfu"}

void   *Math(void *argument);
#define MathOperators                     {"abs", "int", "nint", "sqr", "sqrt", "exp", "ln", "log10", "sin", "cos", "tan", "asin", "acos", "atan", "pow", "rand", "reci", "not", "conj", "re", "im", "arg"}

void   *Merge(void *argument);
#define MergeOperators                    {"merge"}

void   *Mergegrid(void *argument);
#define MergegridOperators                {"mergegrid"}

void   *Mergetime(void *argument);
#define MergetimeOperators                {"mergetime"}

void   *Merstat(void *argument);
#define MerstatOperators                  {"merrange", "mermin", "mermax", "mersum", "mermean", "meravg", "merstd", "merstd1", "mervar", "mervar1", "merskew", "merkurt", "mermedian", "merpctl"}

void   *Monarith(void *argument);
#define MonarithOperators                 {"monadd", "monsub", "monmul", "mondiv"}

void   *Mrotuv(void *argument);
#define MrotuvOperators                   {"mrotuv"}

void   *Mrotuvb(void *argument);
#define MrotuvbOperators                  {"mrotuvb"}

void   *NCL_wind(void *argument);
#define NCL_windOperators                 {"uv2dv_cfd", "uv2vr_cfd"}

void   *Ninfo(void *argument);
#define NinfoOperators                    {"nyear", "nmon", "ndate", "ntime", "ncode", "npar", "nlevel", "ngridpoints", "ngrids"}

void   *Nmldump(void *argument);
#define NmldumpOperators                  {"nmldump", "kvldump"}

void   *Output(void *argument);
#define OutputOperators                   {"output", "outputint", "outputsrv", "outputext", "outputf", "outputts", "outputfld", "outputarr", "outputxyz"}
#define OutputtabOperators                {"outputtab"}
 
void   *Outputgmt(void *argument);
#define OutputgmtOperators                {"gmtxyz", "gmtcells", "outputcenter2", "outputcentercpt", "outputboundscpt", "outputvector", "outputtri", "outputvrml", "outputkml"}

void   *Pack(void *argument);
#define PackOperators                     {"pack"}

void   *Pardup(void *argument);
#define PardupOperators                   {"pardup", "parmul"}

void   *Pinfo(void *argument);
#define PinfoOperators                    {"pinfo", "pinfov"}

void   *Query(void *argument);
#define QueryOperators                    {"query"}

void   *Pressure(void *argument);
#define PressureOperators                 {"pressure_fl", "pressure_hl", "deltap"}

void   *Recttocomplex(void *argument);
#define RecttocomplexOperators            {"recttocomplex"}

void   *Regres(void *argument);
#define RegresOperators                   {"regres"}

void   *Remap(void *argument);
#define RemapOperators                    {"remap"}
#define RemapbilOperators                 {"remapbil", "genbil"}
#define RemapbicOperators                 {"remapbic", "genbic"}
#define RemapnnOperators                  {"remapnn", "gennn"}
#define RemapdisOperators                 {"remapdis", "gendis"}
#define RemapconOperators                 {"remapcon", "gencon"}
#define Remapycon2Operators               {"remapycon2test", "genycon2test"}
#define RemapsconOperators                {"remapscon", "genscon"}
#define Remapcon2Operators                {"remapcon2", "gencon2"}
#define RemaplafOperators                 {"remaplaf", "genlaf"}
#define RemapavgOperators                 {"remapavgtest"}

void   *Remapweights(void *argument);
#define RemapweightsOperators             {"verifyweights", "writeremapscrip"}

void   *Remapeta(void *argument);
#define RemapetaOperators                 {"remapeta", "remapeta_s", "remapeta_z"}

void   *Remapstat(void *argument);
#define RemapstatOperators                {"remaprange", "remapmin", "remapmax", "remapsum", "remapmean", "remapavg", "remapstd", "remapstd1", "remapvar", "remapvar1", "remapskew", "remapkurt", "remapmedian"}

void   *Replace(void *argument);
#define ReplaceOperators                  {"replace"}

void   *Replacevalues(void *argument);
#define ReplacevaluesOperators            {"setvals", "setrtoc", "setrtoc2"}

void   *Rotuv(void *argument);
#define RotuvOperators                    {"rotuvb"}

void   *Rhopot(void *argument);
#define RhopotOperators                   {"rhopot"}

void   *Runpctl(void *argument);
#define RunpctlOperators                  {"runpctl"}

void   *Runstat(void *argument);
#define RunstatOperators                  {"runrange", "runmin", "runmax", "runsum", "runmean", "runavg", "runstd", "runstd1", "runvar", "runvar1"}

void   *Samplegridicon(void *argument);
#define SamplegridiconOperators           {"samplegridicon"}

void   *Seascount(void *argument);
#define SeascountOperators                {"seascount"}

void   *Seaspctl(void *argument);
#define SeaspctlOperators                 {"seaspctl"}

void   *Seasstat(void *argument);
#define SeasstatOperators                 {"seasrange", "seasmin", "seasmax", "seassum", "seasmean", "seasavg", "seasstd", "seasstd1", "seasvar", "seasvar1"}

void   *Seasmonstat(void *argument);
#define SeasmonstatOperators              {"seasmonmean", "seasmonavg"}

void   *Selbox(void *argument);
#define SelboxOperators                   {"sellonlatbox", "selindexbox"}

void   *Selgridcell(void *argument);
#define SelgridcellOperators              {"selgridcell", "delgridcell"}

void   *Select(void *argument);
#define SelectOperators                   {"select", "delete"}

void   *Selvar(void *argument);
#define SelvarOperators                   {"selparam", "selcode", "selname", "selstdname", "sellevel", "sellevidx", "selgrid", "selzaxis", "selzaxisname", "seltabnum", "delparam", "delcode", "delname", "selltype"}

void   *Seloperator(void *argument);
#define SeloperatorOperators              {"seloperator"}

void   *Selrec(void *argument);
#define SelrecOperators                   {"selrec"}

void   *Selregion(void *argument);
#define SelregionOperators                {"selregion", "selcircle"}

void   *Selsurface(void *argument);
#define SelsurfaceOperators               {"isosurface", "bottomvalue", "topvalue"}

void   *Seltime(void *argument);
#define SeltimeOperators                  {"seltimestep", "selyear", "selseason", "selmonth", "selday", "selhour", "seldate", "seltime", "selsmon"}

void   *Selyearidx(void *argument);
#define SelyearidxOperators               {"selyearidx", "seltimeidx"}

void   *Set(void *argument);
#define SetOperators                      {"setcode", "setparam", "setname", "setunit", "setlevel", "setltype", "settabnum", "setmaxsteps"}

void   *Setattribute(void *argument);
#define SetattributeOperators             {"setattribute"}

void   *Setbox(void *argument);
#define SetboxOperators                   {"setclonlatbox", "setcindexbox"}

void   *Setgrid(void *argument);
#define SetgridOperators                  {"setgrid", "setgridtype", "setgridarea", "setgridmask", "unsetgridmask", "setgridnumber", "setgriduri", "usegridnumber", "setprojparams"}

void   *Setgridcell(void *argument);
#define SetgridcellOperators              {"setgridcell"}

void   *Sethalo(void *argument);
#define SethaloOperators                  {"sethalo", "tpnhalo"}

void   *Setmiss(void *argument);
#define SetmisstonnOperators              {"setmisstonn", "setmisstodis"}

#define SetmissOperators                  {"setmissval", "setctomiss", "setmisstoc", "setrtomiss", "setvrange"}
void   *Setpartab(void *argument);

#define SetpartabOperators                {"setpartabc", "setpartabp", "setpartabn"}
#define SetcodetabOperators               {"setcodetab"}

void   *Setrcaname(void *argument);
#define SetrcanameOperators               {"setrcaname"}

void   *Settime(void *argument);
#define SettimeOperators                  {"setyear", "setmon", "setday", "setdate", "settime", "settunits", "settaxis", "settbounds", "setreftime", "setcalendar", "shifttime"}

void   *Setzaxis(void *argument);
#define SetzaxisOperators                 {"setzaxis", "genlevelbounds"}

void   *Shiftxy(void *argument);
#define ShiftxyOperators                  {"shiftx", "shifty"}

void   *Showinfo(void *argument);
#define ShowinfoOperators                 {"showyear", "showmon", "showdate", "showtime", "showtimestamp", "showcode", "showunit", "showparam", "showname", "showstdname", "showlevel", "showltype", "showformat", "showgrid", "showatts", "showattsglob"}

void   *Showattribute(void *argument);
#define ShowattributeOperators            {"showattribute", "showattsvar"}

void   *Sinfo(void *argument);
#define SinfoOperators                    {"sinfo", "sinfop", "sinfon", "sinfoc", "seinfo", "seinfop", "seinfon", "seinfoc"}
#define XSinfoOperators                   {"xsinfo", "xsinfop", "xsinfon", "xsinfoc"}

void   *Smooth(void *argument);
#define SmoothOperators                   {"smooth", "smooth9"}

void   *Sort(void *argument);
#define SortOperators                     {"sortcode", "sortparam", "sortname", "sortlevel"}

void   *Sorttimestamp(void *argument);
#define SorttimestampOperators            {"sorttimestamp", "sorttaxis"}

void   *Specinfo(void *argument);
#define SpecinfoOperators                 {"specinfo"}

void   *Spectral(void *argument);
#define SpectralOperators                 {"gp2sp", "gp2spl", "sp2gp", "sp2gpl"}
#define SpecconvOperators                 {"sp2sp", "spcut"}

void   *Spectrum(void *argument);
#define SpectrumOperators                 {"spectrum"}

void   *Split(void *argument);
#define SplitOperators                    {"splitcode", "splitparam", "splitname", "splitlevel", "splitgrid", "splitzaxis", "splittabnum"}

void   *Splitdate(void *argument);
#define SplitdateOperators                {"splitdate"}

void   *Splitrec(void *argument);
#define SplitrecOperators                 {"splitrec"}

void   *Splitsel(void *argument);
#define SplitselOperators                 {"splitsel"}

void   *Splittime(void *argument);
#define SplittimeOperators                {"splithour", "splitday", "splitmon", "splitseas"}

void   *Splityear(void *argument);
#define SplityearOperators                {"splityear", "splityearmon"}

void   *Tee(void *argument);
#define TeeOperators                      {"tee"}

void   *Template1(void *argument);
#define Template1Operators                {"template1"}

void   *Template2(void *argument);
#define Template2Operators                {"template2"}

void   *Test(void *argument);
#define TestOperators                     {"test"}

void   *Test2(void *argument);
#define Test2Operators                    {"test2"}

void   *Testdata(void *argument);
#define TestdataOperators                 {"testdata"}

void   *Tests(void *argument);
#define TestsOperators                    {"normal", "studentt", "chisquare", "beta", "fisher"}

void   *Timfill(void *argument);
#define TimfillOperators                  {"setmisstonntime"}

void   *Timsort(void *argument);
#define TimsortOperators                  {"timsort"}

void   *Timcount(void *argument);
#define TimcountOperators                 {"timcount"}
#define YearcountOperators                {"yearcount"}
#define MoncountOperators                 {"moncount"}
#define DaycountOperators                 {"daycount"}
#define HourcountOperators                {"hourcount"}

void   *Timcumsum(void *argument);
#define TimcumsumOperators                {"timcumsum"}

void   *Timpctl(void *argument);
#define TimpctlOperators                  {"timpctl"}
#define YearpctlOperators                 {"yearpctl"}
#define MonpctlOperators                  {"monpctl"}
#define DaypctlOperators                  {"daypctl"}
#define HourpctlOperators                 {"hourpctl"}

void   *Timselpctl(void *argument);
#define TimselpctlOperators               {"timselpctl"}

void   *Timselstat(void *argument);
#define TimselstatOperators               {"timselrange", "timselmin", "timselmax", "timselsum", "timselmean", "timselavg", "timselvar", "timselvar1", "timselstd", "timselstd1"}

void   *XTimstat(void *argument);
#define XTimstatOperators                 {"xtimmin", "xtimmax", "xtimsum", "xtimmean", "xtimavg", "xtimvar", "xtimvar1", "xtimstd", "xtimstd1", "xyearmin", "xyearmax", "xyearsum","xyearmean", "xyearavg", "xyearvar", "xyearvar1", "xyearstd", "xyearstd1", "xmonmin", "xmonmax", "xmonsum", "xmonmean", "xmonavg", "xmonvar", "xmonvar1", "xmonstd", "xmonstd1"}

void   *Timstat(void *argument);
#define TimstatOperators                  {"timrange", "timmin", "timmax", "timsum", "timmean", "timavg", "timvar", "timvar1", "timstd", "timstd1", "timminidx", "timmaxidx"}
#define YearstatOperators                 {"yearrange", "yearmin", "yearmax", "yearsum", "yearmean", "yearavg", "yearvar", "yearvar1", "yearstd", "yearstd1", "yearminidx", "yearmaxidx"}
#define MonstatOperators                  {"monrange", "monmin", "monmax", "monsum", "monmean", "monavg", "monvar", "monvar1", "monstd", "monstd1"}
#define DaystatOperators                  {"dayrange", "daymin", "daymax", "daysum", "daymean", "dayavg", "dayvar", "dayvar1", "daystd", "daystd1"}
#define HourstatOperators                 {"hourrange", "hourmin", "hourmax", "hoursum", "hourmean", "houravg", "hourvar", "hourvar1", "hourstd", "hourstd1"}

void   *Timstat2(void *argument);
#define TimcorOperators                   {"timcor"}
#define TimcovarOperators                 {"timcovar"}
#define TimrmsdOperators                  {"timrmsd"}

void   *Timstat3(void *argument);
#define Timstat3Operators                 {"meandiff2test", "varquot2test"}

void   *Tinfo(void *argument);
#define TinfoOperators                    {"tinfo"}

void   *Tocomplex(void *argument);
#define TocomplexOperators                {"retocomplex", "imtocomplex"}

void   *Transpose(void *argument);
#define TransposeOperators                {"transxy"}

void   *Trend(void *argument);
#define TrendOperators                    {"trend"}

void   *Trendarith(void *argument);
#define TrendarithOperators               {"addtrend", "subtrend"}

void   *Tstepcount(void *argument);
#define TstepcountOperators               {"tstepcount"}

void   *Unpack(void *argument);
#define UnpackOperators                   {"unpack"}

void   *Vargen(void *argument);
#define VargenOperators                   {"random", "const", "sincos", "coshill", "testfield", "seq", "topo", "temp", "mask", "stdatm"}

void   *Varrms(void *argument);
#define VarrmsOperators                   {"varrms"}

void   *Varsstat(void *argument);
#define VarsstatOperators                 {"varsrange", "varsmin", "varsmax", "varssum", "varsmean", "varsavg", "varsstd", "varsstd1", "varsvar", "varsvar1"}

void   *Vertintap(void *argument);
#define VertintapOperators                {"ap2pl", "ap2plx", "ap2hl", "ap2hlx"}

void   *Vertintgh(void *argument);
#define VertintghOperators                {"gh2hl", "gh2hlx"}

void   *Vertintml(void *argument);
#define VertintmlOperators                {"ml2pl", "ml2hl", "ml2plx", "ml2hlx", "ml2height", "ml2heightx"}

void   *Vertintzs(void *argument);
#define VertintzsOperators                {"zs2zl", "zs2zlx"}

void   *Vertstat(void *argument);
#define VertstatOperators                 {"vertrange", "vertmin", "vertmax", "vertsum", "vertint", "vertmean", "vertavg", "vertstd", "vertstd1", "vertvar", "vertvar1"}

void   *Vertcum(void *argument);
#define VertcumOperators                  {"vertcum", "vertcumhl"}

void   *Vertwind(void *argument);
#define VertwindOperators                 {"vertwind"}

void   *Verifygrid(void *argument);
#define VerifygridOperators               {"verifygrid"}

void   *Wind(void *argument);
#define WindOperators                     {"uv2dv", "uv2dvl", "dv2uv", "dv2uvl"}

#define Wind2Operators                    {"dv2ps"}
void   *Writegrid(void *argument);

#define WritegridOperators                {"writegrid"}
void   *Writerandom(void *argument);

#define WriterandomOperators              {"writerandom"}
void   *YAR(void *argument);

void   *Yeararith(void *argument);
#define YeararithOperators                {"yearadd", "yearsub", "yearmul", "yeardiv"}

void   *Yearmonstat(void *argument);
#define YearmonstatOperators              {"yearmonmean", "yearmonavg"}

void   *Ydayarith(void *argument);
#define YdayarithOperators                {"ydayadd", "ydaysub", "ydaymul", "ydaydiv"}

void   *Ydaypctl(void *argument);
#define YdaypctlOperators                 {"ydaypctl"}

void   *Ydaystat(void *argument);
#define YdaystatOperators                 {"ydayrange", "ydaymin", "ydaymax", "ydaysum", "ydaymean", "ydayavg", "ydaystd", "ydaystd1", "ydayvar", "ydayvar1"}

void   *Ydrunpctl(void *argument);
#define YdrunpctlOperators                {"ydrunpctl"}

void   *Ydrunstat(void *argument);
#define YdrunstatOperators                {"ydrunmin", "ydrunmax", "ydrunsum", "ydrunmean", "ydrunavg", "ydrunstd", "ydrunstd1", "ydrunvar", "ydrunvar1"}

void   *Yhourarith(void *argument);
#define YhourarithOperators               {"yhouradd", "yhoursub", "yhourmul", "yhourdiv"}

void   *Yhourstat(void *argument);
#define YhourstatOperators                {"yhourrange", "yhourmin", "yhourmax", "yhoursum", "yhourmean", "yhouravg", "yhourstd", "yhourstd1", "yhourvar", "yhourvar1"}
#define DhourstatOperators                {"dhourrange", "dhourmin", "dhourmax", "dhoursum", "dhourmean", "dhouravg", "dhourstd", "dhourstd1", "dhourvar", "dhourvar1"}

void   *Ymonarith(void *argument);
#define YmonarithOperators                {"ymonadd", "ymonsub", "ymonmul", "ymondiv"}
#define YseasarithOperators               {"yseasadd", "yseassub", "yseasmul", "yseasdiv"}

void   *Ymoncomp(void *argument);
#define YmoncompOperators                 {"ymoneq", "ymonne", "ymonle", "ymonlt", "ymonge", "ymongt"}
#define YseascompOperators                {"yseaseq", "yseasne", "yseasle", "yseaslt", "yseasge", "yseasgt"}

void   *Ymonpctl(void *argument);
#define YmonpctlOperators                 {"ymonpctl"}

void   *Ymonstat(void *argument);
#define YmonstatOperators                 {"ymonrange", "ymonmin", "ymonmax", "ymonsum", "ymonmean", "ymonavg", "ymonstd", "ymonstd1", "ymonvar", "ymonvar1"}

void   *Yseaspctl(void *argument);
#define YseaspctlOperators                {"yseaspctl"}

void   *Yseasstat(void *argument);
#define YseasstatOperators                {"yseasrange", "yseasmin", "yseasmax", "yseassum", "yseasmean", "yseasavg", "yseasstd", "yseasstd1", "yseasvar", "yseasvar1"}

void   *Zonstat(void *argument);
#define ZonstatOperators                  {"zonrange", "zonmin", "zonmax", "zonsum", "zonmean", "zonavg", "zonstd", "zonstd1", "zonvar", "zonvar1", "zonskew", "zonkurt", "zonmedian", "zonpctl"}

//------------------------------

void   *EcaCfd(void *argument);
#define EcaCfdOperators                   {"eca_cfd"}
void   *EcaCsu(void *argument);
#define EcaCsuOperators                   {"eca_csu"}
void   *EcaCwdi(void *argument);
#define EcaCwdiOperators                  {"eca_cwdi"}
void   *EcaCwfi(void *argument);
#define EcaCwfiOperators                  {"eca_cwfi", "etccdi_csdi"}
void   *EcaEtr(void *argument);
#define EcaEtrOperators                   {"eca_etr"}
void   *EcaFd(void *argument);
#define EcaFdOperators                    {"eca_fd", "etccdi_fd"}
void   *EcaGsl(void *argument);
#define EcaGslOperators                   {"eca_gsl", "etccdi_gsl"}
void   *EcaHd(void *argument);
#define EcaHdOperators                    {"eca_hd", "etccdi_hd"}
void   *EcaHwdi(void *argument);
#define EcaHwdiOperators                  {"eca_hwdi"}
void   *EcaHwfi(void *argument);
#define EcaHwfiOperators                  {"eca_hwfi", "etccdi_wsdi"}
void   *EcaId(void *argument);
#define EcaIdOperators                    {"eca_id", "etccdi_id"}
void   *EcaSu(void *argument);
#define EcaSuOperators                    {"eca_su", "etccdi_su"}
void   *EcaTr(void *argument);
#define EcaTrOperators                    {"eca_tr", "etccdi_tr"}
void   *EcaTg10p(void *argument);
#define EcaTg10pOperators                 {"eca_tg10p"}
void   *EcaTg90p(void *argument);
#define EcaTg90pOperators                 {"eca_tg90p"}
void   *EcaTn10p(void *argument);
#define EcaTn10pOperators                 {"eca_tn10p"}
void   *EcaTn90p(void *argument);
#define EcaTn90pOperators                 {"eca_tn90p"}
void   *EcaTx10p(void *argument);
#define EcaTx10pOperators                 {"eca_tx10p"}
void   *EcaTx90p(void *argument);
#define EcaTx90pOperators                 {"eca_tx90p"}
void   *EcaCdd(void *argument);
#define EcaCddOperators                   {"eca_cdd", "etccdi_cdd"}
void   *EcaCwd(void *argument);
#define EcaCwdOperators                   {"eca_cwd", "etccdi_cwd"}
void   *EcaRr1(void *argument);
#define EcaRr1Operators                   {"eca_rr1"}
void   *EcaPd(void *argument);
#define EcaPdOperators                    {"eca_pd", "eca_r10mm", "eca_r20mm", "etccdi_r1mm", "etccdi_r10mm", "etccdi_r20mm"}
void   *EcaR75p(void *argument);
#define EcaR75pOperators                  {"eca_r75p"}
void   *EcaR75ptot(void *argument);
#define EcaR75ptotOperators               {"eca_r75ptot"}
void   *EcaR90p(void *argument);
#define EcaR90pOperators                  {"eca_r90p"}
void   *EcaR90ptot(void *argument);
#define EcaR90ptotOperators               {"eca_r90ptot"}
void   *EcaR95p(void *argument);
#define EcaR95pOperators                  {"eca_r95p"}
void   *EcaR95ptot(void *argument);
#define EcaR95ptotOperators               {"eca_r95ptot"}
void   *EcaR99p(void *argument);
#define EcaR99pOperators                  {"eca_r99p"}
void   *EcaR99ptot(void *argument);
#define EcaR99ptotOperators               {"eca_r99ptot"}
void   *EcaRx1day(void *argument);
#define EcaRx1dayOperators                {"eca_rx1day", "etccdi_rx1day", "etccdi_rx1daymon"}
void   *EcaRx5day(void *argument);
#define EcaRx5dayOperators                {"eca_rx5day", "etccdi_rx5day", "etccdi_rx5daymon"}
void   *EcaSdii(void *argument);
#define EcaSdiiOperators                  {"eca_sdii", "etccdi_sdii"}
void   *EcaEtccdi(void *argument);
#define EcaEtccdiOperators                {"etccdi_tx90p","etccdi_tx10p","etccdi_tn90p","etccdi_tn10p","etccdi_r95p","etccdi_r99p","etccdi"}

void   *Fdns(void *argument);
#define FdnsOperators                     {"fdns"}

void   *Strwin(void *argument);
#define StrwinOperators                   {"strwin"}
void   *Strbre(void *argument);
#define StrbreOperators                   {"strbre"}
void   *Strgal(void *argument);
#define StrgalOperators                   {"strgal"}
void   *Hurr(void *argument);
#define HurrOperators                     {"hurr"}

// void   *Hi(void *argument);
//#define HiOperators                     {"hi"}
void   *Wct(void *argument);
#define WctOperators                      {"wct"}

void   *Magplot(void *argument);
#define MagplotOperators                  {"contour", "shaded", "grfill"}
void   *Magvector(void *argument);
#define MagvectorOperators                {"vector"}
void   *Maggraph(void *argument);
#define MaggraphOperators                 {"graph"}

//------------------------------
// HIRLAM_EXTENSIONS
void   *Selmulti(void *argument); // "selmulti", "delmulti"
#define SelmultiOperators                 {"selmulti", "delmulti", "changemulti"}
void   *WindTrans(void *argument); // "uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"
#define WindTransOperators                {"uvDestag", "rotuvN", "rotuvNorth", "projuvLatLon"}
void   *Samplegrid(void *argument); // "samplegrid", "subgrid"
#define SamplegridOperators               {"samplegrid", "subgrid"}

// clang-format on
/* \endcond */
#endif
