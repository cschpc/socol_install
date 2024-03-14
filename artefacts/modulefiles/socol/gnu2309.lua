help(
[[
SOCOL testing version SOCOLv4.0-testing. 

Exports env vars SOCOL_ROOT, PATH and LD_LIBRARY_PATH.

GNU compiler collection version.
]])

local GERACLIS_ROOT="/projappl/project_462000470"
local version = "SOCOLv4.0-testing"
local base = GERACLIS_ROOT

prepend_path("PATH", pathJoin(base, "socol/gnu/usr/bin"))
prepend_path("PATH", pathJoin(base, "socol/gnu/jobscript_toolkit/bin"))
setenv("SOCOL_ROOT", pathJoin(base, "socol/gnu", version))
prepend_path("LD_LIBRARY_PATH", pathJoin(base, "socol/gnu/usr/lib"))

depends_on("LUMI/23.09")
depends_on("partition/C")
depends_on("PrgEnv-gnu/8.4.0")
    
depends_on("cray-mpich/8.1.27")
depends_on("cray-hdf5/1.12.2.1")
depends_on("cray-netcdf/4.9.0.3")
