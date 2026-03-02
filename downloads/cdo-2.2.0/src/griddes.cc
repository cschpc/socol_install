/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <climits>

#include <cdi.h>

#include "cdi_uuid.h"
#include "cdo_options.h"
#include "process_int.h"
#include "griddes.h"
#include "util_wildcards.h"
#include "grid_read_pingo.h"
#include "cdi_lockedIO.h"

int
getoptname(char *optname, const char *optstring, int nopt)
{
  int nerr = 1;
  auto pname = optstring;
  auto pend = optstring;

  for (int i = 0; i < nopt; ++i)
    {
      pend = strchr(pname, ',');
      if (pend == nullptr) break;
      pname = pend + 1;
    }

  if (pend)
    {
      pend = strchr(pname, ',');
      size_t namelen = (pend == nullptr) ? strlen(pname) : pend - pname;
      memcpy(optname, pname, namelen);
      optname[namelen] = '\0';
      nerr = 0;
    }

  return nerr;
}

static void
genXboundsRegular(GridDesciption &grid)
{
  grid.nvertex = 2;
  grid.xbounds.resize(grid.xsize * grid.nvertex);
  for (size_t i = 0; i < grid.xsize - 1; ++i)
    {
      auto value = 0.5 * (grid.xvals[i] + grid.xvals[i + 1]);
      grid.xbounds[2 * i + 1] = value;
      grid.xbounds[2 * (i + 1)] = value;
    }

  if (grid.xsize == 1)
    {
      grid.xbounds[0] = 0.0;
      grid.xbounds[1] = 360.0;
    }
  else
    {
      grid.xbounds[0] = 2 * grid.xvals[0] - grid.xbounds[1];
      grid.xbounds[2 * grid.xsize - 1] = 2 * grid.xvals[grid.xsize - 1] - grid.xbounds[2 * (grid.xsize - 1)];
    }
}

static void
genYboundsRegular(GridDesciption &grid)
{
  auto lrev = (grid.yvals[0] > grid.yvals[grid.ysize - 1]);
  grid.nvertex = 2;
  grid.ybounds.resize(grid.ysize * grid.nvertex);
  for (size_t i = 0; i < grid.ysize - 1; ++i)
    {
      auto value = 0.5 * (grid.yvals[i] + grid.yvals[i + 1]);
      if (lrev)
        {
          grid.ybounds[2 * i] = value;
          grid.ybounds[2 * (i + 1) + 1] = value;
        }
      else
        {
          grid.ybounds[2 * i + 1] = value;
          grid.ybounds[2 * (i + 1)] = value;
        }
    }

  if (lrev)
    {
      grid.ybounds[1] = 90;
      grid.ybounds[2 * grid.ysize - 2] = -90;
    }
  else
    {
      grid.ybounds[0] = -90;
      grid.ybounds[2 * grid.ysize - 1] = 90;
    }
}

int
grid_define(GridDesciption &grid)
{
  int gridID = CDI_UNDEFID;

  if (grid.type == GRID_PROJECTION && grid.projection == "healpix")
    {
      if (grid.healpixNside <= 0) cdo_abort("healpix nside undefined!");
      gridID = gridCreate(grid.type, grid.size);
      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_DIMNAME, "cells");
      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_VARNAME, "healpix");
      cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_GRIDMAP_NAME, grid.projection.c_str());
      cdiDefAttTxt(gridID, CDI_GLOBAL, "grid_mapping_name", (int) grid.projection.size(), grid.projection.c_str());
      cdiDefAttInt(gridID, CDI_GLOBAL, "healpix_nside", CDI_DATATYPE_INT32, 1, &grid.healpixNside);
      cdiDefAttTxt(gridID, CDI_GLOBAL, "healpix_order", (int) grid.healpixOrder.size(), grid.healpixOrder.c_str());
      return gridID;
    }

  switch (grid.type)
    {
    case GRID_GENERIC:
    case GRID_LONLAT:
    case GRID_GAUSSIAN:
    case GRID_GAUSSIAN_REDUCED:
    case GRID_PROJECTION:
    case GRID_CHARXY:
      {
        if (grid.size != 1)
          {
            if (grid.xsize == 0 && grid.type != GRID_PROJECTION && grid.type != GRID_GAUSSIAN_REDUCED)
              cdo_abort("xsize undefined!");
            if (grid.ysize == 0 && grid.type != GRID_PROJECTION) cdo_abort("ysize undefined!");
          }

        if (grid.size == 0) grid.size = grid.xsize * grid.ysize;

        if (grid.type != GRID_PROJECTION && grid.type != GRID_GAUSSIAN_REDUCED && grid.xsize
            && grid.size != grid.xsize * grid.ysize)
          cdo_abort("Inconsistent grid declaration: xsize*ysize!=gridsize (xsize=%zu ysize=%zu gridsize=%zu)", grid.xsize,
                    grid.ysize, grid.size);

        // if ( grid.size < 0 || grid.size > INT_MAX ) cdo_abort("grid size (%ld) out of bounds (0 - %d)!", grid.size, INT_MAX);

        gridID = gridCreate(grid.type, grid.size);

        if (grid.xsize > 0) gridDefXsize(gridID, grid.xsize);
        if (grid.ysize > 0) gridDefYsize(gridID, grid.ysize);
        if (grid.numLPE > 0) gridDefNP(gridID, grid.numLPE);

        if (grid.nvertex) gridDefNvertex(gridID, grid.nvertex);

        auto def_xfirst = IS_NOT_EQUAL(grid.xfirst, undefGridValue);
        auto def_yfirst = IS_NOT_EQUAL(grid.yfirst, undefGridValue);
        auto def_xlast = IS_NOT_EQUAL(grid.xlast, undefGridValue);
        auto def_ylast = IS_NOT_EQUAL(grid.ylast, undefGridValue);
        auto def_xinc = IS_NOT_EQUAL(grid.xinc, undefGridValue);
        auto def_yinc = IS_NOT_EQUAL(grid.yinc, undefGridValue);
        if ((def_xfirst || def_xlast || def_xinc) && grid.xvals.size() == 0)
          {
            auto xfirst = def_xfirst ? grid.xfirst : 0.0;
            auto xlast = def_xlast ? grid.xlast : 0.0;
            auto xinc = def_xinc ? grid.xinc : 0.0;
            grid.xvals.resize(grid.xsize);
            grid_gen_xvals(grid.xsize, xfirst, xlast, xinc, grid.xvals.data());
            if (grid.genBounds && grid.xbounds.size() == 0 /*&& grid.xsize > 1*/) genXboundsRegular(grid);
          }

        if ((def_yfirst || def_ylast || def_yinc) && grid.yvals.size() == 0)
          {
            auto yfirst = def_yfirst ? grid.yfirst : 0.0;
            auto ylast = def_ylast ? grid.ylast : yfirst;
            auto yinc = def_yinc ? grid.yinc : 0.0;
            grid.yvals.resize(grid.ysize);
            if (grid.type == GRID_GAUSSIAN && grid.genBounds)
              {
                grid.ybounds.resize(grid.ysize * 2);
                std::vector<double> lat_bounds(grid.ysize + 1);
                gaussian_latitudes_in_degrees(grid.yvals, lat_bounds, grid.ysize);
                for (size_t i = 0; i < grid.ysize; ++i) grid.ybounds[2 * i] = lat_bounds[i + 1];
                for (size_t i = 0; i < grid.ysize; ++i) grid.ybounds[2 * i + 1] = lat_bounds[i];
              }
            else
              {
                grid_gen_yvals(grid.type, grid.ysize, yfirst, ylast, yinc, grid.yvals.data());
                if (grid.genBounds && grid.ybounds.size() == 0 && grid.ysize > 1) genYboundsRegular(grid);
              }
          }

        if (grid.xcvals)
          {
            // gridDefXCvals(gridID, grid.xcvals);
            for (size_t i = 0; i < grid.xsize; ++i) free(grid.xcvals[i]);
            free(grid.xcvals);
            cdo_warning("CDI function gridDefXCvals() not implemented!");
          }
        if (grid.ycvals)
          {
            // gridDefYCvals(gridID, grid.ycvals);
            for (size_t i = 0; i < grid.ysize; ++i) free(grid.ycvals[i]);
            free(grid.ycvals);
            cdo_warning("CDI function gridDefYCvals() not implemented!");
          }
        if (grid.xvals.size()) gridDefXvals(gridID, grid.xvals.data());
        if (grid.yvals.size()) gridDefYvals(gridID, grid.yvals.data());
        if (grid.xbounds.size()) gridDefXbounds(gridID, grid.xbounds.data());
        if (grid.ybounds.size()) gridDefYbounds(gridID, grid.ybounds.data());
        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());
        if (grid.reducedPoints.size()) gridDefReducedPoints(gridID, grid.ysize, grid.reducedPoints.data());

        break;
      }
    case GRID_CURVILINEAR:
    case GRID_UNSTRUCTURED:
      {
        if (grid.size == 0) grid.size = (grid.type == GRID_CURVILINEAR) ? grid.xsize * grid.ysize : grid.xsize;

        gridID = gridCreate(grid.type, grid.size);

        if (grid.type == GRID_CURVILINEAR)
          {
            if (grid.xsize == 0) cdo_abort("xsize undefined!");
            if (grid.ysize == 0) cdo_abort("ysize undefined!");
            gridDefXsize(gridID, grid.xsize);
            gridDefYsize(gridID, grid.ysize);
          }
        else
          {
            if (grid.nvertex > 0) gridDefNvertex(gridID, grid.nvertex);
            if (grid.number > 0)
              {
                cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDUSED, grid.number);
                if (grid.position >= 0) cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_NUMBEROFGRIDINREFERENCE, grid.position);
              }
            if (grid.path.size()) cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_REFERENCEURI, grid.path.c_str());
          }

        if (grid.xvals.size()) gridDefXvals(gridID, grid.xvals.data());
        if (grid.yvals.size()) gridDefYvals(gridID, grid.yvals.data());
        if (grid.area.size()) gridDefArea(gridID, grid.area.data());
        if (grid.xbounds.size()) gridDefXbounds(gridID, grid.xbounds.data());
        if (grid.ybounds.size()) gridDefYbounds(gridID, grid.ybounds.data());
        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());

        break;
      }
    case GRID_SPECTRAL:
      {
        if (grid.ntr == 0) cdo_abort("truncation undefined!");
        if (grid.size == 0) grid.size = (grid.ntr + 1) * (grid.ntr + 2);

        gridID = gridCreate(grid.type, grid.size);

        gridDefTrunc(gridID, grid.ntr);
        gridDefComplexPacking(gridID, grid.lcomplex);

        break;
      }
    case GRID_GME:
      {
        if (grid.nd == 0) cdo_abort("nd undefined!");
        if (grid.ni == 0) cdo_abort("ni undefined!");
        if (grid.size == 0) cdo_abort("size undefined!");

        gridID = gridCreate(grid.type, grid.size);

        gridDefParamGME(gridID, grid.nd, grid.ni, grid.ni2, grid.ni3);

        if (grid.mask.size()) gridDefMask(gridID, grid.mask.data());

        break;
      }
    default:
      {
        if (grid.type == -1)
          cdo_abort("Undefined grid type!");
        else
          cdo_abort("Unsupported grid type: %s", gridNamePtr(grid.type));

        break;
      }
    }

  if (grid.datatype != CDI_UNDEFID) cdiDefKeyInt(gridID, CDI_GLOBAL, CDI_KEY_DATATYPE, grid.datatype);

  if (!cdiUUIDIsNull(grid.uuid)) cdiDefKeyBytes(gridID, CDI_GLOBAL, CDI_KEY_UUID, grid.uuid, CDI_UUID_SIZE);

  // clang-format off
  if (grid.xname.size())     cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_NAME, grid.xname.c_str());
  if (grid.xlongname.size()) cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_LONGNAME, grid.xlongname.c_str());
  if (grid.xunits.size())    cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_UNITS, grid.xunits.c_str());
  if (grid.yname.size())     cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_NAME, grid.yname.c_str());
  if (grid.ylongname.size()) cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_LONGNAME, grid.ylongname.c_str());
  if (grid.yunits.size())    cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_UNITS, grid.yunits.c_str());
  if (grid.xdimname.size())  cdiDefKeyString(gridID, CDI_XAXIS, CDI_KEY_DIMNAME, grid.xdimname.c_str());
  if (grid.ydimname.size())  cdiDefKeyString(gridID, CDI_YAXIS, CDI_KEY_DIMNAME, grid.ydimname.c_str());
  if (grid.vdimname.size())  cdiDefKeyString(gridID, CDI_GLOBAL, CDI_KEY_VDIMNAME, grid.vdimname.c_str());
  // clang-format on

  return gridID;
}

int
gridFromMPAS(const char *filename)
{
  int mpasGridID = -1;

  open_lock();
  auto streamID = streamOpenRead(filename);
  open_unlock();
  if (streamID >= 0)
    {
      auto vlistID = streamInqVlist(streamID);
      auto nvars = vlistNvars(vlistID);
      int latCellID = -1, lonCellID = -1, areaCellID = -1;
      int latVertexID = -1, lonVertexID = -1, verticesOnCellID = -1;
      for (int varID = 0; varID < nvars; ++varID)
        {
          auto varname = cdo::inq_var_name(vlistID, varID);
          // clang-format off
          if      (varname == "latCell")   latCellID = varID;
          else if (varname == "lonCell")   lonCellID = varID;
          else if (varname == "latVertex") latVertexID = varID;
          else if (varname == "lonVertex") lonVertexID = varID;
          else if (varname == "areaCell")  areaCellID = varID;
          else if (varname == "verticesOnCell") verticesOnCellID = varID;
          // clang-format on
        }

      // printf("%d %d %d %d %d %d\n", latCellID, lonCellID, areaCellID, latVertexID, lonVertexID, verticesOnCellID);
      if (latCellID == -1 || lonCellID == -1 || areaCellID == -1 || latVertexID == -1 || lonVertexID == -1
          || verticesOnCellID == -1)
        {
          cdo_warning("MPAS grid variable missing!");
          return -1;
        }

      auto maxEdges = gridInqXsize(vlistInqVarGrid(vlistID, verticesOnCellID));
      auto ncells = gridInqSize(vlistInqVarGrid(vlistID, latCellID));
      std::vector<double> latCell(ncells), lonCell(ncells), areaCell(ncells);
      auto nvertices = gridInqSize(vlistInqVarGrid(vlistID, latVertexID));
      std::vector<double> latVertex(nvertices), lonVertex(nvertices);
      std::vector<double> verticesOnCell(maxEdges * ncells);
      std::vector<double> xbounds(maxEdges * ncells), ybounds(maxEdges * ncells);

      auto nrecs = streamInqTimestep(streamID, 0);
      for (int recID = 0; recID < nrecs; ++recID)
        {
          size_t nmiss;
          int varID, levelID;
          streamInqRecord(streamID, &varID, &levelID);
          // clang-format off
          if      (varID == latCellID)   streamReadRecord(streamID, latCell.data(), &nmiss);
          else if (varID == lonCellID)   streamReadRecord(streamID, lonCell.data(), &nmiss);
          else if (varID == latVertexID) streamReadRecord(streamID, latVertex.data(), &nmiss);
          else if (varID == lonVertexID) streamReadRecord(streamID, lonVertex.data(), &nmiss);
          else if (varID == areaCellID)  streamReadRecord(streamID, areaCell.data(), &nmiss);
          else if (varID == verticesOnCellID) streamReadRecord(streamID, verticesOnCell.data(), &nmiss);
          // clang-format on
        }
      streamClose(streamID);

      mpasGridID = gridCreate(GRID_UNSTRUCTURED, ncells);
      gridDefXvals(mpasGridID, lonCell.data());
      gridDefYvals(mpasGridID, latCell.data());
      cdiDefKeyString(mpasGridID, CDI_XAXIS, CDI_KEY_UNITS, "radians");
      cdiDefKeyString(mpasGridID, CDI_YAXIS, CDI_KEY_UNITS, "radians");
      gridDefArea(mpasGridID, areaCell.data());

      size_t nv = 0;
      for (size_t i = 0; i < ncells; ++i)
        for (size_t k = 0; k < maxEdges; ++k)
          {
            auto cellIndex = std::lround(verticesOnCell[i * maxEdges + k]);
            if (cellIndex == 0)
              {
                if (k > nv) nv = k;
                break;
              }
          }

      std::vector<double> xc(nv), yc(nv);

      auto withBounds = true;
      for (size_t i = 0; i < ncells; ++i)
        for (size_t k = 0; k < nv; ++k)
          {
            auto cellIndex = std::lround(verticesOnCell[i * maxEdges + k]);
            if (cellIndex < 0 || cellIndex > (long) nvertices) { withBounds = false; }
            else if (cellIndex == 0)
              {
                if (k == (nv - 1))
                  {
                    xc[k] = xc[k - 1];
                    yc[k] = yc[k - 1];
                  }
                else
                  {
                    withBounds = false;
                  }
              }
            else
              {
                xc[k] = lonVertex[cellIndex - 1];
                yc[k] = latVertex[cellIndex - 1];
              }
            xbounds[i * nv + k] = xc[k];
            ybounds[i * nv + k] = yc[k];
          }

      if (withBounds)
        {
          gridDefNvertex(mpasGridID, nv);
          gridDefXbounds(mpasGridID, xbounds.data());
          gridDefYbounds(mpasGridID, ybounds.data());
        }
    }

  return mpasGridID;
}

static int
cdo_define_grid(const char *pgridfile)
{
  int gridID = -1;
  bool isreg = false;
  bool lalloc = false;
  bool lmpas = false;

  auto len = strlen(pgridfile);
  if (len > 5 && strncmp(pgridfile, "mpas:", 5) == 0)
    {
      lmpas = true;
      pgridfile += 5;
      len -= 5;
    }

  char *gridfile = strdup(pgridfile);

  int gridno = 1;
  if (len > 2 && gridfile[len - 2] == ':' && isdigit(gridfile[len - 1]))
    {
      gridno = gridfile[len - 1] - '0';
      gridfile[len - 2] = 0;
    }

  char *filename = expand_filename(gridfile);
  if (filename)
    lalloc = true;
  else
    filename = gridfile;

  auto fileno = open(filename, O_RDONLY);
  if (fileno >= 0)
    {
      struct stat filestat;
      if (fstat(fileno, &filestat) == 0) isreg = S_ISREG(filestat.st_mode);
    }

  if (fileno == -1 || !isreg)
    {
      if (isreg) close(fileno);

      gridID = grid_from_name(gridfile);

      if (gridID == -1) cdo_abort("Open failed on %s!", gridfile);
    }
  else if (lmpas)
    {
      Debug(cdoDebug, "Grid from MPAS file");
      gridID = gridFromMPAS(filename);
      if (gridID == -1) cdo_abort("Unsupported MPAS grid file %s!", filename);
    }
  else
    {
      char buffer[4];
      if (read(fileno, buffer, 4) != 4) cdo_sys_error("Read grid from %s failed!", filename);

      close(fileno);

      if (cdo_cmpstrLenRhs(buffer, "CDF"))
        {
          Debug(cdoDebug, "Grid from NetCDF file");
          gridID = gird_from_nc_file(filename);
        }

      if (gridID == -1)
        {
          if (cdo_cmpstrLenRhs(buffer + 1, "HDF"))
            {
              Debug(cdoDebug, "Grid from HDF5 file");
              gridID = grid_from_h5_file(filename);
            }
        }

      if (gridID == -1)
        {
          if (cdo_cmpstrLenRhs(buffer + 1, "HDF"))
            {
              Debug(cdoDebug, "Grid from NetCDF4 file");
              gridID = gird_from_nc_file(filename);
            }
        }

      if (gridID == -1)
        {
          Debug(cdoDebug, "Grid from CDI file");
          open_lock();
          auto streamID = streamOpenRead(filename);
          open_unlock();
          if (streamID >= 0)
            {
              auto vlistID = streamInqVlist(streamID);
              auto ngrids = vlistNgrids(vlistID);
              if (gridno < 1 || gridno > ngrids) cdo_abort("Grid number %d not available in %s!", gridno, filename);
              gridID = vlistGrid(vlistID, gridno - 1);
              streamClose(streamID);
            }
        }

      if (gridID == -1)
        {
          Debug(cdoDebug, "grid from ASCII file");
          auto gfp = std::fopen(filename, "r");
          gridID = grid_read(gfp, filename);
          std::fclose(gfp);
        }

      if (gridID == -1)
        {
          Debug(cdoDebug, "grid from PINGO file");
          auto gfp = std::fopen(filename, "r");
          gridID = grid_read_pingo(gfp);
          std::fclose(gfp);
        }

      if (gridID == -1) cdo_abort("Invalid grid description file %s!", filename);
    }

  if (lalloc) free(filename);
  free(gridfile);

  return gridID;
}

int
cdo_define_grid(const std::string &gridfile)
{
  return cdo_define_grid(gridfile.c_str());
}

void
cdo_set_grids(const char *gridarg)
{
  char gridfile[4096];
  int nfile = 0;

  while (getoptname(gridfile, gridarg, nfile++) == 0) { (void) cdo_define_grid(gridfile); }
}
