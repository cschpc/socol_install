#ifndef REMAP_VARS_H
#define REMAP_VARS_H

#include <cstdio>  // size_t
#include <array>

#include "field.h"

struct RemapGradients
{
  Varray<double> grad_lat;
  Varray<double> grad_lon;
  Varray<double> grad_latlon;

  void
  init(size_t size)
  {
    grad_lat.resize(size);
    grad_lon.resize(size);
    grad_latlon.resize(size);
  }

  explicit RemapGradients(size_t size) { init(size); }
  RemapGradients() {}
};

enum class RemapMethod
{
  UNDEF,
  BILINEAR,
  BICUBIC,
  DISTWGT,
  CONSERV,
  CONSERV_SCRIP
};

enum class SubmapType
{
  NONE,
  LAF,
  SUM,
  AVG
};

enum class NormOpt
{
  NONE,
  DESTAREA,
  FRACAREA
};

struct RemapLink
{
  bool option;
  size_t maxLinks;
  size_t numBlocks;
  Varray<size_t> numLinks;
  Varray2D<size_t> src_add;
  Varray2D<size_t> tgt_add;
  Varray2D<size_t> w_index;
};

struct RemapSwitches
{
  RemapMethod mapType{ RemapMethod::UNDEF };
  SubmapType submapType{ SubmapType::NONE };
  int numNeighbors{ 0 };
  int remapOrder{ 0 };
};

struct RemapVars
{
  bool sort_add;
  bool pinit;           // true: if the pointers are initialized
  RemapMethod mapType;  // identifier for remapping method
  NormOpt normOpt;      // option for normalization (conserv only)
  long linksPerValue;
  size_t maxLinks;          // current size of link arrays
  size_t numLinks;          // actual number of links for remapping
  size_t num_wts;           // num of weights used in remapping
  size_t resize_increment;  // default amount to increase array size

  Varray<size_t> srcCellIndices;  // source grid indices for each link
  Varray<size_t> tgtCellIndices;  // target grid indices for each link
  Varray<double> wts;             // map weights for each link [maxLinks*num_wts]

  RemapLink links;
};

void remap_field(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1,
                 RemapGradients &gradients);
void remap_laf(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1);
void remap_avg(Field &field2, double missval, size_t gridsize2, const RemapVars &rv, const Field &field1);
void remap_vars_init(RemapMethod mapType, int remapOrder, RemapVars &rv);
void remap_vars_ensure_size(RemapVars &rv, size_t size);
void remap_vars_resize(RemapVars &rv, size_t size);
void remap_vars_reorder(RemapVars &rv);
void remap_vars_free(RemapVars &rv);
void remap_vars_check_weights(const RemapVars &rv);

#endif /* REMAP_VARS_H */
