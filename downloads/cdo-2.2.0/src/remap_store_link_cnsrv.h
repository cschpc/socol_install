/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef REMAP_STORE_LINK_CNSRV_H
#define REMAP_STORE_LINK_CNSRV_H

// used for store_link_cnsrv

#define BLK_SIZE 4096
#define BLK_NUM(x) (x / grid_store.blk_size)
#define BLK_IDX(x) (x % grid_store.blk_size)

// Predeclarations
struct RemapVars;

struct GridLayer
{
  long *grid2_link;
  GridLayer *next;
};

struct GridStore
{
  long blk_size;
  long max_size;
  long nblocks;
  long *blksize;
  long *nlayers;
  GridLayer **layers;
};

void grid_store_init(GridStore &grid_store, long gridsize);
void grid_store_delete(GridStore &grid_store);

void store_link_cnsrv(RemapVars &rv, long add1, long add2, long num_wts, double *weights, GridStore &grid_store);

#endif /* REMAP_STORE_LINK_CNSRV */
