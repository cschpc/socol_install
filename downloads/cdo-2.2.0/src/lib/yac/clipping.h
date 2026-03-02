/**
 * @file clipping.h
 * @brief Structs and interfaces for cell clipping
 *
 * @copyright Copyright  (C)  2013 Moritz Hanke <hanke@dkrz.de>
 *                                 Rene Redler <rene.redler@mpimet.mpg.de>
 *
 * @version 1.0
 * @author Moritz Hanke <hanke@dkrz.de>
 *         Rene Redler <rene.redler@mpimet.mpg.de>
 */
/*
 * Keywords:
 * Maintainer: Moritz Hanke <hanke@dkrz.de>
 *             Rene Redler <rene.redler@mpimet.mpg.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
 *
 * This file is part of YAC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef CLIPPING_H
#define CLIPPING_H

#include "grid_cell.h"

/** \example test_clipping.c
 * This contains some examples on how to use the \ref yac_cell_clipping
 * routine.
 */
/** \example test_lat_clipping.c
 * This contains some examples on how to use the yac_cell_lat_clipping
 * routine.
 */

/**
  * \brief cell clipping to get the cells describing the intersections
  *
  * The routine takes (a list of) source cells and a target cell. It sets the
  * target cell data and does some further initialisation. Thus it needs to be
  * called for each new target cell intersection calculation
  *
  * The vertices of source and target cells can be either provided in a clockwise
  * or anticlockwise sense.
  *
  * @param[in] N              number of source cells
  * @param[in] source_cell    list of source cells
  * @param[in] target_cell    target cell
  * @param[in] overlap_buffer buffer for the overlaps between the target and
  *                           the source cells
  *
  * \remark source cells can be either convex or concave
  * \remark target cell has to be convex
  * \remark cells in overlap_buffer can be concave (even if source cell was convex)
  * \remark overlap_buffer must contain valid grid_cells (have to be initialised
  *         using \ref yac_init_grid_cell; initialisation have to be done only once,
  *         in consecutive calls, the cells can be reused with have to be
  *         reinitialised)
  *
 **/
void yac_cell_clipping (size_t N,
                        struct grid_cell * source_cell,
                        struct grid_cell target_cell,
                        struct grid_cell * overlap_buffer);

/**
  * \brief cell clipping to get the cells describing the intersections
  *
  * The routine takes (a list of) cells and two latitude bounds.
  *
  * @param[in] N              number of cells
  * @param[in] cells          list of cells
  * @param[in] lat_bounds     latitude bounds in radiant
  * @param[in] overlap_buffer buffer for the overlaps between the cells and
  *                           latitude band
  *
  * \remark cells in overlap_buffer can be concave
  * \remark overlap_buffer must contain valid grid_cells (have to be initialised
  *         using \ref yac_init_grid_cell; initialisation have to be done only once,
  *         in consecutive calls, the cells can be reused with have to be
  *         reinitialised)
  * \remark this routine is currently not being used within YAC but potentially
  *         used within the CDOs
  *
 **/
void yac_cell_lat_clipping (size_t N,
                            struct grid_cell * cells,
                            double lat_bounds[2],
                            struct grid_cell * overlap_buffer);

/** \example test_partial_areas.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/** \example test_compute_overlap_area.c
 * This contains examples on how to use \ref yac_compute_overlap_areas.
 */

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with triangular target cells. This is required for
  *        conservative remapping
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a convex target cell. As
  * a triangle is always convex we recommend to use this routine only for
  * triangular target cells. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N             number of source cells
  * @param[in]  source_cell   list of source cells
  * @param[in]  target_cell   target cell
  * @param[out] partial_areas list of N partial weights, one weight for each
  *                           source-target intersection
  *
  * \remark source and target cell have to be convex
  *
 **/
void yac_compute_overlap_areas (size_t N,
                                struct grid_cell * source_cell,
                                struct grid_cell target_cell,
                                double * partial_areas);

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with arbitrary target cells, this is required for conservative
  *        remapping.
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a target cell. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N               number of source cells
  * @param[in]  source_cell     list of source cells
  * @param[in]  target_cell     target cell
  * @param[out] partial_areas   list of N partial weights, one weight for each
  *                             source-target intersection
  *
 **/
void yac_compute_concave_overlap_areas (size_t N,
                                        struct grid_cell * source_cell,
                                        struct grid_cell target_cell,
                                        double * partial_areas);

/**
  * \brief calculates partial areas for all overlapping parts of the source
  *        cells with arbitrary target cells, this is required for conservative
  *        remapping. In addition, the barycenter of each overlap is calculated.
  *
  * Some of the underlying concepts can be found in
  *
  * See e.g. Joseph O'Rourke, Computational Geometry in C, 2nd Ed., 1998
  *          Sec. 7.6 Intersections of Convex Polygons, page 253.
  *
  * The routine takes (a list of) source cells and a target cell. It determines the
  * clipping points of the intersection between a source and the target cells using
  * cell_clipping internally. In a second step areas are calculated for each
  * intersection described in the overlap cells. If a target cell is fully
  * covered by N source cells the N partial areas should add up to the area of
  * the target cell.
  *
  * @param[in]  N                   number of source cells
  * @param[in]  source_cell         list of source cells
  * @param[in]  target_cell         target cell
  * @param[out] overlap_areas       list of N partial weights, one weight for
  *                                 each source-target intersection
  * @param[out] overlap_barycenters coordinates of the barycenters of the
  *                                 overlap cell
  *
 **/
void yac_compute_concave_overlap_info (size_t N,
                                       struct grid_cell * source_cell,
                                       struct grid_cell target_cell,
                                       double * overlap_areas,
                                       double (*overlap_barycenters)[3]);
/**
  * \brief correct interpolation weights
  *
  * Returns weights with a sum close to 1.0
  *
  * @param[in]  N                 number of source cells
  * @param[out] weight            list of N partial weights
  *
 **/
void yac_correct_weights (size_t N, double * weight);

#endif // CLIPPING_H

