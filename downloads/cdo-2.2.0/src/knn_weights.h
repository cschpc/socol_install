/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef KNN_WEIGHTS_H
#define KNN_WEIGHTS_H

#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <vector>

#include "varray.h"
#include "interpol.h"

class knnWeightsType
{
private:
  size_t m_maxNeighbors = 0;

public:
  size_t m_numNeighbors = 0;
  std::vector<uint8_t> m_mask;  // mask at nearest neighbors
  std::vector<size_t> m_addr;   // source address at nearest neighbors
  std::vector<double> m_dist;   // angular distance four nearest neighbors
  std::vector<size_t> m_tmpaddr;
  std::vector<double> m_tmpdist;

  inline void
  init()
  {
    m_mask.resize(m_maxNeighbors);
    m_addr.resize(m_maxNeighbors);
    m_dist.resize(m_maxNeighbors);
  }

  knnWeightsType(size_t maxNeighbors) : m_maxNeighbors(maxNeighbors) { init(); }

  inline size_t
  maxNeighbors() const
  {
    return m_maxNeighbors;
  }

  inline size_t
  numNeighbors() const
  {
    return m_numNeighbors;
  }

  inline void
  initAddr(size_t numNeighbors)
  {
    for (size_t i = 0; i < numNeighbors; ++i) m_addr[i] = SIZE_MAX;
  }

  inline void
  initDist(size_t numNeighbors)
  {
    for (size_t i = 0; i < numNeighbors; ++i) m_dist[i] = DBL_MAX;
  }

  inline void
  initAddr()
  {
    initAddr(m_maxNeighbors);
  }

  inline void
  initDist()
  {
    initDist(m_maxNeighbors);
  }

  inline bool
  distance_is_less(double distance, double distx, size_t addr, size_t addrx)
  {
    constexpr double cmp_tol = 1.e-12;
    // return (distance < distx || (distance <= distx && addr < addrx));
    return (distance + cmp_tol) < distx || (addr < addrx && std::fabs(distance - distx) < cmp_tol);
  }

  inline void
  storeDistance(size_t addr, double distance, size_t numNeighbors)
  {
    assert(numNeighbors <= m_maxNeighbors);
    m_numNeighbors = numNeighbors;

    if (numNeighbors == 1)
      {
        if (distance_is_less(distance, m_dist[0], addr, m_addr[0]))
          {
            m_addr[0] = addr;
            m_dist[0] = distance;
          }
      }
    else
      {
        for (size_t i = 0; i < numNeighbors; ++i)
          {
            if (distance_is_less(distance, m_dist[i], addr, m_addr[i]))
              {
                for (size_t n = numNeighbors - 1; n > i; --n)
                  {
                    m_addr[n] = m_addr[n - 1];
                    m_dist[n] = m_dist[n - 1];
                  }
                m_addr[i] = addr;
                m_dist[i] = distance;
                break;
              }
          }
      }
  }

  inline void
  storeDistance(size_t addr, double distance)
  {
    storeDistance(addr, distance, m_maxNeighbors);
  }

  inline void
  setDistance(const size_t *addr, const double *distance, size_t numNeighbors)
  {
    assert(numNeighbors <= m_maxNeighbors);
    m_numNeighbors = numNeighbors;

    for (size_t i = 0; i < numNeighbors; ++i) m_addr[i] = addr[i];
    for (size_t i = 0; i < numNeighbors; ++i) m_dist[i] = distance[i];
  }

  inline void
  checkDistance()
  {
    constexpr double eps = 1.e-14;
    // If distance is zero, set to small number
    for (size_t i = 0; i < m_numNeighbors; ++i)
      if (m_addr[i] < SIZE_MAX && m_dist[i] <= 0.0) m_dist[i] = eps;
  }

  size_t
  normalizeWeights(double dist_tot, size_t numNeighbors)
  {
    // Normalize weights and store the link
    size_t nadds = 0;

    for (size_t n = 0; n < numNeighbors; ++n)
      {
        if (m_mask[n])
          {
            m_dist[nadds] = m_dist[n] / dist_tot;
            m_addr[nadds] = m_addr[n];
            nadds++;
          }
      }

    m_numNeighbors = nadds;
    return nadds;
  }

  size_t
  computeWeights()
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.0;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_maxNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX)
          {
            m_dist[n] = 1.0 / m_dist[n];
            dist_tot += m_dist[n];
            m_mask[n] = true;
          }
      }

    return normalizeWeights(dist_tot, m_maxNeighbors);
  }

  size_t
  computeWeights(const Varray<short> &gridMask)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.0;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_maxNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX && gridMask[m_addr[n]])
          {
            m_dist[n] = 1.0 / m_dist[n];
            dist_tot += m_dist[n];
            m_mask[n] = true;
          }
      }

    return normalizeWeights(dist_tot, m_maxNeighbors);
  }

  size_t
  computeWeights(const Varray<uint8_t> &grid_mask, double searchRadius, double weight0, double weightR)
  {
    // Compute weights based on inverse distance if mask is false, eliminate those points

    double dist_tot = 0.0;  // sum of neighbor distances (for normalizing)

    for (size_t n = 0; n < m_numNeighbors; ++n)
      {
        m_mask[n] = false;
        if (m_addr[n] < SIZE_MAX && grid_mask[m_addr[n]])
          {
            m_dist[n] = intlin(m_dist[n], weight0, 0, weightR, searchRadius);
            dist_tot += m_dist[n];
            m_mask[n] = true;
          }
      }

    return normalizeWeights(dist_tot, m_numNeighbors);
  }

  template <typename T>
  double
  arrayWeightsSum(const Varray<T> &array) const
  {
    double result = 0.0;
    for (size_t n = 0; n < m_numNeighbors; ++n) result += array[m_addr[n]] * m_dist[n];
    return result;
  }
};

#endif
