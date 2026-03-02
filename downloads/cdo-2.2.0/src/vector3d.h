/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef VECTOR3D_H
#define VECTOR3D_H

#define _USE_MATH_DEFINES

#include <cstdio>
#include <cmath>

// clang-format off
class  // Vector3d
#ifdef WARN_UNUSED
[[gnu::warn_unused]]
#endif
Vector3d
// clang-format on
{
private:
  double X = 0.0, Y = 0.0, Z = 0.0;

public:
  Vector3d() noexcept {}
  Vector3d(const double &x, const double &y, const double &z) noexcept : X(x), Y(y), Z(z) {}

  explicit Vector3d(const double coords[3]) noexcept : X(coords[0]), Y(coords[1]), Z(coords[2]) {}

  double
  operator[](size_t index) const noexcept
  {
    switch (index)
      {
      case 0: return X;
      case 1: return Y;
      case 2: return Z;
      }

    return double();
  }

  Vector3d
  operator+(const Vector3d &other) const noexcept
  {
    return Vector3d(X + other.X, Y + other.Y, Z + other.Z);
  }

  Vector3d
  operator-(const Vector3d &other) const noexcept
  {
    return Vector3d(X - other.X, Y - other.Y, Z - other.Z);
  }

  Vector3d
  operator-(void) const noexcept
  {
    return Vector3d(-X, -Y, -Z);
  }

  // Calculate the cross/outer/vector product
  Vector3d
  operator%(const Vector3d &other) const noexcept
  {
    return Vector3d(Y * other.Z - Z * other.Y, Z * other.X - X * other.Z, X * other.Y - Y * other.X);
  }

  // Division by scalars
  Vector3d
  operator/(const double scalar) const noexcept
  {
    return Vector3d(X / scalar, Y / scalar, Z / scalar);
  }

  Vector3d
  operator/=(const double scalar) noexcept
  {
    return *this = *this / scalar;
  }

  // Calculate the dot/inner/scalar  product
  double
  operator*(const Vector3d &other) const
  {
    return (X * other.X) + (Y * other.Y) + (Z * other.Z);
  }

  double
  magnitude() const noexcept
  {
    return std::sqrt((X * X) + (Y * Y) + (Z * Z));
  }

  Vector3d
  normalised() const noexcept
  {
    return Vector3d(*this) / magnitude();
  }
  /*
  Vector3d d_normalize()
  {
    const double dnorm = std::sqrt(*this * *this);
    return Vector3d(*this) / dnorm;
  }
  */
  double
  longitude() const noexcept
  {
    return std::atan2(Y, X);
  }

  double
  latitude() const noexcept
  {
    return M_PI_2 - std::acos(Z);
  }
};

#endif
