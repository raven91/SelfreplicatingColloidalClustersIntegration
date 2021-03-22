//
// Created by Nikita Kruk on 27.02.20.
//

#include "PeriodicBoundaryConditions.hpp"

PeriodicBoundaryConditions::PeriodicBoundaryConditions(Real x_size, Real y_size, Real z_size) :
    x_size_(x_size),
    y_size_(y_size),
    z_size_(z_size),
    x_rsize_(1.0 / x_size),
    y_rsize_(1.0 / y_size),
    z_rsize_(1.0 / z_size)
{

}

PeriodicBoundaryConditions::PeriodicBoundaryConditions(Real size_of_each_dimension) :
    x_size_(size_of_each_dimension),
    y_size_(size_of_each_dimension),
    z_size_(size_of_each_dimension),
    x_rsize_(1.0 / x_size_),
    y_rsize_(1.0 / y_size_),
    z_rsize_(1.0 / z_size_)
{

}

PeriodicBoundaryConditions::~PeriodicBoundaryConditions() = default;

// periodic signed distance if all interaction sites are in the same simulation box
void PeriodicBoundaryConditions::ClassAEffectiveParticleDistance(Real x_i,
                                                                 Real y_i,
                                                                 Real z_i,
                                                                 Real x_j,
                                                                 Real y_j,
                                                                 Real z_j,
                                                                 Real &dx,
                                                                 Real &dy,
                                                                 Real &dz)
{
  dx = x_j - x_i;
  dx -= static_cast<int>(dx * 2.0 * x_rsize_) * x_size_;

  dy = y_j - y_i;
  dy -= static_cast<int>(dy * 2.0 * y_rsize_) * y_size_;

  dz = z_j - z_i;
  dz -= static_cast<int>(dz * 2.0 * z_rsize_) * z_size_;
}

// periodic signed distance if all interaction sites are in the same simulation box
void PeriodicBoundaryConditions::ClassAEffectiveParticleDistance(const Eigen::Vector3d &r_i,
                                                                 const Eigen::Vector3d &r_j,
                                                                 Eigen::Vector3d &r_ij)
{
  r_ij = r_j - r_i;
  r_ij.x() -= static_cast<int>(r_ij.x() * 2.0 * x_rsize_) * x_size_;
  r_ij.y() -= static_cast<int>(r_ij.y() * 2.0 * y_rsize_) * y_size_;
  r_ij.z() -= static_cast<int>(r_ij.z() * 2.0 * z_rsize_) * z_size_;
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditions::ClassBEffectiveParticleDistanceSigned(Real x_i,
                                                                       Real y_i,
                                                                       Real z_i,
                                                                       Real x_j,
                                                                       Real y_j,
                                                                       Real z_j,
                                                                       Real &dx,
                                                                       Real &dy,
                                                                       Real &dz)
{
  dx = x_j - x_i;
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  } else if (dx < -0.5 * x_size_)
  {
    dx += x_size_;
  }

  dy = y_j - y_i;
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  } else if (dy < -0.5 * y_size_)
  {
    dy += y_size_;
  }

  dz = z_j - z_i;
  if (dz > 0.5 * z_size_)
  {
    dz -= z_size_;
  } else if (dz < -0.5 * z_size_)
  {
    dz += z_size_;
  }
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditions::ClassBEffectiveParticleDistanceUnsigned(Real x_i,
                                                                         Real y_i,
                                                                         Real z_i,
                                                                         Real x_j,
                                                                         Real y_j,
                                                                         Real z_j,
                                                                         Real &dx,
                                                                         Real &dy,
                                                                         Real &dz)
{
  dx = std::fabs(x_j - x_i);
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  }

  dy = std::fabs(y_j - y_i);
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  }

  dz = std::fabs(z_j - z_i);
  if (dz > 0.5 * z_size_)
  {
    dz -= z_size_;
  }
}

// calculate remaindars at any distance
//if the sign of the distance is relevant
void PeriodicBoundaryConditions::ClassCEffectiveParticleDistanceSigned(Real x_i,
                                                                       Real y_i,
                                                                       Real z_i,
                                                                       Real x_j,
                                                                       Real y_j,
                                                                       Real z_j,
                                                                       Real &dx,
                                                                       Real &dy,
                                                                       Real &dz)
{
  dx = x_j - x_i;
  dx -= x_size_ * std::nearbyint(dx * x_rsize_);

  dy = y_j - y_i;
  dy -= y_size_ * std::nearbyint(dy * y_rsize_);

  dz = z_j - z_i;
  dz -= z_size_ * std::nearbyint(dz * z_rsize_);
}

// calculate remaindars at any distance
//if the sign of the distance is not relevant
//#pragma acc routine seq
void PeriodicBoundaryConditions::ClassCEffectiveParticleDistanceUnsigned(Real x_i,
                                                                         Real y_i,
                                                                         Real z_i,
                                                                         Real x_j,
                                                                         Real y_j,
                                                                         Real z_j,
                                                                         Real &dx,
                                                                         Real &dy,
                                                                         Real &dz)
{
  dx = std::fabs(x_j - x_i);
  dx -= static_cast<int>(dx * x_rsize_ + 0.5) * x_size_;

  dy = std::fabs(y_j - y_i);
  dy -= static_cast<int>(dy * y_rsize_ + 0.5) * y_size_;

  dz = std::fabs(z_j - z_i);
  dz -= static_cast<int>(dz * z_rsize_ + 0.5) * z_size_;
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(Real x,
                                                                 Real y,
                                                                 Real z,
                                                                 Real &x_pbc,
                                                                 Real &y_pbc,
                                                                 Real &z_pbc)
{
  x_pbc = x - std::floor(x * x_rsize_) * x_size_;
  y_pbc = y - std::floor(y * y_rsize_) * y_size_;
  z_pbc = z - std::floor(z * z_rsize_) * z_size_;
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(std::vector<Sphere> &system_state)
{
  static Eigen::Vector3d position;
  for (Sphere &sphere : system_state)
  {
    position = sphere.GetPosition();
    position.x() -= std::floor(position.x() * x_rsize_) * x_size_;
    position.y() -= std::floor(position.y() * y_rsize_) * y_size_;
    position.z() -= std::floor(position.z() * z_rsize_) * z_size_;
    sphere.SetPosition(position);
  } // i
}

void PeriodicBoundaryConditions::ApplyPeriodicBoundaryConditions(Sphere *const system_state, long size)
{
  static Eigen::Vector3d position;
  for (int i = 0; i < size; ++i)
  {
    position = system_state[i].GetPosition();
    position.x() -= std::floor(position.x() * x_rsize_) * x_size_;
    position.y() -= std::floor(position.y() * y_rsize_) * y_size_;
    position.z() -= std::floor(position.z() * z_rsize_) * z_size_;
    system_state[i].SetPosition(position);
  } // i
}

Real PeriodicBoundaryConditions::GetXSize()
{
  return x_size_;
}

Real PeriodicBoundaryConditions::GetYSize()
{
  return y_size_;
}

Real PeriodicBoundaryConditions::GetZSize()
{
  return z_size_;
}