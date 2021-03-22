//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP

#include "../Definitions.hpp"
#include "../ParticleTypes/Sphere.hpp"

#include <vector>

#include <eigen3/Eigen/Dense>

class PeriodicBoundaryConditions
{
 public:

  explicit PeriodicBoundaryConditions(Real x_size, Real y_size, Real z_size);
  explicit PeriodicBoundaryConditions(Real size_of_each_dimension);
  ~PeriodicBoundaryConditions();

  void ClassAEffectiveParticleDistance(Real x_i,
                                       Real y_i,
                                       Real z_i,
                                       Real x_j,
                                       Real y_j,
                                       Real z_j,
                                       Real &dx,
                                       Real &dy,
                                       Real &dz);
  void ClassAEffectiveParticleDistance(const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j, Eigen::Vector3d &r_ij);

  void ClassBEffectiveParticleDistanceSigned(Real x_i,
                                             Real y_i,
                                             Real z_i,
                                             Real x_j,
                                             Real y_j,
                                             Real z_j,
                                             Real &dx,
                                             Real &dy,
                                             Real &dz);
  void ClassBEffectiveParticleDistanceUnsigned(Real x_i,
                                               Real y_i,
                                               Real z_i,
                                               Real x_j,
                                               Real y_j,
                                               Real z_j,
                                               Real &dx,
                                               Real &dy,
                                               Real &dz);

  void ClassCEffectiveParticleDistanceSigned(Real x_i,
                                             Real y_i,
                                             Real z_i,
                                             Real x_j,
                                             Real y_j,
                                             Real z_j,
                                             Real &dx,
                                             Real &dy,
                                             Real &dz);
  void ClassCEffectiveParticleDistanceUnsigned(Real x_i,
                                               Real y_i,
                                               Real z_i,
                                               Real x_j,
                                               Real y_j,
                                               Real z_j,
                                               Real &dx,
                                               Real &dy,
                                               Real &dz);

  void ApplyPeriodicBoundaryConditions(Real x, Real y, Real z, Real &x_pbc, Real &y_pbc, Real &z_pbc);
  void ApplyPeriodicBoundaryConditions(std::vector<Sphere> &system_state);
  void ApplyPeriodicBoundaryConditions(Sphere *const system_state, long size);

  Real GetXSize();
  Real GetYSize();
  Real GetZSize();

 private:

  Real x_size_;
  Real y_size_;
  Real z_size_;
  Real x_rsize_;
  Real y_rsize_;
  Real z_rsize_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_PERIODICBOUNDARYCONDITIONS_HPP
