//
// Created by Nikita Kruk on 02.03.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_VELOCITYVERLETSTEPPER_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_VELOCITYVERLETSTEPPER_HPP

#include "../Definitions.hpp"
#include "../ParticleTypes/Sphere.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedron.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedronNondimensionalized.hpp"

#include <vector>

class VelocityVerletStepper
{
 public:

  explicit VelocityVerletStepper(int number_of_all_particles, Real dt);
  ~VelocityVerletStepper();

  void DoStep(SelfreplicationForOctahedron &system, std::vector<Sphere> &system_state);
  void DoStep(SelfreplicationForOctahedronNondimensionalized &system, std::vector<Sphere> &system_state);

 private:

  int number_of_all_particles_;
  Real dt_;
  int order_;
  std::vector<Sphere> k_1_;
  std::vector<Sphere> k_2_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_VELOCITYVERLETSTEPPER_HPP
