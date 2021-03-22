//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_RUNGEKUTTA2STEPPER_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_RUNGEKUTTA2STEPPER_HPP

#include "../Definitions.hpp"
#include "../ParticleTypes/Sphere.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedron.hpp"

#include <vector>

class RungeKutta2Stepper
{
 public:

  explicit RungeKutta2Stepper(int number_of_all_particles, Real dt);
  ~RungeKutta2Stepper();

  void DoStep(SelfreplicationForOctahedron &system, std::vector<Sphere> &system_state);

 private:

  int number_of_all_particles_;
  Real dt_;
  int order_;
  std::vector<Sphere> k_1_;
  std::vector<Sphere> k_2_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_RUNGEKUTTA2STEPPER_HPP
