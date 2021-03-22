//
// Created by Nikita Kruk on 28.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_STOCHASTICEULERSTEPPER_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_STOCHASTICEULERSTEPPER_HPP

#include "../Definitions.hpp"
#include "../ParticleTypes/Sphere.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedron.hpp"

#include <cmath>
#include <algorithm> // std::fill

#include <eigen3/Eigen/Dense>

class StochasticEulerStepper
{
 public:

  explicit StochasticEulerStepper(Real dt) :
      dt_(dt)
  {

  }

  ~StochasticEulerStepper() = default;

  void DoStep(SelfreplicationForOctahedron &particle_system, std::vector<Sphere> &x)
  {
    static std::vector<Sphere> deterministic(x.size(), Sphere()), stochastic(x.size(), Sphere());
    std::fill(deterministic.begin(), deterministic.end(), Sphere());
    std::fill(stochastic.begin(), stochastic.end(), Sphere());

    particle_system.EvaluateRhsWithAllPairs(x, deterministic, deterministic, 0.0, dt_, dt_);
    particle_system.AddNoise(x, stochastic);

    Eigen::Vector3d position(Eigen::Vector3d::Zero()), velocity(Eigen::Vector3d::Zero());
    for (int i = 0; i < x.size(); ++i)
    {
      position = x[i].GetPosition();
      velocity = x[i].GetVelocity();
      position += deterministic[i].GetPosition() * dt_ + stochastic[i].GetPosition() * std::sqrt(dt_);
      velocity += deterministic[i].GetVelocity() * dt_ + stochastic[i].GetVelocity() * std::sqrt(dt_);
      x[i].SetPosition(position);
      x[i].SetVelocity(velocity);
    } // i
  }

 private:

  Real dt_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_STOCHASTICEULERSTEPPER_HPP
