//
// Created by Nikita Kruk on 27.02.20.
//

#include "RungeKutta2Stepper.hpp"

#include <iostream>
#include <eigen3/Eigen/Dense>

RungeKutta2Stepper::RungeKutta2Stepper(int number_of_all_particles, Real dt) :
    number_of_all_particles_(number_of_all_particles),
    dt_(dt),
    order_(2),
    k_1_(number_of_all_particles, Sphere()),
    k_2_(number_of_all_particles, Sphere())
{

}

RungeKutta2Stepper::~RungeKutta2Stepper()
{
  k_1_.clear();
  k_2_.clear();
}

void RungeKutta2Stepper::DoStep(SelfreplicationForOctahedron &system, std::vector<Sphere> &system_state)
{
  system.EvaluateRhsWithAllPairs(system_state, k_1_, k_1_, 0.0, dt_, dt_);
  system.EvaluateRhsWithAllPairs(system_state, k_1_, k_2_, 1.0, dt_, dt_);
  static Eigen::Vector3d position, velocity;
  for (int i = 0; i < system_state.size(); ++i)
  {
    position = system_state[i].GetPosition();
    velocity = system_state[i].GetVelocity();
    position += (k_1_[i].GetPosition() + k_2_[i].GetPosition()) * dt_ / 2.0;
    velocity += (k_1_[i].GetVelocity() + k_2_[i].GetVelocity()) * dt_ / 2.0;
    system_state[i].SetPosition(position);
    system_state[i].SetVelocity(velocity);
//    if (!std::isfinite(system_state[i].GetVelocity().x()) || !std::isfinite(system_state[i].GetVelocity().y()) || !std::isfinite(system_state[i].GetVelocity().z()))
//    {
//      std::cout << "debug\n";
//    }
  } // i
}