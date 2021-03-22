//
// Created by Nikita Kruk on 02.03.20.
//

#include "VelocityVerletStepper.hpp"

#include <eigen3/Eigen/Dense>

VelocityVerletStepper::VelocityVerletStepper(int number_of_all_particles, Real dt) :
    number_of_all_particles_(number_of_all_particles),
    dt_(dt),
    order_(2),
    k_1_(number_of_all_particles, Sphere()),
    k_2_(number_of_all_particles, Sphere())
{

}

VelocityVerletStepper::~VelocityVerletStepper()
{
  k_1_.clear();
  k_2_.clear();
}

void VelocityVerletStepper::DoStep(SelfreplicationForOctahedron &system, std::vector<Sphere> &system_state)
{
  system.EvaluateRhsWithVerletNeighborList(system_state, k_1_, k_1_, 0.0, dt_, dt_);
  static std::vector<Sphere> new_system_state(system_state.size(), Sphere());
  const Real lambda = 0.5;
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    new_system_state[i] = system_state[i]; // to transfer additional variables
    new_system_state[i].SetPosition(
        system_state[i].GetPosition() + dt_ * k_1_[i].GetPosition() + 0.5 * dt_ * dt_ * k_1_[i].GetVelocity());
    // the velocity is updated only half step first
    new_system_state[i].SetVelocity(system_state[i].GetVelocity() + lambda * dt_ * k_1_[i].GetVelocity());
  } // i
  system.EvaluateRhsWithVerletNeighborList(new_system_state, k_2_, k_2_, 0.0, dt_, dt_);
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    new_system_state[i].SetVelocity(
        system_state[i].GetVelocity() + 0.5 * dt_ * (k_1_[i].GetVelocity() + k_2_[i].GetVelocity()));
  } // i

  system_state = new_system_state;
}

void VelocityVerletStepper::DoStep(SelfreplicationForOctahedronNondimensionalized &system,
                                   std::vector<Sphere> &system_state)
{
  system.EvaluateRhsWithVerletNeighborList(system_state, k_1_, k_1_, 0.0, dt_, dt_);
  static std::vector<Sphere> new_system_state(system_state.size(), Sphere());
  const Real lambda = 0.5;
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    new_system_state[i] = system_state[i]; // to transfer additional variables
    new_system_state[i].SetPosition(
        system_state[i].GetPosition() + dt_ * k_1_[i].GetPosition() + 0.5 * dt_ * dt_ * k_1_[i].GetVelocity());
    // the velocity is updated only half step first
    new_system_state[i].SetVelocity(system_state[i].GetVelocity() + lambda * dt_ * k_1_[i].GetVelocity());
  } // i
  system.EvaluateRhsWithVerletNeighborList(new_system_state, k_2_, k_2_, 0.0, dt_, dt_);
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    new_system_state[i].SetVelocity(
        system_state[i].GetVelocity() + 0.5 * dt_ * (k_1_[i].GetVelocity() + k_2_[i].GetVelocity()));
  } // i

  system_state = new_system_state;
}