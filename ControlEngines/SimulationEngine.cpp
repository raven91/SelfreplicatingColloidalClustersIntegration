//
// Created by Nikita Kruk on 27.02.20.
//

#include "SimulationEngine.hpp"
#include "../Steppers/RungeKutta2Stepper.hpp"
#include "../Steppers/StochasticEulerStepper.hpp"
#include "../Steppers/VelocityVerletStepper.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedron.hpp"
#include "../DynamicalSystems/SelfreplicationForOctahedronNondimensionalized.hpp"
#include "../Observers/BinaryObserver.hpp"

#include <random>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_map>

#include <eigen3/Eigen/Dense>

SimulationEngine::SimulationEngine(int number_of_colloidal_particles,
                                   Real colloid_diameter,
                                   Real solvent_diameter,
                                   Real colloid_volume_fraction) :
    number_of_colloidal_particles_(number_of_colloidal_particles),
    system_state_(),
    colloid_diameter_(colloid_diameter),
    solvent_diameter_(solvent_diameter),
    colloid_volume_fraction_(colloid_volume_fraction),
//    pbc_config_(10.0)
    pbc_config_(std::cbrt(
        number_of_colloidal_particles * M_PI * std::pow(colloid_diameter, 3.0) / (6.0 * colloid_volume_fraction)))
{
  number_of_solvent_particles_ = int(3.0 * pbc_config_.GetXSize() * pbc_config_.GetYSize() * pbc_config_.GetZSize());
  number_of_all_particles_ = number_of_colloidal_particles_ + number_of_solvent_particles_;
}
SimulationEngine::SimulationEngine(int number_of_colloidal_particles,
                                   int number_of_solvent_particles,
                                   Real colloid_diameter,
                                   Real solvent_diameter,
                                   Real colloid_volume_fraction) :
    number_of_colloidal_particles_(number_of_colloidal_particles),
    number_of_solvent_particles_(number_of_solvent_particles),
    number_of_all_particles_(number_of_colloidal_particles + number_of_solvent_particles),
    system_state_(),
    colloid_diameter_(colloid_diameter),
    solvent_diameter_(solvent_diameter),
    colloid_volume_fraction_(colloid_volume_fraction),
//    pbc_config_(5.0)
    pbc_config_(std::cbrt(
        ((number_of_colloidal_particles > 0) ? number_of_colloidal_particles : number_of_solvent_particles)
            * M_PI * std::pow(colloid_diameter, 3.0) / (6.0 * colloid_volume_fraction)))
{

}

SimulationEngine::~SimulationEngine()
{
  system_state_.clear();
}

void SimulationEngine::InitializeRandomSystemState()
{
  std::mt19937 mersenne_twister_generator(std::random_device{}());
  std::uniform_real_distribution<Real> unif_real_dist0L(0.0, pbc_config_.GetXSize());
  std::uniform_real_distribution<Real> unif_real_dist0pi(0.0, M_PI);
  std::uniform_real_distribution<Real> unif_real_dist02pi(0.0, 2.0 * M_PI);
#if defined(NONDIMENSIONALIZED)
  std::normal_distribution<Real> maxwell_boltzmann_distribution_for_colloids(0.0, std::sqrt(1.0 / kMassOfPolystyrene));
  std::normal_distribution<Real> maxwell_boltzmann_distribution_for_solvents(0.0, std::sqrt(1.0 / kMassOfPolyNipam));
#elif
  std::normal_distribution<Real> maxwell_boltzmann_distribution_for_colloids(0.0, std::sqrt(kBoltzmannConstant * kTemperature / kMassOfPolystyrene));
  std::normal_distribution<Real> maxwell_boltzmann_distribution_for_solvents(0.0, std::sqrt(kBoltzmannConstant * kTemperature / kMassOfPolyNipam));
#endif

  int number_of_prepared_colloids = 0;//ReadInOctahedronAndCatalyst();

  Eigen::Vector3d position(Eigen::Vector3d::Zero()), velocity(Eigen::Vector3d::Zero());
  const Real maximal_cutoff_radius = 1.5 * colloid_diameter_;
  int i = 0;
  while (i < number_of_colloidal_particles_ - number_of_prepared_colloids)
  {
    position.x() = unif_real_dist0L(mersenne_twister_generator);
    position.y() = unif_real_dist0L(mersenne_twister_generator);
    position.z() = unif_real_dist0L(mersenne_twister_generator);
    bool is_valid = true;
    for (Sphere &already_inserted_particles : system_state_)
    {
      if ((position - already_inserted_particles.GetPosition()).norm() <= maximal_cutoff_radius)
      {
        is_valid = false;
      }
    } // already_inserted_particles
    if (is_valid)
    {
      velocity.x() = maxwell_boltzmann_distribution_for_colloids(mersenne_twister_generator);
      velocity.y() = maxwell_boltzmann_distribution_for_colloids(mersenne_twister_generator);
      velocity.z() = maxwell_boltzmann_distribution_for_colloids(mersenne_twister_generator);
      system_state_.emplace_back(ParticleType::kColloidal,
                                 Sphere::GetNewPeriodicSubtype(),
//                                 ParticleSubtype::k0,
//                                 Sphere::GetNewComplementaryPeriodicSubtype(),
                                 position,
                                 velocity,
                                 kMassOfPolystyrene,
                                 PolymerType::kNone,
                                 -1);
      ++i;
    }
  } // i

  //  Real theta = 0.0, phi = 0.0;
  i = 0;
  while (i < number_of_solvent_particles_)
  {
    position.x() = unif_real_dist0L(mersenne_twister_generator);
    position.y() = unif_real_dist0L(mersenne_twister_generator);
    position.z() = unif_real_dist0L(mersenne_twister_generator);
    bool is_valid = true;
//    for (int j = 0; j < system_state_.size(); ++j)
    for (int j = 0; j < number_of_colloidal_particles_; ++j)
    {
      if ((position - system_state_[j].GetPosition()).norm() <= maximal_cutoff_radius)
      {
        is_valid = false;
      }
    } // already_inserted_particles
    if (is_valid)
    {
//      theta = unif_real_dist0pi(mersenne_twister_generator);
//      phi = unif_real_dist02pi(mersenne_twister_generator);
//      velocity.x() = 0.0 * std::sin(theta) * std::cos(phi);
//      velocity.y() = 0.0 * std::sin(theta) * std::sin(phi);
//      velocity.z() = 0.0 * std::cos(theta);
      velocity.x() = maxwell_boltzmann_distribution_for_solvents(mersenne_twister_generator);
      velocity.y() = maxwell_boltzmann_distribution_for_solvents(mersenne_twister_generator);
      velocity.z() = maxwell_boltzmann_distribution_for_solvents(mersenne_twister_generator);
      system_state_.emplace_back(ParticleType::kSolvent,
                                 ParticleSubtype::k0,
                                 position,
                                 velocity,
                                 kMassOfPolyNipam,
                                 PolymerType::kNone,
                                 -1);
      ++i;
    }
  } // i
//  position = {pbc_config_.GetXSize() * 0.1, pbc_config_.GetYSize() * 0.5, pbc_config_.GetZSize() * 0.5};
//  velocity = {-1.0, 0.0, 0.0};
//  velocity *= maxwell_boltzmann_distribution_for_colloids.stddev();
//  system_state_.emplace_back(ParticleType::kColloidal, position, velocity, kMassOfPolystyrene);
//  position = {pbc_config_.GetXSize() * 0.9, pbc_config_.GetYSize() * 0.5, pbc_config_.GetZSize() * 0.5};
//  velocity = {1.0, 0.0, 0.0};
//  velocity *= maxwell_boltzmann_distribution_for_colloids.stddev();
//  system_state_.emplace_back(ParticleType::kColloidal, position, velocity, kMassOfPolystyrene);
}

int SimulationEngine::ReadInOctahedronAndCatalyst()
{
#if defined(BCS_CLUSTER)
  std::string folder("/home/nkruk/cpp/SelfreplicatingColloidalClustersIntegration/input/");
#else
  std::string folder("/Users/nikita/Documents/Projects/SelfreplicatingColloidalClusters/");
#endif
  std::ifstream parameter_file(folder + "simulation_parameters_octahedron_and_catalyst_assembly.txt", std::ios::in);
  assert(parameter_file.is_open());
  std::string key("");
  Real value = 0.0;
  std::unordered_map<std::string, Real> parameters_dictionary;
  while (parameter_file >> key >> value)
  {
    parameters_dictionary[key] = value;
  }
  parameter_file.close();

  const int number_of_state_variables = 6, n_time_steps = 2e+4;
  std::ifstream file_for_colloidal_particles
      (folder + "colloidal_particles_octahedron_and_catalyst_assembly.bin", std::ios::binary | std::ios::in);
  assert(file_for_colloidal_particles.is_open());
  file_for_colloidal_particles.seekg(
      n_time_steps * (1l + number_of_state_variables * parameters_dictionary["number_of_colloidal_particles"])
          * sizeof(RealOutput), std::ios::beg);
  std::ifstream file_with_particle_subtypes
      (folder + "particle_subtypes_octahedron_and_catalyst_assembly.bin", std::ios::binary | std::ios::in);
  assert(file_with_particle_subtypes.is_open());
  file_with_particle_subtypes.seekg(
      n_time_steps * (sizeof(RealOutput) + parameters_dictionary["number_of_colloidal_particles"] * sizeof(int)),
      std::ios::beg);

  RealOutput time_point = 0.0;
  file_for_colloidal_particles.read((char *) &time_point, sizeof(RealOutput));
  assert(file_for_colloidal_particles.good());
  file_with_particle_subtypes.read((char *) &time_point, sizeof(RealOutput));
  assert(file_with_particle_subtypes.good());

  PeriodicBoundaryConditions old_periodic_boundary_conditions
      (parameters_dictionary["x_size"], parameters_dictionary["y_size"], parameters_dictionary["z_size"]);
  Eigen::Vector3d position(Eigen::Vector3d::Zero()), velocity(Eigen::Vector3d::Zero()), dr(Eigen::Vector3d::Zero());
  std::vector<RealOutput> input_position(position.rows() * position.cols(), 0.0),
      input_velocity(velocity.rows() * velocity.cols());
  int particle_subtype = static_cast<int>(ParticleSubtype::k0);
  PolymerType polymer_type = PolymerType::kNone;
  int polymer_index = 0;
  for (int i = 0; i < parameters_dictionary["number_of_colloidal_particles"]; ++i)
  {
    file_for_colloidal_particles.read((char *) &input_position[0], input_position.size() * sizeof(RealOutput));
    file_for_colloidal_particles.read((char *) &input_velocity[0], input_velocity.size() * sizeof(RealOutput));
    file_with_particle_subtypes.read((char *) &particle_subtype, sizeof(int));
    if (i == 0)
    {
      dr.x() = input_position[0] - old_periodic_boundary_conditions.GetXSize() / 2.0;
      dr.y() = input_position[1] - old_periodic_boundary_conditions.GetYSize() / 2.0;
      dr.z() = input_position[2] - old_periodic_boundary_conditions.GetZSize() / 2.0;
    }

    // take into account that the previous simulation might have been done with different scaling
    position.x() = input_position[0] - dr.x();// / parameters_dictionary["x_size"] * pbc_config_.GetXSize();
    position.y() = input_position[1] - dr.y();// / parameters_dictionary["y_size"] * pbc_config_.GetYSize();
    position.z() = input_position[2] - dr.z();// / parameters_dictionary["z_size"] * pbc_config_.GetZSize();
    old_periodic_boundary_conditions.ApplyPeriodicBoundaryConditions(position.x(), position.y(), position.z(),
                                                                     position.x(), position.y(), position.z());
    position.x() += (pbc_config_.GetXSize() - old_periodic_boundary_conditions.GetXSize()) / 2.0;
    position.y() += (pbc_config_.GetYSize() - old_periodic_boundary_conditions.GetYSize()) / 2.0;
    position.z() += (pbc_config_.GetZSize() - old_periodic_boundary_conditions.GetZSize()) / 2.0;
    pbc_config_.ApplyPeriodicBoundaryConditions(position.x(), position.y(), position.z(),
                                                position.x(), position.y(), position.z());
    // todo: consider the scaling for the dimensional velocity
    velocity.x() = input_velocity[0];
    velocity.y() = input_velocity[1];
    velocity.z() = input_velocity[2];
    if (kParentOctahedronParticles.find(static_cast<ParticleSubtype>(particle_subtype))
        != kParentOctahedronParticles.end())
    {
      polymer_type = PolymerType::kOctahedron;
      polymer_index = 0;
    } else if (kParentCatalystParticles.find(static_cast<ParticleSubtype>(particle_subtype))
        != kParentCatalystParticles.end())
    {
      polymer_type = PolymerType::kCatalyst;
      polymer_index = 1;
    } else
    {
      polymer_type = PolymerType::kNone;
      polymer_index = -1;
    }
    system_state_.emplace_back(ParticleType::kColloidal,
                               static_cast<ParticleSubtype>(particle_subtype),
                               position,
                               velocity,
                               kMassOfPolystyrene,
                               polymer_type,
                               polymer_index);
  } // i

  return system_state_.size();
}

void SimulationEngine::RunSimulation()
{
  const Real t_0 = 0.0;
  const Real t_1 = 1e+7;
  const long n_time_steps = 2e+6;
  Real dt = 0.0005;
  InitializeRandomSystemState();

  VelocityVerletStepper stepper(number_of_all_particles_, dt);
#if defined(NONDIMENSIONALIZED)
  SelfreplicationForOctahedronNondimensionalized system(pbc_config_, number_of_all_particles_, colloid_diameter_);
#elif
  SelfreplicationForOctahedron system(pbc_config_, number_of_all_particles_, colloid_diameter_);
#endif
  BinaryObserver binary_observer(pbc_config_);

  Real t = t_0;
  long time_step = 0;
  binary_observer.SaveSimulationParameters(number_of_colloidal_particles_,
                                           number_of_solvent_particles_,
                                           number_of_all_particles_,
                                           colloid_diameter_,
                                           solvent_diameter_,
                                           colloid_volume_fraction_,
                                           dt);
  binary_observer.SaveSystemState(system_state_, t);
  binary_observer.SaveColloidalParticleSubtype(system_state_, t);
//  while (t <= t_1)
  while (time_step < n_time_steps)
  {
    t += dt;
    ++time_step;
    stepper.DoStep(system, system_state_);
    pbc_config_.ApplyPeriodicBoundaryConditions(system_state_);
    {
      binary_observer.SaveSystemState(system_state_, t);
      binary_observer.SaveColloidalParticleSubtype(system_state_, t);
    }
  } // t
}

void SimulationEngine::DetermineNumberOfBonds()
{

}