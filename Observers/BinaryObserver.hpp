//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_BINARYOBSERVER_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_BINARYOBSERVER_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../ParticleTypes/Sphere.hpp"

#include <string>
#include <chrono>
#include <unordered_map>

class BinaryObserver
{
 public:

  explicit BinaryObserver(PeriodicBoundaryConditions &pbc_config);
  ~BinaryObserver();

  void SaveSimulationParameters(int number_of_colloidal_particles,
                                int number_of_solvent_particles,
                                int number_of_all_particles,
                                Real colloid_diameter,
                                Real solvent_diameter,
                                Real colloid_volume_fraction,
                                Real dt);
  void SaveSystemState(const std::vector<Sphere> &system_state, Real t);
  void SaveColloidalParticleSubtype(const std::vector<Sphere> &system_state, Real t);

 private:

  PeriodicBoundaryConditions &pbc_config_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int output_time_counter_[2];
  int output_time_threshold_[2];

  std::string parameter_file_name_;
  std::unordered_map<ParticleType, std::string> simulation_file_names_;
  std::string particle_subtype_file_name_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_BINARYOBSERVER_HPP
