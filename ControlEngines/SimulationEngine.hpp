//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SIMULATIONENGINE_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../ParticleTypes/Sphere.hpp"

#include <vector>

class SimulationEngine
{
 public:

  explicit SimulationEngine(int number_of_colloidal_particles,
                            Real colloid_diameter,
                            Real solvent_diameter,
                            Real colloid_volume_fraction);
  explicit SimulationEngine(int number_of_colloidal_particles,
                            int number_of_solvent_particles,
                            Real colloid_diameter,
                            Real solvent_diameter,
                            Real colloid_volume_fraction);
  ~SimulationEngine();

  void RunSimulation();

 private:

  int number_of_colloidal_particles_;
  int number_of_solvent_particles_;
  int number_of_all_particles_;
  std::vector<Sphere> system_state_;
  Real colloid_diameter_;
  Real solvent_diameter_;
  Real colloid_volume_fraction_;
  PeriodicBoundaryConditions pbc_config_;

  void InitializeRandomSystemState();
  int ReadInOctahedronAndCatalyst();
  void DetermineNumberOfBonds();

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SIMULATIONENGINE_HPP
