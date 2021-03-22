#include "ControlEngines/SimulationEngine.hpp"

int main()
{
  const int number_of_colloidal_particles = 4 + 0, number_of_solvent_particles = 0;
  const Real colloid_volume_fraction = 1.0 / 100.0;
  SimulationEngine engine(number_of_colloidal_particles,
                          number_of_solvent_particles,
                          kColloidDiameter,
                          kSolventDiameter,
                          colloid_volume_fraction);
//  SimulationEngine engine(number_of_colloidal_particles,
//                          kColloidDiameter,
//                          kSolventDiameter,
//                          colloid_volume_fraction);
  engine.RunSimulation();

  return 0;
}
