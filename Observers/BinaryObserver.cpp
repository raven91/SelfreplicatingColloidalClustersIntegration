//
// Created by Nikita Kruk on 27.02.20.
//

#include "BinaryObserver.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include <eigen3/Eigen/Dense>

BinaryObserver::BinaryObserver(PeriodicBoundaryConditions &pbc_config) :
    pbc_config_(pbc_config),
    output_time_counter_{0},
    output_time_threshold_{100, 100} // mod 1 - save at every dt
{
  integration_step_timer_ = std::chrono::system_clock::now();

#if defined(BCS_CLUSTER)
  std::string folder("/home/nkruk/cpp/SelfreplicatingColloidalClustersIntegration/output/");
#else
  std::string folder("/Users/nikita/Documents/Projects/SelfreplicatingColloidalClusters/");
#endif
  parameter_file_name_ = folder + "simulation_parameters.txt";

  std::ostringstream file_name_buffer_for_colloidal_particles;
  file_name_buffer_for_colloidal_particles << folder << "colloidal_particles.bin";
  simulation_file_names_[ParticleType::kColloidal] = file_name_buffer_for_colloidal_particles.str();
  std::ofstream file_for_colloidal_particles(simulation_file_names_[ParticleType::kColloidal],
                                             std::ios::binary | std::ios::out | std::ios::trunc);
  assert(file_for_colloidal_particles.is_open());
  file_for_colloidal_particles.close();

  std::ostringstream file_name_buffer_for_solvent_particles;
  file_name_buffer_for_solvent_particles << folder << "solvent_particles.bin";
  simulation_file_names_[ParticleType::kSolvent] = file_name_buffer_for_solvent_particles.str();
  std::ofstream file_for_solvent_particles(simulation_file_names_[ParticleType::kSolvent],
                                           std::ios::binary | std::ios::out | std::ios::trunc);
  assert(file_for_solvent_particles.is_open());
  file_for_solvent_particles.close();

  std::ostringstream particle_subtype_file_name_buffer;
  particle_subtype_file_name_buffer << folder << "particle_subtypes.bin";
  particle_subtype_file_name_ = particle_subtype_file_name_buffer.str();
  std::ofstream particle_subtype_file(particle_subtype_file_name_, std::ios::binary | std::ios::out | std::ios::trunc);
  assert(particle_subtype_file.is_open());
  particle_subtype_file.close();
}

BinaryObserver::~BinaryObserver()
{
  simulation_file_names_.clear();
}

void BinaryObserver::SaveSimulationParameters(int number_of_colloidal_particles,
                                              int number_of_solvent_particles,
                                              int number_of_all_particles,
                                              Real colloid_diameter,
                                              Real solvent_diameter,
                                              Real colloid_volume_fraction,
                                              Real dt)
{
  std::ofstream parameter_file(parameter_file_name_, std::ios::out | std::ios::trunc);
  assert(parameter_file.is_open());
  parameter_file << "number_of_colloidal_particles" << '\t' << number_of_colloidal_particles << std::endl;
  parameter_file << "number_of_solvent_particles" << '\t' << number_of_solvent_particles << std::endl;
  parameter_file << "number_of_all_particles" << '\t' << number_of_all_particles << std::endl;
  parameter_file << "colloid_diameter" << '\t' << colloid_diameter << std::endl;
  parameter_file << "solvent_diameter" << '\t' << solvent_diameter << std::endl;
  parameter_file << "colloid_volume_fraction" << '\t' << colloid_volume_fraction << std::endl;
  parameter_file << "x_size" << '\t' << pbc_config_.GetXSize() << std::endl;
  parameter_file << "y_size" << '\t' << pbc_config_.GetYSize() << std::endl;
  parameter_file << "z_size" << '\t' << pbc_config_.GetZSize() << std::endl;
  parameter_file << "friction_coefficient" << '\t' << kGamma << std::endl;
  parameter_file << "temperature" << '\t' << kTemperature << std::endl;
  parameter_file << "bond_strength" << '\t' << kEpsilon << std::endl;
  parameter_file << "noise_strength" << '\t' << kSigma << std::endl;
  parameter_file << "mass_of_colloid" << '\t' << kMassOfPolystyrene << std::endl;
  parameter_file << "mass_of_solvent" << '\t' << kMassOfPolyNipam << std::endl;
  parameter_file << "time_step" << '\t' << dt << std::endl;
  parameter_file.close();
}

void BinaryObserver::SaveSystemState(const std::vector<Sphere> &system_state, Real t)
{
  if (!(output_time_counter_[0] % output_time_threshold_[0]))
  {
    std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
    std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
    integration_step_timer_ = std::chrono::system_clock::now();

    std::ofstream file_for_colloidal_particles(simulation_file_names_[ParticleType::kColloidal],
                                               std::ios::binary | std::ios::out | std::ios::app);
    assert(file_for_colloidal_particles.is_open());
    std::ofstream file_for_solvent_particles(simulation_file_names_[ParticleType::kSolvent],
                                             std::ios::binary | std::ios::out | std::ios::app);
    assert(file_for_solvent_particles.is_open());

    RealOutput time_point = t;
    file_for_colloidal_particles.write((char *) &time_point, sizeof(RealOutput));
    file_for_solvent_particles.write((char *) &time_point, sizeof(RealOutput));

    static Eigen::Vector3d position(Eigen::Vector3d::Zero()), velocity(Eigen::Vector3d::Zero());
    std::vector<RealOutput> output_position, output_velocity;
    for (const Sphere &sphere : system_state)
    {
      position = sphere.GetPosition();
      output_position = std::vector<RealOutput>(position.data(), position.data() + position.rows() * position.cols());
      velocity = sphere.GetVelocity();
      output_velocity = std::vector<RealOutput>(velocity.data(), velocity.data() + velocity.rows() * velocity.cols());
      switch (sphere.GetType())
      {
        case ParticleType::kSolvent :
          file_for_solvent_particles.write((char *) &output_position[0],
                                           output_position.size() * sizeof(RealOutput));
          file_for_solvent_particles.write((char *) &output_velocity[0], output_velocity.size() * sizeof(RealOutput));
          break;
        case ParticleType::kColloidal:
          file_for_colloidal_particles.write((char *) &output_position[0],
                                             output_position.size() * sizeof(RealOutput));
          file_for_colloidal_particles.write((char *) &output_velocity[0], output_velocity.size() * sizeof(RealOutput));
          break;
        default: std::cerr << "binary observer: wrong particle type" << std::endl;
          break;
      }
    } // sphere
    file_for_colloidal_particles.close();
    file_for_solvent_particles.close();
  }
  ++output_time_counter_[0];
}

void BinaryObserver::SaveColloidalParticleSubtype(const std::vector<Sphere> &system_state, Real t)
{
  if (!(output_time_counter_[1] % output_time_threshold_[1]))
  {
    std::ofstream particle_subtype_file(particle_subtype_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(particle_subtype_file.is_open());
    RealOutput time_point = t;
    particle_subtype_file.write((char *) &time_point, sizeof(RealOutput));
    for (const Sphere &sphere : system_state)
    {
      if (ParticleType::kColloidal == sphere.GetType())
      {
        int subtype = static_cast<int>(sphere.GetSubtype());
        particle_subtype_file.write((char *) &subtype, sizeof(int));
      }
    } // sphere
    particle_subtype_file.close();
  }
  ++output_time_counter_[1];
}