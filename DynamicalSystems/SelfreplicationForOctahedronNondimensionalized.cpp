//
// Created by Nikita Kruk on 04.03.20.
//

#include "SelfreplicationForOctahedronNondimensionalized.hpp"

#include <vector>
#include <iostream>

#include <eigen3/Eigen/Dense>

SelfreplicationForOctahedronNondimensionalized::SelfreplicationForOctahedronNondimensionalized(
    PeriodicBoundaryConditions &pbc_config,
    int number_of_all_particles,
    Real colloid_diameter) :
    pbc_config_(pbc_config),
    mersenne_twister_generator_(std::random_device{}()),
    unif01_(0.0, 1.0),
    number_of_all_particles_(number_of_all_particles),
    colloid_diameter_(colloid_diameter)
{
  const Real maximal_cutoff_radius = 1.5 * colloid_diameter_;
  num_subcells_x_ = int(pbc_config_.GetXSize() / maximal_cutoff_radius);
  num_subcells_y_ = int(pbc_config_.GetYSize() / maximal_cutoff_radius);
  num_subcells_z_ = int(pbc_config_.GetZSize() / maximal_cutoff_radius);
//  pre_linked_list_ = std::vector<int>(number_of_all_particles_, 0);
  linked_list_ = std::vector<std::vector<int>>(num_subcells_x_ * num_subcells_y_ * num_subcells_z_, std::vector<int>());

  // all neighbors
  /*for (int ix = -1; ix <= 1; ++ix)
  {
    for (int iy = -1; iy <= 1; ++iy)
    {
      for (int iz = -1; iz <= 1; ++iz)
      {
        neighboring_cells_.push_back({ix, iy, iz});
      } // iz
    } // iy
  } // ix*/
  // half of neighbors
  neighboring_cells_ =
      {
          {-1, -1, 1}, {0, -1, 1}, {1, -1, 1},
          {-1, 0, 1}, {0, 0, 1}, {1, 0, 1},
          {-1, 1, 1}, {0, 1, 1}, {1, 1, 1},
          {-1, 1, 0}, {0, 1, 0}, {1, 1, 0},
          {0, 0, 0}, {1, 0, 0}
      };

  // Verlet neighbor list
  verlet_distance_ = 2.0 * maximal_cutoff_radius;
  should_update_lists_ = true;
  accumulated_displacement_ = 0.0;
  verlet_list_ = std::vector<std::vector<int>>(number_of_all_particles_, std::vector<int>());
}

SelfreplicationForOctahedronNondimensionalized::~SelfreplicationForOctahedronNondimensionalized()
{
  neighboring_cells_.clear();
  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list
  linked_list_.clear();
  for (std::vector<int> &verlet_list : verlet_list_)
  {
    verlet_list.clear();
  } // verlet_list
  verlet_list_.clear();
}

void SelfreplicationForOctahedronNondimensionalized::EvaluateRhsWithAllPairs(std::vector<Sphere> &system_state,
                                                                             const std::vector<Sphere> &k_prev,
                                                                             std::vector<Sphere> &k_next,
                                                                             Real k_coef,
                                                                             Real dt_method_substep,
                                                                             Real dt_numerical_scheme)
{
  std::vector<Eigen::Vector3d> dissipative_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      random_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      conservative_forces(number_of_all_particles_, Eigen::Vector3d::Zero());
  Eigen::Vector3d dissipative_force(Eigen::Vector3d::Zero()), random_force(Eigen::Vector3d::Zero()),
      conservative_force(Eigen::Vector3d::Zero());
  static Eigen::Vector3d position_i, velocity_i;
  static Eigen::Vector3d position_j, velocity_j;
  static Eigen::Vector3d position_ij, velocity_ij, position_ij_normalized;
  Real r_ij = 0.0, r_ij_squared = 0.0, v_ij = 0.0, theta_ij = 0.0;
  Real omega_ij_r = 0.0, omega_ij_d = 0.0;
  ParticleType type_i, type_j;
  ParticleSubtype subtype_i, subtype_j;

  static std::vector<Sphere> rk_system_state(system_state.begin(), system_state.end());
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    rk_system_state[i] = system_state[i]; // to transfer additional variables
    rk_system_state[i].SetPosition(
        system_state[i].GetPosition() + k_coef * dt_method_substep * k_prev[i].GetPosition());
    rk_system_state[i].SetVelocity(
        system_state[i].GetVelocity() + k_coef * dt_method_substep * k_prev[i].GetVelocity());
  } // i

  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    type_i = rk_system_state[i].GetType();
    subtype_i = rk_system_state[i].GetSubtype();
    for (int j = i + 1; j < number_of_all_particles_; ++j)
    {
      position_j = rk_system_state[j].GetPosition();
      velocity_j = rk_system_state[j].GetVelocity();
      type_j = rk_system_state[j].GetType();
      subtype_j = rk_system_state[j].GetSubtype();
      pbc_config_.ClassAEffectiveParticleDistance(position_j, position_i, position_ij);
      r_ij_squared = position_ij.squaredNorm();
      r_ij = std::sqrt(r_ij_squared);
      velocity_ij = velocity_i - velocity_j;

      position_ij_normalized = position_ij.normalized();
      omega_ij_r = SoftWeightingFunction(r_ij, type_i, type_j);
      omega_ij_d = omega_ij_r * omega_ij_r;
      dissipative_force =
          -kGamma * omega_ij_d * position_ij_normalized.dot(velocity_ij) * position_ij_normalized;
      dissipative_forces[i] += dissipative_force;
      dissipative_forces[j] -= dissipative_force;

      theta_ij = std::sqrt(3.0) * (2.0 * unif01_(mersenne_twister_generator_) - 1.0);
      random_force = kSigma * omega_ij_r * theta_ij * position_ij_normalized / std::sqrt(dt_numerical_scheme);
      random_forces[i] += random_force;
      random_forces[j] -= random_force;

      conservative_force =
          ConservativeInteractionForce(r_ij, type_i, type_j, subtype_i, subtype_j, system_state[i], system_state[j])
              * position_ij_normalized;
      conservative_forces[i] += conservative_force;
      conservative_forces[j] -= conservative_force;
    } // j
  } // i

  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    k_next[i].SetPosition(velocity_i);
    k_next[i].SetVelocity(
        (dissipative_forces[i] + random_forces[i] + conservative_forces[i]) / rk_system_state[i].GetMass());
  } // i
}

void SelfreplicationForOctahedronNondimensionalized::EvaluateRhsWithLinkedList(std::vector<Sphere> &system_state,
                                                                               const std::vector<Sphere> &k_prev,
                                                                               std::vector<Sphere> &k_next,
                                                                               Real k_coef,
                                                                               Real dt_method_substep,
                                                                               Real dt_numerical_scheme)
{
  std::vector<Eigen::Vector3d> dissipative_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      random_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      conservative_forces(number_of_all_particles_, Eigen::Vector3d::Zero());
  Eigen::Vector3d dissipative_force(Eigen::Vector3d::Zero()), random_force(Eigen::Vector3d::Zero()),
      conservative_force(Eigen::Vector3d::Zero());
  static Eigen::Vector3d position_i, velocity_i;
  static Eigen::Vector3d position_j, velocity_j;
  static Eigen::Vector3d position_ij, velocity_ij, position_ij_normalized;
  Real r_ij = 0.0, r_ij_squared = 0.0, v_ij = 0.0, theta_ij = 0.0;
  Real omega_ij_r = 0.0, omega_ij_d = 0.0;
  ParticleType type_i, type_j;
  ParticleSubtype subtype_i, subtype_j;

  static std::vector<Sphere> rk_system_state(system_state.begin(), system_state.end());
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    rk_system_state[i] = system_state[i]; // to transfer additional variables
    rk_system_state[i].SetPosition(
        system_state[i].GetPosition() + k_coef * dt_method_substep * k_prev[i].GetPosition());
    rk_system_state[i].SetVelocity(
        system_state[i].GetVelocity() + k_coef * dt_method_substep * k_prev[i].GetVelocity());
  } // i
  // to construct the linked list, particles must strictly remain inside the simulation box
  pbc_config_.ApplyPeriodicBoundaryConditions(rk_system_state);

  for (std::vector<int> &linked_list : linked_list_)
  {
    linked_list.clear();
  } // linked_list

  // construct a linked list
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0, i_cell_z = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0, j_cell_z = 0;
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();

    i_cell_x = int(position_i.x() / pbc_config_.GetXSize() * num_subcells_x_);
    i_cell_y = int(position_i.y() / pbc_config_.GetYSize() * num_subcells_y_);
    i_cell_z = int(position_i.z() / pbc_config_.GetZSize() * num_subcells_z_);
    utilities::ThreeDimIdxToOneDimIdx(i_cell_x,
                                      i_cell_y,
                                      i_cell_z,
                                      i_cell,
                                      num_subcells_x_,
                                      num_subcells_y_,
                                      num_subcells_z_);

//    pre_linked_list_[i] = i_cell;
    linked_list_[i_cell].push_back(i);
  } // i
//  for (int i = 0; i < number_of_all_particles_; ++i)
//  {
//    linked_list_[pre_linked_list_[i]].push_back(i);
//  } // i

// loop using the linked list
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    type_i = rk_system_state[i].GetType();
    subtype_i = rk_system_state[i].GetSubtype();

    i_cell_x = int(position_i.x() / pbc_config_.GetXSize() * num_subcells_x_);
    i_cell_y = int(position_i.y() / pbc_config_.GetYSize() * num_subcells_y_);
    i_cell_z = int(position_i.z() / pbc_config_.GetZSize() * num_subcells_z_);
    utilities::ThreeDimIdxToOneDimIdx(i_cell_x,
                                      i_cell_y,
                                      i_cell_z,
                                      i_cell,
                                      num_subcells_x_,
                                      num_subcells_y_,
                                      num_subcells_z_);
    for (auto &neighboring_cell : neighboring_cells_)
    {
      j_cell_x = i_cell_x + neighboring_cell[0];
      j_cell_y = i_cell_y + neighboring_cell[1];
      j_cell_z = i_cell_z + neighboring_cell[2];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y, j_cell_z);
      utilities::ThreeDimIdxToOneDimIdx(j_cell_x,
                                        j_cell_y,
                                        j_cell_z,
                                        j_cell,
                                        num_subcells_x_,
                                        num_subcells_y_,
                                        num_subcells_z_);

      for (int j : linked_list_[j_cell])
      {
        if ((i_cell != j_cell) || (j > i))
        {
          position_j = rk_system_state[j].GetPosition();
          velocity_j = rk_system_state[j].GetVelocity();
          type_j = rk_system_state[j].GetType();
          subtype_j = rk_system_state[j].GetSubtype();

          pbc_config_.ClassAEffectiveParticleDistance(position_j, position_i, position_ij);
          r_ij_squared = position_ij.squaredNorm();
          r_ij = std::sqrt(r_ij_squared);
          velocity_ij = velocity_i - velocity_j;

          position_ij_normalized = position_ij.normalized();
          omega_ij_r = SoftWeightingFunction(r_ij, type_i, type_j);
          omega_ij_d = omega_ij_r * omega_ij_r;
          dissipative_force =
              -kGamma * omega_ij_d * position_ij_normalized.dot(velocity_ij) * position_ij_normalized;
          dissipative_forces[i] += dissipative_force;
          dissipative_forces[j] -= dissipative_force;

          theta_ij = std::sqrt(3.0) * (2.0 * unif01_(mersenne_twister_generator_) - 1.0);
          random_force = kSigma * omega_ij_r * theta_ij * position_ij_normalized / std::sqrt(dt_numerical_scheme);
          random_forces[i] += random_force;
          random_forces[j] -= random_force;

          conservative_force =
              ConservativeInteractionForce(r_ij, type_i, type_j, subtype_i, subtype_j, system_state[i], system_state[j])
                  * position_ij_normalized;
          conservative_forces[i] += conservative_force;
          conservative_forces[j] -= conservative_force;
        }
      } // j
    } // neighboring_cell
  } // i

  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    k_next[i].SetPosition(velocity_i);
    k_next[i].SetVelocity(
        (dissipative_forces[i] + random_forces[i] + conservative_forces[i]) / rk_system_state[i].GetMass());
  } // i
}

void SelfreplicationForOctahedronNondimensionalized::EvaluateRhsWithVerletNeighborList(std::vector<Sphere> &system_state,
                                                                                       const std::vector<Sphere> &k_prev,
                                                                                       std::vector<Sphere> &k_next,
                                                                                       Real k_coef,
                                                                                       Real dt_method_substep,
                                                                                       Real dt_numerical_scheme)
{
  static Eigen::Vector3d position_i, velocity_i;
  static Eigen::Vector3d position_j, velocity_j;
  static Eigen::Vector3d position_ij, velocity_ij, position_ij_normalized;
  Real r_ij = 0.0, r_ij_squared = 0.0, v_ij = 0.0, theta_ij = 0.0;
  Real omega_ij_r = 0.0, omega_ij_d = 0.0;
  ParticleType type_i, type_j;
  ParticleSubtype subtype_i, subtype_j;

  static std::vector<Sphere> rk_system_state(system_state.begin(), system_state.end());
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    rk_system_state[i] = system_state[i]; // to transfer additional variables
    rk_system_state[i].SetPosition(
        system_state[i].GetPosition() + k_coef * dt_method_substep * k_prev[i].GetPosition());
    rk_system_state[i].SetVelocity(
        system_state[i].GetVelocity() + k_coef * dt_method_substep * k_prev[i].GetVelocity());
  } // i
  // to construct the linked list, particles must strictly remain inside the simulation box
  pbc_config_.ApplyPeriodicBoundaryConditions(rk_system_state);

  // if it is time to update the linked list and the Verlet list
  if (should_update_lists_)
  {
    for (auto &linked_list_entry : linked_list_)
    {
      linked_list_entry.clear();
    } // linked_list_entry

    // construct a linked list
    int i_cell = 0, i_cell_x = 0, i_cell_y = 0, i_cell_z = 0;
    int j_cell = 0, j_cell_x = 0, j_cell_y = 0, j_cell_z = 0;
    for (int i = 0; i < number_of_all_particles_; ++i)
    {
      position_i = rk_system_state[i].GetPosition();
      velocity_i = rk_system_state[i].GetVelocity();

      i_cell_x = int(position_i.x() / pbc_config_.GetXSize() * num_subcells_x_);
      i_cell_y = int(position_i.y() / pbc_config_.GetYSize() * num_subcells_y_);
      i_cell_z = int(position_i.z() / pbc_config_.GetZSize() * num_subcells_z_);
      utilities::ThreeDimIdxToOneDimIdx(i_cell_x,
                                        i_cell_y,
                                        i_cell_z,
                                        i_cell,
                                        num_subcells_x_,
                                        num_subcells_y_,
                                        num_subcells_z_);

//    pre_linked_list_[i] = i_cell;
      linked_list_[i_cell].push_back(i);
    } // i
//  for (int i = 0; i < number_of_all_particles_; ++i)
//  {
//    linked_list_[pre_linked_list_[i]].push_back(i);
//  } // i

    // construct the Verlet list using the linked list
    for (auto &verlet_list_entry : verlet_list_)
    {
      verlet_list_entry.clear();
    } // verlet_list_entry

    for (int i = 0; i < number_of_all_particles_; ++i)
    {
      position_i = rk_system_state[i].GetPosition();
      velocity_i = rk_system_state[i].GetVelocity();

      i_cell_x = int(position_i.x() / pbc_config_.GetXSize() * num_subcells_x_);
      i_cell_y = int(position_i.y() / pbc_config_.GetYSize() * num_subcells_y_);
      i_cell_z = int(position_i.z() / pbc_config_.GetZSize() * num_subcells_z_);
      utilities::ThreeDimIdxToOneDimIdx(i_cell_x,
                                        i_cell_y,
                                        i_cell_z,
                                        i_cell,
                                        num_subcells_x_,
                                        num_subcells_y_,
                                        num_subcells_z_);

      for (auto &neighboring_cell : neighboring_cells_)
      {
        j_cell_x = i_cell_x + neighboring_cell[0];
        j_cell_y = i_cell_y + neighboring_cell[1];
        j_cell_z = i_cell_z + neighboring_cell[2];
        AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y, j_cell_z);
        utilities::ThreeDimIdxToOneDimIdx(j_cell_x,
                                          j_cell_y,
                                          j_cell_z,
                                          j_cell,
                                          num_subcells_x_,
                                          num_subcells_y_,
                                          num_subcells_z_);

        for (int j : linked_list_[j_cell])
        {
          if ((i_cell != j_cell) || (j > i))
          {
            position_j = rk_system_state[j].GetPosition();
            velocity_j = rk_system_state[j].GetVelocity();

            pbc_config_.ClassAEffectiveParticleDistance(position_j, position_i, position_ij);
            r_ij_squared = position_ij.squaredNorm();
            if (r_ij_squared <= verlet_distance_ * verlet_distance_)
            {
              verlet_list_[i].push_back(j);
            }
          }
        } // j
      } // neighboring_cell
    } // i
    should_update_lists_ = false;
  } // if should_update_lists_

  // loop using the Verlet list
  // todo: make static?
  std::vector<Eigen::Vector3d> dissipative_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      random_forces(number_of_all_particles_, Eigen::Vector3d::Zero()),
      conservative_forces(number_of_all_particles_, Eigen::Vector3d::Zero());
  Eigen::Vector3d dissipative_force(Eigen::Vector3d::Zero()), random_force(Eigen::Vector3d::Zero()),
      conservative_force(Eigen::Vector3d::Zero());

  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    type_i = rk_system_state[i].GetType();
    subtype_i = rk_system_state[i].GetSubtype();

    for (int j : verlet_list_[i])
    {
      position_j = rk_system_state[j].GetPosition();
      velocity_j = rk_system_state[j].GetVelocity();
      type_j = rk_system_state[j].GetType();
      subtype_j = rk_system_state[j].GetSubtype();

      pbc_config_.ClassAEffectiveParticleDistance(position_j, position_i, position_ij);
      r_ij_squared = position_ij.squaredNorm();
      r_ij = std::sqrt(r_ij_squared);
      velocity_ij = velocity_i - velocity_j;

      position_ij_normalized = position_ij.normalized();
      omega_ij_r = SoftWeightingFunction(r_ij, type_i, type_j);
      omega_ij_d = omega_ij_r * omega_ij_r;
      dissipative_force =
          -kGamma * omega_ij_d * position_ij_normalized.dot(velocity_ij) * position_ij_normalized;
      dissipative_forces[i] += dissipative_force;
      dissipative_forces[j] -= dissipative_force;

      theta_ij = std::sqrt(3.0) * (2.0 * unif01_(mersenne_twister_generator_) - 1.0);
      random_force = kSigma * omega_ij_r * theta_ij * position_ij_normalized / std::sqrt(dt_numerical_scheme);
      random_forces[i] += random_force;
      random_forces[j] -= random_force;

      conservative_force =
          ConservativeInteractionForce(r_ij, type_i, type_j, subtype_i, subtype_j, system_state[i], system_state[j])
              * position_ij_normalized;
      conservative_forces[i] += conservative_force;
      conservative_forces[j] -= conservative_force;
    } // j
  } // i

  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    position_i = rk_system_state[i].GetPosition();
    velocity_i = rk_system_state[i].GetVelocity();
    k_next[i].SetPosition(velocity_i);
    k_next[i].SetVelocity(
        (dissipative_forces[i] + random_forces[i] + conservative_forces[i]) / rk_system_state[i].GetMass());
  } // i

  CalculateMaxDisplacement(k_next, dt_method_substep);
}

Real SelfreplicationForOctahedronNondimensionalized::SoftWeightingFunction(Real r_ij,
                                                                           ParticleType type_i,
                                                                           ParticleType type_j) const
{
  static const Real r_cut_cc = 1.5 * colloid_diameter_;
  static const Real r_cut_ss = 0.5 * colloid_diameter_;
  static const Real r_cut_cs = 1.0 * colloid_diameter_;

  Real r_cut = 0.0;
  if (type_i == type_j)
  {
    if (type_i == ParticleType::kSolvent)
    {
      r_cut = r_cut_ss;
    } else if (type_i == ParticleType::kColloidal)
    {
      r_cut = r_cut_cc;
    }
  } else
  {
    r_cut = r_cut_cs;
  }

  if (r_ij <= r_cut)
  {
    return 1.0 - r_ij / r_cut;
  } else
  {
    return 0.0;
  }
}

Real SelfreplicationForOctahedronNondimensionalized::ConservativeInteractionForce(Real r_ij,
                                                                                  ParticleType type_i,
                                                                                  ParticleType type_j,
                                                                                  ParticleSubtype subtype_i,
                                                                                  ParticleSubtype subtype_j,
                                                                                  Sphere &sphere_i,
                                                                                  Sphere &sphere_j) const
{
  static const Real r_m = 1.0 * colloid_diameter_;
  static const Real r_cc_favorable = 1.05 * colloid_diameter_;
  static const Real r_cc_unfavorable = 1.0 * colloid_diameter_;
  static const Real r_cut_ss = 0.5 * colloid_diameter_;
  static const Real r_cut_cs = 0.75 * colloid_diameter_;

  Real force = 0.0;
  if ((type_i == type_j) && (type_i == ParticleType::kColloidal))
  {
    if (r_ij <= r_cc_unfavorable)
    {
      if (subtype_i == ParticleSubtype::k0)
      {
        sphere_i.SetSubtype(Sphere::QueryComplementaryParticleSubtype(subtype_j));
        subtype_i = sphere_i.GetSubtype();
//      std::cout << "[" << static_cast<int>(subtype_i) << "][" << static_cast<int>(subtype_j) << "]" << std::endl;
      } else if (subtype_j == ParticleSubtype::k0)
      {
        sphere_j.SetSubtype(Sphere::QueryComplementaryParticleSubtype(subtype_i));
        subtype_j = sphere_j.GetSubtype();
//      std::cout << "[" << static_cast<int>(subtype_i) << "][" << static_cast<int>(subtype_j) << "]" << std::endl;
      }
    }

    Real interaction_coefficient = Sphere::GetInteractionCoefficient(subtype_i, subtype_j);
    if (Sphere::IsFavorableInteraction(subtype_i, subtype_j))
    {
      if ((sphere_i.GetPolymerType() != PolymerType::kNone) || (sphere_j.GetPolymerType() != PolymerType::kNone))
      {
        if (sphere_i.GetPolymerIndex() == sphere_j.GetPolymerIndex())
        { // if both colloids belong to the same polymer
          interaction_coefficient = 5.0;
        } else if ((sphere_i.GetPolymerIndex() == -1) || (sphere_j.GetPolymerIndex() == -1))
        { // if one of the colloids belongs to a polymer but the other does not
          interaction_coefficient = 2.0;
        }
      }

      if (r_ij <= r_cc_favorable)
      {
        Real r_ij_48 = std::pow(r_m / r_ij, 48.0);
        force = 96.0 * interaction_coefficient * 1.0 * r_ij_48 * (r_ij_48 - 1.0)
            / r_ij; // kEpsilon disappears during nondimensionalization
      }
    } else
    {
      if (r_ij <= r_cc_unfavorable)
      {
        Real r_ij_48 = std::pow(r_m / r_ij, 48.0);
        force = 96.0 * interaction_coefficient * 1.0 * r_ij_48 * (r_ij_48 - 1.0)
            / r_ij; // kEpsilon disappears during nondimensionalization
      }
    }
  } else
  {
    if ((type_i == type_j) && (type_i == ParticleType::kSolvent))
    {
      if (r_ij <= r_cut_ss)
      {
        force = kRepulsionForce * (1.0 - r_ij / r_cut_ss);
      }
    } else
    {
      if (r_ij <= r_cut_cs)
      {
        force = kRepulsionForce * (1.0 - r_ij / r_cut_cs);
      }
    }
  }

  return force;
}

void SelfreplicationForOctahedronNondimensionalized::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x,
                                                                                               int &cell_y,
                                                                                               int &cell_z) const
{
  if (cell_x < 0)
  {
    cell_x = num_subcells_x_ - 1;
  } else if (cell_x >= num_subcells_x_)
  {
    cell_x = 0;
  }

  if (cell_y < 0)
  {
    cell_y = num_subcells_y_ - 1;
  } else if (cell_y >= num_subcells_y_)
  {
    cell_y = 0;
  }

  if (cell_z < 0)
  {
    cell_z = num_subcells_z_ - 1;
  } else if (cell_z >= num_subcells_z_)
  {
    cell_z = 0;
  }
}

void SelfreplicationForOctahedronNondimensionalized::CalculateMaxDisplacement(const std::vector<Sphere> &derivative,
                                                                              Real dt)
{
  Real max2 = 0.0, dist2 = 0.0;
  Eigen::Vector3d dr = Eigen::Vector3d::Zero();
  for (int i = 0; i < number_of_all_particles_; ++i)
  {
    dr = derivative[i].GetPosition() * dt;
    dist2 = dr.squaredNorm();
    if (max2 < dist2)
    {
      max2 = dist2;
    }
  } // i
  accumulated_displacement_ += std::sqrt(max2);
  if (accumulated_displacement_ > 0.5 * verlet_distance_)
  {
    accumulated_displacement_ = 0.0;
    should_update_lists_ = true;
  }
}