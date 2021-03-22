//
// Created by Nikita Kruk on 04.03.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SELFREPLICATIONFOROCTAHEDRONNONDIMENSIONALIZED_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SELFREPLICATIONFOROCTAHEDRONNONDIMENSIONALIZED_HPP

#include "../Definitions.hpp"
#include "../BoundaryConfiguration/PeriodicBoundaryConditions.hpp"
#include "../ParticleTypes/Sphere.hpp"

#include <random>
#include <vector>

class SelfreplicationForOctahedronNondimensionalized
{
 public:

  explicit SelfreplicationForOctahedronNondimensionalized(PeriodicBoundaryConditions &pbc_config,
                                                          int number_of_all_particles,
                                                          Real colloid_diameter);
  ~SelfreplicationForOctahedronNondimensionalized();

  void EvaluateRhsWithAllPairs(std::vector<Sphere> &system_state,
                               const std::vector<Sphere> &k_prev,
                               std::vector<Sphere> &k_next,
                               Real k_coef,
                               Real dt_method_substep,
                               Real dt_numerical_scheme);
  void EvaluateRhsWithLinkedList(std::vector<Sphere> &system_state,
                                 const std::vector<Sphere> &k_prev,
                                 std::vector<Sphere> &k_next,
                                 Real k_coef,
                                 Real dt_method_substep,
                                 Real dt_numerical_scheme);
  void EvaluateRhsWithVerletNeighborList(std::vector<Sphere> &system_state,
                                         const std::vector<Sphere> &k_prev,
                                         std::vector<Sphere> &k_next,
                                         Real k_coef,
                                         Real dt_method_substep,
                                         Real dt_numerical_scheme);

 private:

  PeriodicBoundaryConditions &pbc_config_;
  std::mt19937 mersenne_twister_generator_;
  std::uniform_real_distribution<Real> unif01_;

  int number_of_all_particles_;
  Real colloid_diameter_;

  // Cell Subdivision Routine
  int num_subcells_x_;
  int num_subcells_y_;
  int num_subcells_z_;
//  std::vector<int> pre_linked_list_;
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  // Verlet Neighbor List Routine
  Real verlet_distance_;
  std::vector<std::vector<int>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;

  [[nodiscard]] Real SoftWeightingFunction(Real r_ij, ParticleType type_i, ParticleType type_j) const;
  [[nodiscard]] Real ConservativeInteractionForce(Real r_ij,
                                                  ParticleType type_i,
                                                  ParticleType type_j,
                                                  ParticleSubtype subtype_i,
                                                  ParticleSubtype subtype_j,
                                                  Sphere &sphere_i,
                                                  Sphere &sphere_j) const;
  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y, int &cell_z) const;
  void CalculateMaxDisplacement(const std::vector<Sphere> &derivative, Real dt);

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SELFREPLICATIONFOROCTAHEDRONNONDIMENSIONALIZED_HPP
