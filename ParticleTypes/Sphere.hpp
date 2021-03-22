//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SPHERE_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SPHERE_HPP

#include "../Definitions.hpp"

#include <utility> // std::pair
#include <map>
#include <set>
#include <vector>

#include <eigen3/Eigen/Dense>

enum class ParticleType
{
  kColloidal, kSolvent
};

const int kNumberOfPeriodicColloidSubtypes = 2;
enum class ParticleSubtype
{
  // periodic subtypes must be declared first
  kA = 0, kB, kC, kBStar, kCStar,
  kAPrime, kBPrime, kCPrime, kBStarPrime, kCStarPrime,
  k0
};

enum class PolymerType
{
  kOctahedron, kCatalyst, kNone
};

const std::set<ParticleSubtype> kParentOctahedronParticles =
    {ParticleSubtype::kA, ParticleSubtype::kB, ParticleSubtype::kC,
     ParticleSubtype::kAPrime, ParticleSubtype::kBPrime, ParticleSubtype::kCPrime};
const std::set<ParticleSubtype> kParentCatalystParticles =
    {ParticleSubtype::kBStar, ParticleSubtype::kCStar, ParticleSubtype::kBStarPrime, ParticleSubtype::kCStarPrime};

class Sphere
{
 public:

  Sphere();
  Sphere(ParticleType type,
         ParticleSubtype subtype,
         const Eigen::Vector3d &position,
         const Eigen::Vector3d &velocity,
         Real mass,
         PolymerType polymer_type,
         int polymer_index);
  ~Sphere();

  ParticleType GetType() const;
  ParticleSubtype GetSubtype() const;
  void SetSubtype(ParticleSubtype subtype);
  [[nodiscard]] const Eigen::Vector3d &GetPosition() const;
  void SetPosition(const Eigen::Vector3d &position);
  [[nodiscard]] const Eigen::Vector3d &GetVelocity() const;
  void SetVelocity(const Eigen::Vector3d &velocity);
  Real GetMass() const;
  PolymerType GetPolymerType() const;
  void SetPolymerType(PolymerType polymer_type);
  int GetPolymerIndex() const;
  void SetPolymerIndex(int polymer_index);

  static ParticleSubtype GetNewPeriodicSubtype();
  static ParticleSubtype GetNewNoncomplementaryPeriodicSubtype();
  static ParticleSubtype GetNewComplementaryPeriodicSubtype();
  static Real GetInteractionCoefficient(ParticleSubtype subtype_1, ParticleSubtype subtype_2);
  static bool IsFavorableInteraction(ParticleSubtype subtype_1, ParticleSubtype subtype_2);
  static ParticleSubtype QueryComplementaryParticleSubtype(ParticleSubtype subtype);

 private:

  ParticleType type_;
  ParticleSubtype subtype_;
  Eigen::Vector3d position_;
  Eigen::Vector3d velocity_;
  Real mass_;
  PolymerType polymer_type_;
  int polymer_index_;
  std::vector<int> attracted_neighbors_;

  static int new_periodic_subtype_;
  static const std::map<std::pair<ParticleSubtype, ParticleSubtype>, Real> interaction_matrix_;
  static const std::set<std::pair<ParticleSubtype, ParticleSubtype>> favorable_interaction_pairs_;
  static const std::map<ParticleSubtype, ParticleSubtype> complementary_subtypes_;

};

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_SPHERE_HPP
