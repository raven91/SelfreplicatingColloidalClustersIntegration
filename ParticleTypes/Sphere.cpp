//
// Created by Nikita Kruk on 27.02.20.
//

#include "Sphere.hpp"

int Sphere::new_periodic_subtype_ = 0;
const std::map<std::pair<ParticleSubtype, ParticleSubtype>, Real>Sphere::interaction_matrix_ =
    {{{ParticleSubtype::kA, ParticleSubtype::kA}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kB}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kC}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kBStar}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kCStar}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kAPrime}, 1.5},
     {{ParticleSubtype::kA, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kA, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kB}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kC}, 1.0},//1.5
     {{ParticleSubtype::kB, ParticleSubtype::kBStar}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kCStar}, 1.5},
     {{ParticleSubtype::kB, ParticleSubtype::kAPrime}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kBPrime}, 1.5},
     {{ParticleSubtype::kB, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kC}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kBStar}, 1.5},
     {{ParticleSubtype::kC, ParticleSubtype::kCStar}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kAPrime}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kCPrime}, 1.5},
     {{ParticleSubtype::kC, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::kBStar}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::kCStar}, 1.5},
     {{ParticleSubtype::kBStar, ParticleSubtype::kAPrime}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::kBStarPrime}, 1.5},
     {{ParticleSubtype::kBStar, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kCStar}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kAPrime}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::kCStarPrime}, 1.5},
     {{ParticleSubtype::kAPrime, ParticleSubtype::kAPrime}, 1.0},
     {{ParticleSubtype::kAPrime, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kAPrime, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kAPrime, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kAPrime, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kBPrime, ParticleSubtype::kBPrime}, 1.0},
     {{ParticleSubtype::kBPrime, ParticleSubtype::kCPrime}, 1.5},
     {{ParticleSubtype::kBPrime, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kBPrime, ParticleSubtype::kCStarPrime}, 1.5},
     {{ParticleSubtype::kCPrime, ParticleSubtype::kCPrime}, 1.0},
     {{ParticleSubtype::kCPrime, ParticleSubtype::kBStarPrime}, 1.5},
     {{ParticleSubtype::kCPrime, ParticleSubtype::kCStarPrime}, 1.0},
     {{ParticleSubtype::kBStarPrime, ParticleSubtype::kBStarPrime}, 1.0},
     {{ParticleSubtype::kBStarPrime, ParticleSubtype::kCStarPrime}, 1.5},
     {{ParticleSubtype::kCStarPrime, ParticleSubtype::kCStarPrime}, 1.0},

     {{ParticleSubtype::kA, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kB, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kC, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kBStar, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kCStar, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kAPrime, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kBPrime, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kCPrime, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kBStarPrime, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::kCStarPrime, ParticleSubtype::k0}, 1.0},
     {{ParticleSubtype::k0, ParticleSubtype::k0}, 1.0}
    };
const std::set<std::pair<ParticleSubtype, ParticleSubtype>> Sphere::favorable_interaction_pairs_ =
    {
//        {ParticleSubtype::kA, ParticleSubtype::kA},
//        {ParticleSubtype::kB, ParticleSubtype::kB},
        {ParticleSubtype::kA, ParticleSubtype::kB},
        {ParticleSubtype::kA, ParticleSubtype::kC},
        {ParticleSubtype::kA, ParticleSubtype::kBStar},
        {ParticleSubtype::kA, ParticleSubtype::kCStar},
        {ParticleSubtype::kA, ParticleSubtype::kAPrime},
        {ParticleSubtype::kB, ParticleSubtype::kC},
        {ParticleSubtype::kB, ParticleSubtype::kCStar},
        {ParticleSubtype::kB, ParticleSubtype::kBPrime},
        {ParticleSubtype::kC, ParticleSubtype::kBStar},
        {ParticleSubtype::kC, ParticleSubtype::kCPrime},
        {ParticleSubtype::kBStar, ParticleSubtype::kCStar},
        {ParticleSubtype::kBStar, ParticleSubtype::kBStarPrime},
        {ParticleSubtype::kCStar, ParticleSubtype::kCStarPrime},
        {ParticleSubtype::kAPrime, ParticleSubtype::kBPrime},
        {ParticleSubtype::kAPrime, ParticleSubtype::kCPrime},
        {ParticleSubtype::kAPrime, ParticleSubtype::kBStarPrime},
        {ParticleSubtype::kAPrime, ParticleSubtype::kCStarPrime},
        {ParticleSubtype::kBPrime, ParticleSubtype::kCPrime},
        {ParticleSubtype::kBPrime, ParticleSubtype::kCStarPrime},
        {ParticleSubtype::kCPrime, ParticleSubtype::kBStarPrime},
        {ParticleSubtype::kBStarPrime, ParticleSubtype::kCStarPrime}
    };
const std::map<ParticleSubtype, ParticleSubtype> Sphere::complementary_subtypes_ =
    {
        {ParticleSubtype::kA, ParticleSubtype::kAPrime},
        {ParticleSubtype::kB, ParticleSubtype::kBPrime},
        {ParticleSubtype::kC, ParticleSubtype::kCPrime},
        {ParticleSubtype::kBStar, ParticleSubtype::kBStarPrime},
        {ParticleSubtype::kCStar, ParticleSubtype::kCStarPrime},
        {ParticleSubtype::kAPrime, ParticleSubtype::kA},
        {ParticleSubtype::kBPrime, ParticleSubtype::kB},
        {ParticleSubtype::kCPrime, ParticleSubtype::kC},
        {ParticleSubtype::kBStarPrime, ParticleSubtype::kBStar},
        {ParticleSubtype::kCStarPrime, ParticleSubtype::kCStar},
        {ParticleSubtype::k0, ParticleSubtype::k0}
    };

Sphere::Sphere() :
    type_(ParticleType::kSolvent),
    subtype_(ParticleSubtype::k0),
    position_(Eigen::Vector3d::Zero()),
    velocity_(Eigen::Vector3d::Zero()),
    mass_(1.0),
    polymer_type_(PolymerType::kNone),
    polymer_index_(-1),
    attracted_neighbors_()
{

}

Sphere::Sphere(ParticleType type,
               ParticleSubtype subtype,
               const Eigen::Vector3d &position,
               const Eigen::Vector3d &velocity,
               Real mass,
               PolymerType polymer_type,
               int polymer_index) :
    type_(type),
    subtype_(subtype),
    position_(position),
    velocity_(velocity),
    mass_(mass),
    polymer_type_(polymer_type),
    polymer_index_(polymer_index),
    attracted_neighbors_()
{

}

Sphere::~Sphere() = default;

ParticleType Sphere::GetType() const
{
  return type_;
}

ParticleSubtype Sphere::GetSubtype() const
{
  return subtype_;
}

void Sphere::SetSubtype(ParticleSubtype subtype)
{
  subtype_ = subtype;
}

const Eigen::Vector3d &Sphere::GetPosition() const
{
  return position_;
}

void Sphere::SetPosition(const Eigen::Vector3d &position)
{
  position_ = position;
}

const Eigen::Vector3d &Sphere::GetVelocity() const
{
  return velocity_;
}

void Sphere::SetVelocity(const Eigen::Vector3d &velocity)
{
  velocity_ = velocity;
}

Real Sphere::GetMass() const
{
  return mass_;
}

PolymerType Sphere::GetPolymerType() const
{
  return polymer_type_;
}

void Sphere::SetPolymerType(PolymerType polymer_type)
{
  polymer_type_ = polymer_type;
}

int Sphere::GetPolymerIndex() const
{
  return polymer_index_;
}

void Sphere::SetPolymerIndex(int polymer_index)
{
  polymer_index_ = polymer_index;
}

ParticleSubtype Sphere::GetNewPeriodicSubtype()
{
  ParticleSubtype new_subtype = static_cast<ParticleSubtype>(new_periodic_subtype_);
  new_periodic_subtype_ = (new_periodic_subtype_ + 1) % kNumberOfPeriodicColloidSubtypes;
  return new_subtype;
}

ParticleSubtype Sphere::GetNewNoncomplementaryPeriodicSubtype()
{
  ParticleSubtype new_subtype = static_cast<ParticleSubtype>(new_periodic_subtype_);
  new_periodic_subtype_ = (new_periodic_subtype_ + 1) % (kNumberOfPeriodicColloidSubtypes / 2);
  return new_subtype;
}

ParticleSubtype Sphere::GetNewComplementaryPeriodicSubtype()
{
  ParticleSubtype
      new_subtype = static_cast<ParticleSubtype>(new_periodic_subtype_ + kNumberOfPeriodicColloidSubtypes / 2);
  new_periodic_subtype_ = (new_periodic_subtype_ + 1) % (kNumberOfPeriodicColloidSubtypes / 2);
  return new_subtype;
}

Real Sphere::GetInteractionCoefficient(ParticleSubtype subtype_1, ParticleSubtype subtype_2)
{
  if (subtype_1 < subtype_2)
  {
    return interaction_matrix_.at(std::make_pair(subtype_1, subtype_2));
  } else
  {
    return interaction_matrix_.at(std::make_pair(subtype_2, subtype_1));
  }
}

bool Sphere::IsFavorableInteraction(ParticleSubtype subtype_1, ParticleSubtype subtype_2)
{
  return
      ((favorable_interaction_pairs_.find(std::make_pair(subtype_1, subtype_2)) != favorable_interaction_pairs_.end())
          || (favorable_interaction_pairs_.find(std::make_pair(subtype_2, subtype_1))
              != favorable_interaction_pairs_.end()));
}

ParticleSubtype Sphere::QueryComplementaryParticleSubtype(ParticleSubtype subtype)
{
  return complementary_subtypes_.at(subtype);
}