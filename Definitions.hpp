//
// Created by Nikita Kruk on 27.02.20.
//

#ifndef SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_DEFINITIONS_HPP
#define SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_DEFINITIONS_HPP

//#define BCS_CLUSTER
#define NONDIMENSIONALIZED

#include <cmath>

typedef double Real;
typedef float RealOutput;

#if defined(NONDIMENSIONALIZED)
// Independent units
const Real kGamma = 10.0;
const Real kTemperature = 1.0;
const Real kMassOfPolystyrene = 1.0;
const Real kMassOfPolyNipam = 1.0;//5.6731e-4;
const Real kColloidDiameter = 1.0;
const Real kSolventDiameter = 1.0;//0.08;
const Real kRepulsionForce = 25.0;

// Dependent units
const Real kEpsilon = kTemperature / 0.14;
const Real kSigma = std::sqrt(2.0 * kGamma / kEpsilon);
//const Real kGamma = kSigma * kSigma / 2.0;
#else
// Conversion between units
const Real kMicrometerPerMeter = 1e+6;
const Real kGramPerKilogram = 1e+3;

// Independent units
const Real kBoltzmannConstant =
    1.3806503e-23 * (kMicrometerPerMeter * kMicrometerPerMeter * kGramPerKilogram); // [(um^2 g)/(s^2 K)]
const Real kGamma = 10.0 * 0.52 * 1e-15 * kGramPerKilogram; // [g/s]
const Real kT = 22.0 + 273.15;//1.0 / kBoltzmannConstant;//22.0 + 273.15; // [K]
const Real kMassOfPolystyrene = 0.52 * 1e-15 * kGramPerKilogram; // [g]
const Real kMassOfPolyNipam = 2.95 * 1e-19 * kGramPerKilogram; // [g]
const Real kColloidDiameter = 1.0 * 1e-6 * kMicrometerPerMeter; // [um]
const Real kSolventDiameter = 80.0 * 1e-9 * kMicrometerPerMeter; // [um]
const Real kRepulsionForce = 25.0 * kGramPerKilogram * kMicrometerPerMeter; // [g um/s^2]

// Dependent units
const Real kEpsilon = 4.0 * kBoltzmannConstant * kT; // [um^2 g / s^2]
const Real kSigmaColloidal = std::sqrt(2.0 * kGamma * kBoltzmannConstant * kT); // [g um / s^1.5]
const Real kSigmaSolvent = std::sqrt(2.0 * kGamma * kBoltzmannConstant * kT); // [g um / s^3]
#endif

// Characteristic scales for dynamical systems
//const Real kTimeScale = 1.0;//0.5 * kColloidDiameter * std::sqrt(kMassOfPolystyrene / kEpsilon);
//const Real kSpatialScale = 1.0;//kColloidDiameter;
//const Real kVelocityScale = 1.0; // kSpatialScale / kTimeScale;
//const Real kTemperatureScale = 1.0;//kEpsilon / kBoltzmannConstant;
//const Real kEnergyScale = 1.0;//kBoltzmannConstant * 1.0;

namespace utilities
{
  inline void ThreeDimIdxToOneDimIdx(int x, int y, int z, int &idx, int num_cells_x, int num_cells_y, int num_cells_phi)
  {
    // the winding order is x->y->z
    idx = x + num_cells_x * (y + num_cells_y * z);
  }
} // namespace utilities

#endif //SELFREPLICATINGCOLLOIDALCLUSTERSINTEGRATION_DEFINITIONS_HPP
