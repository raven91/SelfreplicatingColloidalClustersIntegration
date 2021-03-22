cmake_minimum_required(VERSION 3.15)
project(SelfreplicatingColloidalClustersIntegration)

set(CMAKE_CXX_STANDARD 17)

include_directories(/Users/nikita/Documents/mpich-install/include/ /usr/local/include/)
link_directories(/Users/nikita/Documents/mpich-install/lib/ /usr/local/lib/)
set(CMAKE_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk)

set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -O3 -flto")
set(GCC_COVERAGE_LINK_FLAGS "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

add_executable(SelfreplicatingColloidalClustersIntegration main.cpp Definitions.hpp ParticleTypes/Sphere.cpp ParticleTypes/Sphere.hpp ControlEngines/SimulationEngine.cpp ControlEngines/SimulationEngine.hpp BoundaryConfiguration/PeriodicBoundaryConditions.cpp BoundaryConfiguration/PeriodicBoundaryConditions.hpp Steppers/RungeKutta2Stepper.cpp Steppers/RungeKutta2Stepper.hpp DynamicalSystems/SelfreplicationForOctahedron.cpp DynamicalSystems/SelfreplicationForOctahedron.hpp Observers/BinaryObserver.cpp Observers/BinaryObserver.hpp Steppers/StochasticEulerStepper.hpp Steppers/VelocityVerletStepper.cpp Steppers/VelocityVerletStepper.hpp DynamicalSystems/SelfreplicationForOctahedronNondimensionalized.cpp DynamicalSystems/SelfreplicationForOctahedronNondimensionalized.hpp Definitions.cpp)