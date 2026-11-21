# USTSolver

## About

USTSolver is a C++ based CFD solver project that includes the core solver library (ust_cfd) and several test cases. The following test cases are currently available, 

- 3D Hypersonic intake simulation using Gas Kinetic Scheme
- 3D Incompressible flow simulation using SIMPLE algorithm
- 2D Heat conduction simulation
- 2D Flow and energy coupling
- 1D Volumetric reaction calculation

## Prerequisites
- C++ compiler with C++17 support
- CMake 3.21+

## Build
Run the following commands from the project root (x64 Native Tools Command Prompt) :

```bash
cmake -S . -B build
cmake --build build --config Release
```

After the build completes, the executables are located in build/applications/.

CMake automatically copies the mesh files under applications/mesh/ to build/applications/, so no manual action is needed.

## Run Test Cases
Enter the applications directory inside the build directory :

```bash
cd build/applications
```

Run any test case from below :

```bash
heat_conduction.exe
flow_energy.exe
incompress_flow.exe
hypersonic_intake.exe
volumetric_reaction.exe
```

Calculation results will be stored at build/applications/mesh/ directory under the same name as the executable.
