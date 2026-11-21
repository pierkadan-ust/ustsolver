// Global Parameters

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <string>
#include <cmath>
#include <climits>
#include <vector>
#include <functional>
#include <unordered_map>

using label = long long int;
using scalar = double;
using string = std::string;
using vector = std::vector<double>;
using tensor = std::vector<std::vector<double>>;

inline string path_main;        // path to the folder that includes setup files for different cases
inline label num_dimension;     // 2D or 3D
inline label num_write;         // output mid-step solutions

inline scalar tresMomt;         // target residual for Momentum Equation
inline scalar tresP;            // target residual for Pressure Correction Equation
inline scalar tresEnergy;       // target residual for Energy Equation

inline scalar sorMomt;          // relaxation factor for Momentum Equation
inline scalar sorP;             // relaxation factor for Pressure Correction Equation
inline scalar sorEnergy;        // relaxation factor for Energy Equation

inline label max_iter_V;        // maximum number of iteration in steady Flow calculation
inline label max_iter_T;        // maximum number of iteration in steady Heat calculation

inline label max_iter_flow;     // maximum number of iteration in transient Flow calculation (inner loop)
inline scalar time_flow;        // end time for transient Flow calculation
inline scalar dt_flow;          // time step for transient Flow calculation

inline label max_iter_energy;   // maximum number of iteration in transient Energy calculation (inner loop)
inline scalar time_energy;      // end time for transient Energy calculation
inline scalar dt_energy;        // time step for transient Energy calculation

inline label if_reference_p;    // set zero pressure in the first cell or not 
inline label if_tune_decimal;   // tune the decimal digits of input mesh or not
inline scalar if_tune_scale;    // tune the scale of input mesh or not
inline scalar if_offset_axis;   // offset the origin of axis or not
inline string if_reverse_axis;  // reverse the direction of axis or not
