// FACE-based Information

#pragma once

#include "base/Global.hpp"

class FaceField {

public:

	label cell;       // index of owner cell  
	label cell_nb;    // index of neighbour cell
	label side;       // index of owner side
	label side_nb;    // index of neighbour side

	scalar area;      // face area
	scalar distance;  // center to center_nb distance
	scalar gc;        // linear interpolation factor of cell side
	scalar gnb;       // linear interpolation factor of cell_nb side
	scalar center[3] = { 0.0 };     // face center

	scalar normal[3] = { 0.0 };     // face normal (unit) 
	scalar direction[3] = { 0.0 };  // center to center_nb vector (unit)

	label if_correction;            // non-orthogonal correction or not
	label if_repeat;                // repeated face or not

public:

	FaceField() : cell(0), cell_nb(0), side(0), side_nb(0), area(0.0), distance(0.0), gc(0.0), gnb(0.0), if_correction(0), if_repeat(0) {}
};
