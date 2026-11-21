// Calculation Domain
// store mesh, boundary location and geometry information for each zone
// interpret mesh information from different input formats

// DATA STRUCTURE:
//-----------------------------------------------------------------------------------------------------------------------------------------------------
//                         C(2)                 
//            n1(1)   n1(1)---f0----n0(3)     * node_local(Global), face_local, cell(Global)
//           /   \      \           /  C(4)         
//   BC    f1     f0     \   C(1)  /          <nodes> : { n(0), n(1), n(2), n(3) }
//         /  C(0) \     f1      f2           <cells> : { C(0)[ n0:n(0), n1:n(1), n2:n(2) ], C(1)[ n0:n(3), n1:n(1), n2:n(0) ] } 
//        /         \      \     /            <faces> : { C(0)[ f0:C(1), f1:BC,   f2:C(3) ], C(1)[ f0:C(2), f1:R,    f2:C(4) ] }
//    n2(2)---f2---n0(0)    n2(0)                                 BC_type = -1,-2, ...                    R = Repeated
//               C(3)                         <types> : { C(0): 1, C(1): 2 }  1 - Quad, 2 - Tri, 3 - Hex(C3D8), 4 - Tetra(C3D4), ... 
//-----------------------------------------------------------------------------------------------------------------------------------------------------

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "base/Tool.hpp"
#include "field/FaceField.hpp"

class Zone {

public:

	// mesh information 
	Array<FaceField> faces;
	Array<scalar> nodes;
	Array<label> cells;
	label num_faces;
	label num_nodes;
	label num_cells;

	Array<label> types;
	label type_mesh;
	label num_face;
	label num_vertex;
	label num_vertex_face;

	label id;
	string name;
	string type_input;
	label num_coordinate;

	// geometry information
	Array<scalar> volume;
	Array<scalar> center;

	// boundary information
	label num_surface;
	Array<string> surface;

public:

	Zone(const label id_zone, const string name_zone, const string path, const string format);

	void implement_bc(const label id_bc, const string path);

	// Interpret Mesh 
	void interpret_inp(const string& path);

	void interpret_tecplot(const string& path);

	void bc_tecplot(const label& id_bc, const string& path);

	// Geometry Calculation
	void geometry_quad(label& i);

	void geometry_triangle(label& i);

	void geometry_hex(label& i);

	void geometry_tetra(label& i);

	// Offset Calculation
	void tune_decimal(Array<scalar>& array);

	void tune_decimal(scalar& value);

	void tune_scale(Array<scalar>& array);

	void reverse_axis(Array<scalar>& array);

	void offset_axis(Array<scalar>& array);
};
