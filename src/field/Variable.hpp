// Variables and Boundary Conditions

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "field/Zone.hpp"

class Variable {

public:

	Zone* ZONE{ nullptr };
	string NAME = "";
	string TYPE = "";             // scalar or vector

	Array<scalar> X;              // cell_center_value
	Array<scalar> X_faces;        // interpolate to face center
	Array<scalar> X_grad;         // gradient in cell (for scalar) 
	Array<scalar> X_divg;         // divergence in cell (for vector)

	vector BC;                    // used in Laplacian
	vector BC_exact;              // used in Convection and Interpolation
	std::vector<string> BC_type;  // store types of BCs in order

	Array<scalar> wallRotation;   // store line speed

public:

	Variable() {}

	Variable(Zone* zone, string type, string name, label ifPrint = 1);

	void setVariable(Zone* zone, string type, string name);

	void initial_value(scalar value_x = 0, scalar value_y = 0, scalar value_z = 0);

	void initial_value(string path, label dim = -1);

	// Boundary Conditions 
	void BC_exact_value(label id, scalar value_x, scalar value_y = 0, scalar value_z = 0);

	void BC_zero_gradient(label id);

	void BC_wall_temperature(label id, scalar temp);

	void BC_wall_heatflux(label id, scalar influx);

	void BC_wall_convection(label id, scalar coeff_convection, scalar temp);

	void BC_wall_velocity(label id, scalar value_x, scalar value_y, scalar value_z = 0);

	void BC_wall_rotation(label id, scalar angularSpeed, string axis, scalar origin1, scalar origin2);

	void BC_velocity_inlet(label id, scalar value_x, scalar value_y, scalar value_z = 0);

	void BC_pressure_inlet(label id, scalar pressure);

	void BC_pressure_outlet(label id, scalar pressure);

	void BC_copy(Variable* variable);

	void BC_all(string type);

	void BC_implement(string type_physics, Array<label>& bc_index, Array<string>& bc_type, Array<string>& bc_variable, Array<vector>& bc_value);

	// Field Operations 
	void get_face(string type_interp = "linear");

	void get_face_PWIM(Variable* p, Variable* A_diag);

	void get_gradient(Array<scalar>* phi_grad = nullptr, string type = "original");

	void get_gradient_face(Array<scalar>& phi_grad_face, string type = "zero grad");

	void get_divergence(Array<scalar>* phi_divg = nullptr, string type = "original");
};
