#include "field/Variable.hpp"

Variable::Variable(Zone* zone, string type, string name, label ifPrint) : ZONE(zone), TYPE(type), NAME(name) {
	if (TYPE == "scalar") {
		X.resize(ZONE->num_cells);
		X_faces.resize(ZONE->num_faces);
		X_grad.resize(ZONE->num_cells, num_dimension);
	}
	else if (TYPE == "vector") {
		X.resize(ZONE->num_cells, num_dimension);
		X_faces.resize(ZONE->num_faces, num_dimension);
		X_divg.resize(ZONE->num_cells);
	}
	if (ifPrint) {
		std::cout << "\nInitialize  : Zone = " << ZONE->name << ", Variable = " << NAME << ", Type = " << TYPE << std::endl;
	}
}

void Variable::setVariable(Zone* zone, string type, string name) {
	ZONE = zone;
	TYPE = type;
	NAME = name;

	if (TYPE == "scalar") {
		X.resize(ZONE->num_cells);
		X_faces.resize(ZONE->num_faces);
		X_grad.resize(ZONE->num_cells, num_dimension);
	}
	else if (TYPE == "vector") {
		X.resize(ZONE->num_cells, num_dimension);
		X_faces.resize(ZONE->num_faces, num_dimension);
		X_divg.resize(ZONE->num_cells);
	}
}

void Variable::initial_value(scalar value_x, scalar value_y, scalar value_z) {
	std::cout << "              " << "Initial Value = ";
	if (TYPE == "scalar") {
		X.set_all(value_x);
		std::cout << value_x << std::endl;
	}
	else if (num_dimension == 2) {
		for (label i = 0; i < X.NX; i++) {
			X(i, 0) = value_x;
			X(i, 1) = value_y;
		}
		std::cout << value_x << ", " << value_y << std::endl;
	}
	else {
		for (label i = 0; i < X.NX; i++) {
			X(i, 0) = value_x;
			X(i, 1) = value_y;
			X(i, 2) = value_z;
		}
		std::cout << value_x << ", " << value_y << ", " << value_z << std::endl;
	}
}

void Variable::initial_value(string path, label dim) {

	Array<scalar> X_read;
	X_read.read(path);

	if (dim < 0) {
		for (label i = 0; i < X_read.size; i++) {
			X(i) = X_read(i);
		}
	}
	else {
		for (label i = 0; i < X_read.size; i++) {
			X(i, dim) = X_read(i);
		}
	}
	std::cout << "              " << "Initialize with File = " << path << "\n";
}

void Variable::BC_exact_value(label id, scalar value_x, scalar value_y, scalar value_z) {  // <scalar> BC{ nan, nan, nan }, BC_exact{ value }  
	BC.push_back(nan("")); BC.push_back(nan("")); BC.push_back(nan(""));                   // <vector> BC{ nan, nan, nan, nan, nan } BC_exact{ value_x, value_y, value_z }
	BC_exact.push_back(value_x);
	std::cout << "              " << "BC_ID = " << id << ", Exact Value = " << value_x;
	if (TYPE == "vector") {
		BC.push_back(nan("")); BC_exact.push_back(value_y); std::cout << ", " << value_y;
		if (num_dimension == 3) { BC.push_back(nan("")); BC_exact.push_back(value_z); std::cout << ", " << value_z; }
	}
	BC_type.push_back("exact value");
	std::cout << "\n";
}

void Variable::BC_zero_gradient(label id) {                               // <scalar> BC{ nan, nan, nan }, BC_exact{ nan } 
	BC.push_back(nan("")); BC.push_back(nan("")); BC.push_back(nan(""));  // <vector> BC{ nan, nan, nan, nan, nan } BC_exact{ nan, nan, nan }
	BC_exact.push_back(nan(""));
	if (TYPE == "vector") {
		BC.push_back(nan("")); BC_exact.push_back(nan(""));
		if (num_dimension == 3) { BC.push_back(nan("")); BC_exact.push_back(nan("")); }
	}
	BC_type.push_back("zero grad");
	std::cout << "              " << "BC_ID = " << id << ", Zero Gradient\n";
}

void Variable::BC_wall_temperature(label id, scalar temp) {  // <scalar> BC{ 1, 1.0, temp }, BC_exact{ temp }
	BC.push_back(1);  // index of bc_type
	BC.push_back(1.0);
	BC.push_back(temp);
	BC_exact.push_back(temp);
	BC_type.push_back("wall temperature");
	std::cout << "              " << "BC_ID = " << id << ", Wall Temperature = " << temp << std::endl;
}

void Variable::BC_wall_heatflux(label id, scalar influx) {  // <scalar> BC{ 2, 0.0, influx }, BC_exact{ nan }
	BC.push_back(2);  // index of bc_type
	BC.push_back(0.0);
	BC.push_back(influx);
	BC_exact.push_back(nan(""));
	BC_type.push_back("wall heatflux");
	std::cout << "              " << "BC_ID = " << id << ", Wall Heat Flux = " << influx << std::endl;
}

void Variable::BC_wall_convection(label id, scalar coeff_convection, scalar temp) {  // <scalar> BC{ 3, coeff, temp }, BC_exact{ nan }
	BC.push_back(3);  // index of bc_type
	BC.push_back(coeff_convection);
	BC.push_back(temp);
	BC_exact.push_back(nan(""));
	BC_type.push_back("wall convection");
	std::cout << "              " << "BC_ID = " << id << ", Wall Convection Coefficient = " << coeff_convection << ", Temperature = " << temp << std::endl;
}

void Variable::BC_wall_velocity(label id, scalar value_x, scalar value_y, scalar value_z) {  // <vector> BC{ 1, 1.0, value_x, value_y, value_z }, BC_exact{ value_x, value_y, value_z }
	BC.push_back(1);  // index of bc_type                                                          
	BC.push_back(1.0);
	BC.push_back(value_x);  BC_exact.push_back(value_x);  std::cout << "              " << "BC_ID = " << id << ", Wall Velocity = " << value_x;
	BC.push_back(value_y);  BC_exact.push_back(value_y);  std::cout << ", " << value_y;
	if (num_dimension == 3) { BC.push_back(value_z);  BC_exact.push_back(value_z);  std::cout << ", " << value_z; }
	BC_type.push_back("wall velocity");
	std::cout << "\n";
}

void Variable::BC_velocity_inlet(label id, scalar value_x, scalar value_y, scalar value_z) {  // <vector> BC{ 2, 1.0, value_x, value_y, value_z }, BC_exact{ value_x, value_y, value_z } 
	BC.push_back(2);  // index of bc_type
	BC.push_back(1.0);
	BC.push_back(value_x);  BC_exact.push_back(value_x);  std::cout << "              " << "BC_ID = " << id << ", Velocity Inlet = " << value_x;
	BC.push_back(value_y);  BC_exact.push_back(value_y);  std::cout << ", " << value_y;
	if (num_dimension == 3) { BC.push_back(value_z);  BC_exact.push_back(value_z);  std::cout << ", " << value_z; }
	BC_type.push_back("velocity inlet");
	std::cout << "\n";
}

void Variable::BC_wall_rotation(label id, scalar angularSpeed, string axis, scalar origin1, scalar origin2) {  // <vector> BC{ 3, 1.0, 0, 0, 0 }, BC_exact{ 0, 0, 0 }
	BC.push_back(3);  // index of bc_type                                                          
	BC.push_back(1.0);
	BC.push_back(0);  BC_exact.push_back(0);
	BC.push_back(0);  BC_exact.push_back(0);
	if (num_dimension == 3) { BC.push_back(0);  BC_exact.push_back(0); }
	BC_type.push_back("wall rotation");

	if (wallRotation.size < ZONE->num_faces) {
		wallRotation.resize(ZONE->num_faces, 3);
	}
	scalar distance1, distance2;

	if (axis == "x") {
		for (label i = 0; i < ZONE->num_faces; i++) {
			if (ZONE->faces(i).cell_nb == -id) {
				distance1 = ZONE->faces(i).center[1] - origin1;
				distance2 = ZONE->faces(i).center[2] - origin2;
				wallRotation(i, 0) = 0;
				wallRotation(i, 1) = angularSpeed * distance1;
				wallRotation(i, 2) = angularSpeed * distance2;
			}
		}
	}
	else if (axis == "y") {
		for (label i = 0; i < ZONE->num_faces; i++) {
			if (ZONE->faces(i).cell_nb == -id) {
				distance1 = ZONE->faces(i).center[0] - origin1;
				distance2 = ZONE->faces(i).center[2] - origin2;
				wallRotation(i, 0) = angularSpeed * distance1;
				wallRotation(i, 1) = 0;
				wallRotation(i, 2) = angularSpeed * distance2;
			}
		}
	}
	else if (axis == "z") { 
		for (label i = 0; i < ZONE->num_faces; i++) {
			if (ZONE->faces(i).cell_nb == -id) {
				distance1 = ZONE->faces(i).center[0] - origin1;
				distance2 = ZONE->faces(i).center[1] - origin2;
				wallRotation(i, 0) = angularSpeed * distance1;
				wallRotation(i, 1) = angularSpeed * distance2;
				wallRotation(i, 2) = 0;
			}
		}
	}
	std::cout << "              " << "BC_ID = " << id << ", Wall Rotation, Angular Velocity = " << angularSpeed << "\n";
}

void Variable::BC_pressure_inlet(label id, scalar pressure) {  // <scalar> BC{ 4, 1.0, pressure }, BC_exact{ pressure }
	BC.push_back(4);  // index of bc_type
	BC.push_back(1.0);
	BC.push_back(pressure);
	BC_exact.push_back(pressure);
	BC_type.push_back("pressure inlet");
	std::cout << "              " << "BC_ID = " << id << ", Pressure Inlet = " << pressure << std::endl;
}

void Variable::BC_pressure_outlet(label id, scalar pressure) {  // <scalar> BC{ 5, 1.0, pressure }, BC_exact{ pressure }
	BC.push_back(5);  // index of bc_type
	BC.push_back(1.0);
	BC.push_back(pressure);
	BC_exact.push_back(pressure);
	BC_type.push_back("pressure outlet");
	std::cout << "              " << "BC_ID = " << id << ", Pressure Outlet = " << pressure << std::endl;
}

void Variable::BC_copy(Variable* variable) {  // copy Boundary Conditions of other variables
	for (size_t i = 0; i < variable->BC.size(); i++) { this->BC.push_back(variable->BC[i]); }
	for (size_t i = 0; i < variable->BC_exact.size(); i++) { this->BC_exact.push_back(variable->BC_exact[i]); }
	for (size_t i = 0; i < variable->BC_type.size(); i++) { this->BC_type.push_back(variable->BC_type[i]); }
}

void Variable::BC_all(string type) {  // set all the boundary conditions as zero-gradient extrapolation
	if (type == "zero grad") {
		for (size_t i = 0; i < BC.size(); i++) { BC[i] = nan(""); }
		for (size_t i = 0; i < BC_exact.size(); i++) { BC_exact[i] = nan(""); }
		for (size_t i = 0; i < BC_type.size(); i++) { BC_type[i] = type; }
	}
}

void Variable::BC_implement(string type_physics, Array<label>& bc_index, Array<string>& bc_type, Array<string>& bc_variable, Array<vector>& bc_value) {  // called from <Setup.hpp> to implement BC
	label index;
	label type = 0;
	if (type_physics == "TEMP") { type = 1; }
	else if (type_physics == "VT") { type = 2; }
	else if (type_physics == "P") { type = 3; }

	for (label i = 0; i < ZONE->num_surface; i++) {
		index = bc_index(i);

		if (bc_type(index) == "BoundaryWall") {

			if (bc_variable(index) == "V1,V2,V3,Temperature") {
				if (type == 1) { BC_wall_temperature(index + 1, bc_value(index)[3]); }
				else if (type == 2) { BC_wall_velocity(index + 1, bc_value(index)[0], bc_value(index)[1], bc_value(index)[2]); }
				else if (type == 3) { BC_zero_gradient(index + 1); }
			}

			else if (bc_variable(index) == "V1,V2,V3,HeatFlux") {
				if (type == 1) { BC_wall_heatflux(index + 1, bc_value(index)[3]); }
				else if (type == 2) { BC_wall_velocity(index + 1, bc_value(index)[0], bc_value(index)[1], bc_value(index)[2]); }
				else if (type == 3) { BC_zero_gradient(index + 1); }
			}

			else if (bc_variable(index) == "V1,V2,V3,Temperature,Convection") {
				if (type == 1) { BC_wall_convection(index + 1, bc_value(index)[4], bc_value(index)[3]); }
				else if (type == 2) { BC_wall_velocity(index + 1, bc_value(index)[0], bc_value(index)[1], bc_value(index)[2]); }
				else if (type == 3) { BC_zero_gradient(index + 1); }
			}
		}

		else if (bc_type(index) == "Velocitylnlet") {

			if (bc_variable(index) == "Temperature,Magnitude") {
				if (type == 1) { BC_wall_temperature(index + 1, bc_value(index)[0]); }
				else if (type == 2) { BC_velocity_inlet(index + 1, bc_value(index)[1], nan(""), nan("")); }
				else if (type == 3) { BC_zero_gradient(index + 1); }
			}

			else if (bc_variable(index) == "Temperature,Components") {
				if (type == 1) { BC_wall_temperature(index + 1, bc_value(index)[0]); }
				else if (type == 2) { BC_velocity_inlet(index + 1, bc_value(index)[1], bc_value(index)[2], bc_value(index)[3]); }
				else if (type == 3) { BC_zero_gradient(index + 1); }
			}
		}

		else if (bc_type(index) == "PressureInlet") {
			if (type == 1) { BC_wall_temperature(index + 1, bc_value(index)[1]); }
			else if (type == 2) { BC_zero_gradient(index + 1); }
			else if (type == 3) { BC_pressure_inlet(index + 1, bc_value(index)[0]); }
		}

		else if (bc_type(index) == "PressureOutlet") {
			if (type == 1) { BC_zero_gradient(index + 1); }
			else if (type == 2) { BC_zero_gradient(index + 1); }
			else if (type == 3) { BC_pressure_outlet(index + 1, bc_value(index)[0]); }
		}
	}

	scalar initial;
	scalar initial_x, initial_y, initial_z;
	if (TYPE == "scalar") {
		initial = BC_exact[0];
		for (size_t i = 1; i < BC_exact.size(); i++) {
			initial = find_min(initial, BC_exact[i]);
		}
		initial_value(initial);
	}
	else if (TYPE == "vector") {
		initial_x = BC_exact[0];
		initial_y = BC_exact[1];
		initial_z = BC_exact[2];
		for (size_t i = 1; i < BC_exact.size() / 3; i++) {
			initial_x = find_min(initial_x, BC_exact[i * 3]);
			initial_y = find_min(initial_y, BC_exact[i * 3 + 1]);
			initial_z = find_min(initial_z, BC_exact[i * 3 + 2]);
		}
		initial_value(initial_x, initial_y, initial_z);
	}
}

void Variable::get_face(string type_interp) {

	auto& faces = ZONE->faces;
	auto& num_faces = ZONE->num_faces;
	label index;

	if (TYPE == "scalar") {
		if (type_interp == "linear") {
			for (label i = 0; i < num_faces; i++) {
				if (faces(i).cell_nb >= 0) {  // inner faces
					X_faces(i) = faces(i).gc * X(faces(i).cell) + faces(i).gnb * X(faces(i).cell_nb);
				}
				else {  // boundary facas
					index = -faces(i).cell_nb - 1;
					if (BC_exact.size() == 0) { X_faces(i) = X(faces(i).cell); }
					else if (!isnan(BC_exact[index])) { X_faces(i) = BC_exact[index]; }
					else {
						if (BC_type[index] == "zero grad") { X_faces(i) = X(faces(i).cell); }
					}
				}
			}
		}
	}
	else if (TYPE == "vector") {
		if (type_interp == "linear") {
			for (label i = 0; i < num_faces; i++) {
				if (faces(i).cell_nb >= 0) {  // inner faces
					for (label j = 0; j < num_dimension; j++) {
						X_faces(i, j) = faces(i).gc * X(faces(i).cell, j) + faces(i).gnb * X(faces(i).cell_nb, j);
					}
				}
				else {  // boundary faces
					index = num_dimension * (-faces(i).cell_nb - 1);
					if (BC_exact.size() == 0) {
						for (label j = 0; j < num_dimension; j++) {
							X_faces(i, j) = X(faces(i).cell, j);
						}
					}
					else if (!isnan(BC_exact[index])) {
						for (label j = 0; j < num_dimension; j++) {
							X_faces(i, j) = BC_exact[index + j];
						}
					}
					else {
						if (BC_type[-faces(i).cell_nb - 1] == "zero grad") {
							for (label j = 0; j < num_dimension; j++) {
								X_faces(i, j) = X(faces(i).cell, j);
							}
						}
					}
				}
			}
		}
	}
}

void Variable::get_face_PWIM(Variable* p, Variable* A_diag) {  // for Vector

	auto& faces = ZONE->faces;
	auto& num_faces = ZONE->num_faces;
	auto& p_grad = p->X_grad;
	auto& diag = A_diag->X;
	label cell, cell_nb, index;
	scalar gc, gnb, correction;
	p->get_gradient();

	for (label i = 0; i < num_faces; i++) {

		cell = faces(i).cell;
		cell_nb = faces(i).cell_nb;
		gc = faces(i).gc;
		gnb = faces(i).gnb;

		if (cell_nb >= 0) {
			for (label j = 0; j < num_dimension; j++) {
				correction = diag(cell) * p_grad(cell, j) + diag(cell_nb) * p_grad(cell_nb, j)
					- (diag(cell) + diag(cell_nb)) * (gc * p_grad(cell, j) + gnb * p_grad(cell_nb, j));
				X_faces(i, j) = gc * X(cell, j) + gnb * X(cell_nb, j) + 0.5 * correction;
			}
		}
		else {
			for (label j = 0; j < num_dimension; j++) {
				index = num_dimension * (-cell_nb - 1) + j;
				if (!isnan(BC_exact[index])) { X_faces(i, j) = BC_exact[index]; }
				else {
					if (BC_type[-cell_nb - 1] == "zero grad") { X_faces(i, j) = X(cell, j); }
				}
			}
		}
	}
}

void Variable::get_gradient(Array<scalar>* phi_grad, string type) {  // for Scalar

	auto& num_cells = ZONE->num_cells;
	auto& num_faces = ZONE->num_faces;
	auto& faces = ZONE->faces;
	auto& volume = ZONE->volume;
	X_grad.set_all(0);

	if (type == "original") {
		for (label i = 0; i < num_faces; i++) {
			if (faces(i).if_repeat == 1) { continue; }

			X_grad(faces(i).cell, 0) += X_faces(i) * faces(i).normal[0] * faces(i).area / volume(faces(i).cell);
			X_grad(faces(i).cell, 1) += X_faces(i) * faces(i).normal[1] * faces(i).area / volume(faces(i).cell);
			if (num_dimension == 3) { X_grad(faces(i).cell, 2) += X_faces(i) * faces(i).normal[2] * faces(i).area / volume(faces(i).cell); }

			if (faces(i).cell_nb >= 0) {
				X_grad(faces(i).cell_nb, 0) -= X_faces(i) * faces(i).normal[0] * faces(i).area / volume(faces(i).cell_nb);
				X_grad(faces(i).cell_nb, 1) -= X_faces(i) * faces(i).normal[1] * faces(i).area / volume(faces(i).cell_nb);
				if (num_dimension == 3) { X_grad(faces(i).cell_nb, 2) -= X_faces(i) * faces(i).normal[2] * faces(i).area / volume(faces(i).cell_nb); }
			}
		}
	}
	else if (type == "SIMPLE") {
		for (label i = 0; i < num_faces; i++) {
			if (faces(i).if_repeat == 1) { continue; }

			X_grad(faces(i).cell, 0) -= X_faces(i) * faces(i).normal[0] * faces(i).area;
			X_grad(faces(i).cell, 1) -= X_faces(i) * faces(i).normal[1] * faces(i).area;
			if (num_dimension == 3) { X_grad(faces(i).cell, 2) -= X_faces(i) * faces(i).normal[2] * faces(i).area; }

			if (faces(i).cell_nb >= 0) {
				X_grad(faces(i).cell_nb, 0) += X_faces(i) * faces(i).normal[0] * faces(i).area;
				X_grad(faces(i).cell_nb, 1) += X_faces(i) * faces(i).normal[1] * faces(i).area;
				if (num_dimension == 3) { X_grad(faces(i).cell_nb, 2) += X_faces(i) * faces(i).normal[2] * faces(i).area; }
			}
		}
	}

	if (phi_grad != nullptr) {
		for (label i = 0; i < num_cells; i++) {
			for (label j = 0; j < num_dimension; j++) {
				(*phi_grad)(i, j) += X_grad(i, j);
			}
		}
	}
}

void Variable::get_gradient_face(Array<scalar>& phi_grad_face, string type) {

	auto& faces = ZONE->faces;
	label cell, cell_nb;
	scalar gc, gnb, correction;

	for (label i = 0; i < ZONE->num_faces; i++) {

		cell = faces(i).cell;
		cell_nb = faces(i).cell_nb;
		gc = faces(i).gc;
		gnb = faces(i).gnb;

		if (cell_nb >= 0) {
			correction = (X(cell_nb) - X(cell)) / faces(i).distance;
			for (label j = 0; j < num_dimension; j++) {
				phi_grad_face(i, j) = gc * X_grad(cell, j) + gnb * X_grad(cell_nb, j);
				correction -= phi_grad_face(i, j) * faces(i).normal[j];
			}
			for (label j = 0; j < num_dimension; j++) {
				phi_grad_face(i, j) += correction * faces(i).normal[j];
			}
		}
		else {
			if (type == "zero grad") {
				for (label j = 0; j < num_dimension; j++) {
					phi_grad_face(i, j) = 0;
				}
			}
		}
	}
}

void Variable::get_divergence(Array<scalar>* phi_divg, string type) {  // for Vector

	auto& num_cells = ZONE->num_cells;
	auto& num_faces = ZONE->num_faces;
	auto& faces = ZONE->faces;
	auto& volume = ZONE->volume;
	scalar value;
	X_divg.set_all(0);

	if (type == "original") {
		for (label i = 0; i < num_faces; i++) {
			value = 0;
			for (label j = 0; j < num_dimension; j++) {
				value += X_faces(i, j) * faces(i).normal[j];
			}
			value = value * faces(i).area;
			X_divg(faces(i).cell) += value / volume(faces(i).cell);
		}
	}
	else if (type == "SIMPLE") {
		for (label i = 0; i < num_faces; i++) {
			value = 0;
			for (label j = 0; j < num_dimension; j++) {
				value += X_faces(i, j) * faces(i).normal[j];
			}
			value = value * faces(i).area;
			X_divg(faces(i).cell) -= value;
		}
	}

	if (phi_divg != nullptr) {
		for (label i = 0; i < num_cells; i++) {
			(*phi_divg)(i) += X_divg(i);
		}
	}
}
