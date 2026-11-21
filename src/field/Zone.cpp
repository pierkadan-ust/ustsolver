#include "field/Zone.hpp"

Zone::Zone(const label id_zone, const string name_zone, const string path, const string format) : type_mesh(0), num_coordinate(0), num_faces(0), num_nodes(0), num_cells(0), num_face(0), num_vertex(0), num_vertex_face(0), num_surface(0) {
	id = id_zone;
	name = name_zone;
	type_input = format;
	std::cout << "\nRead Mesh   : ";
	if (type_input == "Tecplot") { interpret_tecplot(path); }
	else if (type_input == "inp") { interpret_inp(path); }
	std::cout << "Zone Created Successfully, Zone_ID = " << id << ", Zone_Name = " << name << std::endl;
}

void Zone::implement_bc(const label id_bc, const string path) { 
	if (type_input == "Tecplot") { bc_tecplot(id_bc, path); }
}

void Zone::interpret_inp(const string& path) {  

	std::ifstream file;
	file.open(path);

	Array<scalar> node;
	Array<label> element;
	std::istringstream iss;
	size_t pos;
	string line;
	string value;
	label cell;
	label side;
	label value_label;

	// [1] store Node
	while (getline_clean(file, line)) {
		if (line.length() >= 5 && line.substr(0, 5) == "*Node") {
			while (getline_clean(file, line) && line[0] != '*') {
				iss.str(line);
				std::getline(iss, value, ',');
				std::getline(iss, value, ',');  node.push_back(std::stod(value));
				std::getline(iss, value, ',');  node.push_back(std::stod(value));
				std::getline(iss, value);       node.push_back(std::stod(value));
				iss.clear();
			}
			break;
		}
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	num_nodes = node.size / 3;
	num_coordinate = 3;
	nodes.resize(num_nodes, num_coordinate);
	nodes.set_array(node);

	if (if_tune_decimal) { tune_decimal(nodes); }
	if (if_tune_scale) { tune_scale(nodes); }
	if (if_reverse_axis != "") { reverse_axis(nodes); }
	if (if_offset_axis) { offset_axis(nodes); }

	std::cout << "Node Information Found, Number of Nodes = " << num_nodes << std::endl;

	// [2] store Element
	while (getline_clean(file, line)) {
		if (line.length() >= 8 && line.substr(0, 8) == "*Element") {
			iss.str(line);
			std::getline(iss, value, ',');
			std::getline(iss, value, ',');
			iss.clear();
			pos = value.find('=');
			if (value.substr(pos + 1) == "C3D8" || value.substr(pos + 1) == "C3D8R") {
				type_mesh = 3;
				num_face = 6;
				num_vertex = 8;
				num_vertex_face = 4;
			}
			else if (value.substr(pos + 1) == "C3D4") {
				type_mesh = 4;
				num_face = 4;
				num_vertex = 4;
				num_vertex_face = 3;
			}
			while (getline_clean(file, line) && line[0] != '*') {
				iss.str(line);
				std::getline(iss, value, ',');
				for (label i = 0; i < num_vertex - 1; i++) {
					std::getline(iss, value, ',');
					value_label = std::stoi(value);
					value_label--;
					element.push_back(value_label);
				}
				std::getline(iss, value, ',');
				value_label = std::stoi(value);
				value_label--;
				element.push_back(value_label);
				iss.clear();
			}
			break;
		}
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	num_cells = element.size / num_vertex;
	cells.resize(num_cells, num_vertex);  // num_vertex = max_num_vertex, if multiple mesh types
	cells.set_array(element);

	center.resize(num_cells, num_coordinate);
	volume.resize(num_cells);

	types.resize(num_cells, 4);
	types.set_all(type_mesh, 0);
	types.set_all(num_face, 1);
	types.set_all(num_vertex, 2);
	types.set_all(num_vertex_face, 3);

	std::cout << "              " << "Element Information Found, Number of Elements = " << num_cells << "\n";
	std::cout << "              ";  return;

	// [3] find Boundary Surface [NAME ELSET SIDE]
	while (getline_clean(file, line)) {
		if (line.length() >= 8 && line.substr(0, 8) == "*Surface") {
			iss.str(line);
			std::getline(iss, value, ',');
			std::getline(iss, value, ',');
			iss.clear();
			pos = value.find('=');
			surface.push_back(value.substr(pos + 1));  // NAME
			getline_clean(file, line);
			iss.str(line);
			std::getline(iss, value, ',');
			surface.push_back(value);  // ELSET
			std::getline(iss, value);
			surface.push_back(value);  // SIDE
			iss.clear();
		}
	}
	file.clear();
	file.seekg(0, std::ios::beg);

	num_surface = surface.size / 3;
	num_face = find_max(types, 1);  // num_face = max_num_face, if multiple mesh types
	num_faces = num_cells * num_face;
	faces.resize(num_cells, num_face);
	for (label i = 0; i < num_faces; i++) {
		faces(i).cell_nb = LLONG_MIN;  // the initial value of cell_nb
	}

	// [4] store ElementSet related to Boundary Surface
	for (label i = 0; i < num_surface; i++) {
		while (getline_clean(file, line)) {
			if (line.length() >= 6 && line.substr(0, 6) == "*ELSET") {
				iss.str(line);
				std::getline(iss, value, ',');
				std::getline(iss, value, ',');
				iss.clear();
				pos = value.find('=');
				if (value.substr(pos + 1) == surface(i * 3 + 1)) {
					value = surface(i * 3 + 2);
					value = value[1];
					side = std::stoi(value);
					side--;
					while (getline_clean(file, line) && line[0] != '*') {
						std::replace(line.begin(), line.end(), ',', ' ');
						iss.str(line);
						while (iss >> cell) {
							cell--;
							faces(cell, side).cell_nb = -(i + 1);  // all the boundary faces are indexed as {cell_nb < 0}
						}
						iss.clear();
					}
					break;
				}
				else { continue; }
			}
		}
		file.clear();
		file.seekg(0, std::ios::beg);
	}
	file.close();

	// [5] obtain Connectivity by comparing the face center
	label* face = nullptr;
	label index_node = 0;
	label number_face = 0;
	label number_vertex_face = 0;

	label face_C3D8[6][4] = {
		{0, 1, 2, 3},
		{4, 7 ,6, 5},
		{0, 4, 5, 1},
		{1, 5, 6, 2},
		{2, 6, 7, 3},
		{3, 7, 4, 0}
	};

	label face_C3D4[4][3] = {
		{0, 1, 2},
		{0, 3, 1},
		{1, 3, 2},
		{2, 3, 0},
	};

	for (label i = 0; i < num_cells; i++) {
		if (types(i, 0) == 3) { face = face_C3D8[0]; number_face = 6; number_vertex_face = 4; }
		else if (types(i, 0) == 4) { face = face_C3D4[0]; number_face = 4; number_vertex_face = 3; }

		for (label j = 0; j < number_face; j++) {
			for (label k = 0; k < number_vertex_face; k++) {
				index_node = *(face + j * number_vertex_face + k);
				index_node = cells(i, index_node);
				for (label d = 0; d < num_coordinate; d++) {
					faces(i, j).center[d] += nodes(index_node, d);
				}
			}
		}
	}

	bool if_break = 0;
	scalar x = 0, y = 0, z = 0;
	for (label i = 0; i < num_cells; i++) {
		for (label j = 0; j < num_face; j++) {
			faces(i, j).cell = i;
			faces(i, j).side = j;

			if (faces(i, j).cell_nb == LLONG_MIN && faces(i, j).if_repeat == 0) {
				x = faces(i, j).center[0];
				y = faces(i, j).center[1];
				z = faces(i, j).center[2];

				for (label i_nb = 0; i_nb < num_cells; i_nb++) {
					if (i_nb == i) { continue; }
					for (label j_nb = 0; j_nb < num_face; j_nb++) {
						if (x == faces(i_nb, j_nb).center[0] &&
							y == faces(i_nb, j_nb).center[1] &&
							z == faces(i_nb, j_nb).center[2]) {

							faces(i, j).cell_nb = i_nb;
							faces(i, j).side_nb = j_nb;

							faces(i_nb, j_nb).cell_nb = i;
							faces(i_nb, j_nb).side_nb = j;
							faces(i_nb, j_nb).if_repeat = 1;

							if_break = 1;
							break;
						}
					}
					if (if_break) { if_break = 0; break; }
				}
			}
		}
	}

	// [6] obtain Geometry Information
	for (label i = 0; i < num_cells; i++) {
		if (types(i, 0) == 3) { geometry_hex(i); }
	}

	// [7] obtain cell_center to cell_center vector
	scalar norm = 0;
	label cell_nb = 0;

	for (label i = 0; i < num_cells; i++) {
		for (label j = 0; j < num_face; j++) {
			cell_nb = faces(i, j).cell_nb;
			if (cell_nb >= 0) {  // inner faces
				x = center(cell_nb, 0) - center(i, 0);
				y = center(cell_nb, 1) - center(i, 1);
				z = center(cell_nb, 2) - center(i, 2);
				if (if_tune_decimal) { tune_decimal(x); tune_decimal(y); tune_decimal(z); }
				norm = sqrt(x * x + y * y + z * z);
				faces(i, j).direction[0] = x / norm;
				faces(i, j).direction[1] = y / norm;
				faces(i, j).direction[2] = z / norm;
				faces(i, j).distance = norm;
				faces(i, j).gc = volume(i) / (volume(i) + volume(cell_nb));
				faces(i, j).gnb = 1 - faces(i, j).gc;
				if (faces(i, j).direction[0] == faces(i, j).normal[0] && faces(i, j).direction[1] == faces(i, j).normal[1] && faces(i, j).direction[2] == faces(i, j).normal[2]) { continue; }
				else { faces(i, j).if_correction = 1; }
			}
			else {  // boundary faces
				x = faces(i, j).center[0] - center(i, 0);
				y = faces(i, j).center[1] - center(i, 1);
				z = faces(i, j).center[2] - center(i, 2);
				if (if_tune_decimal) { tune_decimal(x); tune_decimal(y); tune_decimal(z); }
				norm = sqrt(x * x + y * y + z * z);
				faces(i, j).direction[0] = x / norm;
				faces(i, j).direction[1] = y / norm;
				faces(i, j).direction[2] = z / norm;
				faces(i, j).distance = norm;
				faces(i, j).gc = 0.5;
				faces(i, j).gnb = 0.5;
				if (faces(i, j).direction[0] == faces(i, j).normal[0] && faces(i, j).direction[1] == faces(i, j).normal[1] && faces(i, j).direction[2] == faces(i, j).normal[2]) { continue; }
				else { faces(i, j).if_correction = 1; }
			}
		}
	}
}

void Zone::interpret_tecplot(const string& path) {

	Array<scalar> mesh;
	mesh.read(path);

	num_coordinate = 3;
	label num_connects;
	label count = 0;

	type_mesh = static_cast<int>(mesh(count)); count++;
	num_nodes = static_cast<int>(mesh(count)); count++;
	num_cells = static_cast<int>(mesh(count)); count++;
	num_connects = static_cast<int>(mesh(count)); count++;
	nodes.resize(num_nodes, num_coordinate);
	center.resize(num_cells, num_coordinate);
	volume.resize(num_cells);

	types.resize(num_cells, 4);
	if (type_mesh == 1 || type_mesh == -1) { num_face = 4; num_vertex = 4; num_vertex_face = 2; }
	else if (type_mesh == 2) { num_face = 3; num_vertex = 3; num_vertex_face = 2; }
	else if (type_mesh == 3) { num_face = 6; num_vertex = 8; num_vertex_face = 4; }
	else if (type_mesh == 4) { num_face = 4; num_vertex = 4; num_vertex_face = 3; }
	types.set_all(abs(type_mesh), 0); types.set_all(num_face, 1); types.set_all(num_vertex, 2); types.set_all(num_vertex_face, 3);

	num_faces = num_cells * num_face;
	faces.resize(num_cells, num_face);
	cells.resize(num_cells, num_vertex);

	for (label i = 0; i < num_nodes * num_coordinate; i++) { nodes(i) = mesh(count); count++; }
	for (label i = 0; i < num_cells * num_vertex; i++) { cells(i) = static_cast<int>(mesh(count)); cells(i) -= 1; count++; }

	if (if_tune_decimal) { tune_decimal(nodes); }
	if (if_tune_scale) { tune_scale(nodes); }
	if (if_reverse_axis != "") { reverse_axis(nodes); }
	if (if_offset_axis) { offset_axis(nodes); }

	// store Connectivity to <faces> ------------------------------------
	label cell, face, cell_nb;
	for (label i = 0; i < num_cells; i++) {
		for (label j = 0; j < num_face; j++) {
			faces(i, j).cell = i;
			faces(i, j).side = j;
			faces(i, j).cell_nb = -100;
		}
	}
	for (label i = 0; i < num_connects; i++) {
		cell = count + 3 * i; face = cell + 1; cell_nb = cell + 2;
		cell = static_cast<int>(mesh(cell)); cell -= 1;
		face = static_cast<int>(mesh(face)); face -= 1;
		cell_nb = static_cast<int>(mesh(cell_nb)); cell_nb -= 1;
		faces(cell, face).cell_nb = cell_nb;  // now all the boundary faces are marked with cell_nb = -100
	}
	for (label i = 0; i < num_cells; i++) {
		for (label j = 0; j < num_face; j++) {
			cell_nb = faces(i, j).cell_nb;
			if (cell_nb == -100) { continue; }
			if (cell_nb < i) { faces(i, j).if_repeat = 1; }  // marked all the repeated faces
			for (label k = 0; k < num_face; k++) {
				if (faces(cell_nb, k).cell_nb == i) { faces(i, j).side_nb = k; }
			}
		}
	}

	// store Geometry Information to <faces> ----------------------------
	for (label i = 0; i < num_cells; i++) {
		if (types(i, 0) == 1) { geometry_quad(i); }
		else if (types(i, 0) == 2) { geometry_triangle(i); }
		else if (types(i, 0) == 3) { geometry_hex(i); }
		else if (types(i, 0) == 4) { geometry_tetra(i); }
	}

	// cell_center to cell_center vector --------------------------------
	scalar x, y, z, norm;
	for (label i = 0; i < num_cells; i++) {
		for (label j = 0; j < num_face; j++) {
			cell_nb = faces(i, j).cell_nb;
			if (cell_nb >= 0) {
				x = center(cell_nb, 0) - center(i, 0);
				y = center(cell_nb, 1) - center(i, 1);
				z = center(cell_nb, 2) - center(i, 2);
				if (if_tune_decimal) { tune_decimal(x); tune_decimal(y); tune_decimal(z); }
				norm = sqrt(x * x + y * y + z * z);
				faces(i, j).direction[0] = x / norm;
				faces(i, j).direction[1] = y / norm;
				faces(i, j).direction[2] = z / norm;
				faces(i, j).distance = norm;
				faces(i, j).gc = volume(i) / (volume(i) + volume(cell_nb));
				faces(i, j).gnb = 1 - faces(i, j).gc;
				if (faces(i, j).direction[0] == faces(i, j).normal[0] && faces(i, j).direction[1] == faces(i, j).normal[1] && faces(i, j).direction[2] == faces(i, j).normal[2]) { continue; }
				else { faces(i, j).if_correction = 1; }
			}
			else {
				x = faces(i, j).center[0] - center(i, 0);
				y = faces(i, j).center[1] - center(i, 1);
				z = faces(i, j).center[2] - center(i, 2);
				if (if_tune_decimal) { tune_decimal(x); tune_decimal(y); tune_decimal(z); }
				norm = sqrt(x * x + y * y + z * z);
				faces(i, j).direction[0] = x / norm;
				faces(i, j).direction[1] = y / norm;
				faces(i, j).direction[2] = z / norm;
				faces(i, j).distance = norm;
				faces(i, j).gc = 0.5;
				faces(i, j).gnb = 0.5;
				if (faces(i, j).direction[0] == faces(i, j).normal[0] && faces(i, j).direction[1] == faces(i, j).normal[1] && faces(i, j).direction[2] == faces(i, j).normal[2]) { continue; }
				else { faces(i, j).if_correction = 1; }
			}
		}
	}
}

void Zone::bc_tecplot(const label& id_bc, const string& path) {

	Array<scalar> bc;
	bc.read(path);

	label count = 0;
	label count_find = 0;
	label num_nodes_bc;
	label num_faces_bc;
	scalar x_bc = 0, y_bc = 0, z_bc = 0;

	num_nodes_bc = static_cast<int>(bc(count)); count++;
	num_faces_bc = static_cast<int>(bc(count)); count++;

	Array<scalar> nodes_bc(num_nodes_bc, num_coordinate);
	Array<label>  faces_bc(num_faces_bc, num_vertex_face);

	for (label i = 0; i < num_nodes_bc * num_coordinate; i++) { nodes_bc(i) = bc(count); count++; }
	for (label i = 0; i < num_faces_bc * num_vertex_face; i++) { faces_bc(i) = static_cast<int>(bc(count)); faces_bc(i) -= 1; count++; }

	if (if_tune_decimal) { tune_decimal(nodes_bc); }
	if (if_tune_scale) { tune_scale(nodes_bc); }
	if (if_reverse_axis != "") { reverse_axis(nodes_bc); }
	if (if_offset_axis) { offset_axis(nodes_bc); }

	// find boundary location by comparing face center ------------------------------
	for (label i = 0; i < num_faces_bc; i++) {
		if (num_vertex_face == 2) {
			x_bc = 0.5 * (nodes_bc(faces_bc(i, 0), 0) + nodes_bc(faces_bc(i, 1), 0));
			y_bc = 0.5 * (nodes_bc(faces_bc(i, 0), 1) + nodes_bc(faces_bc(i, 1), 1));
		}
		else if (num_vertex_face == 3) {
			x_bc = nodes_bc(faces_bc(i, 0), 0) + nodes_bc(faces_bc(i, 1), 0) + nodes_bc(faces_bc(i, 2), 0);
			y_bc = nodes_bc(faces_bc(i, 0), 1) + nodes_bc(faces_bc(i, 1), 1) + nodes_bc(faces_bc(i, 2), 1);
			z_bc = nodes_bc(faces_bc(i, 0), 2) + nodes_bc(faces_bc(i, 1), 2) + nodes_bc(faces_bc(i, 2), 2);
			x_bc = x_bc / 3.0;
			y_bc = y_bc / 3.0;
			z_bc = z_bc / 3.0;
		}
		else if (num_vertex_face == 4) {
			x_bc = nodes_bc(faces_bc(i, 0), 0) + nodes_bc(faces_bc(i, 1), 0) + nodes_bc(faces_bc(i, 2), 0) + nodes_bc(faces_bc(i, 3), 0);
			y_bc = nodes_bc(faces_bc(i, 0), 1) + nodes_bc(faces_bc(i, 1), 1) + nodes_bc(faces_bc(i, 2), 1) + nodes_bc(faces_bc(i, 3), 1);
			z_bc = nodes_bc(faces_bc(i, 0), 2) + nodes_bc(faces_bc(i, 1), 2) + nodes_bc(faces_bc(i, 2), 2) + nodes_bc(faces_bc(i, 3), 2);
			x_bc = x_bc / 4.0;
			y_bc = y_bc / 4.0;
			z_bc = z_bc / 4.0;
		}

		for (label j = 0; j < num_faces; j++) {
			if (abs(x_bc - faces(j).center[0]) < 0.0001 && abs(y_bc - faces(j).center[1]) < 0.0001 && abs(z_bc - faces(j).center[2]) < 0.0001) {
				faces(j).cell_nb = -id_bc;
				count_find++;
			}
		}
	}
	if (count_find == num_faces_bc) { std::cout << "              " << "Boundary Implemented Successfully, BC_ID = " << id_bc << std::endl; }
	else { std::cout << "              " << "ERROR: Boundary Implementation" << std::endl; }
}

void Zone::geometry_quad(label& i) {

	label index;
	label nb = 0;
	scalar x[5] = { 0.0 };
	scalar y[5] = { 0.0 };
	scalar a1, b1, c1, a2, b2, c2, m;
	scalar d0, d1, d2, d3, d4, k1, k2;

	for (label j = 0; j < num_vertex; j++) {
		index = cells(i, j);
		x[j] = nodes(index, 0);
		y[j] = nodes(index, 1);
	}

	a1 = y[2] - y[0];
	b1 = x[0] - x[2];
	c1 = x[2] * y[0] - x[0] * y[2];
	a2 = y[3] - y[1];
	b2 = x[1] - x[3];
	c2 = x[3] * y[1] - x[1] * y[3];
	m = a1 * b2 - a2 * b1;

	// center of each cell --------------------------
	center(i, 0) = (b1 * c2 - b2 * c1) / m;
	center(i, 1) = (a2 * c1 - a1 * c2) / m;

	x[4] = x[0];
	y[4] = y[0];
	for (int j = 0; j < num_face; j++) {
		a1 = y[j + 1] - y[j];
		b1 = x[j] - x[j + 1];
		c1 = sqrt(a1 * a1 + b1 * b1);

		// normal of each face --------------------------
		faces(i, j).normal[0] = a1 / c1;
		faces(i, j).normal[1] = b1 / c1;
		faces(i, j).area = c1;

		// center of each face --------------------------
		faces(i, j).center[0] = 0.5 * (x[j] + x[j + 1]);
		faces(i, j).center[1] = 0.5 * (y[j] + y[j + 1]);
	}

	d0 = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2));
	d1 = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2));
	d2 = sqrt(pow(x[3] - x[2], 2) + pow(y[3] - y[2], 2));
	d3 = sqrt(pow(x[0] - x[3], 2) + pow(y[0] - y[3], 2));
	d4 = sqrt(pow(x[3] - x[1], 2) + pow(y[3] - y[1], 2));
	k1 = 0.5 * (d0 + d3 + d4);
	k2 = 0.5 * (d1 + d2 + d4);

	// volume of each cell --------------------------
	volume(i) = sqrt(k1 * (k1 - d0) * (k1 - d3) * (k1 - d4)) + sqrt(k2 * (k2 - d1) * (k2 - d2) * (k2 - d4));
}

void Zone::geometry_triangle(label& i) {

	label index;
	scalar x[4] = { 0.0 };
	scalar y[4] = { 0.0 };
	scalar d0, d1, d2, dk;

	for (label j = 0; j < num_vertex; j++) {
		index = cells(i, j);
		x[j] = nodes(index, 0);
		y[j] = nodes(index, 1);
	}
	x[3] = x[0];
	y[3] = y[0];

	// center of each cell --------------------------
	center(i, 0) = (x[0] + x[1] + x[2]) / 3.0;
	center(i, 1) = (y[0] + y[1] + y[2]) / 3.0;

	// volume and area of each cell --------------------------
	d0 = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2));
	d1 = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2));
	d2 = sqrt(pow(x[3] - x[2], 2) + pow(y[3] - y[2], 2));
	dk = 0.5 * (d0 + d1 + d2);

	volume(i) = sqrt(dk * (dk - d0) * (dk - d1) * (dk - d2));
	faces(i, 0).area = d0;
	faces(i, 1).area = d1;
	faces(i, 2).area = d2;

	for (int j = 0; j < num_face; j++) {

		// center of each face --------------------------
		faces(i, j).center[0] = 0.5 * (x[j] + x[j + 1]);
		faces(i, j).center[1] = 0.5 * (y[j] + y[j + 1]);

		// normal of each face --------------------------
		d1 = y[j + 1] - y[j];
		d2 = x[j] - x[j + 1];
		dk = sqrt(d1 * d1 + d2 * d2);

		faces(i, j).normal[0] = d1 / dk;
		faces(i, j).normal[1] = d2 / dk;
	}
}

void Zone::geometry_hex(label& i) {

	label index;
	label index_face[6][4] = { 0 };
	scalar x[8] = { 0 }, y[8] = { 0 }, z[8] = { 0 };
	scalar fx[4] = { 0 }, fy[4] = { 0 }, fz[4] = { 0 };
	scalar d0, d1, d2, d3, d4, k1, k2;
	scalar a, b, c, d, e, f, g, h, m, k;
	scalar V = 0;

	if (type_input == "Tecplot") {
		index_face[0][0] = 0;  index_face[0][1] = 4;  index_face[0][2] = 7;  index_face[0][3] = 3;
		index_face[1][0] = 1;  index_face[1][1] = 2;  index_face[1][2] = 6;  index_face[1][3] = 5;
		index_face[2][0] = 4;  index_face[2][1] = 0;  index_face[2][2] = 1;  index_face[2][3] = 5;
		index_face[3][0] = 7;  index_face[3][1] = 6;  index_face[3][2] = 2;  index_face[3][3] = 3;
		index_face[4][0] = 0;  index_face[4][1] = 3;  index_face[4][2] = 2;  index_face[4][3] = 1;
		index_face[5][0] = 4;  index_face[5][1] = 5;  index_face[5][2] = 6;  index_face[5][3] = 7;
	}
	else if (type_input == "inp") {
		index_face[0][0] = 0;  index_face[0][1] = 1;  index_face[0][2] = 2;  index_face[0][3] = 3;         // !!!-----ALL FACES HAVE INNER NORMAL-----!!! //
		index_face[1][0] = 4;  index_face[1][1] = 7;  index_face[1][2] = 6;  index_face[1][3] = 5;
		index_face[2][0] = 0;  index_face[2][1] = 4;  index_face[2][2] = 5;  index_face[2][3] = 1;
		index_face[3][0] = 1;  index_face[3][1] = 5;  index_face[3][2] = 6;  index_face[3][3] = 2;
		index_face[4][0] = 2;  index_face[4][1] = 6;  index_face[4][2] = 7;  index_face[4][3] = 3;
		index_face[5][0] = 3;  index_face[5][1] = 7;  index_face[5][2] = 4;  index_face[5][3] = 0;
	}

	for (label j = 0; j < 8; j++) {
		index = cells(i, j);
		x[j] = nodes(index, 0);
		y[j] = nodes(index, 1);
		z[j] = nodes(index, 2);
	}

	// center of each cell ---------------------------------
	center(i, 0) = (x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]) / 8.0;
	center(i, 1) = (y[0] + y[1] + y[2] + y[3] + y[4] + y[5] + y[6] + y[7]) / 8.0;
	center(i, 2) = (z[0] + z[1] + z[2] + z[3] + z[4] + z[5] + z[6] + z[7]) / 8.0;

	for (label j = 0; j < 6; j++) {
		for (label k = 0; k < 4; k++) {
			index = index_face[j][k];
			fx[k] = x[index];
			fy[k] = y[index];
			fz[k] = z[index];
		}

		// center of each face ---------------------------------
		faces(i, j).center[0] = (fx[0] + fx[1] + fx[2] + fx[3]) / 4.0;
		faces(i, j).center[1] = (fy[0] + fy[1] + fy[2] + fy[3]) / 4.0;
		faces(i, j).center[2] = (fz[0] + fz[1] + fz[2] + fz[3]) / 4.0;

		// area of each face ---------------------------------
		d0 = sqrt(pow(fx[1] - fx[0], 2) + pow(fy[1] - fy[0], 2) + pow(fz[1] - fz[0], 2));
		d1 = sqrt(pow(fx[2] - fx[1], 2) + pow(fy[2] - fy[1], 2) + pow(fz[2] - fz[1], 2));
		d2 = sqrt(pow(fx[3] - fx[2], 2) + pow(fy[3] - fy[2], 2) + pow(fz[3] - fz[2], 2));
		d3 = sqrt(pow(fx[0] - fx[3], 2) + pow(fy[0] - fy[3], 2) + pow(fz[0] - fz[3], 2));
		d4 = sqrt(pow(fx[3] - fx[1], 2) + pow(fy[3] - fy[1], 2) + pow(fz[3] - fz[1], 2));
		k1 = 0.5 * (d0 + d3 + d4);
		k2 = 0.5 * (d1 + d2 + d4);
		faces(i, j).area = sqrt(k1 * (k1 - d0) * (k1 - d3) * (k1 - d4)) + sqrt(k2 * (k2 - d1) * (k2 - d2) * (k2 - d4));

		// normal of each face ---------------------------------
		a = fx[1] - fx[0]; b = fy[1] - fy[0]; c = fz[1] - fz[0];
		d = fx[2] - fx[0]; e = fy[2] - fy[0]; f = fz[2] - fz[0];
		g = b * f - c * e;
		h = c * d - a * f;
		m = a * e - b * d;
		k = sqrt(g * g + h * h + m * m);
		faces(i, j).normal[0] = g / k;
		faces(i, j).normal[1] = h / k;
		faces(i, j).normal[2] = m / k;

		a = faces(i, j).center[0] - center(i, 0);
		b = faces(i, j).center[1] - center(i, 1);
		c = faces(i, j).center[2] - center(i, 2);

		k = a * faces(i, j).normal[0] + b * faces(i, j).normal[1] + c * faces(i, j).normal[2];

		if (k < 0) {
			faces(i, j).normal[0] = -faces(i, j).normal[0];
			faces(i, j).normal[1] = -faces(i, j).normal[1];
			faces(i, j).normal[2] = -faces(i, j).normal[2];
		}

		// volume of each cell ---------------------------------
		a = fx[0] - center(i, 0); b = fy[0] - center(i, 1); c = fz[0] - center(i, 2);
		d = fx[1] - center(i, 0); e = fy[1] - center(i, 1); f = fz[1] - center(i, 2);
		g = fx[3] - center(i, 0); h = fy[3] - center(i, 1); m = fz[3] - center(i, 2);

		V += fabs(a * (e * m - f * h) - b * (d * m - g * f) + c * (d * h - e * g));

		a = fx[1] - center(i, 0); b = fy[1] - center(i, 1); c = fz[1] - center(i, 2);
		d = fx[2] - center(i, 0); e = fy[2] - center(i, 1); f = fz[2] - center(i, 2);
		g = fx[3] - center(i, 0); h = fy[3] - center(i, 1); m = fz[3] - center(i, 2);

		V += fabs(a * (e * m - f * h) - b * (d * m - g * f) + c * (d * h - e * g));
	}
	volume(i) = V / 6.0;
}

void Zone::geometry_tetra(label& i) {

	label index;
	label index_face[4][3] = { 0 };
	scalar x[4] = { 0.0 };
	scalar y[4] = { 0.0 };
	scalar z[4] = { 0.0 };
	scalar a, b, c, d, e, f, g, h, m, k;

	if (type_input == "inp") {
		index_face[0][0] = 1;  index_face[0][1] = 2;  index_face[0][2] = 3;
		index_face[1][0] = 1;  index_face[1][1] = 4;  index_face[1][2] = 2;
		index_face[2][0] = 2;  index_face[2][1] = 4;  index_face[2][2] = 3;
		index_face[3][0] = 3;  index_face[3][1] = 4;  index_face[3][2] = 1;
	}

	for (label j = 0; j < num_vertex; j++) {
		index = cells(i, j);
		x[j] = nodes(index, 0);
		y[j] = nodes(index, 1);
		z[j] = nodes(index, 2);
	}

	// center of each cell ---------------------------------
	center(i, 0) = (x[0] + x[1] + x[2] + x[3]) / 4.0;
	center(i, 1) = (y[0] + y[1] + y[2] + y[3]) / 4.0;
	center(i, 2) = (z[0] + z[1] + z[2] + z[3]) / 4.0;

	// volume of each cell ---------------------------------
	a = x[1] - x[0]; b = y[1] - y[0]; c = z[1] - z[0];
	d = x[2] - x[0]; e = y[2] - y[0]; f = z[2] - z[0];
	g = x[3] - x[0]; h = y[3] - y[0]; m = z[3] - z[0];

	volume(i) = a * (e * m - f * h) - b * (d * m - g * f) + c * (d * h - e * g);
	volume(i) = volume(i) / 6.0;

	// area, center, normal of each face --------------------------------------
	a = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2) + pow(z[1] - z[0], 2));
	b = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2) + pow(z[2] - z[1], 2));
	c = sqrt(pow(x[0] - x[2], 2) + pow(y[0] - y[2], 2) + pow(z[0] - z[2], 2));
	k = 0.5 * (a + b + c);
	faces(i, 0).area = sqrt(k * (k - a) * (k - b) * (k - c));
	faces(i, 0).center[0] = (x[0] + x[1] + x[2]) / 3.0;
	faces(i, 0).center[1] = (y[0] + y[1] + y[2]) / 3.0;
	faces(i, 0).center[2] = (z[0] + z[1] + z[2]) / 3.0;
	a = x[1] - x[0]; b = y[1] - y[0]; c = z[1] - z[0];
	d = x[2] - x[0]; e = y[2] - y[0]; f = z[2] - z[0];
	g = b * f - c * e;
	h = c * d - a * f;
	m = a * e - b * d;
	k = sqrt(g * g + h * h + m * m);
	faces(i, 0).normal[0] = g / k;
	faces(i, 0).normal[1] = h / k;
	faces(i, 0).normal[2] = m / k;

	a = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2) + pow(z[1] - z[0], 2));
	b = sqrt(pow(x[3] - x[1], 2) + pow(y[3] - y[1], 2) + pow(z[3] - z[1], 2));
	c = sqrt(pow(x[0] - x[3], 2) + pow(y[0] - y[3], 2) + pow(z[0] - z[3], 2));
	k = 0.5 * (a + b + c);
	faces(i, 1).area = sqrt(k * (k - a) * (k - b) * (k - c));
	faces(i, 1).center[0] = (x[0] + x[1] + x[3]) / 3.0;
	faces(i, 1).center[1] = (y[0] + y[1] + y[3]) / 3.0;
	faces(i, 1).center[2] = (z[0] + z[1] + z[3]) / 3.0;
	a = x[1] - x[0]; b = y[1] - y[0]; c = z[1] - z[0];
	d = x[3] - x[0]; e = y[3] - y[0]; f = z[3] - z[0];
	g = b * f - c * e;
	h = c * d - a * f;
	m = a * e - b * d;
	k = sqrt(g * g + h * h + m * m);
	faces(i, 1).normal[0] = g / k;
	faces(i, 1).normal[1] = h / k;
	faces(i, 1).normal[2] = m / k;

	a = sqrt(pow(x[2] - x[0], 2) + pow(y[2] - y[0], 2) + pow(z[2] - z[0], 2));
	b = sqrt(pow(x[3] - x[2], 2) + pow(y[3] - y[2], 2) + pow(z[3] - z[2], 2));
	c = sqrt(pow(x[0] - x[3], 2) + pow(y[0] - y[3], 2) + pow(z[0] - z[3], 2));
	k = 0.5 * (a + b + c);
	faces(i, 2).area = sqrt(k * (k - a) * (k - b) * (k - c));
	faces(i, 2).center[0] = (x[0] + x[2] + x[3]) / 3.0;
	faces(i, 2).center[1] = (y[0] + y[2] + y[3]) / 3.0;
	faces(i, 2).center[2] = (z[0] + z[2] + z[3]) / 3.0;
	a = x[2] - x[0]; b = y[2] - y[0]; c = z[2] - z[0];
	d = x[3] - x[0]; e = y[3] - y[0]; f = z[3] - z[0];
	g = b * f - c * e;
	h = c * d - a * f;
	m = a * e - b * d;
	k = sqrt(g * g + h * h + m * m);
	faces(i, 2).normal[0] = g / k;
	faces(i, 2).normal[1] = h / k;
	faces(i, 2).normal[2] = m / k;

	a = sqrt(pow(x[2] - x[1], 2) + pow(y[2] - y[1], 2) + pow(z[2] - z[1], 2));
	b = sqrt(pow(x[3] - x[2], 2) + pow(y[3] - y[2], 2) + pow(z[3] - z[2], 2));
	c = sqrt(pow(x[1] - x[3], 2) + pow(y[1] - y[3], 2) + pow(z[1] - z[3], 2));
	k = 0.5 * (a + b + c);
	faces(i, 3).area = sqrt(k * (k - a) * (k - b) * (k - c));
	faces(i, 3).center[0] = (x[1] + x[2] + x[3]) / 3.0;
	faces(i, 3).center[1] = (y[1] + y[2] + y[3]) / 3.0;
	faces(i, 3).center[2] = (z[1] + z[2] + z[3]) / 3.0;
	a = x[2] - x[1]; b = y[2] - y[1]; c = z[2] - z[1];
	d = x[3] - x[1]; e = y[3] - y[1]; f = z[3] - z[1];
	g = b * f - c * e;
	h = c * d - a * f;
	m = a * e - b * d;
	k = sqrt(g * g + h * h + m * m);
	faces(i, 3).normal[0] = g / k;
	faces(i, 3).normal[1] = h / k;
	faces(i, 3).normal[2] = m / k;

	for (label j = 0; j < num_face; j++) {
		a = faces(i, j).center[0] - center(i, 0);
		b = faces(i, j).center[1] - center(i, 1);
		c = faces(i, j).center[2] - center(i, 2);
		k = a * faces(i, j).normal[0] + b * faces(i, j).normal[1] + c * faces(i, j).normal[2];
		if (k < 0) {
			faces(i, j).normal[0] = -faces(i, j).normal[0];
			faces(i, j).normal[1] = -faces(i, j).normal[1];
			faces(i, j).normal[2] = -faces(i, j).normal[2];
		}
	}
}

void Zone::tune_decimal(Array<scalar>& array) {
	scalar num;
	num = pow(10, if_tune_decimal);
	for (label i = 0; i < array.size; i++) {
		array(i) = round(array(i) * num) / num;
	}
}

void Zone::tune_decimal(scalar& value) {
	scalar num;
	num = pow(10, if_tune_decimal);
	value = round(value * num) / num;
}

void Zone::tune_scale(Array<scalar>& array) {
	for (label i = 0; i < array.size; i++) {
		array(i) = array(i) * if_tune_scale;
	}
}

void Zone::reverse_axis(Array<scalar>& array) {
	if (if_reverse_axis == "all") {
		for (label i = 0; i < array.size; i++) {
			array(i) = -array(i);
		}
	}
}

void Zone::offset_axis(Array<scalar>& array) {
	for (label i = 0; i < array.size; i++) {
		array(i) += if_offset_axis;
	}
}
