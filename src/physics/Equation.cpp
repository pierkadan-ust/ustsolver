#include "physics/Equation.hpp"

Equation::Equation(label ifPrint) {
	if (ifPrint) { std::cout << "\nCalculation : "; }
}

void Equation::resize_matrix(Variable* variable) {
	auto& num_cells = variable->ZONE->num_cells;
	label length;
	A.resize(num_cells);
	for (label row = 0; row < num_cells; row++) {
		length = variable->ZONE->types(row, 1) + 1;
		A.resize(row, length);
	}
	if (variable->TYPE == "scalar") { B.resize(num_cells); }
	else if (variable->TYPE == "vector") { B.resize(num_cells, num_dimension);  A.resize(variable->ZONE->faces); }
}

void Equation::get_diagonal(Variable* variable, string type) {
	if (type == "original") {
		for (label i = 0; i < (*variable).X.size; i++) {
			(*variable).X(i) = this->A(i);
		}
	}
	else if (type == "SIMPLE") {
		for (label i = 0; i < (*variable).X.size; i++) {
			(*variable).X(i) = (*variable).ZONE->volume(i) / this->A(i);
		}
	}
}

void Equation::steadyFlow(Variable* VT, Variable* P, Variable* RHO, Variable* MU, label if_print, scalar sor_momt, scalar sor_p, label iter_momt, label iter_p) {
	std::cout << "<Incompressible Fluid>, solving for [ " << VT->NAME << ", " << P->NAME << " ],";
	label num_cells = VT->ZONE->num_cells;
	label num_faces = VT->ZONE->num_faces;
	resize_matrix(VT);
	VT->get_face();

	Variable diagA(P->ZONE, "scalar", "Matrix Diagonal", 0);   diagA.BC_copy(P);  diagA.BC_all("zero grad");
	Variable PC(P->ZONE, "scalar", "Pressure Correction", 0);  PC.BC_copy(P);     PC.BC_all("zero grad");
	Array<scalar> gradFacePC(num_faces, num_dimension);
	label iter, final_iter_momt, final_iter_p;

	for (iter = 0; iter < max_iter_V; iter++) {
		// [1] Momentum Equation
		laplacian(VT, MU);
		div(VT, VT, RHO);
		P->get_face();
		P->get_gradient(&B, "SIMPLE");
		final_iter_momt = Solve::GaussSeidel(tresMomt, resMomt, A, B, VT, sor_momt, iter_momt, 0, "Momentum", 0);  B.set_all(0);

		// [2] Pressure Correction Equation
		get_diagonal(&diagA, "SIMPLE");
		A.set_all(0);
		diagA.get_face();
		poission(&diagA);
		VT->get_face_PWIM(P, &diagA);
		VT->get_divergence(nullptr, "SIMPLE");
		final_iter_p = Solve::GaussSeidel(tresP, resP, A, VT->X_divg, &PC, 1.0, iter_p, 0, "Pressure", if_reference_p);  A.set_all(0);

		// [3] Update
		PC.get_face();
		PC.get_gradient();
		PC.get_gradient_face(gradFacePC);
		P->X += PC.X * sor_p;
		VT->X -= PC.X_grad * diagA.X;
		VT->X_faces -= gradFacePC * diagA.X_faces;

		if (if_print) { std::cout << "\n              Iteration = " << iter + 1 << ", Res_Momentum = " << resMomt << ", Res_Continuity = " << resP; }
		if (isnan(resP)) { std::cerr << "\n\n              ERROR : NAN value occured under bad setup, please recheck calculation setup\n";  return; }
	}
	if (if_print) { std::cout << "  Finished\n"; }
	else { std::cout << " Final Iteration = " << iter << ", Res_Momentum = " << resMomt << ", Res_Continuity = " << resP << "\n"; }
}

void Equation::steadyEnergy(Variable* TEMP, Variable* VT, Variable* RHO, Variable* K, Variable* CP, label if_print, scalar sor) {
	std::cout << "<Heat Convection>, solving for [ " << TEMP->NAME << ", " << VT->NAME << " ],";
	resize_matrix(TEMP);

	laplacian(TEMP, K);     // [1]
	VT->get_face();
	div(TEMP, VT, RHO, CP); // [2]

	label iter;
	iter = Solve::GaussSeidel(tresEnergy, resEnergy, A, B, TEMP, sor, max_iter_T, if_print);

	std::cout << " Final Iteration = " << iter << ", Residual = " << resEnergy << "\n";
}

void Equation::steadyConduction(Variable* TEMP, Variable* K, label if_print, scalar sor) {
	std::cout << "<Heat Conduction>, solving for [ " << TEMP->NAME << " ],";
	resize_matrix(TEMP);

	laplacian(TEMP, K); // [1]

	label iter;
	iter = Solve::GaussSeidel(tresEnergy, resEnergy, A, B, TEMP, sor, max_iter_T, if_print);

	std::cout << " Final Iteration = " << iter << ", Residual = " << resEnergy << "\n";
}

void Equation::transientFlow(Variable* VT, Variable* P, Variable* RHO, Variable* MU, scalar sor_momt, scalar sor_p, label iter_momt, label iter_p) {

	Variable diagA(P->ZONE, "scalar", "Matrix Diagonal", 0);   diagA.BC_copy(P);  diagA.BC_all("zero grad");
	Variable PC(P->ZONE, "scalar", "Pressure Correction", 0);  PC.BC_copy(P);     PC.BC_all("zero grad");
	Array<scalar> gradFacePC(P->ZONE->num_faces, num_dimension);
	Array<scalar> lastVT(P->ZONE->num_cells, num_dimension);   lastVT.set_array(VT->X);
	VT->get_face();

	for (label iter = 0; iter < max_iter_flow; iter++) {
		// [1] Momentum Equation
		laplacian(VT, MU);
		div(VT, VT, RHO);
		P->get_face();
		P->get_gradient(&B, "SIMPLE");
		flowForwardEulerFirst(&lastVT, RHO);
		Solve::GaussSeidel(tresMomt, resMomt, A, B, VT, sor_momt, iter_momt, 0, "Momentum", 0);  B.set_all(0);

		// [2] Pressure Correction Equation
		get_diagonal(&diagA, "SIMPLE");
		A.set_all(0);
		PC.X.set_all(0);
		diagA.get_face();
		poission(&diagA);
		VT->get_face_PWIM(P, &diagA);
		VT->get_divergence(nullptr, "SIMPLE");
		Solve::GaussSeidel(tresP, resP, A, VT->X_divg, &PC, 1.0, iter_p, 0, "Pressure", if_reference_p);  A.set_all(0);

		// [3] Update
		PC.get_face();
		PC.get_gradient();
		PC.get_gradient_face(gradFacePC);
		P->X += PC.X * sor_p;
		VT->X -= PC.X_grad * diagA.X;
		VT->X_faces -= gradFacePC * diagA.X_faces;

		if (isnan(resP)) { std::cerr << "\n\n              ERROR : NAN value occured under bad setup, please recheck calculation setup\n";  return; }
	}
}

void Equation::transientEnergy(Variable* TEMP, Variable* VT, Variable* RHO, Variable* K, Variable* CP, scalar sor) {

	A.set_all(0);
	B.set_all(0);

	VT->get_face();
	RHO->get_face();
	K->get_face();
	CP->get_face();

	laplacian(TEMP, K);                     // [1]                 
	div(TEMP, VT, RHO, CP);                 // [2]
	energyForwardEulerFirst(TEMP, RHO, CP); // [3]  

	Solve::GaussSeidel(tresEnergy, resEnergy, A, B, TEMP, sor, max_iter_energy, 0);
}

void Equation::flowForwardEulerFirst(Array<scalar>* lastVT, Variable* RHO) {
	scalar value = 0;
	for (label i = 0; i < RHO->ZONE->num_cells; i++) {
		value = RHO->X(i) * RHO->ZONE->volume(i) / dt_flow;
		A(i) += value;
		if (A.if_boundary(i)) {
			A(i, 'y') += value;
			if (num_dimension == 3) { A(i, 'z') += value; }
		}
		for (label j = 0; j < num_dimension; j++) {
			B(i, j) += value * (*lastVT)(i, j);
		}
	}
}

void Equation::energyForwardEulerFirst(Variable* TEMP, Variable* RHO, Variable* CP) {
	scalar value = 0;
	for (label i = 0; i < TEMP->ZONE->num_cells; i++) {
		value = RHO->X(i) * CP->X(i) * TEMP->ZONE->volume(i) / dt_energy;
		A(i) += value;
		B(i) += value * TEMP->X(i);
	}
}

void Equation::laplacian(Variable* phi, Variable* k) {

	auto& type = phi->TYPE;
	auto& faces = phi->ZONE->faces;
	auto& num_faces = phi->ZONE->num_faces;
	auto& k_faces = k->X_faces;
	auto& bc = phi->BC;
	scalar coeff;
	label index;

	if (type == "scalar") {
		for (label i = 0; i < num_faces; i++) {

			if (faces(i).if_repeat == 1) { continue; }

			if (faces(i).cell_nb >= 0) {  // inner faces
				coeff = -k_faces(i) * faces(i).area / faces(i).distance;
				A(faces(i).cell, faces(i).side) += coeff;
				A(faces(i).cell) -= coeff;
				A(faces(i).cell_nb, faces(i).side_nb) += coeff;
				A(faces(i).cell_nb) -= coeff;
			}
			else {  // boundary faces
				index = -3 * faces(i).cell_nb - 3;
				if (bc.size() == 0) {}
				else if (bc[index] == 1) {  // BC_wall_temperature
					coeff = k_faces(i) * faces(i).area / faces(i).distance;
					A(faces(i).cell) += coeff * bc[index + 1];
					B(faces(i).cell) += coeff * bc[index + 2];
				}
				else if (bc[index] == 2) {  // BC_wall_heatflux
					B(faces(i).cell) += bc[index + 2] * faces(i).area;
				}
				else if (bc[index] == 3) {  // BC_wall_convection
					coeff = k_faces(i) / faces(i).distance;
					coeff = faces(i).area * (bc[index + 1] * coeff) / (bc[index + 1] + coeff);
					A(faces(i).cell) += coeff;
					B(faces(i).cell) += coeff * bc[index + 2];
				}
			}
		}
	}

	else if (type == "vector") {
		if (num_dimension == 2) {
			for (label i = 0; i < num_faces; i++) {
				if (faces(i).if_repeat == 1) { continue; }

				if (faces(i).cell_nb >= 0) {  // inner faces
					coeff = -k_faces(i) * faces(i).area / faces(i).distance;
					A(faces(i).cell, faces(i).side) += coeff;
					A(faces(i).cell) -= coeff;
					A(faces(i).cell_nb, faces(i).side_nb) += coeff;
					A(faces(i).cell_nb) -= coeff;

					if (A.if_boundary(faces(i).cell)) {
						A(faces(i).cell, faces(i).side, 'y') += coeff;
						A(faces(i).cell, 'y') -= coeff;
					}
					if (A.if_boundary(faces(i).cell_nb)) {
						A(faces(i).cell_nb, faces(i).side_nb, 'y') += coeff;
						A(faces(i).cell_nb, 'y') -= coeff;
					}
				}
				else {  // boundary faces
					index = -4 * faces(i).cell_nb - 4;
					if (bc.size() == 0) {}
					else if (bc[index] == 1 || bc[index] == 2) {  // BC_wall_velocity or BC_velocity_inlet
						coeff = k_faces(i) * faces(i).area / faces(i).distance;
						A(faces(i).cell) += coeff;
						A(faces(i).cell, 'y') += coeff;
						B(faces(i).cell, 0) += coeff * bc[index + 2];
						B(faces(i).cell, 1) += coeff * bc[index + 3];
					}
					else if (bc[index] == 3) {  // BC_wall_rotation  
						coeff = k_faces(i) * faces(i).area / faces(i).distance;
						A(faces(i).cell) += coeff;
						A(faces(i).cell, 'y') += coeff;
						B(faces(i).cell, 0) += coeff * phi->wallRotation(i, 0);
						B(faces(i).cell, 1) += coeff * phi->wallRotation(i, 1);
					}
				}
			}
		}

		else if (num_dimension == 3) {
			for (label i = 0; i < num_faces; i++) {
				if (faces(i).if_repeat == 1) { continue; }

				if (faces(i).cell_nb >= 0) {  // inner faces
					coeff = -k_faces(i) * faces(i).area / faces(i).distance;
					A(faces(i).cell, faces(i).side) += coeff;
					A(faces(i).cell) -= coeff;
					A(faces(i).cell_nb, faces(i).side_nb) += coeff;
					A(faces(i).cell_nb) -= coeff;

					if (A.if_boundary(faces(i).cell)) {
						A(faces(i).cell, faces(i).side, 'y') += coeff;
						A(faces(i).cell, 'y') -= coeff;
						A(faces(i).cell, faces(i).side, 'z') += coeff;
						A(faces(i).cell, 'z') -= coeff;
					}
					if (A.if_boundary(faces(i).cell_nb)) {
						A(faces(i).cell_nb, faces(i).side_nb, 'y') += coeff;
						A(faces(i).cell_nb, 'y') -= coeff;
						A(faces(i).cell_nb, faces(i).side_nb, 'z') += coeff;
						A(faces(i).cell_nb, 'z') -= coeff;
					}
				}
				else {  // boundary faces
					index = -5 * faces(i).cell_nb - 5;
					if (bc.size() == 0) {}
					else if (bc[index] == 1 || bc[index] == 2) {  // BC_wall_velocity or BC_velocity_inlet
						coeff = k_faces(i) * faces(i).area / faces(i).distance;
						A(faces(i).cell) += coeff;
						A(faces(i).cell, 'y') += coeff;
						A(faces(i).cell, 'z') += coeff;
						B(faces(i).cell, 0) += coeff * bc[index + 2];
						B(faces(i).cell, 1) += coeff * bc[index + 3];
						B(faces(i).cell, 2) += coeff * bc[index + 4];
					}
					else if (bc[index] == 3) {  // BC_wall_rotation
						coeff = k_faces(i) * faces(i).area / faces(i).distance;
						A(faces(i).cell) += coeff;
						A(faces(i).cell, 'y') += coeff;
						A(faces(i).cell, 'z') += coeff;
						B(faces(i).cell, 0) += coeff * phi->wallRotation(i, 0);
						B(faces(i).cell, 1) += coeff * phi->wallRotation(i, 1);
						B(faces(i).cell, 2) += coeff * phi->wallRotation(i, 2);
					}
				}
			}
		}
	}
}

void Equation::div(Variable* phi, Variable* U, Variable* rho, Variable* cp) {

	auto& type = phi->TYPE;
	auto& faces = phi->ZONE->faces;
	auto& num_faces = phi->ZONE->num_faces;
	auto& bc = phi->BC_exact;
	auto& bc_type = phi->BC_type;
	auto& U_faces = U->X_faces;
	auto& rho_faces = rho->X_faces;
	label index;
	scalar mf;

	if (type == "scalar") {
		for (label i = 0; i < num_faces; i++) {

			if (faces(i).if_repeat == 1) { continue; }
			mf = 0;
			for (label j = 0; j < num_dimension; j++) { mf += U_faces(i, j) * faces(i).normal[j]; }
			if (mf == 0) { continue; }
			mf = mf * faces(i).area * rho_faces(i);
			if (cp != nullptr) { mf = mf * cp->X_faces(i); }

			if (faces(i).cell_nb >= 0) {  // inner faces
				if (mf < 0) {
					A(faces(i).cell, faces(i).side) += mf;
					A(faces(i).cell_nb) -= mf;
				}
				else {
					A(faces(i).cell) += mf;
					A(faces(i).cell_nb, faces(i).side_nb) -= mf;
				}
			}
			else {  // boundary faces
				if (mf < 0) {
					index = -faces(i).cell_nb - 1;
					if (bc.size() == 0) { B(faces(i).cell) -= mf * phi->X(faces(i).cell); }
					else if (!isnan(bc[index])) { B(faces(i).cell) -= mf * bc[index]; }  // BC_exact_value
					else if (bc_type[index] == "zero grad") { B(faces(i).cell) -= mf * phi->X(faces(i).cell); }  // zero-gradient interpolation	
				}
				else { A(faces(i).cell) += mf; }
			}
		}
	}

	else if (type == "vector") {
		if (num_dimension == 2) {
			for (label i = 0; i < num_faces; i++) {

				if (faces(i).if_repeat == 1) { continue; }
				mf = 0;
				for (label j = 0; j < num_dimension; j++) { mf += U_faces(i, j) * faces(i).normal[j]; }
				if (mf == 0) { continue; }
				mf = mf * faces(i).area * rho_faces(i);
				if (cp != nullptr) { mf = mf * cp->X_faces(i); }

				if (faces(i).cell_nb >= 0) {  // inner faces
					if (mf < 0) {
						A(faces(i).cell, faces(i).side) += mf;
						A(faces(i).cell_nb) -= mf;

						if (A.if_boundary(faces(i).cell)) {
							A(faces(i).cell, faces(i).side, 'y') += mf;
						}
						if (A.if_boundary(faces(i).cell_nb)) {
							A(faces(i).cell_nb, 'y') -= mf;
						}
					}
					else {
						A(faces(i).cell) += mf;
						A(faces(i).cell_nb, faces(i).side_nb) -= mf;

						if (A.if_boundary(faces(i).cell)) {
							A(faces(i).cell, 'y') += mf;
						}
						if (A.if_boundary(faces(i).cell_nb)) {
							A(faces(i).cell_nb, faces(i).side_nb, 'y') -= mf;
						}
					}
				}
				else {  // boundary faces
					if (mf < 0) {
						index = num_dimension * (-faces(i).cell_nb - 1);
						if (bc.size() == 0) {
							B(faces(i).cell, 0) -= mf * phi->X(faces(i).cell, 0);
							B(faces(i).cell, 1) -= mf * phi->X(faces(i).cell, 1);
						}
						else if (!isnan(bc[index])) {  // BC_exact_value
							B(faces(i).cell, 0) -= mf * bc[index];
							B(faces(i).cell, 1) -= mf * bc[index + 1];
						}
						else if (bc_type[-faces(i).cell_nb - 1] == "zero grad") {  // zero-gradient interpolation
							B(faces(i).cell, 0) -= mf * phi->X(faces(i).cell, 0);
							B(faces(i).cell, 1) -= mf * phi->X(faces(i).cell, 1);
						}
					}
					else {
						A(faces(i).cell) += mf;
						A(faces(i).cell, 'y') += mf;
					}
				}
			}
		}

		else if (num_dimension == 3) {
			for (label i = 0; i < num_faces; i++) {

				if (faces(i).if_repeat == 1) { continue; }
				mf = 0;
				for (label j = 0; j < num_dimension; j++) { mf += U_faces(i, j) * faces(i).normal[j]; }
				if (mf == 0) { continue; }
				mf = mf * faces(i).area * rho_faces(i);
				if (cp != nullptr) { mf = mf * cp->X_faces(i); }

				if (faces(i).cell_nb >= 0) {  // inner faces
					if (mf < 0) {
						A(faces(i).cell, faces(i).side) += mf;
						A(faces(i).cell_nb) -= mf;

						if (A.if_boundary(faces(i).cell)) {
							A(faces(i).cell, faces(i).side, 'y') += mf;
							A(faces(i).cell, faces(i).side, 'z') += mf;
						}
						if (A.if_boundary(faces(i).cell_nb)) {
							A(faces(i).cell_nb, 'y') -= mf;
							A(faces(i).cell_nb, 'z') -= mf;
						}
					}
					else {
						A(faces(i).cell) += mf;
						A(faces(i).cell_nb, faces(i).side_nb) -= mf;

						if (A.if_boundary(faces(i).cell)) {
							A(faces(i).cell, 'y') += mf;
							A(faces(i).cell, 'z') += mf;
						}
						if (A.if_boundary(faces(i).cell_nb)) {
							A(faces(i).cell_nb, faces(i).side_nb, 'y') -= mf;
							A(faces(i).cell_nb, faces(i).side_nb, 'z') -= mf;
						}
					}
				}
				else {  // boundary faces
					if (mf < 0) {
						index = num_dimension * (-faces(i).cell_nb - 1);
						if (bc.size() == 0) {
							B(faces(i).cell, 0) -= mf * phi->X(faces(i).cell, 0);
							B(faces(i).cell, 1) -= mf * phi->X(faces(i).cell, 1);
							B(faces(i).cell, 2) -= mf * phi->X(faces(i).cell, 2);
						}
						else if (!isnan(bc[index])) {  // BC_exact_value
							B(faces(i).cell, 0) -= mf * bc[index];
							B(faces(i).cell, 1) -= mf * bc[index + 1];
							B(faces(i).cell, 2) -= mf * bc[index + 2];
						}
						else if (bc_type[-faces(i).cell_nb - 1] == "zero grad") {  // zero-gradient interpolation
							B(faces(i).cell, 0) -= mf * phi->X(faces(i).cell, 0);
							B(faces(i).cell, 1) -= mf * phi->X(faces(i).cell, 1);
							B(faces(i).cell, 2) -= mf * phi->X(faces(i).cell, 2);
						}
					}
					else {
						A(faces(i).cell) += mf;
						A(faces(i).cell, 'y') += mf;
						A(faces(i).cell, 'z') += mf;
					}
				}
			}
		}
	}
}

void Equation::poission(Variable* diagA) {  // obtain Stiffness Matrix of Pressure Correction Equation in SIMPLE

	auto& faces = diagA->ZONE->faces;
	auto& num_faces = diagA->ZONE->num_faces;
	auto& face_property = diagA->X_faces;
	scalar coeff;

	if (num_dimension == 2) {
		for (label i = 0; i < num_faces; i++) {

			if (faces(i).if_repeat == 1) { continue; }

			if (faces(i).cell_nb >= 0) {
				coeff = -face_property(i) / faces(i).distance;
				A(faces(i).cell, faces(i).side) += coeff;
				A(faces(i).cell) -= coeff;
				A(faces(i).cell_nb, faces(i).side_nb) += coeff;
				A(faces(i).cell_nb) -= coeff;

				if (A.if_boundary(faces(i).cell)) {
					A(faces(i).cell, faces(i).side, 'y') += coeff;
					A(faces(i).cell, 'y') -= coeff;
				}
				if (A.if_boundary(faces(i).cell_nb)) {
					A(faces(i).cell_nb, faces(i).side_nb, 'y') += coeff;
					A(faces(i).cell_nb, 'y') -= coeff;
				}
			}
		}
	}

	else if (num_dimension == 3) {
		for (label i = 0; i < num_faces; i++) {

			if (faces(i).if_repeat == 1) { continue; }

			if (faces(i).cell_nb >= 0) {
				coeff = -face_property(i) / faces(i).distance;
				A(faces(i).cell, faces(i).side) += coeff;
				A(faces(i).cell) -= coeff;
				A(faces(i).cell_nb, faces(i).side_nb) += coeff;
				A(faces(i).cell_nb) -= coeff;

				if (A.if_boundary(faces(i).cell)) {
					A(faces(i).cell, faces(i).side, 'y') += coeff;
					A(faces(i).cell, 'y') -= coeff;
					A(faces(i).cell, faces(i).side, 'z') += coeff;
					A(faces(i).cell, 'z') -= coeff;
				}
				if (A.if_boundary(faces(i).cell_nb)) {
					A(faces(i).cell_nb, faces(i).side_nb, 'y') += coeff;
					A(faces(i).cell_nb, 'y') -= coeff;
					A(faces(i).cell_nb, faces(i).side_nb, 'z') += coeff;
					A(faces(i).cell_nb, 'z') -= coeff;
				}
			}
		}
	}
}
