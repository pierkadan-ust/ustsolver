#include "solver/Solve.hpp"

label Solve::GaussSeidel(scalar& tres, scalar& residual, SparseMat& A, Array<scalar>& B, Variable* variable, scalar SOR, label num_iter, label if_print_res, string name, label if_reference) {

	auto& type = variable->TYPE;
	label iter = 0;
	scalar sum, last;

	if (type == "scalar") {
		if (if_reference) { A(1) = 2.0 * A(1); }  // set first_cell_value as zero
		for (iter = 0; iter < num_iter; iter++) {
			residual = 0.0;
			for (label i = 0; i < A.num_row; i++) {
				sum = 0.0;
				if (iter == 0 && SOR != 1.0) { A(i) = A(i) / SOR; }  // relaxation
				for (label j = 0; j < A.num_column(i) - 1; j++) {
					if (variable->ZONE->faces(i, j).cell_nb >= 0) {
						sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb);
					}
				}
				last = variable->X(i);
				variable->X(i) = (B(i) + last * A(i) * (1 - SOR) - sum) / A(i);
				last -= variable->X(i);
				residual = fabs_max(residual, last);
			}
			if (if_print_res) { std::cout << "\n              " << name << " = " << iter + 1 << ", Residual = " << residual; }
			if (residual <= tres || iter + 1 == num_iter) { if (if_print_res) { std::cout << "\n             "; } break; }  // end condition
		}
	}

	else if (type == "vector") {
		if (num_dimension == 2) {
			if (if_reference) { A(1) = 2.0 * A(1);  A(1, 'y') = 2.0 * A(1, 'y'); }  // set first_cell_value as zero
			for (iter = 0; iter < num_iter; iter++) {
				residual = 0.0;
				for (label i = 0; i < A.num_row; i++) {
					// X-direction
					sum = 0;
					if (iter == 0 && SOR != 1.0) { A(i) = A(i) / SOR; }  // relaxation
					for (label j = 0; j < A.num_column(i) - 1; j++) {
						if (variable->ZONE->faces(i, j).cell_nb >= 0) {
							sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb, 0);
						}
					}
					last = variable->X(i, 0);
					variable->X(i, 0) = (B(i, 0) + last * A(i) * (1 - SOR) - sum) / A(i);
					last -= variable->X(i, 0);
					residual = fabs_max(residual, last);
					// Y-direction
					sum = 0;
					if (A.if_boundary(i)) {
						if (iter == 0 && SOR != 1.0) { A(i, 'y') = A(i, 'y') / SOR; }  // relaxation
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j, 'y') * variable->X(variable->ZONE->faces(i, j).cell_nb, 1);
							}
						}
						last = variable->X(i, 1);
						variable->X(i, 1) = (B(i, 1) + last * A(i, 'y') * (1 - SOR) - sum) / A(i, 'y');
						last -= variable->X(i, 1);
					}
					else {
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb, 1);
							}
						}
						last = variable->X(i, 1);
						variable->X(i, 1) = (B(i, 1) + last * A(i) * (1 - SOR) - sum) / A(i);
						last -= variable->X(i, 1);
					}
					residual = fabs_max(residual, last);
				}
				if (if_print_res) { std::cout << "\n              " << name << " = " << iter + 1 << ", Residual = " << residual; }
				if (residual <= tres || iter + 1 == num_iter) { if (if_print_res) { std::cout << "\n             "; } break; }  // end condition
			}
		}

		else if (num_dimension == 3) {
			if (if_reference) { A(1) = 2.0 * A(1);  A(1, 'y') = 2.0 * A(1, 'y');  A(1, 'z') = 2.0 * A(1, 'z'); }  // set first_cell_value as zero
			for (iter = 0; iter < num_iter; iter++) {
				residual = 0.0;
				for (label i = 0; i < A.num_row; i++) {
					// X-direction
					sum = 0;
					if (iter == 0 && SOR != 1.0) { A(i) = A(i) / SOR; }  // relaxation
					for (label j = 0; j < A.num_column(i) - 1; j++) {
						if (variable->ZONE->faces(i, j).cell_nb >= 0) {
							sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb, 0);
						}
					}
					last = variable->X(i, 0);
					variable->X(i, 0) = (B(i, 0) + last * A(i) * (1 - SOR) - sum) / A(i);
					last -= variable->X(i, 0);
					residual = fabs_max(residual, last);
					// Y-direction
					sum = 0;
					if (A.if_boundary(i)) {
						if (iter == 0 && SOR != 1.0) { A(i, 'y') = A(i, 'y') / SOR; }  // relaxation
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j, 'y') * variable->X(variable->ZONE->faces(i, j).cell_nb, 1);
							}
						}
						last = variable->X(i, 1);
						variable->X(i, 1) = (B(i, 1) + last * A(i, 'y') * (1 - SOR) - sum) / A(i, 'y');
						last -= variable->X(i, 1);
					}
					else {
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb, 1);
							}
						}
						last = variable->X(i, 1);
						variable->X(i, 1) = (B(i, 1) + last * A(i) * (1 - SOR) - sum) / A(i);
						last -= variable->X(i, 1);
					}
					residual = fabs_max(residual, last);
					// Z-direction
					sum = 0;
					if (A.if_boundary(i)) {
						if (iter == 0 && SOR != 1.0) { A(i, 'z') = A(i, 'z') / SOR; }  // relaxation
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j, 'z') * variable->X(variable->ZONE->faces(i, j).cell_nb, 2);
							}
						}
						last = variable->X(i, 2);
						variable->X(i, 2) = (B(i, 2) + last * A(i, 'z') * (1 - SOR) - sum) / A(i, 'z');
						last -= variable->X(i, 2);
					}
					else {
						for (label j = 0; j < A.num_column(i) - 1; j++) {
							if (variable->ZONE->faces(i, j).cell_nb >= 0) {
								sum += A(i, j) * variable->X(variable->ZONE->faces(i, j).cell_nb, 2);
							}
						}
						last = variable->X(i, 2);
						variable->X(i, 2) = (B(i, 2) + last * A(i) * (1 - SOR) - sum) / A(i);
						last -= variable->X(i, 2);
					}
					residual = fabs_max(residual, last);
				}
				if (if_print_res) { std::cout << "\n              " << name << " = " << iter + 1 << ", Residual = " << residual; }
				if (residual <= tres || iter + 1 == num_iter) { if (if_print_res) { std::cout << "\n             "; } break; }  // end condition
			}
		}
	}

	return iter + 1;
}
