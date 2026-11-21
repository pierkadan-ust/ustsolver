// Sparse Matrix Storage for FVM
// 
// store stiffness coefficients as { [ face_11, face_12, cell_1 ]; 
//                                   [ face_21, face_22, cell_2 ]; }
// 
// with the presence of vector boundary condition { [ face_11,  face_12,  cell_1 ];
//                                                  [ face_21x, face_22x, cell_2x, face_21y, face_22y, cell_2y, face_21z, face_22z, cell_2z ]; }

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "field/FaceField.hpp"

class SparseMat {

public:

	Array<vector> matrix;
	Array<label> num_column;
	Array<label> if_boundary;
	label num_row = 0;

public:

	SparseMat() {}

	scalar& operator () (label i) {
		return matrix(i)[num_column(i) - 1];
	}

	scalar& operator () (label i, label j) {
		return matrix(i)[j];
	}

	scalar& operator () (label i, char drc) {
		if (drc == 'y') { return matrix(i)[2 * num_column(i) - 1]; }
		else if (drc == 'z') { return matrix(i)[3 * num_column(i) - 1]; }
		else { return matrix(i)[num_column(i) - 1]; }
	}

	scalar& operator () (label i, label j, char drc) {
		if (drc == 'y') { return matrix(i)[j + num_column(i)]; }
		else if (drc == 'z') { return matrix(i)[j + num_column(i) * 2]; }
		else { return matrix(i)[j]; }
	}

	void resize(label row) {
		matrix.resize(row);
		num_column.resize(row);
		if_boundary.resize(row);
		num_row = row;
	}

	void resize(label row, label column) {
		matrix(row).resize(column);
		num_column(row) = column;
	}

	void resize(Array<FaceField>& faces) {
		for (label i = 0; i < num_row; i++) {
			for (label j = 0; j < num_column(i) - 1; j++) {
				if (faces(i, j).cell_nb < 0) {
					matrix(i).resize(num_column(i) * num_dimension);
					if_boundary(i) = 1;
					break;
				}
			}
		}
	}

	void set_all(scalar value) {
		for (label i = 0; i < num_row; i++) {
			for (size_t j = 0; j < matrix(i).size(); j++) {
				matrix(i)[j] = value;
			}
		}
	}

	void print() {
		std::cout << "\n";
		for (label i = 0; i < num_row; i++) {
			for (size_t j = 0; j < matrix(i).size(); j++) {
				std::cout << matrix(i)[j] << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
};
