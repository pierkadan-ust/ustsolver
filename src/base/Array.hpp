// Array Storage
// store 1D vector as { [a1, a2] } 
// store 2D vector as { [a11, a12], [a21, a22] }
// store 3D vector as { [(a111, a112), (a121, a122)], [(a211, a212), (a221, a222)] }

#pragma once

#include "base/Global.hpp"

template <typename type>
class Array {

public:

    label NX, NY, NZ;
    label size;
    std::vector<type> data;

public:

    Array(label nx = 0, label ny = 1, label nz = 1) : NX(nx), NY(ny), NZ(nz) {
        size = NX * NY * NZ;
        data.resize(size);
    }

    type& operator () (label i) {
        return data[i];
    }

    type& operator () (label i, label j) {
        return data[i * NY + j];
    }

    type& operator () (label i, label j, label k) {
        return data[i * NY * NZ + j * NZ + k];
    }

    Array<scalar>& operator * (const scalar& value) {
        for (label i = 0; i < this->size; i++) {
            this->data[i] = this->data[i] * value;
        }
        return *this;
    }

    Array<scalar>& operator * (const Array<scalar>& other) {
        if (this->NY == 1) {
            for (label i = 0; i < NX; i++) {
                this->data[i] = this->data[i] * other.data[i];
            }
        }
        else if (this->NY == 2) {
            for (label i = 0; i < NX; i++) {
                for (label j = 0; j < NY; j++) {
                    this->data[i * NY + j] = this->data[i * NY + j] * other.data[i];
                }
            }
        }
        return *this;
    }

    Array<scalar>& operator += (const Array<scalar>& other) {
        for (label i = 0; i < this->size; i++) {
            this->data[i] += other.data[i];
        }
        return *this;
    }

    Array<scalar>& operator -= (const Array<scalar>& other) {
        for (label i = 0; i < this->size; i++) {
            this->data[i] -= other.data[i];
        }
        return *this;
    }

    Array<scalar>& operator = (const Array<scalar>& other) {
        for (label i = 0; i < this->size; i++) {
            this->data[i] = other.data[i];
        }
        return *this;
    }

    void resize(label nx, label ny = 1, label nz = 1) {
        NX = nx;
        NY = ny;
        NZ = nz;
        size = NX * NY * NZ;
        data.resize(size);
    }

    void push_back(type value) {
        data.push_back(value);
        NX++;
        size++;
    }

    void set_all(type value) {
        for (int i = 0; i < data.size(); i++) {
            data[i] = value;
        }
    }

    void set_all(type value, label dim) {
        for (int i = 0; i < NX; i++) {
            data[i * NY + dim] = value;
        }
    }

    void set_array(const type* value) {
        for (int i = 0; i < data.size(); i++) {
            data[i] = value[i];
        }
    }

    void set_array(Array<type>& other) {
        for (label i = 0; i < other.size; i++) {
            this->data[i] = other.data[i];
        }
    }

    void times(scalar value) {
        for (size_t i = 0; i < data.size(); i++) {
            data[i] = data[i] * value;
        }
    }

    void print() {
        std::cout << "\n";
        if (NY == 1 && NZ == 1) {
            for (label i = 0; i < size; i++) {
                std::cout << data[i] << " ";
            }
            std::cout << "\n";
        }
        else if (NY != 1 && NZ == 1) {
            for (label i = 0; i < NX; i++) {
                for (label j = 0; j < NY; j++) {
                    std::cout << (*this)(i, j) << " ";
                }
                std::cout << "\n";
            }
        }
    }

    void read(const std::string& filename) {

        scalar value;
        std::ifstream file;
        file.open(filename);
        while (file >> value) {
            data.push_back(value);
        }
        file.close();
        size = data.size();
    }
};
