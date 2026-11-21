// Pre-defined Equations

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "field/Variable.hpp"
#include "solver/SparseMat.hpp"
#include "solver/Solve.hpp"

class Equation {

public:

	SparseMat A;
	Array<scalar> B;

	scalar resMomt = 0;
	scalar resP = 0;
	scalar resEnergy = 0;
	scalar resMass = 0;

public:

	Equation(label ifPrint = 1);

	void resize_matrix(Variable* variable);

	void get_diagonal(Variable* variable, string type = "original");

	// Mathematical Operations 
	void laplacian(Variable* phi, Variable* k);

	void div(Variable* phi, Variable* U, Variable* rho, Variable* cp = nullptr);

	void poission(Variable* diagA);

	// Physical Modules 
	void steadyFlow(Variable* VT, Variable* P, Variable* RHO, Variable* MU, label if_print = 1, scalar sor_momt = 0.7, scalar sor_p = 0.3, label iter_momt = 200, label iter_p = 200);

	void steadyEnergy(Variable* TEMP, Variable* VT, Variable* RHO, Variable* K, Variable* CP, label if_print = 0, scalar sor = 1.0);

	void steadyConduction(Variable* TEMP, Variable* K, label if_print = 0, scalar sor = 1.0);

	void transientFlow(Variable* VT, Variable* P, Variable* RHO, Variable* MU, scalar sor_momt = 0.7, scalar sor_p = 0.3, label iter_momt = 200, label iter_p = 200);

	void transientEnergy(Variable* TEMP, Variable* VT, Variable* RHO, Variable* K, Variable* CP, scalar sor = 1.0);

	// Time Discretization 
	void flowForwardEulerFirst(Array<scalar>* lastVT, Variable* RHO);

	void energyForwardEulerFirst(Variable* TEMP, Variable* RHO, Variable* CP);
};
