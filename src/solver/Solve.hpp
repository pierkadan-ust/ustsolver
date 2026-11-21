// Solver for Linear Equations

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "field/Variable.hpp"
#include "solver/SparseMat.hpp"

namespace Solve {

	label GaussSeidel(scalar& tres, scalar& residual, SparseMat& A, Array<scalar>& B, Variable* variable, scalar SOR, label num_iter, label if_print_res, string name = "Iteration", label if_reference = 0);
}
