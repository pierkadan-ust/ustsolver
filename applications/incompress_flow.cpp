#include "ustsolver.hpp"

int main() {
	// Method Control
	path_main = "mesh/incompress_flow/";
	num_dimension = 3;
	if_reverse_axis = "all";

	max_iter_V = 200;
	tresMomt = 1e-13;
	tresP = 1e-13;

	sorMomt = 0.7;
	sorP = 0.3;
	
	// Domain and Material
	Zone air(1, "Air", path_main + "air.dat", "Tecplot");
	air.implement_bc(1, path_main + "inlet.dat");
	air.implement_bc(2, path_main + "outlet.dat");
	air.implement_bc(3, path_main + "inwall1.dat");
	air.implement_bc(4, path_main + "inwall2.dat");
	air.implement_bc(5, path_main + "inwall3.dat");
	air.implement_bc(6, path_main + "inwall4.dat");
	air.implement_bc(7, path_main + "outwall1.dat");
	air.implement_bc(8, path_main + "outwall2.dat");
	air.implement_bc(9, path_main + "outwall3.dat");
	air.implement_bc(10, path_main + "outwall4.dat");
	air.implement_bc(11, path_main + "outwall5.dat");

	Material mixture(&air);
	mixture.setMaterial("Density", 0.08189);
	mixture.setMaterial("DynamicViscosity", 8.411e-2);
	
	// Initialization
	Variable VT(&air, "vector", "VT");
	VT.BC_velocity_inlet(1, 0, 0.02189, 0);
	VT.BC_zero_gradient(2);
	VT.BC_wall_velocity(3, 0, 0, 0);
	VT.BC_wall_velocity(4, 0, 0, 0);
	VT.BC_wall_velocity(5, 0, 0, 0);
	VT.BC_wall_velocity(6, 0, 0, 0);
	VT.BC_wall_velocity(7, 0, 0, 0);
	VT.BC_wall_velocity(8, 0, 0, 0);
	VT.BC_wall_velocity(9, 0, 0, 0);
	VT.BC_wall_velocity(10, 0, 0, 0);
	VT.BC_wall_velocity(11, 0, 0, 0);
	VT.initial_value(0, 0, 0);

	Variable P(&air, "scalar", "P");
	P.BC_zero_gradient(1);
	P.BC_pressure_outlet(2, 0);
	P.BC_zero_gradient(3);
	P.BC_zero_gradient(4);
	P.BC_zero_gradient(5);
	P.BC_zero_gradient(6);
	P.BC_zero_gradient(7);
	P.BC_zero_gradient(8);
	P.BC_zero_gradient(9);
	P.BC_zero_gradient(10);
	P.BC_zero_gradient(11);
	P.initial_value(0);

	// Calculation and Solution
	Equation FLOW;
	FLOW.steadyFlow(&VT, &P, &mixture.Density, &mixture.DynamicViscosity, 1, sorMomt, sorP);

	Solution::write_tecplot(&VT, path_main + "out_vt.dat", 1, "magnitude");
	Solution::write_tecplot(&P, path_main + "out_p.dat");

	return 0;
}
