#include "ustsolver.hpp"

int main() {
	// Method Control
	path_main = "mesh/heat_conduction/";
	num_dimension = 2;
	max_iter_T = 100000;
	tresEnergy = 1e-6;
	if_tune_scale = 0.1;
	
	// Domain and Material
	Zone fluid(1, "fluid", path_main + "air.dat", "Tecplot");
	fluid.implement_bc(1, path_main + "up.dat");
	fluid.implement_bc(2, path_main + "down.dat");
	fluid.implement_bc(3, path_main + "left.dat");
	fluid.implement_bc(4, path_main + "right.dat");

	Material air(&fluid);
	air.setMaterial("HeatConductivity", 205);
	
	// Initialization
	Variable TEMP(&fluid, "scalar", "TEMP");
	TEMP.BC_wall_convection(1, 1000, 20);
	TEMP.BC_wall_convection(2, 1000, 20);
	TEMP.BC_wall_temperature(3, 100);
	TEMP.BC_wall_heatflux(4, 0);
	TEMP.initial_value(0);
	
	// Calculation and Solution
	Equation Equ;
	Equ.steadyConduction(&TEMP, &air.HeatConductivity);

	Solution::write_tecplot(&TEMP, path_main + "out_temp.dat");

	return 0;
}
