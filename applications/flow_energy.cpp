#include "ustsolver.hpp"

int main() {
	// Method Control
	path_main = "mesh/flow_energy/";
	num_dimension = 2;
	max_iter_V = 1500;
	max_iter_T = 100000;

	tresMomt = 1e-13;
	tresP = 1e-13;
	tresEnergy = 1e-9;

	sorMomt = 0.7;
	sorP = 0.3;

	if_reverse_axis = "all";
	if_offset_axis = 10;

	// Domain and Material
	Zone fluid(1, "Fluid", path_main + "fluid.dat", "Tecplot");
	fluid.implement_bc(1, path_main + "inlet.dat");
	fluid.implement_bc(2, path_main + "outlet.dat");
	fluid.implement_bc(3, path_main + "wall.dat");
	fluid.implement_bc(4, path_main + "inwall1.dat");
	fluid.implement_bc(5, path_main + "inwall2.dat");
	fluid.implement_bc(6, path_main + "inwall3.dat");
	fluid.implement_bc(7, path_main + "inwall4.dat");
	fluid.implement_bc(8, path_main + "outwall1.dat");
	fluid.implement_bc(9, path_main + "outwall2.dat");
	fluid.implement_bc(10, path_main + "outwall3.dat");
	fluid.implement_bc(11, path_main + "outwall4.dat");
	fluid.implement_bc(12, path_main + "outwall5.dat");

	Material mixture(&fluid);
	mixture.setMaterial("Density", 0.07189);
	mixture.setMaterial("DynamicViscosity", 8.411e-4);
	mixture.setMaterial("HeatConductivity", 0.0242);
	mixture.setMaterial("HeatCapacity", 1006.43);

	// Initialization
	Variable VT(&fluid, "vector", "VT");
	VT.BC_velocity_inlet(1, 0, 0.02189);
	VT.BC_zero_gradient(2);
	VT.BC_wall_velocity(3, 0, 0);
	VT.BC_wall_velocity(4, 0, 0);
	VT.BC_wall_velocity(5, 0, 0);
	VT.BC_wall_velocity(6, 0, 0);
	VT.BC_wall_velocity(7, 0, 0);
	VT.BC_wall_velocity(8, 0, 0);
	VT.BC_wall_velocity(9, 0, 0);
	VT.BC_wall_velocity(10, 0, 0);
	VT.BC_wall_velocity(11, 0, 0);
	VT.BC_wall_velocity(12, 0, 0);
	VT.initial_value(0, 0);

	Variable P(&fluid, "scalar", "P");
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
	P.BC_zero_gradient(12);
	P.initial_value(0);

	Variable TEMP(&fluid, "scalar", "TEMP");
	TEMP.BC_wall_temperature(1, 300);
	TEMP.BC_zero_gradient(2);
	TEMP.BC_wall_temperature(3, 300);
	TEMP.BC_wall_temperature(4, 1023);
	TEMP.BC_wall_temperature(5, 1023);
	TEMP.BC_wall_temperature(6, 1023);
	TEMP.BC_wall_temperature(7, 720);
	TEMP.BC_wall_temperature(8, 473);
	TEMP.BC_wall_temperature(9, 473);
	TEMP.BC_wall_temperature(10, 343);
	TEMP.BC_wall_temperature(11, 303);
	TEMP.BC_wall_temperature(12, 343);
	TEMP.initial_value(300);

	// Calculation and Solution
	Equation FLOW;
	FLOW.steadyFlow(&VT, &P, &mixture.Density, &mixture.DynamicViscosity, 1, sorMomt, sorP);

	Equation HEAT;
	HEAT.steadyEnergy(&TEMP, &VT, &mixture.Density, &mixture.HeatConductivity, &mixture.HeatCapacity);

	Solution::write_tecplot(&VT, path_main + "out_vt.dat", 1, "magnitude");
	Solution::write_tecplot(&P, path_main + "out_p.dat");
	Solution::write_tecplot(&TEMP, path_main + "out_temp.dat");

	return 0;
}
