#include "material/Material.hpp"

Material::Material(Zone* domain) : zone(domain) {
	Density.setVariable(zone, "scalar", "Density");
	HeatCapacity.setVariable(zone, "scalar", "HeatCapacity");
	HeatConductivity.setVariable(zone, "scalar", "HeatConductivity");
	DynamicViscosity.setVariable(zone, "scalar", "DynamicViscosity");
	std::cout << "\nMaterial    : Zone = " << zone->name << std::endl;
}

void Material::setMaterial(string name, scalar value) {
	if (name == "Density") { Density.X.set_all(value);  Density.X_faces.set_all(value);  std::cout << "              set Density = " << value << "\n"; }
	else if (name == "HeatCapacity") { HeatCapacity.X.set_all(value);  HeatCapacity.X_faces.set_all(value);  std::cout << "              set Heat Capacity = " << value << "\n"; }
	else if (name == "HeatConductivity") { HeatConductivity.X.set_all(value);  HeatConductivity.X_faces.set_all(value);  std::cout << "              set Heat Conductivity = " << value << "\n"; }
	else if (name == "DynamicViscosity") { DynamicViscosity.X.set_all(value);  DynamicViscosity.X_faces.set_all(value);  std::cout << "              set Dynamic Viscosity = " << value << "\n"; }
}

void Material::setComponent(string name, scalar valueDensity, scalar valueHeatCapacity, scalar valueHeatConductivity, scalar valueDynamicViscosity) {
	Variable density(zone, "scalar", name, 0);
	density.X.set_all(valueDensity);
	density.X_faces.set_all(valueDensity);
	compDensity.push_back(density);

	Variable heatCapcity(zone, "scalar", name, 0);
	heatCapcity.X.set_all(valueHeatCapacity);
	heatCapcity.X_faces.set_all(valueHeatCapacity);
	compHeatCapacity.push_back(heatCapcity);

	Variable heatConductivity(zone, "scalar", name, 0);
	heatConductivity.X.set_all(valueHeatConductivity);
	heatConductivity.X_faces.set_all(valueHeatConductivity);
	compHeatConductivity.push_back(heatConductivity);

	Variable dynamicViscosity(zone, "scalar", name, 0);
	dynamicViscosity.X.set_all(valueDynamicViscosity);
	dynamicViscosity.X_faces.set_all(valueDynamicViscosity);
	compDynamicViscosity.push_back(dynamicViscosity);

	Variable massFraction(zone, "scalar", name, 0);
	compMassFraction.push_back(massFraction);

	compNumber++;
	compName.push_back(name);

	std::cout << "              set Component<" << name << ">"
		<< ", Density = " << valueDensity
		<< ", Heat Capacity = " << valueHeatCapacity
		<< ", Heat Conductivity = " << valueHeatConductivity
		<< ", Dynamic Viscosity = " << valueDynamicViscosity
		<< "\n";
}
