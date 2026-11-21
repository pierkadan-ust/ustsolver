// Material Properties

#pragma once

#include "base/Global.hpp"
#include "base/Array.hpp"
#include "field/Zone.hpp"
#include "field/Variable.hpp"

class Material {

public:

	Zone* zone{ nullptr };

	// cell property
	Variable Density;
	Variable HeatCapacity;
	Variable HeatConductivity;
	Variable DynamicViscosity;

	// component property
	label compNumber = 0;
	Array<string> compName;
	Array<Variable> compMassFraction;
	Array<Variable> compDensity;
	Array<Variable> compHeatCapacity;
	Array<Variable> compHeatConductivity;
	Array<Variable> compDynamicViscosity;

public:

	Material(Zone* domain);

	void setMaterial(string name, scalar value);

	void setComponent(string name, scalar valueDensity, scalar valueHeatCapacity, scalar valueHeatConductivity, scalar valueDynamicViscosity);
};
