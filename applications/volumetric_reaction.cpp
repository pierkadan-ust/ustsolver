#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

namespace Reaction {

    struct Condition {
        double R{ 8314.47 };  // universal gas constant [J/kmol-K]
        double e{ 2.71828 };  // natural constant [-]

        double p{};  // operating pressure [Pa]
        double V{};  // mesh volume [m3]

        double dt{};        // [s]
        double end_time{};  // [s]
        size_t num_iter{};  // [-]
    };

    struct VolumeReaction {
        double A{};  // pre-exponential factor
        double B{};  // temperature exponent
        double E{};  // activation energy [J/kmol]

        double Kf{};    // forward rate constant
        double R_t1{};  // kinetic reaction rate at time t+1 [kgmol/m3-s]
    };

    struct ReactionCoeff {
        double v1{};  // stoichiometric coefficient for reactant
        double v2{};  // stoichiometric coefficient for product
        double n1{};  // rate exponent for reactant
        double n2{};  // rate exponent for product

        double Cf_t1{};  // forward concentration multiplication without itself at time t+1
    };

    enum class SpeciesEquation {
        TaylorExpansion = 0,
        ExponentReduction = 1
    };

    struct Species {
        std::string formula;
        std::vector<ReactionCoeff> react_coeff;

        double Mw{};  // molecular weight [kg/kmol]
        double Cp{};  // specific heat [J/kg-K]

        double Y_t0{};  // mass fraction at time t
        double Y_t1{};  // mass fraction at time t+1

        double LHS{};  // left hand side of the species equation
        double RHS{};  // right hand side of the species equation

        double LHS_time{};
        double LHS_chem{};

        double RHS_time{};
        double RHS_chem{};
    };

    struct Energy {
        double T_t0{};  // temperature at time t
        double T_t1{};  // temperature at time t+1

        double LHS{};  // left hand side of the energy equation
        double RHS{};  // right hand side of the energy equation

        double LHS_time{};
        double LHS_chem{};

        double RHS_time{};
        double RHS_chem{};
    };

    struct Mixture {
        Condition condition;
        std::vector<Species> species;
        Energy energy;

        double rho_t0{};  // mixture density [kg/m3]
        double rho_t1{};

        double cp_t0{};  // mixture specific heat capacity [J/kg-K]
        double cp_t1{};

        void calc_mass_fraction_of_last_species();
        void calc_mixture_density_by_incomp_ideal_gas();

        bool check_mass_fraction();

        void initialization();
        void update_transient_values();
    };

    void Mixture::calc_mass_fraction_of_last_species() {
        species.back().Y_t1 = 1.0;
        for (size_t i = 0; i < species.size() - 1; ++i) {
            species.back().Y_t1 -= species[i].Y_t1;
        }
    }

    void Mixture::calc_mixture_density_by_incomp_ideal_gas() {
        rho_t1 = 0.0;
        for (size_t i = 0; i < species.size(); ++i) {
            rho_t1 += species[i].Y_t1 / species[i].Mw;
        }
        rho_t1 = condition.p / condition.R / energy.T_t1 / rho_t1;
    }

    void Mixture::update_transient_values() {
        rho_t0 = rho_t1;
        cp_t0 = cp_t1;
        energy.T_t0 = energy.T_t1;
        for (size_t i = 0; i < species.size(); ++i) {
            species[i].Y_t0 = species[i].Y_t1;
        }
    }

    bool Mixture::check_mass_fraction() {
        bool restart = false;
        for (size_t solute = 0; solute < species.size() - 1; ++solute) {
            if (species[solute].Y_t1 > 1.0 || species[solute].Y_t1 < 0.0) {
                std::cerr << "correct the mass fraction of " << species[solute].formula << " from " << species[solute].Y_t1
                    << " to " << species[solute].Y_t0 << std::endl;
                species[solute].Y_t1 = species[solute].Y_t0;
                restart = true;
            }
        }
        return restart;
    }

    void Mixture::initialization() {
        energy.T_t1 = energy.T_t0;
        for (size_t i = 0; i < species.size() - 1; ++i) {
            species[i].Y_t1 = species[i].Y_t0;
        }
        calc_mass_fraction_of_last_species();
        calc_mixture_density_by_incomp_ideal_gas();
        species.back().Y_t0 = species.back().Y_t1;
        rho_t0 = rho_t1;
    }

    struct Solver {
        Mixture mixture;
        std::vector<VolumeReaction> reactions;

        size_t num_species{};
        size_t num_reactions{};

        SpeciesEquation m_species_equation{};

        void calc_forward_rate_constant();
        void calc_kinetic_reaction_rate();
        void calc_forward_concentration_multiplication_without_itself(size_t solute_index);

        void assemble_transient_term_for_species_equation(size_t solute_index);
        void assemble_reaction_source_for_species_equation(size_t solute_index);
        void solve_species_equation(size_t solute_index);

        void assemble_transient_term_for_energy_equation();
        void assemble_reaction_source_for_energy_equation();
        void solve_energy_equation();

        void initialization();
        void print(double time, size_t iter);
    };

    void Solver::calc_forward_rate_constant() {
        for (size_t r = 0; r < num_reactions; ++r) {
            double exponent = -reactions[r].E / mixture.condition.R / mixture.energy.T_t1;
            reactions[r].Kf = reactions[r].A * std::pow(mixture.energy.T_t1, reactions[r].B) * std::pow(mixture.condition.e, exponent);
        }
    }

    void Solver::calc_kinetic_reaction_rate() {
        for (size_t r = 0; r < num_reactions; ++r) {
            reactions[r].R_t1 = reactions[r].Kf;
            for (size_t j = 0; j < num_species; ++j) {
                double Cj_t1 = mixture.rho_t1 * mixture.species[j].Y_t1 / mixture.species[j].Mw;
                reactions[r].R_t1 *= std::pow(Cj_t1, mixture.species[j].react_coeff[r].n1 + mixture.species[j].react_coeff[r].n2);
            }
        }
    }

    void Solver::calc_forward_concentration_multiplication_without_itself(size_t solute_index) {
        for (size_t r = 0; r < num_reactions; ++r) {
            auto& Cf_t1 = mixture.species[solute_index].react_coeff[r].Cf_t1;
            Cf_t1 = 1.0;
            for (size_t j = 0; j < num_species; ++j) {
                if (j != solute_index) {
                    double Cj_t1 = mixture.rho_t1 * mixture.species[j].Y_t1 / mixture.species[j].Mw;
                    Cf_t1 *= std::pow(Cj_t1, mixture.species[j].react_coeff[r].n1 + mixture.species[j].react_coeff[r].n2);
                }
            }
        }
    }

    void Solver::assemble_transient_term_for_species_equation(size_t solute_index) {
        auto& species = mixture.species[solute_index];
        species.RHS_time = mixture.rho_t0 * species.Y_t0 / mixture.condition.dt * mixture.condition.V;
        species.LHS_time = mixture.rho_t0 / mixture.condition.dt * mixture.condition.V;
    }

    void Solver::assemble_reaction_source_for_species_equation(size_t solute_index) {
        auto& species = mixture.species[solute_index];

        for (size_t r = 0; r < num_reactions; ++r) {
            double Rj = species.Mw * (species.react_coeff[r].v2 - species.react_coeff[r].v1)
                * reactions[r].Kf * species.react_coeff[r].Cf_t1;
            double Ci = mixture.rho_t1 * species.Y_t1 / species.Mw;
            double exponent = species.react_coeff[r].n1 + species.react_coeff[r].n2;

            double dRdC = exponent * Rj * std::pow(Ci, exponent - 1.0) * mixture.rho_t1 / species.Mw;

            if (m_species_equation == SpeciesEquation::TaylorExpansion) {
                if (Rj < 0) {
                    species.RHS_chem += (Rj * std::pow(Ci, exponent) - dRdC * species.Y_t1) * mixture.condition.V;
                    species.LHS_chem -= dRdC * mixture.condition.V;
                }
                else {
                    species.RHS_chem += Rj * std::pow(Ci, exponent) * mixture.condition.V;
                }
            }

            if (m_species_equation == SpeciesEquation::ExponentReduction) {
                if (Rj < 0) {
                    species.LHS_chem -= Rj * std::pow(Ci, exponent - 1.0) * mixture.rho_t0 / species.Mw * mixture.condition.V;
                }
                else {
                    species.RHS_chem += Rj * std::pow(Ci, exponent) * mixture.condition.V;
                }
            }
        }
    }

    void Solver::solve_species_equation(size_t solute_index) {
        auto& species = mixture.species[solute_index];

        species.RHS = species.RHS_time + species.RHS_chem;
        species.LHS = species.LHS_time + species.LHS_chem;
        species.Y_t1 = species.RHS / species.LHS;

        species.RHS_chem = 0.0;
        species.LHS_chem = 0.0;
    }

    void Solver::assemble_transient_term_for_energy_equation() {
        mixture.energy.RHS_time = mixture.rho_t0 * mixture.cp_t0 * mixture.energy.T_t0 / mixture.condition.dt * mixture.condition.V;
        mixture.energy.LHS_time = mixture.rho_t0 * mixture.cp_t0 / mixture.condition.dt * mixture.condition.V;
    }

    void Solver::assemble_reaction_source_for_energy_equation() {}

    void Solver::solve_energy_equation() {
        mixture.energy.RHS = mixture.energy.RHS_time + mixture.energy.RHS_chem;
        mixture.energy.LHS = mixture.energy.LHS_time + mixture.energy.LHS_chem;
        mixture.energy.T_t1 = mixture.energy.RHS / mixture.energy.LHS;

        mixture.energy.RHS_chem = 0.0;
        mixture.energy.LHS_chem = 0.0;
    }

    void Solver::initialization() {
        reactions.resize(num_reactions);
        mixture.species.resize(num_species);
        for (size_t i = 0; i < num_species; ++i) {
            mixture.species[i].react_coeff.resize(num_reactions);
        }
    }

    void Solver::print(double time, size_t iter) {
        std::cout << "time = " << time << std::endl;
        std::cout << "iteration = " << iter << std::endl;
        for (size_t i = 0; i < num_species; ++i) {
            std::cout << "mass fraction of " << mixture.species[i].formula << " = "
                << std::fixed << std::setprecision(8) << mixture.species[i].Y_t1 << std::endl;
        }
        for (size_t r = 0; r < num_reactions; ++r) {
            std::cout << "kinetic reaction rate " << r + 1 << " = "
                << std::fixed << std::setprecision(8) << reactions[r].R_t1 << std::endl;
        }
        std::cout << "mixture density = "
            << std::fixed << std::setprecision(9) << mixture.rho_t1 << std::endl;
        std::cout << std::endl;
    }
}

using namespace Reaction;

int main() {
    Solver solver;
    solver.m_species_equation = SpeciesEquation::TaylorExpansion;
    solver.num_species = 14;
    solver.num_reactions = 9;
    solver.initialization();

    auto& mixture = solver.mixture;
    auto& condition = solver.mixture.condition;
    auto& reactions = solver.reactions;

    condition.p = 10000;
    condition.V = 1e-6;

    condition.dt = 1e-4;
    condition.end_time = 1e-3;
    condition.num_iter = 100;

    mixture.energy.T_t0 = 873;

    auto& species_0 = mixture.species[0];
    auto& species_1 = mixture.species[1];
    auto& species_2 = mixture.species[2];
    auto& species_3 = mixture.species[3];
    auto& species_4 = mixture.species[4];
    auto& species_5 = mixture.species[5];
    auto& species_6 = mixture.species[6];
    auto& species_7 = mixture.species[7];
    auto& species_8 = mixture.species[8];
    auto& species_9 = mixture.species[9];
    auto& species_10 = mixture.species[10];
    auto& species_11 = mixture.species[11];
    auto& species_12 = mixture.species[12];
    auto& species_13 = mixture.species[13];

    species_0.formula = "ch3";
    species_0.Mw = 15;
    species_0.Cp = 3582;
    species_0.react_coeff[0].v2 = 1.0;
    species_0.react_coeff[0].n2 = 0.0;
    species_0.react_coeff[1].v2 = 1.0;
    species_0.react_coeff[1].n2 = 0.0;
    species_0.react_coeff[2].v2 = 1.0;
    species_0.react_coeff[2].n2 = 0.0;
    species_0.react_coeff[3].v2 = 1.0;
    species_0.react_coeff[3].n2 = 0.0;
    species_0.react_coeff[4].v1 = 1.0;
    species_0.react_coeff[4].n1 = 1.0;
    species_0.react_coeff[5].v1 = 1.0;
    species_0.react_coeff[5].n1 = 1.0;
    species_0.react_coeff[7].v1 = 1.0;
    species_0.react_coeff[7].n1 = 1.0;
    species_0.react_coeff[8].v1 = 1.0;
    species_0.react_coeff[8].n1 = 1.0;
    species_0.Y_t0 = 0.0;

    species_1.formula = "ch4";
    species_1.Mw = 16;
    species_1.Cp = 2222;
    species_1.react_coeff[4].v2 = 1.0;
    species_1.react_coeff[4].n2 = 0.0;
    species_1.react_coeff[5].v2 = 1.0;
    species_1.react_coeff[5].n2 = 0.0;
    species_1.react_coeff[7].v2 = 1.0;
    species_1.react_coeff[7].n2 = 0.0;
    species_1.react_coeff[8].v2 = 1.0;
    species_1.react_coeff[8].n2 = 0.0;
    species_1.Y_t0 = 0.0;

    species_2.formula = "h";
    species_2.Mw = 1;
    species_2.Cp = 20621;
    species_2.react_coeff[6].v1 = 2.0;
    species_2.react_coeff[6].n1 = 1.0;
    species_2.react_coeff[7].v1 = 1.0;
    species_2.react_coeff[7].n1 = 1.0;
    species_2.react_coeff[8].v2 = 1.0;
    species_2.react_coeff[8].n2 = 0.0;
    species_2.Y_t0 = 0.0;

    species_3.formula = "ph2";
    species_3.Mw = 33;
    species_3.Cp = 1338;
    species_3.react_coeff[5].v2 = 1.0;
    species_3.react_coeff[5].n2 = 0.0;
    species_3.Y_t0 = 0.0;

    species_4.formula = "ph3";
    species_4.Mw = 34;
    species_4.Cp = 1720;
    species_4.react_coeff[5].v1 = 1.0;
    species_4.react_coeff[5].n1 = 1.0;
    species_4.Y_t0 = 0.2;

    species_5.formula = "inme";
    species_5.Mw = 130;
    species_5.Cp = 943;
    species_5.react_coeff[3].v2 = 1.0;
    species_5.react_coeff[3].n2 = 0.0;
    species_5.Y_t0 = 0.0;

    species_6.formula = "inme2";
    species_6.Mw = 145;
    species_6.Cp = 1420;
    species_6.react_coeff[2].v2 = 1.0;
    species_6.react_coeff[2].n2 = 0.0;
    species_6.react_coeff[3].v1 = 1.0;
    species_6.react_coeff[3].n1 = 1.0;
    species_6.Y_t0 = 0.0;

    species_7.formula = "inme3";
    species_7.Mw = 160;
    species_7.Cp = 1767;
    species_7.react_coeff[2].v1 = 1.0;
    species_7.react_coeff[2].n1 = 1.0;
    species_7.Y_t0 = 0.2;

    species_8.formula = "ash2";
    species_8.Mw = 77;
    species_8.Cp = 619;
    species_8.react_coeff[4].v2 = 1.0;
    species_8.react_coeff[4].n2 = 0.0;
    species_8.Y_t0 = 0.0;

    species_9.formula = "ash3";
    species_9.Mw = 78;
    species_9.Cp = 788;
    species_9.react_coeff[4].v1 = 1.0;
    species_9.react_coeff[4].n1 = 1.0;
    species_9.Y_t0 = 0.2;

    species_10.formula = "game2";
    species_10.Mw = 100;
    species_10.Cp = 1420;
    species_10.react_coeff[0].v2 = 1.0;
    species_10.react_coeff[0].n2 = 0.0;
    species_10.react_coeff[1].v1 = 1.0;
    species_10.react_coeff[1].n1 = 1.0;
    species_10.Y_t0 = 0.0;

    species_11.formula = "game";
    species_11.Mw = 85;
    species_11.Cp = 943;
    species_11.react_coeff[1].v2 = 1.0;
    species_11.react_coeff[1].n2 = 0.0;
    species_11.Y_t0 = 0.0;

    species_12.formula = "game3";
    species_12.Mw = 115;
    species_12.Cp = 1767;
    species_12.react_coeff[0].v1 = 1.0;
    species_12.react_coeff[0].n1 = 1.0;
    species_12.Y_t0 = 0.2;

    species_13.formula = "h2";
    species_13.Mw = 2;
    species_13.Cp = 14283;
    species_13.react_coeff[6].v2 = 1.0;
    species_13.react_coeff[6].n2 = 0.0;
    species_13.react_coeff[8].v1 = 1.0;
    species_13.react_coeff[8].n1 = 1.0;

    reactions[0].A = 3.5e15;
    reactions[0].E = 2.47e8;
    reactions[0].B = 0.0;

    reactions[1].A = 8.7e7;
    reactions[1].E = 1.48e8;
    reactions[1].B = 0.0;

    reactions[2].A = 5.25e15;
    reactions[2].E = 1.979e8;
    reactions[2].B = 0.0;

    reactions[3].A = 3.38e8;
    reactions[3].E = 5.5e7;
    reactions[3].B = 0.0;

    reactions[4].A = 3.9e7;
    reactions[4].E = 6.9e6;
    reactions[4].B = 0.0;

    reactions[5].A = 3.9e7;
    reactions[5].E = 6.9e6;
    reactions[5].B = 0.0;

    reactions[6].A = 1e13;
    reactions[6].E = 0.0;
    reactions[6].B = 0.0;

    reactions[7].A = 2.4e19;
    reactions[7].E = 0.0;
    reactions[7].B = -1.0;

    reactions[8].A = 2.9e9;
    reactions[8].E = 3.645e7;
    reactions[8].B = 3.1;

    mixture.initialization();

    for (double current_time = condition.dt; current_time <= condition.end_time + 1.0e-16; current_time += condition.dt) {

        for (size_t iter = 0; iter < condition.num_iter; ++iter) {

            solver.calc_forward_rate_constant();

            bool restart = false;
            do {
                for (size_t solute = 0; solute < solver.num_species - 1; ++solute) {
                    solver.calc_forward_concentration_multiplication_without_itself(solute);
                    solver.assemble_reaction_source_for_species_equation(solute);
                    solver.assemble_transient_term_for_species_equation(solute);
                    solver.solve_species_equation(solute);
                }
                restart = mixture.check_mass_fraction();
            } while (restart);

            mixture.calc_mass_fraction_of_last_species();
            mixture.calc_mixture_density_by_incomp_ideal_gas();
            solver.calc_kinetic_reaction_rate();
        }

        mixture.update_transient_values();
        solver.print(current_time, condition.num_iter);
    }
}
