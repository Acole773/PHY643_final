#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "reaction_solver.hpp"

const double R = 8.314; // J/mol·K
const double Tol_c = 1e-4; // Conservation tolerance

// ------------------ STRUCTURES ------------------
struct Species {
    std::string name;
    double abundance;   // in mol
    int C, H, O;
    double molar_mass;  // in g/mol
    double enthalpy;  // in J/mol
    Species(std::string n, double a, int c, int h, int o, double m, double Hf)
        : name(n), abundance(a), C(c), H(h), O(o), molar_mass(m), enthalpy(Hf) {}

};

struct Reaction {
    std::vector<int> reactants;
    std::vector<int> products;
    double A_forward, n_forward, Ea_forward; // SI units
    double A_backward = 0.0, n_backward = 0.0, Ea_backward = 0.0;
    double delta_H = 0.0; // J/mol
    double k_fwd, k_rev;
};


// ------------------ FILE READING ------------------
bool read_reactions(const std::string& filename,
                    std::vector<Reaction>& reactions,
                    const std::vector<Species>& species) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }

    std::string line;
while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    Reaction r;

    int num_reactants, num_products;
    bool valid = true;

    if (!(iss >> num_reactants)) continue;
    r.reactants.resize(num_reactants);
    for (int i = 0; i < num_reactants; ++i) {
        if (!(iss >> r.reactants[i])) { valid = false; break; }
        if (r.reactants[i] < 0 || r.reactants[i] >= species.size()) valid = false;
    }

    if (!(iss >> num_products)) continue;
    r.products.resize(num_products);
    for (int i = 0; i < num_products; ++i) {
        if (!(iss >> r.products[i])) { valid = false; break; }
        if (r.products[i] < 0 || r.products[i] >= species.size()) valid = false;
    }

    if (!(iss >> r.A_forward >> r.n_forward >> r.Ea_forward)) continue;

    if (valid) {
        reactions.push_back(r);
    } else {
        std::cerr << "Warning: Skipping invalid reaction due to bad indices.\n";
    }
}

    return true;
}


// ------------------ KINETIC CORE ------------------
void impose_detailed_balance(std::vector<Reaction>& reactions,
                             const std::vector<Species>& species,
                             double T) {
	std::cout << "stop here";
    for (auto& r : reactions) {
	double delta_H = 0.0;
	for (int idx : r.products)
	    if (idx >= 0 && idx < static_cast<int>(species.size())) {
		delta_H += species[idx].enthalpy;
	    } else {
		std::cerr << "Error: Product species index " << idx << " out of range.\n";
	    }
	for (int idx : r.reactants) {
	    if (idx >= 0 && idx < static_cast<int>(species.size())) {
		delta_H -= species[idx].enthalpy;
	    } else {
	        std::cerr << "Error: Reactant species index " << idx << " out of range.\n";
	    }
	}
        r.delta_H = delta_H;

        double k_f = r.A_forward * std::pow(T, r.n_forward) * std::exp(-r.Ea_forward / (R * T));
        double k_b = k_f * std::exp(r.delta_H / (R * T));

	r.k_fwd = k_f;
	r.k_rev = k_b;

    }
}

double compute_flux(const Reaction& r,
                    const std::vector<Species>& species,
                    double k_f, double k_b) {
    double forward = k_f;
    double backward = k_b;

    for (int idx : r.reactants)
        forward *= species[idx].abundance;
    for (int idx : r.products)
        backward *= species[idx].abundance;


    return forward - backward;
}
void update_reactions_FE(std::vector<Species>& species,
                         const std::vector<Reaction>& reactions,
                         double dt, double T,
                         double& net_energy);

void update_reactions_RK45(std::vector<Species>& species,
                           const std::vector<Reaction>& reactions,
                           double& dt, double T,
                           double& net_energy,
                           double abs_tol,
                           double rel_tol, 
			   double dt_min,
			   double dt_max);
//------------Update_reaction-------------------
IntegrationMethod parse_method(std::string s) {
    if (s == "FE") return FORWARD_EULER;
    if (s == "RKF45") return RKF45;
    // ...
    throw std::invalid_argument("Unknown method: " + s);
}

void update_reactions(std::vector<Species>& species,
                      const std::vector<Reaction>& reactions,
                      double& dt, double T,
                      double& net_energy,
                      IntegrationMethod method,
                      double abs_tol,
                      double rel_tol,
		      double dt_min,
		      double dt_max) {
    switch (method) {
        case FORWARD_EULER:
            update_reactions_FE(species, reactions, dt, T, net_energy);
            break;

        case RKF45:
            update_reactions_RK45(species, reactions, dt, T, net_energy, abs_tol, rel_tol, dt_min, dt_max);
            break;

        // Future methods:
        // case ADAPTIVE_IMEX:
        //     update_reactions_IMEX(species, reactions, dt, T, net_energy, ...);
        //     break;

        default:
            std::cerr << "Unknown integration method!" << std::endl;
            std::exit(1);
    }
}

//------------FE update-----------------
void update_reactions_FE(std::vector<Species>& species,
                      const std::vector<Reaction>& reactions,
                      double dt, double T,
                      double& net_energy) {
    net_energy = 0.0;

    for (const Reaction& r : reactions) {
        double k_f = r.A_forward * std::pow(T, r.n_forward) * std::exp(-r.Ea_forward / (R * T));
        double k_b = r.A_backward * std::pow(T, r.n_backward) * std::exp(-r.Ea_backward / (R * T));

        double delta = compute_flux(r, species, k_f, k_b) * dt;

        for (int idx : r.reactants){
	       	species[idx].abundance -= delta;
		if (species[idx].abundance < 0) species[idx].abundance = 1e-20;
	}
	for (int idx : r.products){
	      	species[idx].abundance += delta;
		if (species[idx].abundance < 0) species[idx].abundance = 1e-20;
	}
        net_energy += delta * r.delta_H;
    }
}
//---------------RK45 update ---------------------

void update_reactions_RK45(std::vector<Species>& species,
                           const std::vector<Reaction>& reactions,
                           double& dt, double T,
                           double& net_energy,
                           double abs_tol = 1e-9,
                           double rel_tol = 1e-6,
                           double dt_min = 1e-12,
                           double dt_max = 1e2) {
    // Backup initial state
    std::vector<double> y0(species.size());
    for (size_t i = 0; i < species.size(); ++i) {
        y0[i] = species[i].abundance;
    }

    net_energy = 0.0;
    bool step_accepted = false;

    while (!step_accepted) {
        // Coefficients for RKF45 (Butcher tableau)
        const double a2 = 1.0/4.0;
        const double a3 = 3.0/8.0;
        const double a4 = 12.0/13.0;
        const double a5 = 1.0;
        const double a6 = 1.0/2.0;

        const double b21 = 1.0/4.0;

        const double b31 = 3.0/32.0;
        const double b32 = 9.0/32.0;

        const double b41 = 1932.0/2197.0;
        const double b42 = -7200.0/2197.0;
        const double b43 = 7296.0/2197.0;

        const double b51 = 439.0/216.0;
        const double b52 = -8.0;
        const double b53 = 3680.0/513.0;
        const double b54 = -845.0/4104.0;

        const double b61 = -8.0/27.0;
        const double b62 = 2.0;
        const double b63 = -3544.0/2565.0;
        const double b64 = 1859.0/4104.0;
        const double b65 = -11.0/40.0;

        const double c1 = 16.0/135.0;
        const double c3 = 6656.0/12825.0;
        const double c4 = 28561.0/56430.0;
        const double c5 = -9.0/50.0;
        const double c6 = 2.0/55.0;

        const double dc1 = c1 - 25.0/216.0;
        const double dc3 = c3 - 1408.0/2565.0;
        const double dc4 = c4 - 2197.0/4104.0;
        const double dc5 = c5 + 1.0/5.0;
        const double dc6 = c6;

        std::vector<double> k1(species.size(), 0.0);
        std::vector<double> k2(species.size(), 0.0);
        std::vector<double> k3(species.size(), 0.0);
        std::vector<double> k4(species.size(), 0.0);
        std::vector<double> k5(species.size(), 0.0);
        std::vector<double> k6(species.size(), 0.0);
        std::vector<double> y_temp(species.size());

        // Helper lambda to compute dY/dt
        auto compute_dydt = [&](const std::vector<double>& y) {
            std::vector<double> dydt(species.size(), 0.0);
            for (const auto& r : reactions) {
                double k_f = r.A_forward * std::pow(T, r.n_forward) * std::exp(-r.Ea_forward / (R * T));
                double k_b = r.A_backward * std::pow(T, r.n_backward) * std::exp(-r.Ea_backward / (R * T));
		std::vector<Species> y_temp = species;
		for (size_t i = 0; i < y_temp.size(); ++i) {
		    y_temp[i].abundance = y[i];
		}
		double flux = compute_flux(r, y_temp, k_f, k_b); // mol/m³/s

                for (int idx : r.reactants) dydt[idx] -= flux;
                for (int idx : r.products)  dydt[idx] += flux;
            }
            return dydt;
        };

        // k1 = f(t, y)
        k1 = compute_dydt(y0);

        // k2 = f(t + a2*dt, y + dt*b21*k1)
        for (size_t i = 0; i < species.size(); ++i)
            y_temp[i] = y0[i] + dt * b21 * k1[i];
        k2 = compute_dydt(y_temp);

        for (size_t i = 0; i < species.size(); ++i)
            y_temp[i] = y0[i] + dt * (b31*k1[i] + b32*k2[i]);
        k3 = compute_dydt(y_temp);

        for (size_t i = 0; i < species.size(); ++i)
            y_temp[i] = y0[i] + dt * (b41*k1[i] + b42*k2[i] + b43*k3[i]);
        k4 = compute_dydt(y_temp);

        for (size_t i = 0; i < species.size(); ++i)
            y_temp[i] = y0[i] + dt * (b51*k1[i] + b52*k2[i] + b53*k3[i] + b54*k4[i]);
        k5 = compute_dydt(y_temp);

        for (size_t i = 0; i < species.size(); ++i)
            y_temp[i] = y0[i] + dt * (b61*k1[i] + b62*k2[i] + b63*k3[i] + b64*k4[i] + b65*k5[i]);
        k6 = compute_dydt(y_temp);

        // Compute 5th-order solution and error estimate
        std::vector<double> y5(species.size());
        std::vector<double> error(species.size());
        double max_err_ratio = 0.0;

        for (size_t i = 0; i < species.size(); ++i) {
            y5[i] = y0[i] + dt * (c1*k1[i] + c3*k3[i] + c4*k4[i] + c5*k5[i] + c6*k6[i]);
            error[i] = dt * (dc1*k1[i] + dc3*k3[i] + dc4*k4[i] + dc5*k5[i] + dc6*k6[i]);

            double tol = abs_tol + rel_tol * std::max(std::abs(y0[i]), std::abs(y5[i]));
            max_err_ratio = std::max(max_err_ratio, std::abs(error[i]) / tol);
        }

        // Accept or reject step
        if (max_err_ratio <= 1.0) {
            step_accepted = true;
            net_energy = 0.0;

            for (size_t i = 0; i < species.size(); ++i) {
                species[i].abundance = std::max(1e-20, y5[i]);  // prevent negative abundances
            }

            // Estimate energy release
            for (const auto& r : reactions) {
                double k_f = r.A_forward * std::pow(T, r.n_forward) * std::exp(-r.Ea_forward / (R * T));
                double k_b = r.A_backward * std::pow(T, r.n_backward) * std::exp(-r.Ea_backward / (R * T));
                double delta = compute_flux(r, species, k_f, k_b) * dt;
                net_energy += delta * r.delta_H;
            }

            // Adapt next timestep
            double factor = std::min(2.0, std::max(0.1, 0.9 * std::pow(1.0 / max_err_ratio, 0.2)));
            double clamp = 0.0;
	    if (dt * factor < dt_min)    
	    	clamp = dt_min;
	    else if (dt * factor > dt_max)
		clamp = dt_max;
	    else clamp = dt*factor;
	    dt = clamp;
        } else {
            // Reject step, reduce dt and retry
            double factor = std::max(0.1, 0.9 * std::pow(1.0 / max_err_ratio, 0.2));
            dt = std::max(dt * factor, dt_min);
            // You could add a break or warning here if dt is too small
        }
    }
}



// ------------------ CONSERVATION ------------------
void check_conservation(const std::vector<Species>& species) {
    double total_C = 0, total_H = 0, total_O = 0;

    for (const auto& sp : species) {
        total_C += sp.C * sp.abundance;
        total_H += sp.H * sp.abundance;
        total_O += sp.O * sp.abundance;
    }

    std::cout << "Atoms -> C: " << total_C << " H: " << total_H << " O: " << total_O << "\n";
}

// ------------------ MAIN DRIVER ------------------
int main() {
    std::vector<Species> species = {
        {"CH4",   1000000000.0, 1, 4, 0, 16.04,  -74870},	//0
        {"O2",    2000000000.0, 0, 0, 2, 32.00,   1e-20},	//1
        {"CO2",	 1e-20, 1, 0, 2, 44.01, -393520},	//2
        {"H2O",  1e-20, 0, 2, 1, 18.02, -241820},	//3
        {"CO",   1e-20, 1, 0, 1, 28.01, -110530},	//4
        {"H2",   1e-20, 0, 2, 0,  2.02,   1e-20},	//5
        {"OH",   1e-20, 0, 1, 1, 17.01,  +39450},	//6
        {"H",    1e-20, 0, 1, 0, 1.008, +218000},	//7
        {"O",    1e-20, 0, 0, 1, 16.00, +249170},	//8
        {"HO2",  1e-20, 0, 1, 2, 33.00, +120300},  	//9
	{"CH3",  1e-20, 1, 3, 0, 15.03, +145690},  	//10
	{"CH3OO", 1e-20, 1, 3, 2, 47.05,  +12000},  	//11 //ESTIMATED
	{"CH3OOH",1e-20, 1, 4, 2, 48.06, -131000},  	//12
	{"CH3OH", 1e-20, 1, 4, 1, 32.04, -205000},  	//13
	{"CH2",  1e-20, 1, 2, 0, 14.03, +386390},  	//14
	{"CH3O", 1e-20, 1, 3, 1, 31.03,  +17200},  	//15
    	{"CH2O", 1e-20, 1, 2, 1, 30.03, -115900},	//16
    	{"CH",   1e-20, 1, 1, 0, 13.02, +594130},	//17
    	{"CHO",  1e-20, 1, 1, 1, 29.02,  +43510},	//18
    	{"CH2O2",1e-20, 1, 2, 2, 46.03, -378600}	//19
};

    std::vector<Reaction> reactions;
    if (!read_reactions("rates32.txt", reactions, species)) return 1;

    double T = 1000.0;

    IntegrationMethod method = RKF45; // or FORWARD_EULER

    double t = 0.0;
    double t_end = 1.0e-3;
    double dt = 1.0e-12; // Initial timestep
    double net_energy = 0.0;
    double dt_min = 1e-16, dt_max = 1.0;
    double total_energy = 0.0;
    int step = 0;

    double abs_tol = 1e-6;
    double rel_tol = 1e-4;

    //bool restep = 0;
    //int sub_step = 0; 
    //int max_restep = 20;  
    //std::vector<Species> trial_species;

    //-------------
std::cout << "Initial CH4: " << species[0].abundance << ", O2: " << species[1].abundance << "before impose detailed balance" <<std::endl;
    //-------------

    impose_detailed_balance(reactions, species, T);

    std::cout << "Initial CH4: " << species[0].abundance << ", O2: " << species[1].abundance << "after impose detailed balance" <<std::endl;

    while (t < t_end) {
	double dt_actual = std::min(dt, t_end - t); // avoid overshooting
    	update_reactions(species, reactions, dt_actual, T, net_energy, method, abs_tol, rel_tol, dt_min, dt_max);
    	t += dt_actual;

    	double net_energy = 0.0;

        total_energy += net_energy;

        double total_C = 0, total_H = 0, total_O = 0;
        for (const auto& sp : species) {
            total_C += sp.C * sp.abundance;
            total_H += sp.H * sp.abundance;
            total_O += sp.O * sp.abundance;
        }

        double max_dev = std::max({std::abs(total_C - 10), std::abs(total_H - 40), std::abs(total_O - 40)});
        if (max_dev > Tol_c) {
            dt = std::max(dt * 0.5, dt_min);
        } else {
            dt = std::min(dt * 1.1, dt_max);
        }

        std::cout << "Step " << step << ": ";
        for (const auto& sp : species)
            std::cout << sp.name << "=" << sp.abundance << " ";
        std::cout << "| ΔE = " << net_energy << " J | ";
        check_conservation(species);

	step++;
    }

    std::cout << "Total energy released: " << total_energy << " J\n";
    return 0;
}


