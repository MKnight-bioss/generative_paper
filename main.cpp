#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iomanip>

using namespace std;

struct farm {
    int disease_state = 0; //disease status of farm: 1 - Infected, 0 - Suceptible
    int first_inf = 0;  //index identifying first farm infected. If equal to 1, farm is the initial case
    int recovered = 0;  //index indicating whether farm has recovered from initial infection, assuming farm is
                        //initial case
    int infections = 0; //number of infections during first infection for initial case, i.e. while recovered = 0.

    vector<int> traders;    //list of trading partners
    vector<int> not_traders;    //list of farms that aren't trading partners

    long double tot_add_rate = 0.0; //total partnership formation rate
    long double tot_delete_rate = 0.0;  //total partnership cessation rate
    long double tot_trade_rate = 0.0;   //total trade rate
};

vector<farm> farms; //each element of structure contained in vector

int main() {
    long double sim_start_time = time(NULL);

    //output file initialisation
    string path("FILE PATH HERE");
    string name("FILE NAME HERE");
    string type("FILE TYPE HERE");

    ofstream out_file;
    out_file.open(path + name + type);
    out_file << "inf_prob,inf_period,theta,ptnr_duration,r0_theory,r0_sim,equi_inf" << endl;

    //GLOBAL PARAMS
    const double N = 200.0; //number of farms
    const double t_max = 500.0; //time length
    const double equi_begin = 400.0; //time to introduce disease
    const double time_equi_inf = t_max - 50.0;  //time period in which to calculate equilibrium prevalence

    double inf_prob = 0.4;  //per-animal infection probability (lambda in paper)
    double gamma = 0.2; //recovery rate

    //INITIALISE RNG
    unsigned long int seed = time (NULL);
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    //INITIALISE STRUCTURE
    for (int i = 0; i < N; ++i) {
        farm f;
        farms.push_back(f);
    }

    //DISTRIBUTION QUANTITIES
    double eta = 2.0;   //average in-flow
    double zeta = 2.0;  //average out-flow
    double k = 2.0;     //expected number of trading partners
    double trades = 2.0;    //expected unit-time number of trades
    double d = eta*zeta;    //constant in partnership cessation rate (initial value decided so average
                            //partnership duration is one unit of time)
    double delta = d / (eta * zeta);    //partnership cessation rate
    double current_time = 1 / delta;    //average partnership duration
    double a = delta * k * (N - 1) / (eta * zeta * (N - 1 - k));    //constant in partnership formation rate. Value
                                                                    //set so desired number of trading partners is
                                                                    //obtained
    double alpha = a * eta * zeta / (N - 1);    //partnership formation rate
    double theta = eta / trades;    //batch size. value chosen so initial theta = 1
    double b = eta / (k * theta * min(eta, zeta));  //constant in trade rate.

    //BEGINNING SIM
    int num_params = 20;    //number of parameter values to run sims for
    int num_reps = (int) N; //number of runs per parameter value. Set to N as every farm is chosen as the initial
                            //infected
    int num_reps_per_farm = 10; //number of reps per farm initially infected

    double avg_r0_theory;   //value of R0 for param value as predicted by theoretical R0 (see paper)
    double avg_r0_sim;      //R0 from simulation (number of infections initial farm causes before recovery)
    double avg_equi_inf;    //average system prevalence when infection has equilibriated

    double epsilon_trade;   //scaling constant for batch size and num trades
    double epsilon_ptnr;    //scaling constant for partnership formation and cessation rates

    for (int np = 0; np < num_params; ++np) {
        long double param_start_time = time(NULL);

        auto increment_size = 1.0;
        auto vector_size = int(t_max / increment_size + 1);

        //unit-time average vectors, averaged over num_reps and num_reps_per_farm
        vector<double> avg_avg_in_vol_unit_time(vector_size, 0);
        vector<double> avg_avg_traders_unit_time(vector_size, 0);
        vector<double> avg_avg_trades_unit_time(vector_size, 0);
        vector<double> avg_avg_infected_unit_time(vector_size, 0);
        vector<double> avg_avg_infections_unit_time(vector_size, 0);

        vector<double> r0_values(num_reps, 0);
        vector<double> equi_inf_values(num_reps, 0);

        avg_r0_theory = 0.0;
        avg_r0_sim = 0.0;
        avg_equi_inf = 0.0;

        double new_time = 1.0 + (0.2 * np);
        epsilon_ptnr = current_time / new_time;
        double new_theta = 1.0 + np;
        epsilon_trade = theta / new_theta;

        b *= epsilon_trade;
        theta /= epsilon_trade;
        alpha *= epsilon_ptnr;
        delta *= epsilon_ptnr;
        double phi = b * min(eta, zeta);    //trade rate
        double B = 1 - pow(1 - inf_prob, theta);  //total infection probability
        double beta = phi * B;  //rate of infection
        double T = beta / (beta + delta + gamma); //disease transmissibility

        for (int nr = 0; nr < num_reps; ++nr) {
            cout << np << '\t' << nr << endl;
            for (int nrpf = 0; nrpf < num_reps_per_farm; ++nrpf) {

                //vectors to calculate unit-time averages
                vector<double> traders_unit_time(vector_size, 0);
                vector<double> traders_change_unit_time(vector_size, 0);
                vector<double> avg_traders_unit_time(vector_size, 0);

                vector<double> infected_unit_time(vector_size, 0);
                vector<double> infected_change_unit_time(vector_size, 0);
                vector<double> avg_infected_unit_time(vector_size, 0);

                vector<double> trades_unit_time(vector_size, 0);
                vector<double> infections_unit_time(vector_size, 0);

                vector<double> in_vol_unit_time(vector_size, 0);
                vector<double> avg_in_vol_unit_time(vector_size, 0);

                //SET TRADERS
                double total_traders = 0.0; //total number of partnerships in the system
                long double total_add_rate = 0.0;   //system-wide total partnership formation rate
                long double total_delete_rate = 0.0;    //system-wide total partnership cessation rate

                //system begins in a disconnected state
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        farms[i].not_traders.push_back(j);
                        if (i != j) farms[i].tot_add_rate += alpha;
                    }
                }

                total_add_rate = N * (N - 1) * alpha;

                //SET DISEASE STATES
                double S = N;   //number of susceptible farms
                double I = N - S;   //number of infected farms
                int I_init = 1; //number of farms seeded when disease is introduced
                int first_I;

                //SET RATES
                long double total_rate = 0.0;   //sum of event rates
                long double total_trade_rate = 0.0; //system-wide total trade rate
                long double total_rec_rate = 0.0;   //system-wide total recovery rate

                total_rate = total_add_rate;

                //CALC R0
                //since system is homogeneous, only need to calc theory R0 once
                if (nr == 0 && nrpf == 0) {
                    avg_r0_theory =
                            (alpha * (N - 1) + gamma * k) * beta * T /
                            (beta * (alpha + gamma) - alpha * delta * T);
                }


                //SIMULATION
                double dt = 0.0;    //current time in simulation
                long double event_rate; //next event rate
                long double partial_rate;
                int break_check;

                int disease_intro_check = 0;    //index to identify whether it is time to introduce disease

                int event_check;

                while (dt < t_max) {
                    event_rate = 0.0;
                    partial_rate = 0.0;
                    break_check = 0;
                    event_check = 0;

                    //check to see if it's time to introduce disease
                    if (disease_intro_check != 1) {
                        if (dt >= equi_begin) { //make sure only to introduce disease if it's time
                            disease_intro_check = 1;    //set to 1 so disease isn't introduced multiple times
                            while (I < I_init) {    //only introduce the desired number of initial infecteds
                                auto rand_I = nr;   //each farm is eventually chosen to be the initial infected
                                if (farms[rand_I].disease_state == 0) { //this is kept in case the initial infecteds
                                                                        //is larger than 1
                                    farms[rand_I].disease_state = 1;
                                    farms[rand_I].first_inf = 1;
                                    ++I;    //update number of infecteds
                                    --S;    //update number of susceptibles
                                    total_rec_rate = gamma; //update total recovery rate
                                    total_rate += gamma;    //update total rate
                                    first_I = rand_I;

                                    //update unit time vectors
                                    auto time_index = int(floor(dt / increment_size)) + 1;
                                    if (time_index <= t_max / increment_size) {
                                        infected_unit_time[time_index] += I / N;
                                        ++infected_change_unit_time[time_index];
                                    }
                                }
                            }
                        }
                    }

                    //update time
                    dt += -log(gsl_rng_uniform_pos(r)) / total_rate;

                    //next event rate
                    event_rate = gsl_rng_uniform(r) * total_rate;

                    if (event_rate < total_add_rate && event_check != 1) {
                        //trader addition
                        event_check = 1;

                        for (int i = 0; i < N && break_check != 1; ++i) {
                            partial_rate += farms[i].tot_add_rate;  //decide what farm makes an addition
                            if (partial_rate >= event_rate) {
                                partial_rate -= farms[i].tot_add_rate;
                                for (int j = 0; j < farms[i].not_traders.size() && break_check != 1; ++j) {
                                    if (i != farms[i].not_traders[j]) {
                                        auto trader_index = farms[i].not_traders[j];
                                        auto add_rate = alpha;
                                        partial_rate += add_rate;   //decide what farm is added
                                        if (partial_rate >= event_rate) {
                                            break_check = 1;
                                            farms[i].traders.push_back(trader_index);   //add new partner to list
                                            //remove new partner from list of farms that aren't trade partners
                                            farms[i].not_traders.erase(farms[i].not_traders.begin() + j);
                                            ++total_traders;    //update total number of trade partners

                                            //update rates that have changed
                                            auto delete_rate = delta;
                                            auto trade_rate = phi;
                                            farms[i].tot_add_rate -= add_rate;
                                            farms[i].tot_delete_rate += delete_rate;
                                            farms[i].tot_trade_rate += trade_rate;

                                            total_add_rate -= add_rate;
                                            total_delete_rate += delete_rate;
                                            total_trade_rate += trade_rate;
                                            total_rate += trade_rate + delete_rate - add_rate;

                                            //update unit time vectors
                                            auto time_index = int(floor(dt / increment_size)) + 1;
                                            if (time_index <= t_max / increment_size) {
                                                traders_unit_time[time_index] += total_traders / N;
                                                ++traders_change_unit_time[time_index];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else if (event_rate >= total_add_rate && event_rate < total_add_rate + total_delete_rate &&
                               event_check != 1) {
                        //trader deletion
                        event_check = 1;
                        partial_rate = total_add_rate;

                        for (int i = 0; i < N && break_check != 1; ++i) {
                            partial_rate += farms[i].tot_delete_rate;   //decide what farm removes a trade partner
                            if (partial_rate >= event_rate) {
                                partial_rate -= farms[i].tot_delete_rate;
                                for (int j = 0; j < farms[i].traders.size() && break_check != 1; ++j) {
                                    auto trader_index = farms[i].traders[j];
                                    auto delete_rate = delta;
                                    partial_rate += delete_rate;    //decide what trade partner is removed
                                    if (partial_rate >= event_rate) {
                                        break_check = 1;
                                        farms[i].not_traders.push_back(trader_index);   //add partner to list of farms
                                                                                        //that aren't trade partners
                                        farms[i].traders.erase(farms[i].traders.begin() + j); //remove partner from list
                                        --total_traders;    //update total number of trade partners

                                        //update rates that have changed
                                        auto add_rate = alpha;
                                        auto trade_rate = phi;
                                        farms[i].tot_add_rate += add_rate;
                                        farms[i].tot_delete_rate -= delete_rate;
                                        farms[i].tot_trade_rate -= trade_rate;

                                        total_add_rate += add_rate;
                                        total_delete_rate -= delete_rate;
                                        total_trade_rate -= trade_rate;
                                        total_rate += add_rate - delete_rate - trade_rate;

                                        //update unit time vectors
                                        auto time_index = int(floor(dt / increment_size)) + 1;
                                        if (time_index <= t_max / increment_size) {
                                            traders_unit_time[time_index] += total_traders / N;
                                            ++traders_change_unit_time[time_index];
                                        }
                                    }
                                }
                            }
                        }
                    } else if (event_rate >= total_add_rate + total_delete_rate &&
                               event_rate < total_add_rate + total_delete_rate + total_trade_rate && event_check != 1) {
                        //trade
                        event_check = 1;
                        partial_rate = total_add_rate + total_delete_rate;

                        for (int i = 0; i < N && break_check != 1; ++i) {
                            partial_rate += farms[i].tot_trade_rate;    //decide what farm makes a trade
                            if (partial_rate >= event_rate) {
                                partial_rate -= farms[i].tot_trade_rate;
                                for (int j = 0; j < farms[i].traders.size() && break_check != 1; ++j) {
                                    auto trader_index = farms[i].traders[j];
                                    auto trade_rate = phi;
                                    partial_rate += trade_rate; //decide what trade partner is traded with
                                    if (partial_rate >= event_rate) {
                                        break_check = 1;

                                        //check for infection
                                        if (farms[trader_index].disease_state == 1) {   //make sure partner is infected
                                            if (farms[i].disease_state == 0) {  //make sure buying farm is susceptible
                                                if (gsl_rng_uniform(r) < B) {   //check whether infection occurs
                                                    farms[i].disease_state = 1; //update buying farm's disease status
                                                    ++I;    //update number of infecteds
                                                    --S;    //update number of susceptibles
                                                    total_rec_rate += gamma;    //update total recovery rate
                                                    total_rate += gamma;    //update total rate

                                                    //if selling farm is initial infected and hasn't recovered from
                                                    //initial infection, then update number of infections
                                                    if (farms[trader_index].first_inf == 1) {
                                                        if (farms[trader_index].recovered == 0) {
                                                            ++farms[trader_index].infections;
                                                        }
                                                    }

                                                    //update unit time vectors
                                                    auto time_index = int(floor(dt / increment_size)) + 1;
                                                    if (time_index <= t_max / increment_size) {
                                                        infected_unit_time[time_index] += I / N;
                                                        ++infected_change_unit_time[time_index];

                                                        ++infections_unit_time[time_index];
                                                    }
                                                }
                                            }
                                        }
                                        //update unit time vectors
                                        auto time_index = int(floor(dt / increment_size)) + 1;
                                        if (time_index <= t_max / increment_size) {
                                            in_vol_unit_time[time_index] += theta;
                                            ++trades_unit_time[time_index];
                                        }
                                    }
                                }
                            }
                        }
                    } else if (event_rate >= total_add_rate + total_delete_rate + total_trade_rate &&
                               event_rate < total_rate && event_check != 1) {
                        //recovery
                        event_check = 1;
                        partial_rate = total_add_rate + total_delete_rate + total_trade_rate;
                        for (int i = 0; i < N && break_check != 1; ++i) {
                            if (farms[i].disease_state == 1) partial_rate += gamma; //decide what infected farm recovers
                            if (partial_rate >= event_rate) {
                                break_check = 1;
                                farms[i].disease_state = 0; //update farm's disease status
                                ++S;    //update number of susceptibles
                                --I;    //update number of infecteds
                                total_rec_rate -= gamma;    //update total recovery rate
                                total_rate -= gamma;    //update total rate

                                //if farm is initial infected and hasn't recovered from initial infection, then
                                //R0 is the number of infections this farm caused
                                if (farms[i].first_inf == 1) {
                                    if (farms[i].recovered == 0) {
                                        farms[i].recovered = 1;
                                        avg_r0_sim += (double) farms[i].infections / (num_reps*num_reps_per_farm);
                                        r0_values[nr] = (double) farms[i].infections / num_reps_per_farm;
                                    }
                                }

                                //update unit time vectors
                                auto time_index = int(floor(dt / increment_size)) + 1;
                                if (time_index <= t_max / increment_size) {
                                    infected_unit_time[time_index] += I / N;
                                    ++infected_change_unit_time[time_index];
                                }
                            }
                        }
                    } else {    //if rates don't add up, simulation is halted
                        cout << "no event; rates don't match" << endl;
                        cout << "add_rate " << total_add_rate << endl;
                        cout << "delete_rate " << total_delete_rate << endl;
                        cout << "trade_rate " << total_trade_rate << endl;
                        cout << "tot_rate " << total_rate << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                //calculate unit-time averages and place into respective vectors
                transform(traders_unit_time.begin(), traders_unit_time.end(), traders_change_unit_time.begin(),
                          avg_traders_unit_time.begin(), divides<double>());
                transform(infected_unit_time.begin(), infected_unit_time.end(), infected_change_unit_time.begin(),
                          avg_infected_unit_time.begin(), divides<double>());
                transform(in_vol_unit_time.begin(), in_vol_unit_time.end(), trades_unit_time.begin(),
                          avg_in_vol_unit_time.begin(), divides<double>());

                //set unit-time average vector values at t = 0
                avg_traders_unit_time[0] = 0.0;
                avg_infected_unit_time[0] = 0.0;
                avg_in_vol_unit_time[0] = 0.0;

                for (int i = 0; i < vector_size; ++i) {
                    //check if there were any divisions by 0. If so, set x[t] = x[t-1].
                    if (i > 0) {
                        if (traders_change_unit_time[i] == 0) avg_traders_unit_time[i] = avg_traders_unit_time[i - 1];
                        if (infected_change_unit_time[i] == 0)
                            avg_infected_unit_time[i] = avg_infected_unit_time[i - 1];
                        if (trades_unit_time[i] == 0) avg_in_vol_unit_time[i] = avg_in_vol_unit_time[i - 1];
                    }
                    //calculate average disease prevalence for initial infected farm
                    if (i * increment_size >= time_equi_inf) {
                        equi_inf_values[nr] += avg_infected_unit_time[i] * increment_size / ((t_max - time_equi_inf) * num_reps_per_farm);
                    }
                    //update average average vectors
                    avg_avg_traders_unit_time[i] += avg_traders_unit_time[i] / num_reps;
                    avg_avg_infected_unit_time[i] += avg_infected_unit_time[i] / num_reps;
                    avg_avg_trades_unit_time[i] += trades_unit_time[i] / (N * num_reps);
                    avg_avg_in_vol_unit_time[i] += in_vol_unit_time[i] / (N * num_reps);
                    avg_avg_infections_unit_time[i] += infections_unit_time[i] / num_reps;
                }

                //RESET STRUCTURE
                for (int i = 0; i < N; ++i) {
                    farms[i].traders.clear();
                    farms[i].not_traders.clear();
                    farms[i].tot_add_rate = 0.0;
                    farms[i].tot_delete_rate = 0.0;
                    farms[i].tot_trade_rate = 0.0;
                    farms[i].disease_state = 0;
                    farms[i].first_inf = 0;
                    farms[i].recovered = 0;
                    farms[i].infections = 0;
                }
            }
        }

        //reset trade and partnership params
        b /= epsilon_trade;
        theta *= epsilon_trade;
        alpha /= epsilon_ptnr;
        delta /= epsilon_ptnr;

        for (int i = 0; i < vector_size; ++i) {
            //calculate average disease prevalence averaged over all initially infected farms
            if (i * increment_size >= time_equi_inf) {
                avg_equi_inf += avg_avg_infected_unit_time[i] * increment_size / (t_max - time_equi_inf);
            }
        }
        out_file << inf_prob << "," << 1 / gamma << "," << new_theta << "," << new_time << "," << avg_r0_theory << ","
                 << avg_r0_sim  << "," << avg_equi_inf / num_reps_per_farm << endl;

        long double param_end_time = time(NULL);
        cout << "PARAM TIME TAKEN " << param_end_time - param_start_time << endl;
    }

    out_file.close();   //close output file
    long double sim_end_time = time(NULL);
    cout << "SIM TIME TAKEN " << sim_end_time - sim_start_time << endl;
}