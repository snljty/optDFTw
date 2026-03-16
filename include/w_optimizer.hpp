#pragma once
#ifndef __W_OPTIMIZER_HPP__
#define __W_OPTIMIZER_HPP__

constexpr int ncoords = 3;
enum {coord_x, coord_y, coord_z};

#include <iostream>
#include <fstream>
#include <sstream>
// #include <filesystem>
#include <vector>
#include <string>
#include <regex>
#include <utility>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include <Eigen/Core>

#include <fmt/format.h>
#include <fmt/os.h>

#include <nlopt.hpp>

class ORCA_inp_file {
private:
    int natoms;
    std::vector<std::string> elements;
    Eigen::MatrixXd coordinates;
    std::vector<std::string> keywords;
    int charge;
    int multiplicity;
public:
    ORCA_inp_file(const std::string& template_prefix="template");
    void read_template(const std::string& template_prefix);
    void write_inp(const std::string& prefix, const std::string& extra_keyword, 
        int extra_charge=0, int extra_multiplicity=0, bool stable=false) const;
};

class Omega_optimizer {
public:
    void optimize();
    Omega_optimizer(double w_lower=0.05, double w_upper=0.6, int stable_rounds=1);
    double w() const;
private:
    ORCA_inp_file inp;
    int num_steps;
    double w_lower, w_upper, w_;
    const int stable_rounds;
    std::chrono::steady_clock::time_point time_start, time_end;
    mutable std::string orca_command;

    void generate_input(const std::string& prefix, int extra_charge=0, int extra_multiplicity=0, bool stable=false) const;
    void run(const std::string& prefix) const;
    void format_mkl(const std::string& prefix) const;
    std::pair<double, double> read_info(const std::string& prefix) const;
    double calc_J_squared() const;
    double step(double w);
    bool check_stable(const std::string& prefix, int extra_charge=0, int extra_multiplicity=0) const;
    void remove_garbage(const std::string& prefix) const;
    void get_orca_command() const;
    static double minimize_J_squared_callback(unsigned int n, const double* x, double* grad, void* f_data);
};

#endif // __W_OPTIMIZER_HPP__

