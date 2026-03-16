#include "w_optimizer.hpp"

double Omega_optimizer::w() const {
    return w_;
}

Omega_optimizer::Omega_optimizer(double w_lower, double w_upper, int stable_rounds) : 
    inp(), num_steps(0), w_lower(w_lower), w_upper(w_upper), w_((w_lower + w_upper) / 2.), stable_rounds(stable_rounds) {}

void Omega_optimizer::generate_input(const std::string& prefix, int extra_charge, int extra_multiplicity, bool stable) const {
    inp.write_inp(prefix, fmt::format("%Method\n    RangeSepMu {:.4f}\nEnd", w_), 
        extra_charge, extra_multiplicity, stable);
}

void Omega_optimizer::run(const std::string& prefix) const {
    if (orca_command.empty()) get_orca_command();
    #ifdef _WIN32
    #define popen _popen
    #define pclose _pclose
    FILE* pipe = popen(fmt::format("{:s} {:s}.inp", orca_command, prefix).c_str(), "r");
    #else
    FILE* pipe = popen(fmt::format("'{:s}' '{:s}.inp'", orca_command, prefix).c_str(), "r");
    #endif
    if (!pipe) throw std::ios_base::failure("Error: cannot execute ORCA.");

    char buf[BUFSIZ];
    fmt::ostream ofile(fmt::output_file(prefix + ".out"));
    while (fgets(buf, sizeof(buf), pipe)) {
        ofile.print("{:s}", buf);
    }

    int status = pclose(pipe);
    if (status == -1) throw std::ios_base::failure("Error: fclose failed.");
    #ifdef __linux__
    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
            throw std::ios_base::failure(fmt::format("Error: ORCA exit status: {:d}.", WEXITSTATUS(status)));
        }
    } else if (WIFSIGNALED(status)) {
        throw std::ios_base::failure(fmt::format("Error: ORCA killed by signal: {:d}.", WTERMSIG(status)));
    }
    #else
    if (status) {
        throw std::ios_base::failure(fmt::format("Error: ORCA exit status: {:d}.", status));
    }
    #endif
    ofile.close();
}

void Omega_optimizer::format_mkl(const std::string& prefix) const {
    #ifdef _WIN32
    #define popen _popen
    #define pclose _pclose
    FILE* pipe = popen(fmt::format("{:s} {:s} -molden", "orca_2mkl", prefix).c_str(), "r");
    #else
    FILE* pipe = popen(fmt::format("'{:s}' '{:s}' -molden", "orca_2mkl", prefix).c_str(), "r");
    #endif
    if (!pipe) throw std::ios_base::failure("Error: cannot execute ORCA orca_2mkl.");

    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), pipe)) {;}

    int status = pclose(pipe);
    if (status == -1) throw std::ios_base::failure("Error: fclose failed.");
    #ifdef __linux__
    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
            throw std::ios_base::failure(fmt::format("Error: ORCA orca_2mkl exit status: {:d}.", WEXITSTATUS(status)));
        }
    } else if (WIFSIGNALED(status)) {
        throw std::ios_base::failure(fmt::format("Error: ORCA orca_2mkl by signal: {:d}.", WTERMSIG(status)));
    }
    #else
    if (status) {
        throw std::ios_base::failure(fmt::format("Error: ORCA orca_2mkl exit status: {:d}.", status));
    }
    #endif
    // std::filesystem::rename(prefix + ".molden.input", prefix + ".molden");
    // std::filesystem::remove(prefix + ".chk");
    std::rename((prefix + ".molden.input").c_str(), (prefix + ".molden").c_str());
    std::remove((prefix + ".gbw").c_str());
}

double Omega_optimizer::calc_J_squared() const {
    auto [E_N, HOMO_N] = read_info("N");
    auto [E_Nm1, HOMO_Nm1] = read_info("N-1");
    auto [E_Np1, HOMO_Np1] = read_info("N+1");
    double J_N = HOMO_N + E_Nm1 - E_N;
    double J_Np1 = HOMO_Np1 + E_N - E_Np1;
    double J_squared = J_N * J_N + J_Np1 * J_Np1;
    return J_squared;
}

bool Omega_optimizer::check_stable(const std::string& prefix, int extra_charge, int extra_multiplicity) const {
    inp.write_inp(prefix, fmt::format("%SCF\n    StabPerform True\n    StabDTol 1.E-5\n    StabRTol 1.E-5\nEnd\n%Method\n    RangeSepMu {:.4f}\nEnd", w_), 
        extra_charge, extra_multiplicity);

    run(prefix);

    std::ifstream ifile(prefix + ".out");
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: Cannot open {:s}.out for reading.", prefix));

    bool stable = false;
    std::string line;
    while (std::getline(ifile, line)) {
        if (line.find("The stability analysis shows that the wavefunction is stable") != std::string::npos) {
            stable = true;
            break;
        }
    }

    ifile.close();
    return stable;
}

void Omega_optimizer::remove_garbage(const std::string& prefix) const {
    for (const std::string suffix : {".bibtex", ".densities", ".densitiesinfo", ".property.txt", ".ges"}) {
        std::ifstream ifile(prefix + suffix);
        if (ifile) {
            ifile.close();
            std::remove((prefix + suffix).c_str());
        }
    }
}

void Omega_optimizer::get_orca_command() const {
    #ifdef _WIN32
    const char* orca_variable = std::getenv("ORCA");
    orca_command = orca_variable && *orca_variable ? std::string(orca_variable) : "orca";
    #else
    FILE* orca_pipe = popen("which orca", "r");
    if (!orca_pipe) throw std::ios_base::failure("Error: cannot get ORCA executable.");
    
    char buf[BUFSIZ];
    if (!fgets(buf, sizeof(buf), orca_pipe)) throw std::ios_base::failure("Error: cannot get ORCA executable.");
    size_t buf_len = std::strlen(buf);
    while (std::isspace(buf[buf_len - 1])) {
        buf[-- buf_len] = '\0';
    }
    orca_command = std::string(buf);
    while (fgets(buf, sizeof(buf), orca_pipe)) {;}

    int status = pclose(orca_pipe);
    if (status == -1) throw std::ios_base::failure("Error: fclose failed.");
    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
            throw std::ios_base::failure(fmt::format("Error: fclose exit status: {:d}.", WEXITSTATUS(status)));
        }
    } else if (WIFSIGNALED(status)) {
        throw std::ios_base::failure(fmt::format("Error: fclose killed by signal: {:d}.", WTERMSIG(status)));
    }
    #endif
}

