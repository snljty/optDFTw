#include "w_optimizer.hpp"

double Omega_optimizer::w() const {
    return w_;
}

Omega_optimizer::Omega_optimizer(double w_lower, double w_upper, int stable_rounds) : 
    gjf(), num_steps(0), w_lower(w_lower), w_upper(w_upper), w_((w_lower + w_upper) / 2.), stable_rounds(stable_rounds) {}

void Omega_optimizer::generate_input(const std::string& prefix, int extra_charge, int extra_multiplicity, bool stable) const {
    gjf.write_gjf(prefix, fmt::format("IOp(3/107={0:05d}00000,3/108={0:05d}00000)", static_cast<int>(std::round(w_ * 1.E4))), 
        extra_charge, extra_multiplicity, stable);
}

void Omega_optimizer::run(const std::string& prefix) const {
    #ifdef _WIN32
    #define popen _popen
    #define pclose _pclose
    #endif
    FILE *pipe = popen(fmt::format("{0:s} \"{1:s}.gjf\" \"{1:s}.out\"", "g16", prefix).c_str(), "r");
    if (!pipe) throw std::ios_base::failure("Error: cannot execute Gaussian.");

    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), pipe)) {};

    int status = pclose(pipe);
    if (status == -1) throw std::ios_base::failure("Error: fclose failed.");
    #ifdef __linux__
    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
            throw std::ios_base::failure(fmt::format("Error: Gaussian exit status: {:d}.", WEXITSTATUS(status)));
        }
    } else if (WIFSIGNALED(status)) {
        throw std::ios_base::failure(fmt::format("Error: Gaussian by signal: {:d}.", WTERMSIG(status)));
    }
    #else
    if (status) {
        throw std::ios_base::failure(fmt::format("Error: Gaussian exit status: {:d}.", status));
    }
    #endif
}

void Omega_optimizer::formchk(const std::string& prefix) const {
    #ifdef _WIN32
    #define popen _popen
    #define pclose _pclose
    #endif
    FILE *pipe = popen(fmt::format("{0:s} \"{1:s}.chk\" \"{1:s}.fch\"", "formchk", prefix).c_str(), "r");
    if (!pipe) throw std::ios_base::failure("Error: cannot execute Gaussian formchk.");

    char buf[BUFSIZ];
    while (fgets(buf, sizeof(buf), pipe)) {};

    int status = pclose(pipe);
    if (status == -1) throw std::ios_base::failure("Error: fclose failed.");
    #ifdef __linux__
    if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
            throw std::ios_base::failure(fmt::format("Error: Gaussian formchk exit status: {:d}.", WEXITSTATUS(status)));
        }
    } else if (WIFSIGNALED(status)) {
        throw std::ios_base::failure(fmt::format("Error: Gaussian formchk by signal: {:d}.", WTERMSIG(status)));
    }
    #else
    if (status) {
        throw std::ios_base::failure(fmt::format("Error: Gaussian formchk exit status: {:d}.", status));
    }
    #endif
    std::filesystem::remove(prefix + ".chk");
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
    gjf.write_gjf(prefix, fmt::format("IOp(3/107={0:05d}00000,3/108={0:05d}00000) Stable", static_cast<int>(std::round(w_ * 1.E4))), 
        extra_charge, extra_multiplicity);

    run(prefix);

    std::ifstream ifile(prefix + ".out");
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: Cannot open {:s}.out for reading.", prefix));

    bool stable = false;
    std::string line;
    while (std::getline(ifile, line)) {
        if (line.find("The wavefunction is stable under the perturbations considered") != std::string::npos) {
            stable = true;
            break;
        }
    }

    ifile.close();
    return stable;
}

