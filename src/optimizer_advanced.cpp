#include "w_optimizer.hpp"

double Omega_optimizer::step(double w) {
    std::chrono::steady_clock::time_point current_time_start(std::chrono::steady_clock::now()), current_time_end;
    w_ = w;
    ++ num_steps;
    fmt::print("Step {:d}, w = {:.4f}\n", num_steps, w_);
    fmt::print("Generating input files ...\n");
    generate_input("N");
    generate_input("N-1",  1, 1, num_steps <= stable_rounds);
    generate_input("N+1", -1, 1, num_steps <= stable_rounds);
    fmt::print("Running ORCA for state N   ...\n");
    std::fflush(stdout);
    run("N");
    fmt::print("Running ORCA for state N-1 ...\n");
    std::fflush(stdout);
    run("N-1");
    fmt::print("Running ORCA for state N+1 ...\n");
    std::fflush(stdout);
    run("N+1");
    fmt::print("Calculating J**2 ...\n");
    double ret = calc_J_squared();
    current_time_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span(std::chrono::duration_cast<std::chrono::duration<double> >(current_time_end - current_time_start));
    fmt::print("w = {:.4f}, J**2 = {:.10f}\n", w, ret);
    fmt::print("Time elapsed for this step: {:.1f} s\n", time_span.count());
    std::fflush(stdout);
    fmt::print("\n");
    return ret;
}

double Omega_optimizer::minimize_J_squared_callback(unsigned int n, const double* x, double* grad, void* f_data) {
    Omega_optimizer* opt = static_cast<Omega_optimizer*>(f_data);
    return opt->step(x[0]);
}

void Omega_optimizer::optimize() {
    time_start = std::chrono::steady_clock::now();

    nlopt::opt opt(nlopt::LN_NELDERMEAD, 1);
    opt.set_min_objective(&Omega_optimizer::minimize_J_squared_callback, this);

    opt.set_lower_bounds(std::vector<double>{w_lower});
    opt.set_upper_bounds(std::vector<double>{w_upper});

    opt.set_xtol_rel(1.E-4);

    double min_J_squared;
    std::vector<double> w_init = {w_};
    nlopt::result result = opt.optimize(w_init, min_J_squared);

    if (result < 0) throw std::runtime_error("Error: Brent's optimization failed!");

    std::chrono::steady_clock::time_point stable_time_start(std::chrono::steady_clock::now()), stable_time_end;
    if (!check_stable("N-1") || !check_stable("N+1")) throw std::runtime_error("Error: some of the converged wavefunctions are not stable.");
    stable_time_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> stable_time_span(std::chrono::duration_cast<std::chrono::duration<double> >(stable_time_end - stable_time_start));
    fmt::print("Time elapsed for final wavefunction stability analysys: {:.1f} s\n", stable_time_span.count());

    fmt::print("Converged. The optimal w is {:.4f}.\n", w_);
    fmt::print("you can use \n\"\"\"\n%Method\n    RangeSepMu {:.4f}\nEnd\n\"\"\"\nin your .inp file to use this w value.", w_);

    format_mkl("N");
    format_mkl("N-1");
    format_mkl("N+1");

    time_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span(std::chrono::duration_cast<std::chrono::duration<double> >(time_end - time_start));
    fmt::print("Total time elapsed: {:.1f} s\n", time_span.count());

    fmt::ostream ofile(fmt::output_file("result_option_w.txt"));
    ofile.print("%Method\n    RangeSepMu {:.4f}\nEnd\n", w_);
    ofile.close();
    fmt::print("Result has been saved to \"{:s}\".\n", "result_option_w.txt");

    remove_garbage("N");
    remove_garbage("N-1");
    remove_garbage("N+1");
}

