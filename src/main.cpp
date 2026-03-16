#include "w_optimizer.hpp"

int main(int argc, char const *argv[]) {
    constexpr double w_lower_default = 0.05, w_upper_default = 0.3;
    double w_lower = w_lower_default, w_upper = w_upper_default;

    constexpr int stable_rounds_default = 1;
    int stable_rounds = stable_rounds_default;

    if (argc - 1 == 2 || argc - 1 == 3) {
        w_lower = std::stod(argv[1]);
        w_upper = std::stod(argv[2]);
        if (argc - 1 == 3) stable_rounds = std::stoi(argv[3]);
    } else if (argc - 1 != 0) {
        throw std::invalid_argument(fmt::format("Usage: {:s} [w_lower] [w_upper] [stable_rounds]", argv[0]));
    }

    Omega_optimizer w_optimizer(w_lower, w_upper, stable_rounds);
    w_optimizer.optimize();

    return 0;
}
