#include "w_optimizer.hpp"

void Gjf_file::write_gjf(const std::string& prefix, const std::string& extra_keyword, 
    int extra_charge, int extra_multiplicity, bool stable) const {
    std::ifstream ifile(prefix + ".chk", std::ios::binary);
    bool read_last = static_cast<bool>(ifile);
    if (ifile) ifile.close();

    fmt::ostream ofile(fmt::output_file(prefix + ".gjf"));

    for (const std::string& line : link0) {
        ofile.print("{:s}\n", line);
    }

    ofile.print("%chk={:s}.chk\n", prefix);

    for (size_t i = 0; i < route.size(); ++ i) {
        if (i + 1 != route.size()) {
            ofile.print("{:s}\n", route[i]);
        } else {
            ofile.print("{:s}", route[i]);
            if (!extra_keyword.empty()) ofile.print(" {:s}", extra_keyword);
            if (stable) ofile.print(" Stable=Opt");
            if (read_last) ofile.print(" Guess=Read");
            ofile.print("\n");
        }
    }

    ofile.print("\nfile {:s}\n\n", prefix);
    ofile.print("{:d} {:d}\n", charge + extra_charge, multiplicity + extra_multiplicity);
    for (size_t i = 0; i < natoms; ++ i) {
        ofile.print(" {:<2s}    {:13.8f}    {:13.8f}    {:13.8f}\n", elements[i], 
            coordinates(coord_x, i), coordinates(coord_y, i), coordinates(coord_z, i));
    }
    ofile.print("\n");
    if (!extra.empty()) {
        for (const std::string& line : extra) {
            ofile.print("{:s}\n", line);
        }
        ofile.print("\n");
    }

    ofile.close();
}

