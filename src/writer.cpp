#include "w_optimizer.hpp"

void ORCA_inp_file::write_inp(const std::string& prefix, const std::string& extra_keyword, 
    int extra_charge, int extra_multiplicity, bool stable) const {
    std::ifstream ifile(prefix + ".chk", std::ios::binary);
    bool read_last = static_cast<bool>(ifile);
    if (ifile) ifile.close();

    fmt::ostream ofile(fmt::output_file(prefix + ".inp"));

    for (const std::string& line : keywords) {
        ofile.print("{:s}\n", line);
    }

    ofile.print("{:s}\n", extra_keyword);

    if (stable) {
        ofile.print("{:s}\n", "%SCF\n    StabPerform True\n    StabRestartUHFIfUnstable True\n    StabMaxIter 512\n    StabDTol 1.E-5\n    StabRTol 1.E-5\nEnd");
    }

    ofile.print("* xyz {:d} {:d}\n", charge + extra_charge, multiplicity + extra_multiplicity);
    for (size_t i = 0; i < natoms; ++ i) {
        ofile.print(" {:<2s}    {:13.8f}    {:13.8f}    {:13.8f}\n", elements[i], 
            coordinates(coord_x, i), coordinates(coord_y, i), coordinates(coord_z, i));
    }
    ofile.print("*\n");
}

