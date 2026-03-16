#include "w_optimizer.hpp"

ORCA_inp_file::ORCA_inp_file(const std::string& template_prefix) {
    read_template(template_prefix);
}

void ORCA_inp_file::read_template(const std::string& template_prefix) {
    const std::string ifilename = template_prefix + ".inp";

    std::ifstream ifile(ifilename, std::ios_base::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: cannot open {:s} for reading.", ifilename));

    std::string line;
    std::regex coords_re(R"(\s*\*\s*xyz\s+[-0-9]+\s+[0-9]+\s*)");

    keywords.clear();

    constexpr bool (*is_space)(unsigned char) = [](unsigned char c) { return static_cast<bool>(std::isspace(static_cast<int>(c))); };

    bool found_coords = false;
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, coords_re)) {
            found_coords = true;
            break;
        }
        keywords.emplace_back(line.begin(), std::find_if_not(line.rbegin(), line.rend(), is_space).base());
    }

    if (!found_coords) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");
    std::istringstream iss0(line.substr(line.find("xyz") + std::strlen("xyz")));;
    if (!(iss0 >> charge >> multiplicity)) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");

    std::streampos pos = ifile.tellg();

    natoms = 0;
    while (std::getline(ifile, line)) {
        if (line.find('*') != std::string::npos) break;
        ++ natoms;
    }

    ifile.seekg(pos);

    elements.resize(natoms);
    coordinates.resize(ncoords, natoms);
    for (int i = 0; i < natoms; ++ i) {
        std::getline(ifile, line);
        std::istringstream iss(line);
        if (!(iss >> elements[i] >> coordinates(coord_x, i) >> coordinates(coord_y, i) >> coordinates(coord_z, i))) {
            throw std::ios_base::failure(fmt::format("Error: cannot read atom {:d}.", i + 1));
        }
    }

    ifile.close();
}

std::pair<double, double> Omega_optimizer::read_info(const std::string& prefix) const {
    std::ifstream ifile(prefix + ".out", std::ios_base::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: Cannot open {:s}.out for reading.", prefix));

    bool found_HOMO = false;
    std::string line;
    std::regex blankline_re(R"(\s*)");
    while (std::getline(ifile, line)) {
        if (line.find("ORBITAL ENERGIES") != std::string::npos) {
            found_HOMO = true;
            break;
        }
    }
    if (!found_HOMO) throw std::ios_base::failure("Error: cannot read orbital energy.");
    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy."); // ---
    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy.");

    int orb_index;
    double occ;
    bool found_alpha_HOMO = false, found_beta_HOMO = false;
    double HOMO_ene, alpha_HOMO_ene, beta_HOMO_ene;
    found_HOMO = false;

    if (std::regex_match(line, blankline_re)) {
        // space orbitals
        if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy."); // NO OCC
        while (std::getline(ifile, line)) {
            std::istringstream iss(line);
            if (!(iss >> orb_index)) break;
            if (!(iss >> occ)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            if (occ) {
                found_HOMO = true;
                if (!(iss >> HOMO_ene)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            }
        }
        if (!found_HOMO) throw std::ios_base::failure("Error: cannot read HOMO energy.");
    } else {
        // spin orbitals
        if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy.");
        while (std::getline(ifile, line)) {
            std::istringstream iss(line);
            if (!(iss >> orb_index)) break;
            if (!(iss >> occ)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            if (occ) {
                found_alpha_HOMO = true;
                if (!(iss >> alpha_HOMO_ene)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            }
        }
        if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy."); // blank
        if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read orbital energy."); // NO OCC
        while (std::getline(ifile, line)) {
            std::istringstream iss(line);
            if (!(iss >> orb_index)) break;
            if (!(iss >> occ)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            if (occ) {
                found_beta_HOMO = true;
                if (!(iss >> beta_HOMO_ene)) throw std::ios_base::failure("Error: cannot read orbital energy.");
            }
        }
        if (found_alpha_HOMO) {
            found_HOMO = true;
            HOMO_ene = alpha_HOMO_ene;
            if (found_beta_HOMO && beta_HOMO_ene > alpha_HOMO_ene) HOMO_ene = beta_HOMO_ene;
        } else {
            if (!found_beta_HOMO) throw std::ios_base::failure("Error: cannot read HOMO energy.");
            found_HOMO = true;
            HOMO_ene = beta_HOMO_ene;
        }
    }

    double ene;
    std::string::size_type pos;
    bool found_ene = false;
    while (std::getline(ifile, line)) {
        if ((pos = line.find("FINAL SINGLE POINT ENERGY")) != std::string::npos) {
            std::istringstream iss(line.substr(pos + std::strlen("FINAL SINGLE POINT ENERGY")));
            if (!(iss >> ene)) throw std::ios_base::failure("Error: cannot read single-point energy.");
            found_ene = true;
            break;
        }
    }
    if (!found_ene) throw std::ios_base::failure("Error: cannot read single-point energy.");

    ifile.close();

    return {ene, HOMO_ene};
}

