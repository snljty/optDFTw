#include "w_optimizer.hpp"

Gjf_file::Gjf_file(const std::string& template_prefix) {
    read_template(template_prefix);
}

void Gjf_file::read_template(const std::string& template_prefix) {
    const std::string ifilename = template_prefix + ".gjf";

    std::ifstream ifile(ifilename, std::ios::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: cannot open {:s} for reading.", ifilename));

    std::string line;
    std::regex chk_re(R"(%[Cc][Hh][Kk]=[\s\S]*)");
    std::regex blankline_re(R"(\s*)");

    link0.clear();
    route.clear();
    extra.clear();

    constexpr bool (*is_space)(unsigned char) = [](unsigned char c) { return static_cast<bool>(std::isspace(static_cast<int>(c))); };

    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blankline_re)) break;
        if (line[0] == '%') {
            if (!std::regex_match(line, chk_re)) link0.emplace_back(line.begin(), std::find_if_not(line.rbegin(), line.rend(), is_space).base());
        } else {
            route.emplace_back(line.begin(), std::find_if_not(line.rbegin(), line.rend(), is_space).base());
        }
    }

    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blankline_re)) break;
    }

    if (!std::getline(ifile, line)) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");

    std::istringstream iss0(line);
    if (!(iss0 >> charge >> multiplicity)) throw std::ios_base::failure("Error: cannot read charge and multiplicity.");

    std::streampos pos = ifile.tellg();

    natoms = 0;
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blankline_re)) break;
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

    std::getline(ifile, line);
    while (std::getline(ifile, line)) {
        if (std::regex_match(line, blankline_re)) break;
        extra.emplace_back(line.begin(), std::find_if_not(line.rbegin(), line.rend(), is_space).base());
    }

    ifile.close();
}

std::pair<double, double> Omega_optimizer::read_info(const std::string& prefix) const {
    std::ifstream ifile(prefix + ".out", std::ios_base::binary);
    if (!ifile) throw std::ios_base::failure(fmt::format("Error: Cannot open {:s}.out for reading.", prefix));

    std::string line;
    std::string::size_type pos;
    double ene, alpha_HOMO_ene, beta_HOMO_ene, HOMO_ene;
    bool found_ene = false, found_alpha_HOMO = false, found_beta_HOMO = false;
    std::streampos after_final_ene_pos;

    while (std::getline(ifile, line)) {
        if (line.find("SCF Done:") != std::string::npos) {
            after_final_ene_pos = ifile.tellg();
            if ((pos = line.find('=')) == std::string::npos) {
                throw std::ios_base::failure("Error: cannot read single-point energy.");
            }
            std::istringstream iss(line.substr(pos + std::strlen("-")));
            if (! (iss >> ene)) throw std::ios_base::failure("Error: cannot read single-point energy.");
            found_ene = true;
        }
    }

    if (!found_ene) throw std::ios_base::failure("Error: cannot read single-point energy.");
    ifile.clear();
    ifile.seekg(after_final_ene_pos);

    while (std::getline(ifile, line)) {
        if (line.find("Condensed to atoms") != std::string::npos) break;
        if (line.find("Alpha  occ. eigenvalues") != std::string::npos) {
            found_alpha_HOMO = true;
            if ((pos = line.find("--")) == std::string::npos) {
                throw std::ios_base::failure("Error: cannot read alpha orbital energy.");
            }
            std::istringstream iss(line.substr(pos + std::strlen("--")));
            while (iss >> alpha_HOMO_ene) {;}
        } else if (line.find("Beta  occ. eigenvalues") != std::string::npos) {
            found_beta_HOMO = true;
            if ((pos = line.find("--")) == std::string::npos) {
                throw std::ios_base::failure("Error: cannot read beta orbital energy.");
            }
            std::istringstream iss(line.substr(pos + std::strlen("--")));
            while (iss >> beta_HOMO_ene) {;}
        }
    }

    if (!found_alpha_HOMO && !found_beta_HOMO) throw std::ios_base::failure("Error: cannot read HOMO energy.");

    if (found_alpha_HOMO) {
        HOMO_ene = alpha_HOMO_ene;
        if (found_beta_HOMO && beta_HOMO_ene > alpha_HOMO_ene) HOMO_ene = beta_HOMO_ene;
    } else {
        HOMO_ene = beta_HOMO_ene;
    }

    ifile.close();

    return {ene, HOMO_ene};
}

