//
// Created by dwilde3m on 04.02.21.
//

#ifndef NATRIUM_COMPRESSIBLETURBULENCESTATS_H
#define NATRIUM_COMPRESSIBLETURBULENCESTATS_H

#include "../advection/AdvectionOperator.h"
#include "../utilities/BasicNames.h"
#include "../dataprocessors/DataProcessor.h"
#include "boost/filesystem.hpp"
//#include "../solver/CFDSolver.h"

namespace natrium{
// forward declaration
template <size_t dim>
class CompressibleCFDSolver;

template<size_t dim>
class CompressibleTurbulenceStats {



private:
    CompressibleCFDSolver<dim> & m_solver;
    double m_dilatation;
    double m_solenoidal;
    double m_maxMach;
    double m_gamma;
    double m_totalEnergy;
    vector<string> m_names;
    std::string m_filename;
    std::string m_legendFilename;
    boost::shared_ptr<std::fstream> m_tableFile;
    bool m_outputOff;



std::string outfile(std::string dir) {
    boost::filesystem::path out_dir(dir);
    boost::filesystem::path out_file = out_dir
                                       / "compressible_turbulence_table.txt";
    return out_file.string();
}

std::string legendfile(std::string dir) {
    boost::filesystem::path out_dir(dir);
    boost::filesystem::path out_file = out_dir
                                       / "compressible_turbulence_table.legend";
    return out_file.string();
}

    void printHeaderLine();

    void writeToFile();



public:
    /// constructor
    CompressibleTurbulenceStats(CompressibleCFDSolver<dim> & solver);

    /// destructor
    virtual ~CompressibleTurbulenceStats() {

    }

    void apply();

    void calculate();

};
}
#endif //NATRIUM_COMPRESSIBLETURBULENCESTATS_H
