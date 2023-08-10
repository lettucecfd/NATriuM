/*
 * ShearLayerStats.h
 *
 *  Created on: 06.10.2022
 *      Author: dwilde3m
 */

#ifndef EXAMPLES_STEP_TURBULENTCHANNEL_SHEARLAYERSTATS_H_
#define EXAMPLES_STEP_TURBULENTCHANNEL_SHEARLAYERSTATS_H_

#include "natrium/solver/TurbulenceStats.h"
#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/Math.h"
#include "boost/filesystem.hpp"
//#include "MixingLayerStatistics.h"
#include "natrium/solver/CompressibleCFDSolver.h"

namespace natrium {

class ShearLayerStats: public DataProcessor<3>{
private:

    // Input parameters
    double m_Re0;
    const vector<distributed_vector>& m_u;
    const distributed_vector& m_rho;

    // Output parameters
	boost::filesystem::path m_outDir;
//    vector<string> m_names;
    string m_filename;
    boost::shared_ptr<std::fstream> m_tableFile;
    string m_vectorfilename;
    boost::shared_ptr<std::fstream> m_vectorFile;

    // Y Coordinates
    vector<double> m_yCoordinates;
    int nround;
    std::map<double, size_t, own_double_less> m_yCoordinateToIndex;
    size_t m_nofCoordinates;
    bool m_yCoordsUpToDate;

    // Data
    double m_currentDeltaTheta;
    double m_currentDeltaOmega;
    double m_ReOmega;
    double m_deltaThetaGrowth;
    double m_b11;
    double m_b22;
    double m_b12;
//    // Data stored across output steps
//    double m_currentTime;
    // Data stored across y
    vector<double> m_R11, m_R22, m_R33, m_R12;
    vector<double> m_K;

    void write();
    void calculateRhoU();

    static string scalaroutfile(string dir) {
        boost::filesystem::path out_dir(dir);
        boost::filesystem::path out_file = out_dir / "shearlayer_scalars.txt";
        return out_file.string();
    }
    static string vectoroutfile(string dir) {
        boost::filesystem::path out_dir(dir);
        boost::filesystem::path out_file = out_dir / "shearlayer_vectors.txt";
        return out_file.string();
    }

public:
	ShearLayerStats(CompressibleCFDSolver<3> & solver, string outdir, double starting_delta_theta, double starting_Re);
	void apply() override;
	~ShearLayerStats() override;
    bool isMYCoordsUpToDate() const;
    void updateYValues();
    };

} /* namespace natrium */

#endif /* EXAMPLES_STEP_TURBULENTCHANNEL_SHEARLAYERSTATS_H_ */
