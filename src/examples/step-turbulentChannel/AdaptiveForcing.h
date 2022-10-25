/*
 * AdaptiveForcing.h
 *
 *  Created on: 06.10.2022
 *      Author: dwilde3m
 */

#ifndef EXAMPLES_STEP_TURBULENTCHANNEL_ADAPTIVEFORCING_H_
#define EXAMPLES_STEP_TURBULENTCHANNEL_ADAPTIVEFORCING_H_

#include "natrium/solver/TurbulenceStats.h"
#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/Math.h"
#include "boost/filesystem.hpp"
#include "FinalChannelStatistics.h"
#include "natrium/solver/CompressibleCFDSolver.h"




namespace natrium {


class AdaptiveForcing: public FinalChannelStatistics{
private:

	boost::filesystem::path m_outDir;

    const distributed_vector& m_T;

    double m_targetRhoU;
    double m_lastRhoU;
    double m_currentValueRhoU;
    double m_currentRho;
    double m_force;
    double m_starting_force;
    bool m_restart;

    CompressibleCFDSolver<3> & m_compressibleSolver;


    std::string m_filename;

    boost::shared_ptr<std::fstream> m_tableFile;


    void write();
    void calculateRhoU();
    void calculateForce();
    void changeTemperatureProfile();
    void rescaleDensity();
    void setBoundaryTemperature();


    std::string outfile(std::string dir) {
        boost::filesystem::path out_dir(dir);
        boost::filesystem::path out_file = out_dir
                                           / "adaptive_forcing.txt";
        return out_file.string();
    }

public:
	AdaptiveForcing(CompressibleCFDSolver<3> & solver, std::string outdir, double target, bool restart=false);
	virtual void apply();
	virtual ~AdaptiveForcing();


    void updateCompressibleAverages();

    double m_heatFactor;
    double m_integral;
};

} /* namespace natrium */

#endif /* EXAMPLES_STEP_TURBULENTCHANNEL_ADAPTIVEFORCING_H_ */
