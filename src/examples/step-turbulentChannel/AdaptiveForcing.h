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



namespace natrium {


class AdaptiveForcing: public DataProcessor<3>{
private:

	boost::filesystem::path m_outDir;
	const vector<distributed_vector>& m_u;
	const distributed_vector& m_rho;

    double m_targetRhoU;
    double m_lastRhoU;
    double m_currentValue;
    double m_force;

    std::string m_filename;

    boost::shared_ptr<std::fstream> m_tableFile;


    void write();
    void getRhoU();

    std::string outfile(std::string dir) {
        boost::filesystem::path out_dir(dir);
        boost::filesystem::path out_file = out_dir
                                           / "adaptive_forcing.txt";
        return out_file.string();
    }

public:
	AdaptiveForcing(CFDSolver<3> & solver, std::string outdir, double target);
	virtual void apply();
	virtual ~AdaptiveForcing();




};

} /* namespace natrium */

#endif /* EXAMPLES_STEP_TURBULENTCHANNEL_ADAPTIVEFORCING_H_ */
