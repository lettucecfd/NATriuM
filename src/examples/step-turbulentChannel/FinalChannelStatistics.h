/*
 * FinalChannelStatistics.h
 *
 *  Created on: 03.03.2016
 *      Author: akraem3m
 */

#ifndef EXAMPLES_STEP_TURBULENTCHANNEL_FINALCHANNELSTATISTICS_H_
#define EXAMPLES_STEP_TURBULENTCHANNEL_FINALCHANNELSTATISTICS_H_

#include "natrium/solver/TurbulenceStats.h"
#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/solver/CFDSolver.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/Math.h"
#include "boost/filesystem.hpp"



namespace natrium {


class FinalChannelStatistics: public DataProcessor<3>{
private:
	boost::shared_ptr<std::fstream> m_tableFile;
	boost::filesystem::path m_outDir;
	const vector<distributed_vector>& m_u;
	const distributed_vector& m_rho;
	bool m_yCoordsUpToDate;

	// Data
	std::vector<double> m_yCoordinates;
	std::map<double, size_t, own_double_less> m_yCoordinateToIndex;
	// One point statistics
	size_t m_nofObservables;
	vector<string> m_names;
	vector<size_t> m_number;
	size_t m_nofCoordinates;
	vector<vector<double> > m_averages;
	vector<vector<vector<double> > > m_correlations;
	vector<vector<double> > m_EX3; // for skewness
	vector<vector<double> > m_EX4; // for kurtosis

	// temporal averages
	size_t n_steps;
	vector<vector<double> > m_averages_time;
	vector<vector<vector<double> > > m_correlations_time;
	vector<vector<double> > m_EX3_time; // for skewness
	vector<vector<double> > m_EX4_time; // for kurtosis


public:
	FinalChannelStatistics(const CFDSolver<3> & solver, std::string outdir);
	virtual void apply();
	void update();
	void updateYValues();
	void updateAverages();
	void addToTemporalAverages();
	void write_to_file();
	virtual ~FinalChannelStatistics();
};

} /* namespace natrium */

#endif /* EXAMPLES_STEP_TURBULENTCHANNEL_FINALCHANNELSTATISTICS_H_ */
