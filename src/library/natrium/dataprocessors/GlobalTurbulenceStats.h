/**
 * @file SolverConfiguration.h
 * @short Calculation of physical properties
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef GlobalTurbulenceStats_H_
#define GlobalTurbulenceStats_H_

#include "../advection/AdvectionOperator.h"
#include "../utilities/BasicNames.h"
#include "../dataprocessors/DataProcessor.h"
#include "boost/filesystem.hpp"

namespace natrium {

template<size_t dim>
class GlobalTurbulenceStats: public DataProcessor<dim> {
private:
	std::string m_filename;
	std::string m_legendFilename;
	boost::shared_ptr<std::fstream> m_tableFile;
	bool m_outputOff;

	// Data
	// One point statistics
	size_t m_nofObservables;
	vector<string> m_names;
	vector<size_t> m_number;
	size_t m_nofCoordinates;
	vector<double> m_averages;
	vector<vector<double> > m_correlations;
	vector<double> m_EX3; // for skewness
	vector<double> m_EX4; // for kurtosis
	double m_energy;
	double m_enstrophy;
	double m_energySquared;
	double m_enstrophySquared;

	std::string outfile(std::string dir) {
		boost::filesystem::path out_dir(dir);
		boost::filesystem::path out_file = out_dir
				/ "global_turbulence_table.txt";
		return out_file.string();
	}

	std::string legendfile(std::string dir) {
		boost::filesystem::path out_dir(dir);
		boost::filesystem::path out_file = out_dir
				/ "global_turbulence_table.legend";
		return out_file.string();
	}

	void printHeaderLine();

	void writeToFile();

public:
	/// constructor
	GlobalTurbulenceStats(const CFDSolver<dim> & solver);

	/// destructor
	virtual ~GlobalTurbulenceStats() {

	}

	virtual void apply();

	void calculate();

	const vector<double>& getAverages() const {
		return m_averages;
	}

	const vector<vector<double> >& getCorrelations() const {
		return m_correlations;
	}

	const vector<double>& getEx3() const {
		return m_EX3;
	}

	const vector<double>& getEx4() const {
		return m_EX4;
	}

	const vector<string>& getNames() const {
		return m_names;
	}
};

} /* namespace natrium */

#endif /* GlobalTurbulenceStats_H_ */
