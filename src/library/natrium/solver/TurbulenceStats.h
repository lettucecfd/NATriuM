/*
 * TurbulenceStats.h
 *
 *  Created on: 12.02.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_SOLVER_TURBULENCESTATS_H_
#define LIBRARY_NATRIUM_SOLVER_TURBULENCESTATS_H_

#include <math.h>

#include "../utilities/BasicNames.h"
#include "deal.II/numerics/data_out.h"

namespace natrium {

// forward declaration
template<size_t dim> class CFDSolver;

/**
 * @short This class calculates the turbulence statistics in planes parallel to the wall
 */
template<size_t dim> class StatsInPlane {
private:
	size_t m_wallNormalDirection;
	double m_wallNormalCoordinate;

	dealii::IndexSet m_planeIndices;
	const vector<distributed_vector>& m_u;
public:
	/**
	 * @short constructor
	 * @param wall_normal_direction 0,1,2 for x,y,z; respectively
	 * @param wall_normal_coordinate the value of the wall-normal coordinate which defines the plane (e.g. y=0.01).
	 * @param u reference to the velocity vector
	 */
	StatsInPlane(size_t wall_normal_direction, double wall_normal_coordinate,
			const vector<distributed_vector>& u) :
			m_wallNormalDirection(wall_normal_direction), m_wallNormalCoordinate(
					wall_normal_coordinate), m_u(u) {
	}
	/**
	 * @short Update the index set that stores all indices belonging to the plane.
	 * @param support_points Reference to the map of support points that connects dof indices and grid points.
	 * @param locally_owned index set that contains the indices belonging to locally owned degrees of freedom
	 * @param tolerance Tolerance to compare prescribed wall normal coordinate with grid point coordinate
	 */
	void updatePlaneIndices(
			const map<dealii::types::global_dof_index, dealii::Point<dim> >& support_points,
			const dealii::IndexSet& locally_owned, double tolerance);
	/**
	 * @short Calculate turbulence statistics over the plane.
	 * @param[out] u_average Vector that stores U,V (and W). Is automatically resized to dim.
	 * @param[out] u_average Vector that stores u',v' (and w'). Is automatically resized to dim.
	 */
	void calculate(vector<double>& u_average, vector<double>& u_rms);
};

/**
 * @short A class that calculates and puts out turbulent statistics.
 * It provides the functionality to evaluate the statistics over wall-parallel planes:
 * For many relevant flows (e.g. channels, boundary layers, ...),
 * the turbulence is homogenous in one direction, which allows
 * to calculate time-independent statistics.
 */
template<size_t dim>
class TurbulenceStats {
private:
	CFDSolver<dim> * m_solver;
	boost::shared_ptr<std::fstream> m_tableFile;
	std::string m_filename;
	const bool m_outputOff;
	size_t m_iterationNumber;

	// Reynolds statistics
	size_t m_statSize;
    vector<distributed_vector> m_uAverage;

	// Convergence statistics
	size_t m_wallNormalDirection;
	vector<double> m_wallNormalCoordinates;
	vector<vector<double> > m_averages;

	vector<vector<double> > m_rms;
	vector<StatsInPlane<dim> > m_planes;

public:
	TurbulenceStats(CFDSolver<dim> * solver, size_t wall_normal_direction,
			const vector<double>& wall_normal_coordinates,
			const std::string table_file_name = "");
	virtual ~TurbulenceStats();

	void printHeaderLine();

	void update();

	void printNewLine();

	bool isUpToDate() const {
		return (m_iterationNumber == m_solver->getIteration());
	}

    void addToReynoldsStatistics(const vector<distributed_vector>& u);
    void addReynoldsStatisticsToOutput(dealii::DataOut<dim>& data_out);
};

/**
 * @short A class that calculates and puts out turbulent statistics.
 * It provides the functionality to evaluate the statistics over wall-parallel planes:
 * For many relevant flows (e.g. channels, boundary layers, ...),
 * the turbulence is homogenous in one direction, which allows
 * to calculate time-independent statistics.
 */
template<size_t dim>
class ShearLayerStats {
private:
    CFDSolver<dim> * m_solver;
    boost::shared_ptr<std::fstream> m_tableFile;
    std::string m_filename;
    const bool m_outputOff;
    size_t m_iterationNumber;

    // Y-Axis stats
    size_t m_statSize;
    size_t m_nofCoordinates;
    std::map<double, size_t, own_double_less> m_yCoordinateToIndex;
    std::vector<double> m_yCoordinates;
    bool m_yCoordsUpToDate;
    // One point statistics
    size_t m_nofObservables;
    vector<string> m_names;
    vector<size_t> m_number;
    // temporal averages
//    size_t n_steps;
//    vector<vector<double> > m_averages_time;
//    vector<vector<vector<double> > > m_correlations_time;
//    vector<vector<double> > m_EX3_time; // for skewness
//    vector<vector<double> > m_EX4_time; // for kurtosis

    // Reynolds statistics
    distributed_vector m_RhoUx;
    distributed_vector m_Rho;
    distributed_vector m_UxFavre;
    distributed_vector m_DeltaTheta;
//    vector<distributed_vector> m_uAverage;
//    distributed_vector m_rhoXZAverage;
//    vector<vector<distributed_vector>> m_uXZAverage;

    //vector<distributed_vector> u_;

    // Convergence statistics
    vector<vector<double> > m_averages;

    vector<vector<double> > m_rms;
    vector<StatsInPlane<dim> > m_planes;

public:
    ShearLayerStats(CFDSolver<dim> * solver, const std::string table_file_name = "");
    virtual ~ShearLayerStats();

    void printHeaderLine();

    void update();

    void printNewLine();

    bool isUpToDate() const {
        return (m_iterationNumber == m_solver->getIteration());
    }

    bool isMYCoordsUpToDate() const {
        return m_yCoordsUpToDate;
    }

    void addToReynoldsAveragesXZ(
            const vector<distributed_vector>& u,
            const distributed_vector& rho);
    void addReynoldsAveragesXZToOutput(dealii::DataOut<dim>& data_out);

        void updateYValues();
    };

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_SOLVER_TURBULENCESTATS_H_ */
