/*
 * TurbulenceStats.cpp
 *
 *  Created on: 12.02.2016
 *      Author: akraem3m
 */

#include "TurbulenceStats.h"

#include "deal.II/base/index_set.h"

#include "../utilities/CFDSolverUtilities.h"

namespace natrium {

template<size_t dim>
void StatsInPlane<dim>::updatePlaneIndices(
		const map<dealii::types::global_dof_index, dealii::Point<dim> >& support_points,
		const dealii::IndexSet& locally_owned,
        double tolerance) {
	m_planeIndices.clear();
	// the size of the set of which this index set is a subset
	// has to be specified before adding indices
	m_planeIndices.set_size(locally_owned.size());
	//for all degrees of freedom on current processor
	dealii::IndexSet::ElementIterator it(locally_owned.begin());
	dealii::IndexSet::ElementIterator end(locally_owned.end());
	for (; it != end; it++) {
		size_t i = *it;
		// check if support point is in plane
		const dealii::Point<dim>& p = support_points.at(i);
		if (fabs(p(m_wallNormalDirection) - m_wallNormalCoordinate)
				< tolerance) {
			m_planeIndices.add_index(i);
		}

	}
}
template<size_t dim>
void StatsInPlane<dim>::calculate(vector<double>& u_average,
		vector<double>& u_rms) {
	size_t n = m_planeIndices.n_elements();
	u_average.resize(dim);
	u_rms.resize(dim);
	u_average.at(0) = 0;
	u_rms.at(0) = 0;
	u_average.at(1) = 0;
	u_rms.at(1) = 0;
	if (dim == 3) {
		u_average.at(2) = 0;
		u_rms.at(2) = 0;
	}

	// Average
	dealii::IndexSet::ElementIterator it(m_planeIndices.begin());
	dealii::IndexSet::ElementIterator end(m_planeIndices.end());
	for (; it != end; it++) {
		size_t i = *it;
		u_average[0] += m_u[0][i];
		u_average[1] += m_u[1][i];
		if (dim == 3)
			u_average[2] += m_u[2][i];
	}
	// MPI sync
	vector<double> u_average_mpi(dim);
	dealii::Utilities::MPI::sum(u_average, MPI_COMM_WORLD, u_average_mpi);
	size_t n_mpi = dealii::Utilities::MPI::sum(n, MPI_COMM_WORLD);
	u_average[0] = u_average_mpi[0] / n_mpi;
	u_average[1] = u_average_mpi[1] / n_mpi;
	if (dim == 3)
		u_average[2] = u_average_mpi[2] / n_mpi;

	// RMS value
	it = m_planeIndices.begin();
	for (; it != end; it++) {
		size_t i = *it;
		u_rms[0] += pow(m_u[0][i] - u_average[0], 2);
		u_rms[1] += pow(m_u[1][i] - u_average[1], 2);
		if (dim == 3)
			u_rms[2] += pow(m_u[2][i] - u_average[2], 2);
	}
	// MPI sync
	vector<double> u_rms_mpi(dim);
	dealii::Utilities::MPI::sum(u_rms, MPI_COMM_WORLD, u_rms_mpi);
	u_rms[0] = sqrt(u_rms_mpi[0] / n_mpi);
	u_rms[1] = sqrt(u_rms_mpi[1] / n_mpi);
	if (dim == 3)
		u_rms[2] = sqrt(u_rms_mpi[2] / n_mpi);
}

/*
 const std::string tableFileName = "") :
 m_solver(cfdsolver), m_filename(tableFileName), m_outputOff(
 tableFileName == "")
 */
template<size_t dim>
TurbulenceStats<dim>::TurbulenceStats(CFDSolver<dim> * solver,
		size_t wall_normal_direction,
		const vector<double>& wall_normal_coordinates,
		const std::string table_file_name) :
        m_solver(solver),
        m_filename(table_file_name),
        m_outputOff(table_file_name == ""),
        m_statSize(0),
        m_wallNormalDirection(wall_normal_direction),
        m_wallNormalCoordinates(wall_normal_coordinates) {

	if (not m_outputOff) {
		// create file (if necessary)
		if (m_solver->getIterationStart() > 0) {
			if (is_MPI_rank_0()) {
				m_tableFile = boost::make_shared<std::fstream>(table_file_name,
						std::fstream::out | std::fstream::app);
			}
		} else {
			if (is_MPI_rank_0()) {
				m_tableFile = boost::make_shared<std::fstream>(table_file_name,
						std::fstream::out);
			}
			printHeaderLine();
		}
	}

	double tol = 0.5 * CFDSolverUtilities::getMinimumDoFDistanceGLC<dim>(*m_solver->m_problemDescription->getMesh(),
					m_solver->m_configuration->getSedgOrderOfFiniteElement());
	size_t n_planes = wall_normal_coordinates.size();
	for (size_t i = 0; i < n_planes; i++) {
		StatsInPlane<dim> plane(wall_normal_direction, wall_normal_coordinates.at(i), solver->m_velocity);
		plane.updatePlaneIndices(solver->m_supportPoints, m_solver->m_advectionOperator->getLocallyOwnedDofs(), tol);
		m_planes.push_back(plane);
	}
	m_iterationNumber = -1000;

    // create average vector
    m_uAverage.clear();
    for (size_t i = 0; i < dim; i++) {
        distributed_vector ux_av;
        ux_av.reinit(m_solver->m_velocity.at(i));
        m_uAverage.push_back(ux_av);
    }



	// sync all MPI processes (barrier)
	MPI_sync();
}
template<size_t dim>
TurbulenceStats<dim>::~TurbulenceStats() {
	// TODO Auto-generated destructor stub
}

template<size_t dim>
void TurbulenceStats<dim>::printHeaderLine() {
	assert(not m_outputOff);
	if (is_MPI_rank_0()) {
		string coord;
		if (m_wallNormalDirection == 0) {
			coord = "x";
		}
		if (m_wallNormalDirection == 1) {
			coord = "y";
		}
		if (m_wallNormalDirection == 2) {
			coord = "z";
		}

		(*m_tableFile) << "#  i      t      ";
		for (size_t i = 0; i < m_wallNormalCoordinates.size(); i++) {
			for (size_t j = 0; j < dim; j++) {
				(*m_tableFile) << "u" << j << "_av(" << coord << "="
						<< m_wallNormalCoordinates.at(i) << ")  ";
			}
		}
		for (size_t i = 0; i < m_wallNormalCoordinates.size(); i++) {
			for (size_t j = 0; j < dim; j++) {
				(*m_tableFile) << "u" << j << "_rms(" << coord << "="
						<< m_wallNormalCoordinates.at(i) << ")  ";
			}
		}
		(*m_tableFile) << endl;
	}
}

template<size_t dim>
void TurbulenceStats<dim>::printNewLine() {
	assert(not m_outputOff);
	if (not isUpToDate()) {
		update();
	}
	if (is_MPI_rank_0()) {
		(*m_tableFile) << m_iterationNumber << " " << m_solver->getTime()
				<< " ";
		for (size_t i = 0; i < m_wallNormalCoordinates.size(); i++) {
			for (size_t j = 0; j < dim; j++) {
				(*m_tableFile) << m_averages.at(i).at(j) << "  ";
			}
		}
		for (size_t i = 0; i < m_wallNormalCoordinates.size(); i++) {
			for (size_t j = 0; j < dim; j++) {
				(*m_tableFile) << m_rms.at(i).at(j) << "  ";
			}
		}
		(*m_tableFile) << endl;
	}
}

template<size_t dim>
void TurbulenceStats<dim>::update() {
	m_iterationNumber = m_solver->getIteration();
	size_t n_planes = m_planes.size();
	m_rms.clear();
	m_averages.clear();
	for (size_t i = 0; i < n_planes; i++) {
		vector<double> av_i;
		vector<double> rms_i;
		m_planes.at(i).calculate(av_i, rms_i);
		m_averages.push_back(av_i);
		m_rms.push_back(rms_i);
	}
}

template<size_t dim>
void TurbulenceStats<dim>::addToReynoldsStatistics(
        const vector<distributed_vector>& u) {
    m_uAverage.at(0) *= m_statSize;
    m_uAverage.at(1) *= m_statSize;
    if (3 == dim)
        m_uAverage.at(2) *= m_statSize;
    m_uAverage.at(0).add(u.at(0));
    m_uAverage.at(1).add(u.at(1));
    if (3 == dim)
        m_uAverage.at(2).add(u.at(2));
    m_statSize++;
    m_uAverage.at(0) *= (1.0/m_statSize);
    m_uAverage.at(1) *= (1.0/m_statSize);
    if (3 == dim)
        m_uAverage.at(2) *= (1.0/m_statSize);
}

template<size_t dim>
void TurbulenceStats<dim>::addReynoldsStatisticsToOutput(dealii::DataOut<dim>& data_out) {
    data_out.add_data_vector(m_uAverage.at(0), "ux_average");
    data_out.add_data_vector(m_uAverage.at(1), "uy_average");
    if (3 == dim)
        data_out.add_data_vector(m_uAverage.at(2), "uz_average");

}

// explicit instantiations
template class StatsInPlane<2> ;
template class StatsInPlane<3> ;
template class TurbulenceStats<2> ;
template class TurbulenceStats<3> ;

} /* namespace natrium */
