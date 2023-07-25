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

	double tol = 0.5 * CFDSolverUtilities::getMinimumDoFDistanceGLL<dim>(
					*m_solver->m_problemDescription->getMesh(),
					m_solver->m_configuration->getSedgOrderOfFiniteElement());
	size_t n_planes = wall_normal_coordinates.size();
	for (size_t i = 0; i < n_planes; i++) {
		StatsInPlane<dim> plane(wall_normal_direction,
				wall_normal_coordinates.at(i), solver->m_velocity);
		plane.updatePlaneIndices(solver->m_supportPoints,
				m_solver->m_advectionOperator->getLocallyOwnedDofs(), tol);
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

template<size_t dim>
ShearLayerStats<dim>::ShearLayerStats(CFDSolver<dim> * solver, const std::string table_file_name) :
        m_solver(solver),
        m_filename(table_file_name),
        m_outputOff(table_file_name == ""),
        m_statSize(0) {

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
    m_iterationNumber = -1000;

    updateYValues();

    // sync all MPI processes (barrier)
    MPI_sync();
}
template<size_t dim>
ShearLayerStats<dim>::~ShearLayerStats() {
    // TODO Auto-generated destructor stub
}

template<size_t dim>
void ShearLayerStats<dim>::updateYValues() {
    boost::shared_ptr<AdvectionOperator<dim> > advection = m_solver->getAdvectionOperator();
    m_yCoordsUpToDate = true;

    //////////////////////////
    // Calculate y values ////
    //////////////////////////
    std::set<double, own_double_less> y_coords;

    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points;
    const dealii::DoFHandler<dim> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<dim> fe_values(advection->getMapping(),
                                  *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    // loop
    typename dealii::DoFHandler<dim>::active_cell_iterator cell =
            dof_handler.begin_active(), endc = dof_handler.end();
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {

            // get y coordinates
            fe_values.reinit(cell);
            const std::vector<dealii::Point<dim> >& quad_points =
                    fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y_coords.insert(quad_points.at(i)(1));
            }
        }
    }

    //communicate list of y values
    size_t n_ycoords = y_coords.size();
    size_t max_nycoords = dealii::Utilities::MPI::min_max_avg(n_ycoords,
                                                              MPI_COMM_WORLD).max;
    double *sendbuf = (double*) malloc(max_nycoords * sizeof(double));

    // fill send buffer with y coordinates
    size_t i = 0;
    for (std::set<double>::iterator it = y_coords.begin(); it != y_coords.end();
         ++it) {
        sendbuf[i] = *it;
        i++;
    }
    for (; i < max_nycoords; i++) {
        sendbuf[i] = *y_coords.begin();
    }
    // allocate read buffer
    double *recvbuf = (double*) malloc(dealii::Utilities::MPI::n_mpi_processes(
            MPI_COMM_WORLD) * max_nycoords * sizeof(double));
    // allgather: distribute all ycoords to all processors
    MPI_Allgather(sendbuf, max_nycoords,
                  MPI_DOUBLE, recvbuf, max_nycoords,
                  MPI_DOUBLE,
                  MPI_COMM_WORLD);

    // remove double entries
    std::set<double> y_coords_gathered;
    for (size_t i = 0;
         i < max_nycoords * dealii::Utilities::MPI::n_mpi_processes(
                 MPI_COMM_WORLD); i++) {
        y_coords_gathered.insert(recvbuf[i]);
    }
    // fill member variables
    i = 0;
    for (std::set<double>::iterator it = y_coords_gathered.begin();
         it != y_coords_gathered.end(); ++it) {
        m_yCoordinates.push_back(*it);
        m_yCoordinateToIndex.insert(std::make_pair(*it, i));
        i++;
    }
    m_nofCoordinates = i;
    // free
    free(sendbuf);
    free(recvbuf);

    // resize all other vectors
    for (size_t i = 0; i < m_nofObservables; i++) {
        m_averages.at(i).resize(m_nofCoordinates);
//        m_EX3.at(i).resize(m_nofCoordinates);
//        m_EX4.at(i).resize(m_nofCoordinates);
//        for (size_t j = 0; j <= i; j++) {
//            m_correlations.at(i).at(j).resize(m_nofCoordinates);
//        }
    }
    m_number.resize(m_nofCoordinates);

//    // temporal averages
//    for (size_t i = 0; i < m_nofObservables; i++) {
//        m_averages_time.at(i).resize(m_nofCoordinates);
//        m_EX3_time.at(i).resize(m_nofCoordinates);
//        m_EX4_time.at(i).resize(m_nofCoordinates);
//        for (size_t j = 0; j <= i; j++) {
//            m_correlations_time.at(i).at(j).resize(m_nofCoordinates);
//        }
//    }

}

template<size_t dim>
void ShearLayerStats<dim>::printHeaderLine() {
    assert(not m_outputOff);
    if (is_MPI_rank_0()) {
        (*m_tableFile) << "#  i      t      deltaTheta" << endl;
    }
}

template<size_t dim>
void ShearLayerStats<dim>::printNewLine() {
    assert(not m_outputOff);
    if (not isUpToDate()) {
        update();
    }
    if (is_MPI_rank_0()) {
        (*m_tableFile) << m_iterationNumber << " " << m_solver->getTime()
                       << " " << m_lastDeltaTheta;
        (*m_tableFile) << endl;
    }
}

template<size_t dim>
void ShearLayerStats<dim>::update() {
    m_iterationNumber = m_solver->getIteration();
    m_lastDeltaTheta =
//    size_t n_planes = m_planes.size();
//    m_rms.clear();
//    m_averages.clear();
//    for (size_t i = 0; i < n_planes; i++) {
//        vector<double> av_i;
//        vector<double> rms_i;
//        m_planes.at(i).calculate(av_i, rms_i);
//        m_averages.push_back(av_i);
//        m_rms.push_back(rms_i);
//    }
}

template<size_t dim>
void ShearLayerStats<dim>::addToReynoldsAveragesXZ(
        const vector<distributed_vector>& u,
        const distributed_vector &rho) {
    //////////////////////////
    // Calculate averages ////
    //////////////////////////

    vector<double> rhoux_average;
    vector<double> rho_average;
    vector<double> ux_favre;
    vector<double> integrand;
    rhoux_average.resize(m_nofCoordinates);
    ux_favre.resize(m_nofCoordinates);
    rho_average.resize(m_nofCoordinates);
    integrand.resize(m_nofCoordinates);

    vector<size_t> number;
    number.resize(m_nofCoordinates);
    boost::shared_ptr<AdvectionOperator<dim> > advection = m_solver->getAdvectionOperator();

    const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
                                             | dealii::update_gradients;
    const dealii::DoFHandler<dim> & dof_handler = *(advection->getDoFHandler());
    dealii::FEValues<dim> fe_values(advection->getMapping(),
                                    *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
    size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
    std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);

    // loop
    typename dealii::DoFHandler<dim>::active_cell_iterator cell =
            dof_handler.begin_active(), endc = dof_handler.end();
    double y;
    size_t y_ind;
    size_t dof_ind;
    for (; cell != endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_indices);
            // get averages
            fe_values.reinit(cell);
            const std::vector<dealii::Point<dim> >& quad_points =
                    fe_values.get_quadrature_points();
            for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                y = quad_points.at(i)(1);
                assert(m_yCoordinateToIndex.find(y) != m_yCoordinateToIndex.end());
                y_ind = m_yCoordinateToIndex.at(y);

                dof_ind = local_indices.at(i);
                number.at(y_ind) += 1;
                // fill value vector
                // add to averages:
                rhoux_average.at(y_ind) += rho(dof_ind) * u.at(0)(dof_ind);				// rho u
                rho_average.at(y_ind) += rho(dof_ind);
            }
        } /* if is locally owned */
    }

    // communicate
    for (size_t i = 0; i < m_nofCoordinates; i++) {
        number.at(i) = dealii::Utilities::MPI::sum(number.at(i), MPI_COMM_WORLD);
        ux_favre.at(i) = rhoux_average.at(i) / rho_average.at(i);
        integrand.at(i) = rho_average.at(i) * (1 /* dU / 2 */ - ux_favre.at(i) * (1 /* dU / 2 */ + ux_favre.at(i)));
        ux_favre.at(i) = dealii::Utilities::MPI::sum(ux_favre.at(i), MPI_COMM_WORLD);
        rhoux_average.at(i) = dealii::Utilities::MPI::sum(rhoux_average.at(i), MPI_COMM_WORLD);
        rho_average.at(i) = dealii::Utilities::MPI::sum(rho_average.at(i), MPI_COMM_WORLD);
        ux_favre.at(i) /= number.at(i);
        rhoux_average.at(i) /= number.at(i);
        rho_average.at(i) /= number.at(i);
    }

    // integrate along y (?)
    double integral = 0;
    double rho_int = 0;
    double rhoux_int = 0;
    double ux_favre_int = 0;
    double interval_length = 0;
    for (size_t i = 0; i < m_nofCoordinates-1; i++) {
        double window_size = std::abs( m_yCoordinates.at(i+1) -m_yCoordinates.at(i));
        integral += window_size*0.5*(integrand.at(i)+integrand.at(i+1));
        rho_int += window_size*0.5*(rho_average.at(i)+rho_average.at(i+1));
        rhoux_int += window_size*0.5*(rhoux_average.at(i)+rhoux_average.at(i+1));
        ux_favre_int += window_size*0.5*(ux_favre.at(i)+ux_favre.at(i+1));
        interval_length += window_size;
    }
    m_Rho.add(rho_int / interval_length);
    m_RhoUx.add(rhoux_int / interval_length);
    m_UxFavre.add(ux_favre_int / interval_length);
    m_DeltaTheta.add(integral * (1. / 4. /* rho0 * dU^2 = 1 * 2*2 */));
}

template<size_t dim>
void ShearLayerStats<dim>::addReynoldsAveragesXZToOutput(dealii::DataOut<dim>& data_out) {
    data_out.add_data_vector(m_Rho, "m_Rho");
    data_out.add_data_vector(m_RhoUx, "m_RhoUx");
    data_out.add_data_vector(m_UxFavre, "m_UxFavre");
    data_out.add_data_vector(m_DeltaTheta, "m_DeltaTheta");
}

// explicit instantiations
template class StatsInPlane<2> ;
template class StatsInPlane<3> ;
template class TurbulenceStats<2> ;
template class TurbulenceStats<3> ;
template class ShearLayerStats<3> ;

} /* namespace natrium */
