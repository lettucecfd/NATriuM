/*
 * FinalChannelStatistics.cpp
 *
 *  Created on: 03.03.2016
 *      Author: akraem3m
 */

#include "FinalChannelStatistics.h"

#include "mpi.h"
#include <utility>
#include "deal.II/base/mpi.h"

namespace natrium {

FinalChannelStatistics::FinalChannelStatistics(CFDSolver<3> & solver,
		std::string outdir) :
		DataProcessor<3>(solver), m_outDir(outdir), m_u(solver.getVelocity()), m_rho(
				solver.getDensity()) {

	m_names.push_back("rho");
	m_names.push_back("drho/dx");
	m_names.push_back("drho/dy");
	m_names.push_back("drho/dz");

	m_names.push_back("ux");
	m_names.push_back("uy");
	m_names.push_back("uz");
	m_names.push_back("dux/dx");
	m_names.push_back("dux/dy");
	m_names.push_back("dux/dz");
	m_names.push_back("duy/dx");
	m_names.push_back("duy/dy");
	m_names.push_back("duy/dz");
	m_names.push_back("duz/dx");
	m_names.push_back("duz/dy");
	m_names.push_back("duz/dz");

	m_yCoordsUpToDate = false;

	m_nofObservables = m_names.size();
	m_averages.resize(m_nofObservables);
	m_EX3.resize(m_nofObservables);
	m_EX4.resize(m_nofObservables);
	m_correlations.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		m_correlations.at(i).resize(i + 1);
	}

	m_averages_time.resize(m_nofObservables);
	m_EX3_time.resize(m_nofObservables);
	m_EX4_time.resize(m_nofObservables);
	m_correlations_time.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		m_correlations_time.at(i).resize(i + 1);
	}
	n_steps = 0;
	m_nofCoordinates = 0;
}

void FinalChannelStatistics::update() {
	if (not m_yCoordsUpToDate) {
		updateYValues();
	}
	updateAverages();

}
void FinalChannelStatistics::updateYValues() {
	boost::shared_ptr<AdvectionOperator<3> > advection =
			m_solver.getAdvectionOperator();
	m_yCoordsUpToDate = true;

	//////////////////////////
	// Calculate y values ////
	//////////////////////////
	std::set<double, own_double_less> y_coords;

	const dealii::UpdateFlags update_flags = dealii::update_quadrature_points;
	const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<3> fe_values(advection->getMapping(),
			*(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
	// loop
	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {

			// get y coordinates
			fe_values.reinit(cell);
			const std::vector<dealii::Point<3> >& quad_points =
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
		m_EX3.at(i).resize(m_nofCoordinates);
		m_EX4.at(i).resize(m_nofCoordinates);
		for (size_t j = 0; j <= i; j++) {
			m_correlations.at(i).at(j).resize(m_nofCoordinates);
		}
	}
	m_number.resize(m_nofCoordinates);

	for (size_t i = 0; i < m_nofObservables; i++) {
		m_averages_time.at(i).resize(m_nofCoordinates);
		m_EX3_time.at(i).resize(m_nofCoordinates);
		m_EX4_time.at(i).resize(m_nofCoordinates);
		for (size_t j = 0; j <= i; j++) {
			m_correlations_time.at(i).at(j).resize(m_nofCoordinates);
		}
	}

}

    bool FinalChannelStatistics::isMYCoordsUpToDate() const {
        return m_yCoordsUpToDate;
    }

    void FinalChannelStatistics::updateAverages() {
	boost::shared_ptr<AdvectionOperator<3> > advection =
			m_solver.getAdvectionOperator();

	// prepare local vectors
	vector<size_t> l_number;
	vector<double> l_values;
	vector<vector<double> > l_averages;
	vector<vector<vector<double> > > l_correlations;
	vector<vector<double> > l_EX3; // for skewness
	vector<vector<double> > l_EX4; // for kurtosis
	l_averages.resize(m_nofObservables);
	l_EX3.resize(m_nofObservables);
	l_EX4.resize(m_nofObservables);
	l_correlations.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		l_correlations.at(i).resize(i + 1);
	}
	l_values.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		l_averages.at(i).resize(m_nofCoordinates);
		l_EX3.at(i).resize(m_nofCoordinates);
		l_EX4.at(i).resize(m_nofCoordinates);
		for (size_t j = 0; j <= i; j++) {
			l_correlations.at(i).at(j).resize(m_nofCoordinates);
		}
	}
	l_number.resize(m_nofCoordinates);

	//////////////////////////
	// Calculate averages ////
	//////////////////////////
	const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
			| dealii::update_gradients;
	const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
	dealii::FEValues<3> fe_values(advection->getMapping(),
			*(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
	size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);
	std::vector<dealii::Tensor<1, 3, double> > ux_gradients;
	std::vector<dealii::Tensor<1, 3, double> > uy_gradients;
	std::vector<dealii::Tensor<1, 3, double> > uz_gradients;
	std::vector<dealii::Tensor<1, 3, double> > rho_gradients;
	ux_gradients.resize(dofs_per_cell);
	uy_gradients.resize(dofs_per_cell);
	uz_gradients.resize(dofs_per_cell);
	rho_gradients.resize(dofs_per_cell);
	// loop
	typename dealii::DoFHandler<3>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	double y;
	size_t y_ind;
	size_t dof_ind;
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {

			cell->get_dof_indices(local_indices);

			// get averages
			fe_values.reinit(cell);
			const std::vector<dealii::Point<3> >& quad_points =
					fe_values.get_quadrature_points();

			// calculate gradients (for w and strain rate)
			fe_values.get_function_gradients(m_u.at(0), ux_gradients);
			fe_values.get_function_gradients(m_u.at(1), uy_gradients);
			fe_values.get_function_gradients(m_u.at(2), uz_gradients);
			fe_values.get_function_gradients(m_rho, rho_gradients);

			for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
				y = quad_points.at(i)(1);
				assert(
						m_yCoordinateToIndex.find(y)
								!= m_yCoordinateToIndex.end());
				y_ind = m_yCoordinateToIndex.at(y);
				dof_ind = local_indices.at(i);
				l_number.at(y_ind) += 1;

				// fill value vector
				l_values.at(0) = m_rho(dof_ind);				// rho
				l_values.at(1) = rho_gradients.at(i)[0];		// drho/dx
				l_values.at(2) = rho_gradients.at(i)[1]; 		// drho/dy
				l_values.at(3) = rho_gradients.at(i)[2];		// drho/dz

				l_values.at(4) = m_u.at(0)(dof_ind); 			// ux
				l_values.at(5) = m_u.at(1)(dof_ind);			// uy
				l_values.at(6) = m_u.at(2)(dof_ind);			// uz
				l_values.at(7) = ux_gradients.at(i)[0];			// dux/dx
				l_values.at(8) = ux_gradients.at(i)[1];			// dux/dy
				l_values.at(9) = ux_gradients.at(i)[2];			// dux/dz
				l_values.at(10) = uy_gradients.at(i)[0];		// duy/dx
				l_values.at(11) = uy_gradients.at(i)[1];		// duy/dy
				l_values.at(12) = uy_gradients.at(i)[2];		// duy/dz
				l_values.at(13) = uz_gradients.at(i)[0];		// duz/dx
				l_values.at(14) = uz_gradients.at(i)[1];		// duz/dy
				l_values.at(15) = uz_gradients.at(i)[2];		// duz/dz

				// add to averages:
				for (size_t j = 0; j < m_nofObservables; j++) {
					l_averages.at(j).at(y_ind) += l_values.at(j);
				}
				// add to correlations
				for (size_t j = 0; j < m_nofObservables; j++) {
					for (size_t k = 0; k < j + 1; k++) {
						l_correlations.at(j).at(k).at(y_ind) += (l_values.at(j)
								* l_values.at(k));
					}
				}
				// add to third moment
				for (size_t j = 0; j < m_nofObservables; j++) {
					l_EX3.at(j).at(y_ind) += pow(l_values.at(j), 3);
					l_EX4.at(j).at(y_ind) += pow(l_values.at(j), 4);
				}

			} /* for all quadrature points */
		} /* if locally owned */
	} /* for all cells */

	// communicate
	dealii::Utilities::MPI::sum(l_number, MPI_COMM_WORLD, m_number);
	for (size_t i = 0; i < m_nofObservables; i++) {
		dealii::Utilities::MPI::sum(l_averages.at(i), MPI_COMM_WORLD, m_averages.at(i));
		dealii::Utilities::MPI::sum(l_EX3.at(i), MPI_COMM_WORLD, m_EX3.at(i));
		dealii::Utilities::MPI::sum(l_EX4.at(i), MPI_COMM_WORLD, m_EX4.at(i));
		for (size_t j = 0; j < i + 1; j++) {
			dealii::Utilities::MPI::sum(l_correlations.at(i).at(j),
			MPI_COMM_WORLD, m_correlations.at(i).at(j));
		}
	}

	// divide by number of replicates
	for (size_t i = 0; i < m_nofCoordinates; i++) {
		for (size_t j = 0; j < m_nofObservables; j++) {
			m_averages.at(j).at(i) /= m_number.at(i);
			m_EX3.at(j).at(i) /= m_number.at(i);
			m_EX4.at(j).at(i) /= m_number.at(i);
			for (size_t k = 0; k < j + 1; k++) {
				m_correlations.at(j).at(k).at(i) /= m_number.at(i);
			}
		}
	}
}

void FinalChannelStatistics::write_to_file() {

	if (is_MPI_rank_0()) {
		// averages
		std::stringstream av_name;
		av_name << "statistics_av_i" << m_solver.getIteration() << "_n"
				<< n_steps << ".txt";
		boost::filesystem::path average_file = m_outDir / av_name.str();
		std::ofstream file1(average_file.string(), std::fstream::out);
		file1 << "#";
		for (size_t i = 0; i < m_names.size(); i++) {
			file1 << " " << m_names.at(i);
		}
		file1 << endl;
		for (size_t i = 0; i < m_nofCoordinates; i++) {
			file1 << m_yCoordinates.at(i);
			for (size_t j = 0; j < m_nofObservables; j++) {
				file1 << " " << m_averages_time.at(j).at(i);
			}
			file1 << endl;
		}
		file1.close();

		// EX3
		std::stringstream ex3_name;
		ex3_name << "statistics_EX3_i" << m_solver.getIteration() << "_n"
				<< n_steps << ".txt";
		boost::filesystem::path EX3_file = m_outDir / ex3_name.str();
		std::ofstream file2(EX3_file.string(), std::fstream::out);
		file2 << "#";
		for (size_t i = 0; i < m_names.size(); i++) {
			file2 << " " << m_names.at(i);
		}
		file2 << endl;
		for (size_t i = 0; i < m_nofCoordinates; i++) {
			file2 << m_yCoordinates.at(i);
			for (size_t j = 0; j < m_nofObservables; j++) {
				file2 << " " << m_EX3_time.at(j).at(i);
			}
			file2 << endl;
		}
		file2.close();

		// EX4
		std::stringstream ex4_name;
		ex4_name << "statistics_EX4_i" << m_solver.getIteration() << "_n"
				<< n_steps << ".txt";
		boost::filesystem::path EX4_file = m_outDir / ex4_name.str();
		std::ofstream file3(EX4_file.string(), std::fstream::out);
		file3 << "#";
		for (size_t i = 0; i < m_names.size(); i++) {
			file3 << " " << m_names.at(i);
		}
		file3 << endl;
		for (size_t i = 0; i < m_nofCoordinates; i++) {
			file3 << m_yCoordinates.at(i);
			for (size_t j = 0; j < m_nofObservables; j++) {
				file3 << " " << m_EX4_time.at(j).at(i);
			}
			file3 << endl;
		}
		file3.close();

		// correlation file
		std::stringstream corr_name;
		corr_name << "statistics_corr_i" << m_solver.getIteration() << "_n"
				<< n_steps << ".txt";
		boost::filesystem::path corr_file = m_outDir / corr_name.str();
		std::ofstream file4(corr_file.string(), std::fstream::out);
		file4 << "#";
		for (size_t i = 0; i < m_names.size(); i++) {
			for (size_t j = 0; j < i + 1; j++) {
				file4 << "  <" << m_names.at(i) << "." << m_names.at(j) << ">";
			}
		}
		file4 << endl;
		for (size_t i = 0; i < m_nofCoordinates; i++) {
			file4 << m_yCoordinates.at(i);
			for (size_t j = 0; j < m_nofObservables; j++) {
				for (size_t k = 0; k < j + 1; k++) {
					file4 << " " << m_correlations_time.at(j).at(k).at(i);
				}
			}
			file4 << endl;
		}
		file4.close();

	} /* is mpi rank 0 */
}

void FinalChannelStatistics::apply() {
	if (m_solver.getIteration() % 10 == 0) {
		// add to statistics
		update();
		addToTemporalAverages();
	}
	if (m_solver.getIteration() % 5000 == 0) {
		write_to_file();
	}
}

void FinalChannelStatistics::addToTemporalAverages() {
	for (size_t i = 0; i < m_nofObservables; i++) {
		for (size_t j = 0; j < m_nofCoordinates; j++) {
			m_averages_time.at(i).at(j) *= n_steps;
			m_averages_time.at(i).at(j) += m_averages.at(i).at(j);
			m_averages_time.at(i).at(j) /= (n_steps + 1);

			m_EX3_time.at(i).at(j) *= n_steps;
			m_EX3_time.at(i).at(j) += m_EX3.at(i).at(j);
			m_EX3_time.at(i).at(j) /= (n_steps + 1);

			m_EX4_time.at(i).at(j) *= n_steps;
			m_EX4_time.at(i).at(j) += m_EX4.at(i).at(j);
			m_EX4_time.at(i).at(j) /= (n_steps + 1);
		}
		for (size_t j = 0; j <= i; j++) {
			for (size_t k = 0; k < m_nofCoordinates; k++) {
				m_correlations_time.at(i).at(j).at(k) *= n_steps;
				m_correlations_time.at(i).at(j).at(k) +=
						m_correlations.at(i).at(j).at(k);
				m_correlations_time.at(i).at(j).at(k) /= (n_steps + 1);
			}
		}
	}

	n_steps++;
}

FinalChannelStatistics::~FinalChannelStatistics() {
	// TODO Auto-generated destructor stub
}

} /* namespace natrium */
