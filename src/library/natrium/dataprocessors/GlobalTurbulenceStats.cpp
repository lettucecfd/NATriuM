/**
 * @file PhysicalProperties.cpp
 * @short
 * @date 06.06.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"

#include "GlobalTurbulenceStats.h"

namespace natrium {

template<size_t dim>
GlobalTurbulenceStats<dim>::GlobalTurbulenceStats(const CFDSolver<dim> & solver) :
		DataProcessor<dim>(solver), m_filename(
				outfile(solver.getConfiguration()->getOutputDirectory())), m_legendFilename(
				legendfile(solver.getConfiguration()->getOutputDirectory())), m_outputOff(
				solver.getConfiguration()->isSwitchOutputOff()) {

	// assign names
	m_names.push_back("rho");
	m_names.push_back("drho/dx");
	m_names.push_back("drho/dy");
	if (3 == dim)
		m_names.push_back("drho/dz");

	m_names.push_back("ux");
	m_names.push_back("uy");
	if (3 == dim)
		m_names.push_back("uz");
	m_names.push_back("dux/dx");
	m_names.push_back("dux/dy");
	if (3 == dim)
		m_names.push_back("dux/dz");
	m_names.push_back("duy/dx");
	m_names.push_back("duy/dy");
	if (3 == dim) {
		m_names.push_back("duy/dz");
		m_names.push_back("duz/dx");
		m_names.push_back("duz/dy");
		m_names.push_back("duz/dz");
	}

	// resize data containers
	m_nofObservables = m_names.size();
	m_averages.resize(m_nofObservables);
	m_EX3.resize(m_nofObservables);
	m_EX4.resize(m_nofObservables);
	m_correlations.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		m_correlations.at(i).resize(i + 1);
	}

	// make table file
	if (not m_outputOff) {
		// create file (if necessary)
		if (solver.getIterationStart() > 0) {
			if (is_MPI_rank_0()) {
				m_tableFile = boost::make_shared<std::fstream>(m_filename,
						std::fstream::out | std::fstream::app);
			}
		} else {
			if (is_MPI_rank_0()) {
				m_tableFile = boost::make_shared<std::fstream>(m_filename,
						std::fstream::out);
			}
			printHeaderLine();
		}
	}

}

template<size_t dim>
void GlobalTurbulenceStats<dim>::printHeaderLine() {
	assert(not m_outputOff);
	if (is_MPI_rank_0()) {
		string coord;

		(*m_tableFile) << "# A description of this data is found in "
				<< m_legendFilename << endl;

		std::ofstream legend_file(m_legendFilename, std::fstream::out);

		// write legend file
		legend_file << "1  iteration" << endl;
		legend_file << "3  time" << endl;
		size_t k = 3;
		for (size_t i = 0; i < m_nofObservables; i++) {
			legend_file << k << "  < " << m_names.at(i) << " >" << endl;
			k++;
		}
		for (size_t i = 0; i < m_correlations.size(); i++) {
			for (size_t j = 0; j < m_correlations.at(i).size(); j++) {
				legend_file << k << "  < " << m_names.at(i) << " "
						<< m_names.at(j) << " >" << endl;
				k++;
			}
		}
		for (size_t i = 0; i < m_nofObservables; i++) {
			legend_file << k << "  < " << m_names.at(i) << "^3 >" << endl;
			k++;
		}
		for (size_t i = 0; i < m_nofObservables; i++) {
			legend_file << k << "  < " << m_names.at(i) << "^4 >" << endl;
			k++;
		}
		legend_file << k << "  <energy>" << endl;
		k++;
		legend_file << k << "  <enstrophy>" << endl;
		k++;
		legend_file << k << "  <energy^2>" << endl;
		k++;
		legend_file << k << "  <enstrophy^2>" << endl;
		legend_file.close();

	} /* is_MPI_rank 0 */
}

template<size_t dim>
void GlobalTurbulenceStats<dim>::writeToFile() {
	if ((is_MPI_rank_0()) and (not m_outputOff)) {

		*m_tableFile << this->m_solver.getIteration() << " ";
		*m_tableFile << this->m_solver.getTime() << " ";

		// averages
		for (size_t i = 0; i < m_nofObservables; i++) {
			*m_tableFile << m_averages.at(i) << " ";
		}
		// correlations
		for (size_t i = 0; i < m_nofObservables; i++) {
			for (size_t j = 0; j < m_correlations.at(i).size(); j++) {
				*m_tableFile << m_correlations.at(i).at(j) << " ";
			}
		}
		for (size_t i = 0; i < m_nofObservables; i++) {
			*m_tableFile << m_EX3.at(i) << " ";
		}
		for (size_t i = 0; i < m_nofObservables; i++) {
			*m_tableFile << m_EX4.at(i) << " ";
		}
		*m_tableFile << m_energy << " " ;
		*m_tableFile << m_enstrophy << " " ;
		*m_tableFile << m_energySquared << " " ;
		*m_tableFile << m_enstrophySquared << " " ;
		*m_tableFile << endl;

	} /* is mpi rank 0 */
}

template<size_t dim>
void GlobalTurbulenceStats<dim>::calculate() {

	const vector<distributed_vector> & u(this->m_solver.getVelocity());
	const distributed_vector & rho(this->m_solver.getDensity());
	// prepare local vectors
	vector<double> l_values;
	vector<double> l_averages;
	vector<vector<double> > l_correlations;
	vector<double> l_EX3; // for skewness
	vector<double> l_EX4; // for kurtosis
	l_values.resize(m_nofObservables);
	l_averages.resize(m_nofObservables);
	l_EX3.resize(m_nofObservables);
	l_EX4.resize(m_nofObservables);
	l_correlations.resize(m_nofObservables);
	for (size_t i = 0; i < m_nofObservables; i++) {
		l_correlations.at(i).resize(i + 1);
	}

	//////////////////////////
	// Calculate averages ////
	//////////////////////////
	const dealii::UpdateFlags update_flags = dealii::update_gradients
			| dealii::update_JxW_values;
	const dealii::DoFHandler<dim> & dof_handler =
			*(this->m_solver.getAdvectionOperator()->getDoFHandler());
	dealii::FEValues<dim> fe_values(
			this->m_solver.getAdvectionOperator()->getMapping(),
			*(this->m_solver.getAdvectionOperator()->getFe()),
			*(this->m_solver.getAdvectionOperator()->getQuadrature()), update_flags);
	size_t dofs_per_cell =
			this->m_solver.getAdvectionOperator()->getFe()->dofs_per_cell;
	std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);
	std::vector<dealii::Tensor<1, dim, double> > ux_gradients;
	std::vector<dealii::Tensor<1, dim, double> > uy_gradients;
	std::vector<dealii::Tensor<1, dim, double> > uz_gradients;
	std::vector<dealii::Tensor<1, dim, double> > rho_gradients;
	double ener = 0.0;
	double enst = 0.0;
	double ener_sq = 0.0;
	double enst_sq = 0.0;
	ux_gradients.resize(dofs_per_cell);
	uy_gradients.resize(dofs_per_cell);
	uz_gradients.resize(dofs_per_cell);
	rho_gradients.resize(dofs_per_cell);
	// loop
	typename dealii::DoFHandler<dim>::active_cell_iterator cell =
			dof_handler.begin_active(), endc = dof_handler.end();
	size_t dof_ind;
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {

			cell->get_dof_indices(local_indices);

			// get averages
			fe_values.reinit(cell);
			const std::vector<double>& weights = fe_values.get_JxW_values();

			// calculate gradients (for w and strain rate)
			fe_values.get_function_gradients(u.at(0), ux_gradients);
			fe_values.get_function_gradients(u.at(1), uy_gradients);
			if (3 == dim)
				fe_values.get_function_gradients(u.at(2), uz_gradients);
			fe_values.get_function_gradients(rho, rho_gradients);

			for (size_t i = 0; i < dofs_per_cell; i++) {
				size_t q =
						this->m_solver.getAdvectionOperator()->getCelldofToQIndex().at(
								i);
				dof_ind = local_indices.at(i);

				// fill value vector
				size_t k = 0;
				l_values.at(k) = rho(dof_ind);				// rho
				k++;
				l_values.at(k) = rho_gradients.at(i)[0];	// drho/dx
				k++;
				l_values.at(k) = rho_gradients.at(i)[1]; 	// drho/dy
				k++;
				if (3 == dim) {
					l_values.at(k) = rho_gradients.at(i)[2]; 	// drho/dz
					k++;
				}
				l_values.at(k) = u.at(0)(dof_ind); 			// ux
				k++;
				l_values.at(k) = u.at(1)(dof_ind);			// uy
				k++;
				if (3 == dim) {
					l_values.at(k) = u.at(2)(dof_ind);		// uz
					k++;
				}
				l_values.at(k) = ux_gradients.at(i)[0];		// dux/dx
				k++;
				l_values.at(k) = ux_gradients.at(i)[1];		// dux/dy
				k++;
				if (3 == dim) {
					l_values.at(k) = ux_gradients.at(i)[2];	// dux/dz
					k++;
				}
				l_values.at(k) = uy_gradients.at(i)[0];		// duy/dx
				k++;
				l_values.at(k) = uy_gradients.at(i)[1];		// duy/dy
				k++;
				if (3 == dim) {
					l_values.at(k) = uy_gradients.at(i)[2];	// duy/dz
					k++;
					l_values.at(k) = uz_gradients.at(i)[0];	// duz/dx
					k++;
					l_values.at(k) = uz_gradients.at(i)[1];	// duz/dy
					k++;
					l_values.at(k) = uz_gradients.at(i)[2];	// duz/dz
				}

				// add to averages:
				for (size_t j = 0; j < m_nofObservables; j++) {
					l_averages.at(j) += l_values.at(j) * weights.at(q);
				}

				// add to correlations
				for (size_t j = 0; j < m_nofObservables; j++) {
					for (size_t r = 0; r < j + 1; r++) {
						l_correlations.at(j).at(r) += (l_values.at(j)
								* l_values.at(r)) * weights.at(q);
					}
				}
				// add to third moment
				for (size_t j = 0; j < m_nofObservables; j++) {
					l_EX3.at(j) += pow(l_values.at(j), 3) * weights.at(q);
					l_EX4.at(j) += pow(l_values.at(j), 4) * weights.at(q);
				}
				// add to energies and enstrophies
				double e1 = u.at(0)(dof_ind) * u.at(0)(dof_ind) + u.at(1)(dof_ind) * u.at(1)(dof_ind);
				double e2 = pow(uy_gradients.at(i)[0] - ux_gradients.at(i)[1],2);
				if (3==dim){
					e1 += u.at(2)(dof_ind) * u.at(2)(dof_ind);
					e2 += pow(uz_gradients.at(i)[1] - uy_gradients.at(i)[2],2);
					e2 += pow(ux_gradients.at(i)[2] - uz_gradients.at(i)[0],2);
				}
				ener += e1 * weights.at(q);
				ener_sq += e1 * e1 * weights.at(q);
				enst += e2 * weights.at(q);
				enst_sq += e2*e2*weights.at(q);

			} /* for all quadrature points */
		} /* if locally owned */
	} /* for all cells */

	// communicate
	//for (size_t i = 0; i < m_nofObservables; i++) {
	dealii::Utilities::MPI::sum(l_averages, MPI_COMM_WORLD, m_averages);
	dealii::Utilities::MPI::sum(l_EX3, MPI_COMM_WORLD, m_EX3);
	dealii::Utilities::MPI::sum(l_EX4, MPI_COMM_WORLD, m_EX4);
	for (size_t i = 0; i < m_nofObservables; i++) {
		dealii::Utilities::MPI::sum(l_correlations.at(i),
		MPI_COMM_WORLD, m_correlations.at(i));
	}
	dealii::Utilities::MPI::sum(ener, MPI_COMM_WORLD, m_energy);
	dealii::Utilities::MPI::sum(enst, MPI_COMM_WORLD, m_enstrophy);
	dealii::Utilities::MPI::sum(ener_sq, MPI_COMM_WORLD, m_energySquared);
	dealii::Utilities::MPI::sum(enst_sq, MPI_COMM_WORLD, m_enstrophySquared);
	//}
}

template<size_t dim>
void GlobalTurbulenceStats<dim>::apply() {
	{
		if ((this->m_solver.getIteration()
				% this->m_solver.getConfiguration()->getOutputTableInterval() == 0)
				and (not m_outputOff)) {
			calculate();
			writeToFile();
		}
	}
}

// Explicit instantiation
template class GlobalTurbulenceStats<2> ;
template class GlobalTurbulenceStats<3> ;

} /* namespace natrium */
