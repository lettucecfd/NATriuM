/*
 * Checkpoint.cpp
 *
 *  Created on: 10.03.2016
 *      Author: akraem3m
 */

#include "Checkpoint.h"
#include "deal.II/distributed/solution_transfer.h"

namespace natrium {

template<size_t dim>
Checkpoint<dim>::Checkpoint(size_t iteration,
		boost::filesystem::path checkpoint_dir) {
	std::stringstream status_name;
	std::stringstream data_name;
	status_name << "checkpoint_" << iteration << ".stat";
	data_name << "checkpoint_" << iteration << ".data";
	m_statusFile = checkpoint_dir / status_name.str();
	m_dataFile = checkpoint_dir / data_name.str();
}

template<size_t dim>
void Checkpoint<dim>::write(const Mesh<dim>& mesh,
		const DistributionFunctions& f,
		const dealii::DoFHandler<dim>& dof_handler,
		const CheckpointStatus& status) {

	// write status
	if (is_MPI_rank_0()) {
		std::ofstream outfile(m_statusFile.string());
		outfile << status.iterationNumber << endl;
		// time
		outfile << status.time << endl;
		// stencil scaling
		outfile << status.stencilScaling << endl;
		outfile.close();
	} /*if rank 0*/

	// write mesh and solution
	dealii::parallel::distributed::SolutionTransfer < dim, distributed_vector
			> sol_trans(dof_handler);
	std::vector<const distributed_vector*> to_save;
	size_t Q = f.getQ();
	for (size_t i = 0; i < Q; i++) {
		const distributed_vector* p = &f.at(i);
		to_save.push_back(p);
	}
	sol_trans.prepare_serialization(to_save);
	try {
		cout << "Save: " << m_dataFile.string() << endl;
		mesh.save(m_dataFile.c_str());
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while saving the mesh and distribution functions to checkpoint file."
						"That's bad news. Make sure that you have writing permissions and that the disk quota is not exceeded.");
	}

} /* save */

template<size_t dim>
void Checkpoint<dim>::load(DistributionFunctions& f,
		ProblemDescription<dim>& problem, AdvectionOperator<dim>& advection,
		CheckpointStatus& status) {

	dealii::DoFHandler<dim>& dof_handler = *advection.getDoFHandler();
	boost::shared_ptr<Stencil> new_stencil = advection.getStencil();
	Mesh<dim>& mesh = *advection.getMesh();

	// load status
	// read iteration number and time from file
	double phys_time = 0.0;
	size_t iteration_start = 0;
	double stencil_scaling = 0.0;
	if (is_MPI_rank_0()) {
		std::ifstream ifile(m_statusFile.string());
		ifile >> iteration_start;
		ifile >> phys_time;
		ifile >> stencil_scaling;
		ifile.close();
	}
	// transfer iteration start and time to all mpi processes
	status.time =
			dealii::Utilities::MPI::min_max_avg(phys_time, MPI_COMM_WORLD).max;
	status.iterationNumber = dealii::Utilities::MPI::min_max_avg(
			iteration_start,
			MPI_COMM_WORLD).max;
	status.stencilScaling = dealii::Utilities::MPI::min_max_avg(stencil_scaling,
	MPI_COMM_WORLD).max;

	// load mesh and solution
	try {
		// load triangulation (must not be done with refined grid)
		mesh.load(m_dataFile.c_str());
		problem.refineAndTransform();
		// apply refinement
		// distribute dofs
		advection.setupDoFs();
		f.reinit(new_stencil->getQ(), advection.getLocallyOwnedDofs(),
		MPI_COMM_WORLD);
		// read old solution
		dealii::parallel::distributed::SolutionTransfer < dim, distributed_vector
				> sol_trans(dof_handler);
		std::vector<distributed_vector*> to_load;
		for (size_t i = 0; i < new_stencil->getQ(); i++) {
			distributed_vector* p = &f.at(i);
			to_load.push_back(p);
		}
		sol_trans.deserialize(to_load);
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while reading the mesh and distribution functions from checkpoint file: "
						"Please switch off the restart option to start the simulation from the beginning.");
	}

	// transfer to current stencil scaling, if required
	boost::shared_ptr<Stencil> old_stencil = CFDSolverUtilities::make_stencil(
			new_stencil->getD(), new_stencil->getQ(), status.stencilScaling);
	f.transferFromOtherScaling(*old_stencil, *new_stencil,
			dof_handler.locally_owned_dofs());

} /* load */

template<size_t dim>
void Checkpoint<dim>::loadFromDeprecatedCheckpointVersion(Mesh<dim>& mesh,
		DistributionFunctions& f, const dealii::DoFHandler<dim>& dof_handler,
		CheckpointStatus& status, string& directory) {
/*
	//read file
	std::stringstream filename;
	filename << directory << "/checkpoint_status."
			<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
			<< ".dat";
	std::ifstream infile(filename.str().c_str());

	// check if status file exists
	if (not infile) {
		message = "No checkpoint found. Please disable restart option.";
		return false;
	}
	//number of cells
	size_t tmp;
	infile >> tmp;
	if (tmp != m_tria->n_cells()) {
		message = "Number of cells not equal.";
		return false;
	}
	// order of fe
	infile >> tmp;
	if (tmp != m_fe->get_degree()) {
		message = "Order of finite element not equal.";
		return false;
	}
	// number of dofs
	infile >> tmp;
	if (tmp != this->getNumberOfDoFs()) {
		message = "Number of degrees of freedom not equal.";
		return false;
	}
	// D
	infile >> tmp;
	if (tmp != m_stencil->getD()) {
		message = "Dimension not equal.";
		return false;
	}
	// Q
	infile >> tmp;
	if (tmp != m_stencil->getQ()) {
		message = "Number of particle velocities not equal.";
		return false;
	}
	// magic number of cell geometry
	double dtmp;
	infile >> dtmp;
	if (fabs(dtmp - calcMagicNumber()) > 1e-1) {
		message = "Mesh (or at least its magic number) not equal.";
		return false;
	}
	// dqScaling1
	infile >> dtmp;
	if (fabs(dtmp - m_stencil->getDirection(1)(0)) / dtmp > 1e-2) {
		message = "Scaling of Stencil (1st coordinate) not equal.";
		return false;
	}
	// dqScaling2
	infile >> dtmp;
	if (fabs(dtmp - m_stencil->getDirection(1)(1)) / dtmp > 1e-2) {
		message = "Scaling of Stencil (2nd) not equal.";
		return false;
	}
	// fluxType
	infile >> tmp;
	if (tmp != m_useCentralFlux) {
		message = "Flux not equal.";
		return false;
	}
	// advectionType
	string stmp;
	infile >> stmp;
	if (stmp != "SEDGMinLee") {
		message = "AdvectionOperator Type not equal.";
		return false;
	}

// PRECONDITION: vectors already created with the right sizes
// read the distribution functions from file
	try {
		for (size_t i = 0; i < m_stencil->getQ(); i++) {
			// filename
			std::stringstream filename;
			filename << directory << "/checkpoint_f_" << i << "."
					<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
					<< ".dat";
			std::ifstream file(filename.str().c_str());

			// TODO Write and read functions for Trilinos vectors. This here is really bad.
			numeric_vector tmp(f.at(i));
			tmp.block_read(file);
			f.at(i) = tmp;

		}
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while reading the distribution functions from file: Please switch off the restart option to start the simulation from the beginning.");
	}
*/
}
/* load from deprecated checkpoint version */

template<size_t dim>
Checkpoint<dim>::~Checkpoint() {
}

template class Checkpoint<2> ;
template class Checkpoint<3> ;

} /* namespace natrium */
