/*
 * Checkpoint.cpp
 *
 *  Created on: 10.03.2016
 *      Author: akraem3m
 */

#include "Checkpoint.h"

#include "math.h"
#include <iostream>

#include "boost/filesystem.hpp"
#include "boost/algorithm/string/replace.hpp"

#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/distributed/solution_transfer.h"
#include "deal.II/numerics/vector_tools.h"

#include "../utilities/CFDSolverUtilities.h"


namespace natrium {

template<size_t dim>
Checkpoint<dim>::Checkpoint(size_t iteration,
		boost::filesystem::path checkpoint_dir) {

	size_t it = iteration;

	if (iteration == 1) {
		// find most recent checkpoint automatically
		for (auto& entry : boost::make_iterator_range(
				boost::filesystem::directory_iterator(checkpoint_dir), { })) {
			std::string this_filename = entry.path().filename().string();
			if (boost::filesystem::extension(this_filename)
					!= ".stat") {
				continue;
			}
			boost::replace_all(this_filename, "checkpoint_", "");
			boost::replace_all(this_filename, ".stat", "");

			size_t i = std::stoi(this_filename);
			if (i > it){
				it = i;
			}
		}
	}
	LOG(BASIC) << "Reading checkpoint " << it << " auomatically." <<  endl;
	std::stringstream status_name;
	std::stringstream data_name;
	status_name << "checkpoint_" << it << ".stat";
	data_name << "checkpoint_" << it << ".data";
	m_statusFile = checkpoint_dir / status_name.str();
	m_dataFile = checkpoint_dir / data_name.str();

}

template<size_t dim>
void Checkpoint<dim>::write(const Mesh<dim>& mesh, DistributionFunctions& f,
		const dealii::DoFHandler<dim>& dof_handler,
		const CheckpointStatus& status) {

	f.updateGhosted();
	// write status
	if (is_MPI_rank_0()) {
		std::ofstream outfile(m_statusFile.string());
		outfile << status.iterationNumber << endl;
		// time
		outfile << status.time << endl;
		// stencil scaling
		outfile << status.stencilScaling << endl;
		outfile << status.feOrder << endl;
		outfile.close();
	} /*if rank 0*/

	// write mesh and solution
	dealii::parallel::distributed::SolutionTransfer<dim, distributed_vector> sol_trans(
			dof_handler);
	std::vector<const distributed_vector*> to_save;
	size_t Q = f.getQ();
	for (size_t i = 0; i < Q; i++) {
		const distributed_vector* p = &(f.atGhosted(i));
		to_save.push_back(p);
	}
	sol_trans.prepare_for_serialization(to_save);
	try {
		LOG(BASIC) << "Save: " << m_dataFile.string() << endl;
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
	size_t fe_order = 0;
	if (is_MPI_rank_0()) {
		std::ifstream ifile(m_statusFile.string());
		ifile >> iteration_start;
		ifile >> phys_time;
		ifile >> stencil_scaling;
		ifile >> fe_order;
		ifile.close();
	}
	MPI_sync();

	// container for locally relevant dofs
	dealii::IndexSet locally_relevant_dofs;

	// transfer iteration start and time to all mpi processes
	status.time =
			dealii::Utilities::MPI::min_max_avg(phys_time, MPI_COMM_WORLD).max;
	status.iterationNumber = dealii::Utilities::MPI::min_max_avg(
			iteration_start,
			MPI_COMM_WORLD).max;
	status.stencilScaling = dealii::Utilities::MPI::min_max_avg(stencil_scaling,
	MPI_COMM_WORLD).max;
	status.feOrder = dealii::Utilities::MPI::min_max_avg(fe_order,
	MPI_COMM_WORLD).max;

	if (status.feOrder != advection.getOrderOfFiniteElement()) {
		std::stringstream msg;
		msg
				<< "You are not allowed to restart the simulation with other order of finite element."
						"Previously: " << status.feOrder << "; Now: "
				<< advection.getOrderOfFiniteElement() << endl;
		natrium_errorexit(msg.str().c_str());
	}

	// load mesh and solution
	try {
		// copy triangulation
		// create future mesh just to get difference in refinement level
		Mesh<dim> future_mesh(MPI_COMM_WORLD);
		// copy mesh
		future_mesh.copy_triangulation(mesh);
		// Refine and transform tmp mesh to get the desired refinement level
		problem.refineAndTransform(future_mesh);
		size_t nlevels_new = future_mesh.n_global_levels();

		LOG(DETAILED) << "Read old solution" << endl;
		// Prepare read old solution
		// load mesh (must not be done with refined grid)
		mesh.load(m_dataFile.c_str());
		// on calling load(), the old mesh has been refined, as before saving

		// setup dofs on old mesh
		dealii::parallel::distributed::SolutionTransfer<dim, distributed_vector> sol_trans(
				dof_handler);
		dof_handler.distribute_dofs(*advection.getFe());
		DistributionFunctions tmp_f(f);
		dealii::DoFTools::extract_locally_relevant_dofs(dof_handler,
				locally_relevant_dofs);
		tmp_f.reinit(new_stencil->getQ(), dof_handler.locally_owned_dofs(),
				locally_relevant_dofs,
				MPI_COMM_WORLD, advection.isDG());

		// read old solution
		std::vector<distributed_vector*> all_read;
		for (size_t i = 0; i < new_stencil->getQ(); i++) {
			distributed_vector* ptr = &tmp_f.at(i);
			all_read.push_back(ptr);
		}
		sol_trans.deserialize(all_read);
		tmp_f.updateGhosted();

		LOG(DETAILED) << "Interpolate to refined grid" << endl;
		LOG(DETAILED) << "... from refinement level "
				<< mesh.n_global_levels() - 1 << " to " << nlevels_new - 1
				<< endl;

		// assumption: future_mesh is a globally refined version of mesh
		if (mesh.n_global_levels() > nlevels_new) {
			throw CheckpointException(
					"Restarting from coarser grid is not implemented, yet.");
		}
		if (mesh.n_global_levels() == nlevels_new) {
			advection.setupDoFs();
			f.reinit(new_stencil->getQ(), dof_handler.locally_owned_dofs(),
					locally_relevant_dofs,
					MPI_COMM_WORLD, advection.isDG());
			f = tmp_f;
		}
		while (mesh.n_global_levels() < nlevels_new) {
			// do one refinemenent step
			dealii::parallel::distributed::SolutionTransfer<dim,
					distributed_vector> soltrans_refine(dof_handler);
			mesh.set_all_refine_flags();
			mesh.prepare_coarsening_and_refinement();
			// prepare all in
			std::vector<const distributed_vector*> all_in;
			all_in.clear();
			for (size_t i = 0; i < new_stencil->getQ(); i++) {
				const distributed_vector* ptr = &tmp_f.atGhosted(i);
				all_in.push_back(ptr);
			}
			soltrans_refine.prepare_for_coarsening_and_refinement(all_in);
			mesh.execute_coarsening_and_refinement();
			advection.setupDoFs();
			// after refinement, locally relevant dofs have changed
			dealii::DoFTools::extract_locally_relevant_dofs(dof_handler,
					locally_relevant_dofs);
			f.reinit(new_stencil->getQ(), dof_handler.locally_owned_dofs(),
					locally_relevant_dofs,
					MPI_COMM_WORLD, advection.isDG());
			// write f pointers into std::vector
			std::vector<distributed_vector*> all_out;
			all_out.clear();
			for (size_t i = 0; i < new_stencil->getQ(); i++) {
				distributed_vector* p = &f.at(i);
				all_out.push_back(p);
			}
			// interpolate
			soltrans_refine.interpolate(all_out);
			f.updateGhosted();
			// prepare next cycle
			tmp_f.reinit(new_stencil->getQ(), dof_handler.locally_owned_dofs(),
					locally_relevant_dofs,
					MPI_COMM_WORLD, advection.isDG());
			tmp_f = f;
		}

		// transform mesh
		problem.transform(mesh);

		// clear old dof handler to enable deletion of automatic variable old_fe
		//old_dof_handler.clear();
	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while reading the mesh and distribution functions from checkpoint file: "
						"Please switch off the restart option to start the simulation from the beginning.");
	} catch (dealii::StandardExceptions::ExcMessage& excM) {
		std::stringstream msg;
		msg << "Deal error while restarting:" << excM.what()
				<< "If the error has to do with parallel partitioning, try to restart with a refinement and number of processes"
						"that does not differ much from the previous simulation.";
		natrium_errorexit(msg.str().c_str());

	}

	// TODO Enable transfer to new scaling
	LOG(DETAILED) << "Transfer to new scaling" << endl;
	// transfer to current stencil scaling, if required
	boost::shared_ptr<Stencil> old_stencil = CFDSolverUtilities::make_stencil(
			new_stencil->getD(), new_stencil->getQ(), status.stencilScaling);
	f.transferFromOtherScaling(*old_stencil, *new_stencil,
			dof_handler.locally_owned_dofs());
	f.compress(dealii::VectorOperation::insert);

	LOG(DETAILED) << "Restart successful" << endl;

} /* load */

template<size_t dim>
void Checkpoint<dim>::loadFromDeprecatedCheckpointVersion(
		DistributionFunctions& f, AdvectionOperator<dim>& advection,
		string directory, CheckpointStatus& status) {

	//dealii::DoFHandler<dim>& dof_handler = *advection.getDoFHandler();
	boost::shared_ptr<Stencil> new_stencil = advection.getStencil();
	Mesh<dim>& mesh = *advection.getMesh();
	std::string message;

	//read file
	std::stringstream filename;
	filename << directory << "/checkpoint_status."
			<< dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
			<< ".dat";
	std::ifstream infile(filename.str().c_str());

	bool status_ok = true;
	// while loop is executed only once; only used here to support "break;"
	while (true) {
		// check if status file exists
		if (not infile) {
			message = "No checkpoint found. Please disable restart option.";
			status_ok = false;
			break;
		}
		//number of cells
		size_t tmp;
		infile >> tmp;
		if (tmp != mesh.n_cells()) {
			message = "Number of cells not equal.";
			status_ok = false;
			break;
		}
		// order of fe
		infile >> tmp;
		if (tmp != advection.getFe()->degree) {
			message = "Order of finite element not equal.";
			status_ok = false;
			break;
		}
		// number of dofs
		infile >> tmp;
		if (tmp != advection.getNumberOfDoFs()) {
			message = "Number of degrees of freedom not equal.";
			status_ok = false;
			break;
		}
		// D
		infile >> tmp;
		if (tmp != new_stencil->getD()) {
			message = "Dimension not equal.";
			status_ok = false;
			break;
		}
		// Q
		infile >> tmp;
		if (tmp != new_stencil->getQ()) {
			message = "Number of particle velocities not equal.";
			status_ok = false;
			break;
		}
		// magic number of cell geometry
		double dtmp;
		infile >> dtmp;
		if (fabs(dtmp - 0.0) > 1e-1) {
			message = "Mesh (or at least its magic number) not equal.";
			status_ok = false;
			break;
		}
		// dqScaling1
		infile >> dtmp;
		if (std::fabs(dtmp - new_stencil->getDirection(1)(0)) / dtmp > 1e-2) {
			message = "Scaling of Stencil (1st coordinate) not equal.";
			status_ok = false;
			break;
		}
		// dqScaling2
		infile >> dtmp;
		if (std::fabs(dtmp - new_stencil->getDirection(1)(1)) / dtmp > 1e-2) {
			message = "Scaling of Stencil (2nd) not equal.";
			status_ok = false;
			break;
		}
		// fluxType
		infile >> tmp;
		if (tmp != 0) {
			message = "Flux not equal.";
			status_ok = false;
			break;
		}
		// advectionType
		string stmp;
		infile >> stmp;
		if (stmp != "SEDGMinLee") {
			message = "AdvectionOperator Type not equal.";
			status_ok = false;
			break;
		}
		break;
	}
	if (!status_ok) {
		throw CheckpointException(message);
	}

	// load status
	// read iteration number and time from file
	boost::filesystem::path odir(directory);
	boost::filesystem::path datfile = odir / "checkpoint.dat";
	double phys_time = 0.0;
	size_t iteration_start = 0;
	if (is_MPI_rank_0()) {
		std::ifstream ifile(datfile.string());
		ifile >> iteration_start;
		ifile >> phys_time;
		ifile.close();
	}
	// transfer iteration start and time to all mpi processes
	status.time = dealii::Utilities::MPI::min_max_avg(phys_time,
	MPI_COMM_WORLD).max;
	status.iterationNumber = dealii::Utilities::MPI::min_max_avg(
			iteration_start, MPI_COMM_WORLD).max;
	status.stencilScaling = new_stencil->getScaling();
	status.feOrder = advection.getOrderOfFiniteElement();

// PRECONDITION: vectors already created with the right sizes
// read the distribution functions from file
	try {
		for (size_t i = 0; i < new_stencil->getQ(); i++) {
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

}
/* load from deprecated checkpoint version */

template<size_t dim>
Checkpoint<dim>::~Checkpoint() {
}

template class Checkpoint<2> ;
template class Checkpoint<3> ;

} /* namespace natrium */
