/*
 * Checkpoint.cpp
 *
 *  Created on: 10.03.2016
 *      Author: akraem3m
 */

#include "Checkpoint.h"

#include "math.h"

#include "deal.II/fe/fe_dgq.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/base/quadrature_lib.h"
#include "deal.II/distributed/solution_transfer.h"
#include "deal.II/numerics/vector_tools.h"

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
		outfile << status.feOrder << endl;
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

	// load mesh and solution
	try {
		// copy triangulation
		Mesh<dim> old_mesh(MPI_COMM_WORLD);
		// copy mesh
		old_mesh.copy_triangulation(mesh);

		LOG(DETAILED) << "Read old solution" << endl;
		// Prepare read old solution
		// load mesh (must not be done with refined grid)
		old_mesh.load(m_dataFile.c_str());
		// transform old mesh
		problem.transform(old_mesh);
		cout << "Make old dof handler" << endl;
		// on calling load(), the old mesh has been refined, as before saving
		dealii::DoFHandler<dim> old_dof_handler(old_mesh);
		cout << "Make solution transfer" << endl;
		dealii::parallel::distributed::SolutionTransfer < dim, distributed_vector
				> sol_trans(old_dof_handler);
		cout << "Distribute dofs" << endl;
		old_dof_handler.distribute_dofs(*advection.getFe());
		cout << "Make distribution functions" << endl;
		DistributionFunctions old_f;
		cout << "Reinit" << endl;
		old_f.reinit(new_stencil->getQ(), old_dof_handler.locally_owned_dofs(),
		MPI_COMM_WORLD);
		// read old solution
		std::vector<distributed_vector*> to_load;
		for (size_t i = 0; i < new_stencil->getQ(); i++) {
			distributed_vector* p = &old_f.at(i);
			to_load.push_back(p);
		}
		cout << "deserialize" << endl;
		sol_trans.deserialize(to_load);

		// Refine and transform new mesh
		cout << "Refine" << endl;
		problem.refineAndTransform();
		LOG(DETAILED) << "Interpolate to new grid" << endl;
		// dof handler with old fe on new mesh
		/*cout << mesh.n_levels() << " " << old_mesh.n_levels() << endl;
		cout << mesh.n_cells(0) << " " << old_mesh.n_cells(0) << endl;
		typename Mesh<dim>::cell_iterator cell_1 = mesh.begin(0), cell_2 =
				old_mesh.begin(0), endc = mesh.end(0);
		for (; cell_1 != endc; ++cell_1, ++cell_2)
			for (unsigned int v = 0;
					v < dealii::GeometryInfo<dim>::vertices_per_cell; ++v)
				if (cell_1->vertex(v) != cell_2->vertex(v))
					cout << "vertex " << v << "does not agree on cells" << cell_1->id() << ", " << cell_2->id();
		*/
		//dealii::DoFHandler<dim> dof_handler_old_fe(mesh);
		advection.setupDoFs();
		f.reinit(new_stencil->getQ(), dof_handler.locally_owned_dofs(),
		MPI_COMM_WORLD);
		// interpolate from old to new mesh
		for (size_t i = 0; i < new_stencil->getQ(); i++) {
			cout << "interpolate to new mesh " << i << endl;
			dealii::VectorTools::interpolate_to_different_mesh(old_dof_handler,
					old_f.at(i), dof_handler, f.at(i));
		}

		// clear old dof handler to enable deletion of automatic variable old_fe
		old_dof_handler.clear();

	} catch (dealii::StandardExceptions::ExcIO& excIO) {
		natrium_errorexit(
				"An error occurred while reading the mesh and distribution functions from checkpoint file: "
						"Please switch off the restart option to start the simulation from the beginning.");
	}


	// TODO Enable transfer to new scaling
	LOG(DETAILED) << "Transfer to new scaling" << endl;
	// transfer to current stencil scaling, if required
	boost::shared_ptr<Stencil> old_stencil = CFDSolverUtilities::make_stencil(
			new_stencil->getD(), new_stencil->getQ(), status.stencilScaling);
	f.transferFromOtherScaling(*old_stencil, *new_stencil,
			dof_handler.locally_owned_dofs());

	LOG(DETAILED) << "Restart successful" << endl;

} /* load */

template<size_t dim>
void Checkpoint<dim>::loadFromDeprecatedCheckpointVersion(
		DistributionFunctions& f, AdvectionOperator<dim>& advection,
		string directory, CheckpointStatus& status) {

	dealii::DoFHandler<dim>& dof_handler = *advection.getDoFHandler();
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
		if (tmp != advection.getFe()->get_degree()) {
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
			iteration_start,
			MPI_COMM_WORLD).max;
	status.stencilScaling = new_stencil->getScaling();

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
