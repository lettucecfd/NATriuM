/*
 * Checkpoint.h
 * @short Class for creating checkpoint files.
 *  Created on: 10.03.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_SOLVER_CHECKPOINT_H_
#define LIBRARY_NATRIUM_SOLVER_CHECKPOINT_H_

#include "boost/filesystem.hpp"

#include "deal.II/dofs/dof_handler.h"

#include "DistributionFunctions.h"
#include "../utilities/BasicNames.h"
#include "../stencils/Stencil.h"

#include "../utilities/NATriuMException.h"
#include "../problemdescription/ProblemDescription.h"
#include "../advection/AdvectionOperator.h"

namespace natrium {

/**
 * @short A struct describing the status of a checkpoint (iteration number, time, stencil scaling)
 */
struct CheckpointStatus {
	double time;
	size_t iterationNumber;
	double stencilScaling;
	size_t feOrder;
};

/**
 * @short Exception class for CFDSolver
 */
class CheckpointException: public NATriuMException {
private:
	std::string message;
public:
	CheckpointException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	CheckpointException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~CheckpointException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short A class describing the checkpoint.
 */
template<size_t dim>
class Checkpoint {
private:

	/// file that contains the checkpoint status
	boost::filesystem::path m_statusFile;

	/// file that contains the data (mesh + distribution functions)
	boost::filesystem::path m_dataFile;

	bool m_isG;
    static int m_numberOfRefinements;

public:

	/**
	 * @short Constructor
	 * @param[in] iteration The iteration number of the checkpoint (save or load).
	 * @param[in] checkpoint_dir The directory in which the checkpoint was saved.
	 */
	Checkpoint(size_t iteration, boost::filesystem::path checkpoint_dir, bool isG=false);

	/**
	 * @short Function that determines whether the checkpoint files exist.
	 * @return boolean
	 */
	bool exists() {
		if (not boost::filesystem::exists(m_statusFile)) {
			return false;
		}
		if (not boost::filesystem::exists(m_dataFile)) {
			return false;
		}
		return true;
	}

	/**
	 * @short Write checkpoint to checkpoint directory. A checkpoint contains of three files:
	 *        -# a .stat file that contains some information about the status of the simulations
	 *        (iteration number, time, stencil scaling, ...).
	 *        -# a .data file that contains the serialized data (mesh and distribution functions)
	 *        -# a .data.info file that contains some deal.II information about the mesh, refinement etc.
	 */
	void write(const Mesh<dim>& mesh, DistributionFunctions& f,
			const dealii::DoFHandler<dim>& dof_handler,
			const CheckpointStatus& status);

	/**
	 * @short Load checkpoint from file. A simulation can be resumed from a checkpoint if
	 * 		  -# the mesh is exactly the same as in the previous simulation
	 * 		  -# the mesh is a globally refined version of the previous simulation
	 * 	      Restarts are also possible with arbitrary Mach number.
	 * 	      Varying the order of finite elements is not supported, so far.
	 */
	void load(DistributionFunctions& f, ProblemDescription<dim>& problem,
			AdvectionOperator<dim>& advection, CheckpointStatus& status);
	static void loadFromDeprecatedCheckpointVersion(DistributionFunctions& f,
			AdvectionOperator<dim>& advection, string directory, CheckpointStatus& status);
	virtual ~Checkpoint();

	const boost::filesystem::path& getDataFile() const {
		return m_dataFile;
	}

	const boost::filesystem::path& getStatusFile() const {
		return m_statusFile;
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_SOLVER_CHECKPOINT_H_ */
