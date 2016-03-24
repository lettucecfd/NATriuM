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

namespace natrium {

/**
 * @short A struct describing the status of a checkpoint (iteration number, time, stencil scaling)
 */
struct CheckpointStatus {
	double time;
	size_t iterationNumber;
	double stencilScaling;
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

public:

	/**
	 * @short Constructor
	 * @param[in] iteration The iteration number of the checkpoint (save or load).
	 * @param[in] checkpoint_dir The directory in which the checkpoint was saved.
	 */
	Checkpoint(size_t iteration, boost::filesystem::path checkpoint_dir);

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
	 * @short Write checkpoint int
	 */
	void write(const Mesh<dim>& mesh, const DistributionFunctions& f,
			const dealii::DoFHandler<dim>& dof_handler,
			const CheckpointStatus& status);

	/**
	 * @short load
	 */
	void load(DistributionFunctions& f, AdvectionOperator<dim>& advection,
			CheckpointStatus& status);
	void loadFromDeprecatedCheckpointVersion(Mesh<dim>& mesh,
			DistributionFunctions& f,
			const dealii::DoFHandler<dim>& dof_handler,
			CheckpointStatus& status, string& directory);
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
