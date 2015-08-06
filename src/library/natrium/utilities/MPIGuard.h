/*
 * MPIGuard.h
 *
 *  Created on: Jan 29, 2015
 *      Author: kraemer
 */

#ifndef MPIGUARD_H_
#define MPIGUARD_H_

#include <mpi.h>

#include <deal.II/base/mpi.h>

//#include "BasicNames.h"

namespace natrium {

/**
 * @short Singleton that encapsulates the MPI init and finalize calls so that they are called at most once per run.
 *        This class wraps the function dealii::Utilities::MPI::MPI_InitFinalize(argc, argv).
 */
class MPIGuard {
private:
	static MPIGuard* m_privateInstance;
	dealii::Utilities::MPI::MPI_InitFinalize* m_mpi_initialization;
	MPI_Comm m_mpi_communicator;
protected:
	MPIGuard(int argc, char** argv) :
			m_mpi_communicator(MPI_COMM_WORLD) {
		m_mpi_initialization = new dealii::Utilities::MPI::MPI_InitFinalize(
				argc, argv);
	}
public:
	/**
	 * @short Static constructor.
	 * @argc Command line argument. Can be used to determine the number of parallel processes.
	 * See documentation of  dealii::Utilities::MPI::MPI_InitFinalize(argc, argv) for details.
	 * @argv Command line argument.Can be used to determine the number of parallel processes.
	 * See documentation of  dealii::Utilities::MPI::MPI_InitFinalize(argc, argv) for details.
	 */
	static MPIGuard* getInstance(int argc = 0, char** argv = NULL);
	/**
	 * @short return dealii's MPI_InitFinalize object
	 */
	dealii::Utilities::MPI::MPI_InitFinalize* getMPI_InitFinalize() {
		return m_mpi_initialization;
	}
	static MPI_Comm& getMPICommunicator() {
		return MPIGuard::getInstance()->m_mpi_communicator;
	}
	virtual ~MPIGuard() {
		;
	}

};

} /* namespace natrium */

#endif /* MPIGUARD_H_ */
