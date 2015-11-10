/*
 * MPIGuard.h
 *
 *  Created on: Jan 29, 2015
 *      Author: kraemer
 */

#ifndef MPIGUARD_H_
#define MPIGUARD_H_

#include <mpi.h>

#include "deal.II/base/mpi.h"
#include "deal.II/base/conditional_ostream.h"

#include "Logging.h"
#include "BasicNames.h"

namespace natrium {


// defined in cpp-file
extern dealii::ConditionalOStream perr;
extern dealii::ConditionalOStream pout;

/**
 * @short Singleton that encapsulates the MPI init and finalize calls so that they are called at most once per run.
 *        This class wraps the function dealii::Utilities::MPI::MPI_InitFinalize(argc, argv).
 */
class MPIGuard {
private:
	static MPIGuard* m_privateInstance;
	dealii::Utilities::MPI::MPI_InitFinalize* m_mpi_initialization;

	bool is_rank_0() {
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		return mpi_rank == 0;
	}

protected:
	MPIGuard(int& argc, char**& argv);
public:
	/**
	 * @short Static constructor.
	 * @argc Command line argument. Can be used to determine the number of parallel processes.
	 * See documentation of  dealii::Utilities::MPI::MPI_InitFinalize(argc, argv) for details.
	 * @argv Command line argument.Can be used to determine the number of parallel processes.
	 * See documentation of  dealii::Utilities::MPI::MPI_InitFinalize(argc, argv) for details.
	 */
	static MPIGuard* getInstance(int& argc, char**& argv);

	static MPIGuard* getInstance();
	/**
	 * @short return dealii's MPI_InitFinalize object
	 */
	dealii::Utilities::MPI::MPI_InitFinalize* getMPI_InitFinalize() {
		return m_mpi_initialization;
	}

	virtual ~MPIGuard() {
		;
	}

};


} /* namespace natrium */

#endif /* MPIGUARD_H_ */
