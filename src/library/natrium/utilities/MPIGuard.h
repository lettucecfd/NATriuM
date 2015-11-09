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

#include "Logging.h"
//#include "BasicNames.h"

namespace natrium {

/**
 * @short Singleton that encapsulates the MPI init and finalize calls so that they are called at most once per run.
 *        This class wraps the function dealii::Utilities::MPI::MPI_InitFinalize(argc, argv).
 */
class MPIGuard {
private:
	friend class Rank0Stream;
	static MPIGuard* m_privateInstance;
	dealii::Utilities::MPI::MPI_InitFinalize* m_mpi_initialization;

protected:
	MPIGuard(int& argc, char**& argv) {
		m_mpi_initialization = new dealii::Utilities::MPI::MPI_InitFinalize(
				argc, argv, 1);
		// 1 is the max number of TBB threads
		LOG(BASIC) << "NATriuM runs on "
				<< dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
				<< " MPI process";
		if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) != 1)
			LOG(BASIC) << "es";
		LOG(BASIC) << endl;
	}
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

/**
 * @short class for parallel cout. Before using pout, you have to initialize the MPIGuard
 */
class Rank0Stream {
private:
	static bool is_rank_0() {
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		return mpi_rank == 0;
	}
	shared_ptr<dealii::ConditionalOStream> m_instance;
	std::ostream& m_stream;
public:
	Rank0Stream(std::ostream &stream) :
			m_stream(stream) {
	}
	template<typename T>
	const dealii::ConditionalOStream& operator<<(T &t) {
		if (MPIGuard::m_privateInstance == NULL) {
			LOG(ERROR)
					<< "MPIGuard.getMPI_InitFinalize() has to be called before using pout/perr."
					<< endl;
		}
		if (m_instance == NULL) {
			m_instance = make_shared<dealii::ConditionalOStream>(m_stream,
					is_rank_0());
		}
		return ((*m_instance) << t);
	}

};

static Rank0Stream perr(std::cerr);
static Rank0Stream pout(std::cout);

} /* namespace natrium */

#endif /* MPIGUARD_H_ */
