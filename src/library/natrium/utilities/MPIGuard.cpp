#include "MPIGuard.h"
#include "Logging.h"

namespace natrium {

boost::shared_ptr<MPIGuard> MPIGuard::m_privateInstance = NULL;

/// Parallel cout and cerr
/// Their activity is changed, when MPIGuard is initialized
dealii::ConditionalOStream perr(std::cerr, true);
dealii::ConditionalOStream pout(std::cout, true);

MPIGuard::MPIGuard(int& argc, char**& argv) {
	m_mpi_initialization = boost::make_shared<dealii::Utilities::MPI::MPI_InitFinalize>(argc,
			argv, 1);
	// 1 is the max number of TBB threads

	perr.set_condition(is_rank_0());
	pout.set_condition(is_rank_0());

	// to logfile
	LOG(ALL) << "NATriuM runs on "
			<< dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< " MPI process";
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) != 1)
		LOG(ALL) << "es";
	LOG(ALL) << endl;
	// to pout
	pout << "NATriuM runs on "
			<< dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
			<< " MPI process";
	if (dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) != 1)
		pout << "es";
	pout << endl;

}

boost::shared_ptr<MPIGuard> MPIGuard::getInstance(int& argc, char**& argv) {
	if (m_privateInstance == NULL) {
		LOG(DETAILED) << "MPIGuard  was created." << endl;
		m_privateInstance = boost::make_shared<MPIGuard>(argc, argv);
	} else {
		LOG(DETAILED)
				<< "Double Construction of MPIGuard caught. NATriuM will continue."
				<< endl;
	}
	return m_privateInstance;
}

boost::shared_ptr<MPIGuard> MPIGuard::getInstance() {
	int argc = 0;
	char ** argv = NULL;
	return getInstance(argc, argv);
}

}
