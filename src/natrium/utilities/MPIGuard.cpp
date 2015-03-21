
#include "MPIGuard.h"
#include "Logging.h"

namespace natrium  {

MPIGuard* MPIGuard::m_privateInstance = NULL;


MPIGuard* MPIGuard::getInstance(int argc, char** argv){
	if (m_privateInstance == NULL){
		LOG(DETAILED) << "MPIGuard  was created." << endl;
		m_privateInstance = new MPIGuard(argc, argv);
	} else {
		LOG(WARNING) << "Double Construction of MPIGuard caught. NATriuM will continue." << endl;
	}
	return m_privateInstance;
}

}
