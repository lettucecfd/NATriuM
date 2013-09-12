/**
 * @file DataMinLee2011.h
 * @short Global data which is used by Min and Lee (2011): A spectral-elemennt discontinuous
 *        Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DATAMINLEE2011_H_
#define DATAMINLEE2011_H_

#include "StreamingData.h"

namespace natrium {


/** @short Global data which is used, e.g., by Min and Lee (2011): A spectral-element discontinuous
 *         Galerkin lattice Boltzmann method for nearly incompressible flows, JCP 230 pp. 245-259.
 *         including particle distributions f, system matrix L, diagonal mass matrix M,
 *         gradient matrices Dx, Dy, (Dz) and boundary matrix R
 * @tparam dim The dimension of the flow (2 or 3).
 */
template <int dim> class DataMinLee2011: public  StreamingData<dim>{
public:

	/// constructor
	DataMinLee2011(){};

	/// destructor
	virtual ~DataMinLee2011(){};
};

} /* namespace natrium */
#endif /* DATAMINLEE2011_H_ */
