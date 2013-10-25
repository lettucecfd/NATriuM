/**
 * @file StreamingData.h
 * @short Abstract class to store global streaming data like the particle distributions and assemble the matrices.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef STREAMINGDATA_H_
#define STREAMINGDATA_H_

#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class to store global streaming data, like e.g. the particle distributions.
 *  @note  The data to store differs for different approaches of solving the DBE.
 *         The only common data which will appear in every method
 *         is the particle distribution functions, moments and some global system matrix.
 *  @tparam dim The dimension of the flow (2 or 3).
 *  @tparam Q The number of directions of the finite difference stencil (e.g. D2Q9)
 */
template <int dim> class StreamingData {

public:
	/// particle distribution functions
	vector<distributed_vector> m_f;

	/// macroscopic density
	distributed_vector m_density;

	/// macroscopic velocity
	vector<distributed_vector> m_velocity;

	/// constructor
	StreamingData(){};

	/// destructor
	virtual ~StreamingData(){};

	/// function to (re-)assemble linear system
	virtual void reassemble() = 0;
};

} /* namespace natrium */
#endif /* STREAMINGDATA_H_ */
