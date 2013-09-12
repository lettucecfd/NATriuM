/**
 * @file StreamingData.h
 * @short Abstract class to store global streaming data like the particle distributions and assemble the matrices.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef STREAMINGDATA_H_
#define STREAMINGDATA_H_

namespace natrium {

/** @short Abstract class to store global streaming data, like e.g. the particle distributions.
 *  @note  The data to store differs for different approaches of solving the DBE.
 *         The only common data which will appear in every method
 *         is the particle distribution functions and some global system matrix.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template <int dim> class StreamingData {
public:

	/// constructor
	StreamingData(){};

	/// destructor
	virtual ~StreamingData(){};
};

} /* namespace natrium */
#endif /* STREAMINGDATA_H_ */
