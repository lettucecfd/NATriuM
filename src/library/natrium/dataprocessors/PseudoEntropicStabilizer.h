/**
 * @file BGKStandard.h
 * @short D2Q9 model description for incompressible flow.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PSEUDOENTROPICSTABILIZER_H_
#define PSEUDOENTROPICSTABILIZER_H_

#include "DataProcessor.h"

#include <cassert>
#include <array>

#include "deal.II/base/index_set.h"

#include "../utilities/BasicNames.h"
#include "../solver/CFDSolver.h"

#include "../utilities/Math.h"

using std::array;

namespace natrium {

/**
 * @short Simple BGK model
 */
template<size_t dim>
class PseudoEntropicStabilizer: public DataProcessor<dim> {
private:
	const bool m_withE;
	void apply_d2q9();
	void apply_d2q9_with_e();
	void apply_d3q19();
public:

	/// constructor
	PseudoEntropicStabilizer(CFDSolver<dim> & solver, bool with_e=false):
		DataProcessor<dim>(solver),
		m_withE(with_e) {
		LOG(BASIC) << "Starting pseudo-entropic stabilizer." << endl;

	}


	/// destructor
	virtual ~PseudoEntropicStabilizer(){

	}



	virtual void apply();


};


	template <size_t Q>
	void applyStabilizer(const array<double,Q>& in, array<double,Q>& out);

	// forward declare explicit specializations
	template<> void applyStabilizer<9>(const array<double,9>& in, array<double,9>& out);
	template<> void applyStabilizer<19>(const array<double,19>& in, array<double,19>& out);
	template<> void applyStabilizer<27>(const array<double,27>& in, array<double,27>& out);


} /* namespace natrium */
#endif /* BGKSTANDARD_H_ */
