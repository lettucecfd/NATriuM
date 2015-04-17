/*
 * D2Q9PseudopotentialModel.h
 *
 *  Created on: Nov 11, 2014
 *      Author: kk
 */

#ifndef D2Q9PSEUDOPOTENTIALMODEL_H_
#define D2Q9PSEUDOPOTENTIALMODEL_H_

#include "D2Q9Model.h"

#include "deal.II/grid/tria.h"
#include "deal.II/numerics/vector_tools.h"

//#include "../advection/SEDGMinLee.h"

#include "../utilities/BasicNames.h"

namespace natrium {

// forward declaration
template<size_t dim>
class SEDGMinLee;

class D2Q9PseudopotentialModel: public D2Q9Model {
private:
	shared_ptr<SEDGMinLee<2> > m_advectionOperator;
	const double m_dt;

public:

	/// constructor
	D2Q9PseudopotentialModel(double scaling, const double dt);

	/// destructor
	virtual ~D2Q9PseudopotentialModel();

	virtual double getEquilibriumDistribution(size_t i, const numeric_vector& u,
			const double rho = 1) const;

	//virtual vector<double> getDensityGradient() ;

	void getInteractionForce(const vector<double>& distributions,
			numeric_vector & interactionForce, const double rho = 1);

	const shared_ptr<SEDGMinLee<2> > getAdvectionOperator() const {
		return m_advectionOperator;
	}

	void setAdvectionOperator(
			const shared_ptr<SEDGMinLee<2> >& advectionOperator) {
		m_advectionOperator = advectionOperator;
	}
};

} /* namespace natrium */

#endif /* D2Q9PSEUDOPOTENTIALMODEL_H_ */
