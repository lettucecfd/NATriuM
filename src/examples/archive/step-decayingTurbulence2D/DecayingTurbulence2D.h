/**
 * @file DecayingTurbulence2D.h
 * @short
 * @date 11.01.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DECAYINGTURBULENCE2D_H_
#define DECAYINGTURBULENCE2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {

/** @short Description of a turbulent Channel Flow
 *  The domain is [0,5]x[0,1].
 */
class DecayingTurbulence2D: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		DecayingTurbulence2D* m_flow;
	public:
		InitialVelocity(DecayingTurbulence2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	/// constructor
	DecayingTurbulence2D(double viscosity, size_t refinementLevel);

	/// destructor
	virtual ~DecayingTurbulence2D();


	virtual double getCharacteristicVelocity() const {
		return 1;
	}

	virtual void refine(Mesh<2>& mesh) {
		// refine global
		mesh.refine_global(m_refinementLevel);
	}

	virtual void transform(Mesh<2>& mesh){

	}

	virtual bool isCartesian(){
		return true;
	}
private:

	const size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();


};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
