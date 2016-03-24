/**
 * @file LidDrivenCavity2D.h
 * @short Lid-driven cavity with three static walls and one moving wall
 * @date 31.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef LIDDRIVENCAVIT2D_H_
#define LIDDRIVENCAVIT2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {


/**
 * @short Description of a lid-driven cavity flow.
 */
class LidDrivenCavity2D: public ProblemDescription<2> {
private:
	double topPlateVelocity;
public:

	/// constructor
	LidDrivenCavity2D(double velocity, double viscosity, size_t refinementLevel);

	/// destructor
	virtual ~LidDrivenCavity2D();


	virtual double getCharacteristicVelocity() const {
		return topPlateVelocity;
	}

private:

	/**
	 * @short create triangulation for lid-driven cavity flow.
	 * @param refinementLevel (denoted as N) The grid will have 2^n*2^n even-sized square cells
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid();

	/**
	 * @short create boundaries for lid-driven cavity flow
	 * @return shared pointer to BoundaryCollection. BoundaryCollection is a container class for
	 *         the specified boundary conditions.
	 * @note All boundary types are inherited of the class Boundary
	 *       Here, we have four wall boundaries. One of them (the upper one)
	 *       has a tangential speed.
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();


	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(0.5*(pow(sin(M_PI*(in(0)-0.5)),1)+1), 0.5*(pow(sin(M_PI*(in(1)-0.5)),1)+1));
		}
		UnstructuredGridFunc() {
		}
	};

};

} /* namespace natrium */
#endif /* LIDDRIVENCAVIT2D_H_ */
