/**
 * @file DiamondObstacle2D.h
 * @short Flow around an obstacle
 * @date 14.12.2015
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef DIAMONDOBSTACLE2D_H_
#define DIAMONDOBSTACLE2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"

namespace natrium {

/**
 * @short Description of a flow around a diamond-shaped obstacle
 */
class DiamondObstacle2D: public ProblemDescription<2> {
private:
	double m_meanInflowVelocity;
	size_t m_refinementLevel;
public:

	class InflowVelocity: public dealii::Function<2> {
	private:
		double m_averageU;
	public:
		InflowVelocity(double av_U) :
				m_averageU(av_U) {

		}
		virtual double value(const dealii::Point<2> &x,
				const unsigned int component) const {
			double H = 4.1;
			if (component == 0) {
				return 4 * m_averageU * x(1) * (H - x(1)) / (H * H);
			}
			return 0.0;
		}
		virtual void vector_value(const dealii::Point<2> &x,
				dealii::Vector<double> &return_value) const {
			return_value(0) = value(x, 0);
			return_value(1) = 0.0;
		}
	};
	/// constructor
	DiamondObstacle2D(double velocity, double viscosity,
			size_t refinementLevel);

	/// destructor
	virtual ~DiamondObstacle2D();

	virtual double getCharacteristicVelocity() const {
		return m_meanInflowVelocity;
	}

	virtual void refineAndTransform() {
		// Refine grid to 8 x 8 = 64 cells; boundary indicators are inherited from parent cell
		getMesh()->refine_global(m_refinementLevel);
	}
private:

	/**
	 * @short create triangulation for lid-driven cavity flow.
	 * @param refinementLevel (denoted as N) The grid will have 2^n*2^n even-sized square cells
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(size_t refinementLevel);

	/**
	 * @short create boundaries for lid-driven cavity flow
	 * @return shared pointer to BoundaryCollection. BoundaryCollection is a container class for
	 *         the specified boundary conditions.
	 * @note All boundary types are inherited of the class Boundary
	 *       Here, we have four wall boundaries. One of them (the upper one)
	 *       has a tangential speed.
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();


};

} /* namespace natrium */
#endif /* LIDDRIVENCAVIT2D_H_ */
