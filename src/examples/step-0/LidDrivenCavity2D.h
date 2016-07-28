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
	LidDrivenCavity2D(double velocity, double viscosity,
			size_t refinementLevel);

	/// destructor
	virtual ~LidDrivenCavity2D();

	virtual double getCharacteristicVelocity() const {
		return topPlateVelocity;
	}

	virtual void refine(Mesh<2>& mesh);
	virtual void transform(Mesh<2>& mesh);

private:

	const size_t m_refinementLevel;

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
	private:
		double m_ax;
		double m_bx;
		double m_ay;
		double m_by;

		double shifted_sine(double x) const {
			return 0.5 * (1 + sin(M_PI * x - 0.5 * M_PI));
		}
		double trafo(double x, double a, double b) const {
			assert(a > 0);
			assert(b < 1);
			assert(x >= 0);
			assert(x <= 1);
			return (shifted_sine((b - a) * x + a) - shifted_sine(a))
					/ (shifted_sine(b) - shifted_sine(a));

		}
	public:
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>( trafo(in(0), m_ax, m_bx) , trafo(in(1), m_ay,m_by) );
		}
		UnstructuredGridFunc(double ax, double bx, double ay, double by) :
				m_ax(ax), m_bx(bx), m_ay(ay), m_by(by){
		}
	};

};

} /* namespace natrium */
#endif /* LIDDRIVENCAVIT2D_H_ */
