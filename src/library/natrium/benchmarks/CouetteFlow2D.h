/**
 * @file CouetteFlow2D.h
 * @short Description of a simple Couette Flow (regular channel flow in rectangular domain).
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef COUETTEFLOW2D_H_
#define COUETTEFLOW2D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"
#include "../problemdescription/Benchmark.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Description of a simple Couette Flow (regular channel flow in square domain).
 *  The domain is [0,1]^2. The top plate is moved with constant velocity. The domain consists of
 *  8 x 8 = 64 Elements (contrast to Min and Lee, who have 6 x 6).
 *  @note The analytic solution is obtained by a formula stated in Min and Lee (2011).
 */
class CouetteFlow2D: public Benchmark<2> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class AnalyticVelocity: public dealii::Function<2> {
	private:
		CouetteFlow2D* m_benchmark;
	public:
		AnalyticVelocity(CouetteFlow2D* couette) :
				m_benchmark(couette) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	/// constructor
	CouetteFlow2D(double viscosity, double topPlateVelocity,
			size_t refinementLevel, double L = 1.0, double startTime = 0.0,
			bool isUnstructured = false);

	/// destructor
	virtual ~CouetteFlow2D();

	virtual double getCharacteristicVelocity() const {
		return m_topPlateVelocity;
	}

	virtual void refine(Mesh<2>& mesh);
	virtual void transform(Mesh<2>& mesh);

private:

	const double m_topPlateVelocity;

	const double m_startTime;

	const size_t m_refinementLevel;

	const bool m_isUnstructured;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double L);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries(double topPlateVelocity);

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
		double trans(const double y) const {
			return std::tanh(4 * (y - 0.5)) / tanh(2) / 2 + 0.5;
		}
		dealii::Point<2> operator()(const dealii::Point<2> &in) const {
			return dealii::Point<2>(in(0), trans(in(1)));
		}
	};

	double getStartTime() const {
		return m_startTime;
	}

};

} /* namespace natrium */
#endif /* COUETTEFLOW2D_H_ */
