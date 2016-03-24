/*
 * Droplet2D.h
 *
 *  Created on: Dec 7, 2015
 *      Author: kraemer
 */

#ifndef DROPLET2D_H_
#define DROPLET2D_H_

#include "deal.II/grid/tria.h"

#include "../problemdescription/ProblemDescription.h"
#include "../utilities/BasicNames.h"



namespace natrium {

class Droplet2D: public ProblemDescription<2> {
public:

	/**
	 * @short class to describe the inital density
	 */
	class InitialDensity: public dealii::Function<2> {
	private:
		Droplet2D* m_flow;
	public:
		InitialDensity(Droplet2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component = 0) const;
	};

	/**
	 * @short class to describe the x-component of the analytic solution
	 */
	class InitialVelocity: public dealii::Function<2> {
	private:
		Droplet2D* m_flow;
	public:
		InitialVelocity(Droplet2D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<2>& x, const unsigned int component=0) const;
	};

	/// constructor
	Droplet2D(double viscosity, size_t refinementLevel, double length, double height, double rho_l, double rho_g, double W, double R0);

	/// destructor
	virtual ~Droplet2D();

	virtual double getCharacteristicVelocity() const {
		return 0.0;
	}

	double getLength() const {
		return m_length;
	}

	double getHeight() const {
		return m_height;
	}

	double getR0() const {
		return m_R0;
	}

	double getRhoG() const {
		return m_rhoG;
	}

	double getRhoL() const {
		return m_rhoL;
	}

	double getW() const {
		return m_W;
	}

	virtual void refineAndTransform(){
		// Refine grid
		getMesh()->refine_global(m_refinementLevel);
	}
private:

	const double m_length;
	const double m_height;
	const double m_rhoL;
	const double m_rhoG;
	const double m_W;
	const double m_R0;
	const size_t m_refinementLevel;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<2> > makeGrid(double length, double height);

	/**
	 * @short create periodic boundaries
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<2> > makeBoundaries();

};
} /* namespace natrium */

#endif /* DROPLET2D_H_ */
