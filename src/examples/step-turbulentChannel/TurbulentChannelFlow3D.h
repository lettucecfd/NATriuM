/**
 * @file TurbulentChannelFlow3D.h
 * @short Channel flow in 3D
 * @date 06.01.2016
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef TURBULENTCHANNELFLOW3D_H_
#define TURBULENTCHANNELFLOW3D_H_

#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "natrium/problemdescription/ProblemDescription.h"
#include "natrium/utilities/BasicNames.h"



namespace natrium {

/** @short Description of a turbulent Channel Flow
 *  The domain is [0,5]x[0,1].
 */
class TurbulentChannelFlow3D: public ProblemDescription<3> {
public:

	/**
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0=w0=0, rho0=1)
	 */
	class InitialVelocity: public dealii::Function<3> {
	private:
		TurbulentChannelFlow3D* m_flow;
	public:
		InitialVelocity(TurbulentChannelFlow3D* flow) :
				m_flow(flow) {
		}
		virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;

	};

	/// constructor
	TurbulentChannelFlow3D(double viscosity, size_t refinementLevel, double u_bulk = 1.0, double U_in = 1.0,
			double height = 1.0, double length = 10.0, double width = 3.0, bool is_periodic = true);

	/// destructor
	virtual ~TurbulentChannelFlow3D();


	virtual double getCharacteristicVelocity() const {
		return m_uBulk;
	}
	//TODO: do we really need this?!
	double getRefinementLevel() const {
		return m_refinementLevel;
	}
	double getInletVelocity() const {
		return m_Uin;
	}

	// Functions used by synthetic turbulence generator
	inline void mean(vector<vector<double> > &matrix, double &avg);
	inline void angles(int nmodes, int &iseed, int &iy, vector<int> &iv, vector<double> &fi,
		vector<double> &psi, vector<double> &alfa, vector<double> &teta);
	inline void randf_1(int nmd, int alow, double ahigh, int idum, int &iy, vector<int> &iv,
		vector<double> &out, int &iseed);
	inline void randf_2(int idum, int &iy, vector<int> &iv, double &ran1, int &iseed);

private:

	double m_refinementLevel;
	double m_uBulk;
	double m_Uin;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<3> > makeGrid(double height, double length, double width);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<3> > makeBoundaries(bool is_periodic);

};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
