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
	 * TODO: edit description
	 * @short class to describe the x-component of the analytic solution
	 * @note other are default (v0 = w0 = 0, rho0 = 1)
	 */
	class InitialVelocity: public dealii::Function<3> {
	private:
		TurbulentChannelFlow3D *m_flow;

	public:
		//double *m_maxUtrp;
		InitialVelocity(TurbulentChannelFlow3D *flow) :
				m_flow(flow) {
			//m_maxUtrp = new double (0.0);
		}
//		~InitialVelocity(){
//			delete m_maxUtrp;
//		}

		virtual double value(const dealii::Point<3>& x,
				const unsigned int component = 0) const;
	};

    class InitialDensity: public dealii::Function<3> {
    private:
        TurbulentChannelFlow3D* m_flow;
    public:
        InitialDensity(TurbulentChannelFlow3D* flow) :
                m_flow(flow) {
        }
        virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;


    };

    class InitialTemperature: public dealii::Function<3> {
    private:
        TurbulentChannelFlow3D* m_flow;
    public:
        InitialTemperature(TurbulentChannelFlow3D* flow) :
                m_flow(flow) {
        }
        virtual double value(const dealii::Point<3>& x, const unsigned int component=0) const;
    };

	class IncompressibleU: public dealii::Function<3> {
	private:
		TurbulentChannelFlow3D *m_flow;
		InitialVelocity m_initialU;
		double m_increment;
	public:
		//double *m_maxIncUtrp;
		IncompressibleU(TurbulentChannelFlow3D *flow) :
				m_flow(flow), m_initialU(flow), m_increment(1e-6) {
			//m_maxIncUtrp = new double (0.0);
		}
//		~IncompressibleU(){
//		        delete m_maxIncUtrp;
//		}

		virtual double value(const dealii::Point<3>& x,
				const unsigned int component = 0) const;
	};

	class MeanVelocityProfile: public dealii::Function<3> {
	private:
		TurbulentChannelFlow3D *m_flow;
		IncompressibleU m_initialIncompressibleU;
	public:
		MeanVelocityProfile(TurbulentChannelFlow3D *flow) :
				m_flow(flow), m_initialIncompressibleU(flow) {
		}
		virtual double value(const dealii::Point<3>& x,
				const unsigned int component = 0) const;
	};

	/**
	 * @short function to generate the unstructured mesh grid
	 */
	struct UnstructuredGridFunc {
		double m_length;
		double m_height;
		double m_width;
		double m_gridDensity;
		UnstructuredGridFunc(double length, double height, double width, double gridDensity = 0.8) :
				m_length(length), m_height(height), m_width(width), m_gridDensity(gridDensity) {
		}
		double trans(const double y) const {

            double new_y = (2*y);
            new_y *= M_PI;
            new_y = -m_gridDensity*sin(new_y)/(2*M_PI);
            new_y += y;
            return new_y*m_height;
		}
		dealii::Point<3> operator()(const dealii::Point<3> &in) const {
			return dealii::Point<3>(m_length * in(0), trans(in(1)),
					m_width * in(2));
		}
	};

	/////////////////////////////////
	// CONSTRUCTION // DESTRUCTION //
	/////////////////////////////////

	/// constructor
	TurbulentChannelFlow3D(double viscosity, size_t refinementLevel,
			std::vector<unsigned int> repetitions, double ReTau = 180.0,
			double u_cl = 10, double height = 1.0, double length = 6.0,
			double width = 3.0,
			bool is_periodic = true, double gridDensity=0.8);

	/// destructor
	virtual ~TurbulentChannelFlow3D();

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////

	double getRefinementLevel() const {
		return m_refinementLevel;
	}

	std::vector<unsigned int> getRepetitions() const {
		return m_repetitions;
	}

	double getFrictionReNumber() const {
		return m_ReTau;
	}

	double getCenterLineVelocity() const {
		return m_uCl;
	}

	double getHeight() const {
		return m_height;
	}

	double getLength() const {
		return m_length;
	}

	double getWidth() const {
		return m_width;
	}

	double getMaxUtrp() const {
		return m_maxUtrp;
	}
	double getMaxIncUtrp() const {
		return m_maxIncUtrp;
	}

	// Functions used by synthetic turbulence generator
	inline void mean(vector<vector<double> > &matrix, double &avg);
	inline void angles(int nmodes, int &iseed, int &iy, vector<int> &iv,
			vector<double> &fi, vector<double> &psi, vector<double> &alfa,
			vector<double> &teta);
	inline void randf_1(int nmd, int alow, double ahigh, int idum, int &iy,
			vector<int> &iv, vector<double> &out, int &iseed);
	inline void randf_2(int idum, int &iy, vector<int> &iv, double &ran1,
			int &iseed);

	virtual void refine(Mesh<3>& mesh) {
		// refine global
		mesh.refine_global(m_refinementLevel);
	}
	virtual void transform(Mesh<3>& mesh) {
		// transform grid to unstructured grid
		dealii::GridTools::transform(
				UnstructuredGridFunc(m_length, m_height, m_width, m_gridDensity), mesh);
	}

	virtual bool isCartesian(){
		return true;
	}
private:

	double m_refinementLevel;
	std::vector<unsigned int> m_repetitions;
	double m_ReTau;
	double m_uCl;
	double m_height;
	double m_length;
	double m_width;
	double m_maxUtrp;
	double m_maxIncUtrp;
	double m_gridDensity;

	/**
	 * @short create triangulation for couette flow
	 * @return shared pointer to a triangulation instance
	 */
	boost::shared_ptr<Mesh<3> > makeGrid(std::vector<unsigned int> repetitions);

	/**
	 * @short create boundaries for couette flow
	 * @return shared pointer to a vector of boundaries
	 * @note All boundary types are inherited of BoundaryDescription; e.g. PeriodicBoundary
	 */
	boost::shared_ptr<BoundaryCollection<3> > makeBoundaries(bool is_periodic);

};

} /* namespace natrium */
#endif /* PERIODICFLOW2D_H_ */
