/**
 * @file BoundaryTools.h
 * @short 
 * @date 14.11.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef BOUNDARYTOOLS_H_
#define BOUNDARYTOOLS_H_

#include <string>

#include "deal.II/base/point.h"
#include "deal.II/base/function.h"
#include "deal.II/base/tensor_function.h"

#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../solver/DistributionFunctions.h"
#include "../stencils/Stencil.h"

namespace natrium {

namespace BoundaryTools {

/**
 * @short Exception class for Boundaries
 */
class BoundaryException: public NATriuMException {
private:
	std::string message;
public:
	BoundaryException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	BoundaryException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~BoundaryException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

template<size_t dim>
class BoundaryDensity: public dealii::Function<dim> {
private:
	double m_density;
public:
	BoundaryDensity(double rho = 1) {
		m_density = rho;
	}
	;
	virtual ~BoundaryDensity() {
	}
	;
	virtual double value(const dealii::Point<dim> &,
			const unsigned int  = 0) const {
		return m_density;
	}
};
template<size_t dim>
class BoundaryPressure: public dealii::Function<dim> {
private:
	double m_pressure;
public:
	BoundaryPressure(double p = 0) {
		m_pressure = p;
	}
	;
	virtual ~BoundaryPressure() {
	}
	;
	virtual double value(const dealii::Point<dim> &,
			const unsigned int  = 0) const {
		return m_pressure;
	}
};

    template<size_t dim>
    class BoundaryTemperature: public dealii::Function<dim> {
    private:
        double m_temperature;
    public:
        BoundaryTemperature(double T = 1.0) {
            m_temperature = T;
        }
        ;
        virtual ~BoundaryTemperature() {
        }
        ;
        virtual double value(const dealii::Point<dim> &,
                             const unsigned int  = 0) const {
            return m_temperature;
        }
    };
template<size_t dim>
class BoundaryVelocity: public dealii::Function<dim> {
private:
	dealii::Vector<double> m_Velocity;
public:
	BoundaryVelocity(const dealii::Vector<double>& velocity) :
		dealii::Function<dim>(dim),
			m_Velocity(velocity) {
		assert (velocity.size() == dim);
	}
	BoundaryVelocity(const dealii::Tensor<1,dim>& velocity):
			dealii::Function<dim>(dim){
		m_Velocity.reinit(dim);
		velocity.unroll(m_Velocity);
	}
	virtual ~BoundaryVelocity() {
	}
	;
	virtual double value(const dealii::Point<dim> &,  const unsigned int component=0) const {
		assert (component < dim);
		return m_Velocity(component);
	}
	virtual void vector_value(const dealii::Point<dim> &,
			dealii::Vector<double> &values) const {
		values = m_Velocity;
	}
};

template<size_t dim>
class BoundaryVelocityGradient: public dealii::TensorFunction<2, dim> {
private:
	dealii::Tensor<2, dim> m_VelocityGradient;
public:
	BoundaryVelocityGradient(const dealii::Tensor<2, dim>& velocity_gradient) :
		m_VelocityGradient(velocity_gradient) {
	}
	virtual ~BoundaryVelocityGradient() {
	}
	;
	virtual dealii::Tensor<2, dim> value(const dealii::Point<dim> & ) const {
		return m_VelocityGradient;
	}
};


/**
 * @short function to compare points as map keys;
 *        the points with smaller x components fall before others. in case of equal x,
 *        the points with smaller y components will be placed first.
 */
class point_sorter : public std::binary_function<dealii::Point<2> , dealii::Point<2> , bool > {
public:

  bool operator()(const dealii::Point<2> &left, const dealii::Point<2> &right) const {
    if (left[0] < right[0]){return true;}
        else if ((left[0] == right[0])&&(left[1] < right[1])) {return true;}
    else  return false;
        }
};

/**
 * @short enum to describe couplings at the boundary
 */
enum DistributionCouplingAtBoundary{
	COUPLE_ONLY_OPPOSITE_DISTRIBUTIONS,
	COUPLE_ALL_DISTRIBUTIONS
};

/**
 * @short enum to describe couplings at the boundary
 */
enum PointCouplingAtBoundary{
	COUPLE_ONLY_SINGLE_POINTS,
	COUPLE_WHOLE_FACE
};


/**
 * @short Check if two lines in a 2D plane are parallel and not equal to each other.
 *        If they are antiparallel, swap begin and end of the second line.
 *  @param[in] beginLine1 start point of line 1
 *  @param[in] endLine1 end point of line 1
 *  @param[in/out] beginLine2 start point of line 2
 *  @param[in/out] endLine2 end point of line 2
 *  @param[out] errorMessage contains the reason of failure when lines are not parallel
 *  @return true: if the lines are parallel or antiparallel, else: false
 */
bool checkParallelLines(const dealii::Point<2>& beginLine1,
		const dealii::Point<2>& endLine1, dealii::Point<2>& beginLine2,
		dealii::Point<2>& endLine2, std::string& errorMessage);

/**
 * @short Get the positions of the lines that define interfaces from the
 *        defining Boundary indicators.
 *
 *  @param boundaryIndicator1 boundary indicator of interface line 1
 *  @param boundaryIndicator2 boundary indicator of interface line 2
 *  @param triangulation A (shared ptr to a) triangulation object (the mesh)
 *  @param[out] beginLine1 start point of line 1
 *  @param[out] endLine1 end point of line 1
 *  @param[out] beginLine2 start point of line 2
 *  @param[out] endLine2 end point of line 2
 *  @param[out] errorMessage contains the reason of failure when interfaces are no lines
 *
 *  @return true: if the interfaces are lines, false: else
 */
bool getInterfacialLinesByBoundaryIndicator(size_t boundaryIndicator1,
		size_t boundaryIndicator2,
		boost::shared_ptr<Mesh<2> > triangulation,
		dealii::Point<2>& beginLine1, dealii::Point<2>& endLine1,
		dealii::Point<2>& beginLine2, dealii::Point<2>& endLine2,
		std::string& errorMessage);



/**
 * @short This functions adds elements to the sparsity pattern to couple different distribution functions
 * at the boundary.
 * @param[in/out] cSparse The sparsity pattern to add to.
 * @param[in] doFHandler The doFHandler
 * @param[in] stencil The DQ stencil is required to identify the opposite distributions.
 * @param[in] coupling Describes the point coupling that we apply at the boundary. The possibilities are:
 *				 -# COUPLE_ONLY_SINGLE_POINTS: Couples the dofs that belong to the same integration point.
 *				 -# COUPLE_WHOLE_FACE: Couples all dofs at the face with each other
 *				    (is required when gradients are calculated at the boundary)
 */
template<size_t dim>
void CoupleDoFsAtBoundary(
		dealii::TrilinosWrappers::SparsityPattern& cSparse,
		const dealii::DoFHandler<dim>& doFHandler, size_t boundary_id, PointCouplingAtBoundary coupling);


} /* namespace BoundaryTools */

/**
 * @short Container for global Data that is required by the FEBoundaryValues class
 */
struct GlobalBoundaryData {

	/// old distribution functions (at t-dt)
	const DistributionFunctions& m_fold;

	/// new distribution functions (at t)
	DistributionFunctions& m_fnew;

    /// new distribution functions (at t)
    DistributionFunctions& m_g;

	/// LBM stencil
	const Stencil& m_stencil;

	/// viscosity
	double m_viscosity;

	/// time step
	double m_dt;

	/// number of discrete particle velocities in stencil
	size_t m_Q;

	/// speed of sound
	double m_cs2;

	/***
	 * @short Constructor
	 * @param f_old old distribution functions (at t-dt, read only)
	 * @param f_new new distribution functions (at t, write only)
	 * @param stencil LBM stencil (e.g. D2Q9)
	 * @param viscosity kinematic viscosity of the fluid
	 * @param dt time step
	 */
	GlobalBoundaryData(const DistributionFunctions& f_old,
			DistributionFunctions& f_new, const Stencil& stencil,
			double viscosity, double dt) :
			m_fold(f_old), m_fnew(f_new), m_g(f_new), m_stencil(stencil) {
		m_viscosity = viscosity;
		m_dt = dt;
		m_cs2 = stencil.getSpeedOfSoundSquare();
		m_Q = stencil.getQ();

	}
    GlobalBoundaryData(const DistributionFunctions& f_old,
                       DistributionFunctions& f_new, DistributionFunctions& g, const Stencil& stencil,
                       double viscosity, double dt) :
            m_fold(f_old), m_fnew(f_new), m_g(g), m_stencil(stencil) {
        m_viscosity = viscosity;
        m_dt = dt;
        m_cs2 = stencil.getSpeedOfSoundSquare();
        m_Q = stencil.getQ();
    }



	virtual ~GlobalBoundaryData() {

	}
};

} /* namespace natrium */

#endif /* BOUNDARYTOOLS_H_ */
