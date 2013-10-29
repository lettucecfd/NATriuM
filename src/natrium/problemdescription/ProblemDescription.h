/**
 * @file ProblemDescription.h
 * @short Abstract class for the description of a CFD problem. The description includes the computational mesh,
 *        boundary description, viscosity and initial values.
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef PROBLEMDESCRIPTION_H_
#define PROBLEMDESCRIPTION_H_

#include "deal.II/grid/tria.h"

#include "BoundaryCollection.h"
#include "../utilities/BasicNames.h"

namespace natrium {

/** @short Abstract class for the description of a CFD problem. The description includes the computational mesh,
 *         boundary description, viscosity and initial values.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class ProblemDescription {
private:

	/// computational grid
	shared_ptr<dealii::Triangulation<dim> > m_triangulation;

	/// boundary description
	shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// initial densities
	shared_ptr<distributed_vector> m_initialDensities;

	/// initial velocities
	shared_ptr<vector<shared_ptr<distributed_vector> > > m_initialVelocities;

	/// relaxation parameter
	double m_relaxationParameter;

public:

	/////////////////////////////////
	// CONSTRUCTION // DESTRUCTION //
	/////////////////////////////////

	/// constructor
	ProblemDescription(shared_ptr<dealii::Triangulation<dim> > triangulation,
			double relaxationParameter);

	///  destructor
	virtual ~ProblemDescription() {
	}

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////

	const shared_ptr<distributed_vector>& getInitialDensities() const {
		return m_initialDensities;
	}

	const shared_ptr<vector<distributed_vector> >& getInitialVelocities() const {
		return m_initialVelocities;
	}

	double getRelaxationParameter() const {
		return m_relaxationParameter;
	}

	const shared_ptr<dealii::Triangulation<dim> >& getTriangulation() const {
		return m_triangulation;
	}

	const shared_ptr<BoundaryCollection<dim> >& getBoundaries() const {
		return m_boundaries;
	}

	/** @short set initial density
	 */
	void setInitialDensities(shared_ptr<distributed_vector> initialDensities) {
		m_initialDensities = initialDensities;
	}

	/** @short set constant initial density
	 */
	void setConstantInitialDensity(double initialDensity) {
		// TODO implementation
		// TODO (MPI) different constructor for petsc wrapper vectors
	}

	/** @short set initial density
	 */
	void setInitialVelocities(
			shared_ptr<vector<distributed_vector> > initialVelocities) {
		m_initialVelocities = initialVelocities;
	}

	/** @short set constant initial density
	 */
	void setConstantInitialVelocity(const numeric_vector& initialVelocity) {
		// TODO implementation
		// TODO (MPI) different constructor for petsc wrapper vectors
	}

	void setRelaxationParameter(double relaxationParameter) {
		m_relaxationParameter = relaxationParameter;
	}

	void setTriangulation(
			const shared_ptr<dealii::Triangulation<dim> >& triangulation) {
		m_triangulation = triangulation;
	}

	void setBoundaries(const shared_ptr<BoundaryCollection<dim> >& boundaries) {
		m_boundaries = boundaries;
	}
};
/* class ProblemDescription */


template<size_t dim>
inline ProblemDescription<dim>::ProblemDescription(
		shared_ptr<dealii::Triangulation<dim> > triangulation,
		double relaxationParameter) :
		m_triangulation(triangulation), m_relaxationParameter(
				relaxationParameter) {
}

} /* namespace natrium */

#endif /* PROBLEMDESCRIPTION_H_ */
