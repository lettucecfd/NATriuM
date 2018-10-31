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
#include "deal.II/base/function.h"

#include "BoundaryCollection.h"
#include "ConstantExternalForce.h"

#include "../utilities/Logging.h"
#include "../utilities/BasicNames.h"
#include "../utilities/MPIGuard.h"

namespace natrium {

/** @short Abstract class for the description of a CFD problem. The description includes the computational mesh,
 *         boundary description, viscosity and initial values.
 *  @tparam dim The dimension of the flow (2 or 3).
 */
template<size_t dim> class ProblemDescription {
private:

	/// computational grid
	boost::shared_ptr<Mesh<dim> > m_triangulation;

	/// boundary description
	boost::shared_ptr<BoundaryCollection<dim> > m_boundaries;

	/// constant external force
	boost::shared_ptr<ConstantExternalForce<dim> > m_externalForce;

	/// kinematic viscosity
	double m_viscosity;

	/// characteristic length
	double m_characteristicLength;

	/// function to define initial densities
	boost::shared_ptr<dealii::Function<dim> > m_initialRho;

	/// function to define initial densities
	boost::shared_ptr<dealii::Function<dim> > m_initialT;

	/// function to define initial velocities
	boost::shared_ptr<dealii::Function<dim> > m_initialU;

protected:
	void setInitialRho(boost::shared_ptr<dealii::Function<dim> > ini_rho) {
		m_initialRho = ini_rho;
	}

	void setInitialT(boost::shared_ptr<dealii::Function<dim> > ini_T) {
		m_initialRho = ini_T;
	}

	void setInitialU(boost::shared_ptr<dealii::Function<dim> > ini_u) {
		m_initialU = ini_u;
	}

public:

	/////////////////////////////////
	// CONSTRUCTION // DESTRUCTION //
	/////////////////////////////////

	/// constructor
	ProblemDescription(boost::shared_ptr<Mesh<dim> > triangulation,
			double viscosity, double characteristicLength);

	///  destructor
	virtual ~ProblemDescription() {
	}

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////

	const boost::shared_ptr<Mesh<dim> >& getMesh() const {
		return m_triangulation;
	}

	const boost::shared_ptr<BoundaryCollection<dim> >& getBoundaries() const {
		return m_boundaries;
	}

	virtual const boost::shared_ptr<dealii::Function<dim> >& getInitialRhoFunction() const {
		return m_initialRho;
	}

	virtual const boost::shared_ptr<dealii::Function<dim> >& getInitialUFunction() const {
		return m_initialU;
	}

	virtual const boost::shared_ptr<dealii::Function<dim> >& getInitialTFunction() const {
		return m_initialT;
	}

	virtual void refine(Mesh<dim>& mesh) = 0;
	virtual void transform(Mesh<dim>& mesh) = 0;
	virtual bool isCartesian() = 0;

	void refineAndTransform() {
		refine(*getMesh());
		transform(*getMesh());
	}

	void refineAndTransform(Mesh<dim>& mesh) {
		refine(mesh);
		transform(mesh);
	}

	/**
	 * @short check if boundary conditions are uniquely assigned to boundary indicator
	 * @return true, if boundaries OK
	 */
	bool checkBoundaryConditions() {
		// read boundary ids from triangulation
		bool result = true;
		std::vector<dealii::types::boundary_id> tria_boundary_ids(
				m_triangulation->get_boundary_ids());

		// read boundary ids from boundary collection
		std::vector<dealii::types::boundary_id> collection_boundary_ids;
		typename BoundaryCollection<dim>::ConstIterator boundary;
		typename BoundaryCollection<dim>::ConstIterator end;
		end = m_boundaries->getBoundaries().end();
		for (boundary = m_boundaries->getBoundaries().begin(); boundary != end;
				++boundary) {
			collection_boundary_ids.push_back(boundary->first);
		}
		// check uniqueness in one direction
		std::vector<dealii::types::boundary_id>::iterator it;
		for (size_t i = 0; i < tria_boundary_ids.size(); i++) {
			it = std::find(collection_boundary_ids.begin(),
					collection_boundary_ids.end(), tria_boundary_ids.at(i));
			if (it == collection_boundary_ids.end()) {
				LOG(ERROR) << "Found boundary ID "
						<< size_t(tria_boundary_ids.at(i))
						<< " in mesh, but not in boundaries." << endl;
				result = false;
			} else {
				collection_boundary_ids.erase(it);
			}
		}
		// check uniqueness in other directions
		for (size_t i = 0; i < collection_boundary_ids.size(); i++) {
			LOG(ERROR) << "Found boundary ID "
					<< size_t(collection_boundary_ids.at(i))
					<< " in boundaries, but not in mesh." << endl;
			result = false;
		}
		return result;
	}

	void setMesh(const boost::shared_ptr<Mesh<dim> >& triangulation) {
		m_triangulation = triangulation;
	}

	void setBoundaries(
			const boost::shared_ptr<BoundaryCollection<dim> >& boundaries) {
		m_boundaries = boundaries;
	}

	void setExternalForce(
			const boost::shared_ptr<ConstantExternalForce<dim> >& force) {
		m_externalForce = force;
	}

	boost::shared_ptr<ConstantExternalForce<dim> > getExternalForce() const {
		return m_externalForce;
	}

	bool hasExternalForce() const {
		if (m_externalForce == NULL) {
			return false;
		}
		if (m_externalForce->getForce()[0] != 0) {
			return true;
		}
		if (m_externalForce->getForce()[1] != 0) {
			return true;
		}
		if (3 == dim) {
			if (m_externalForce->getForce()[2] != 0) {
				return true;
			}
		}
		return false;
	}

	double getViscosity() const {
		return m_viscosity;
	}

	void setViscosity(double viscosity) {
		assert(viscosity > 0.0);
		m_viscosity = viscosity;
	}

	double getCharacteristicLength() const {
		return m_characteristicLength;
	}

	void setCharacteristicLength(double characteristicLength) {
		m_characteristicLength = characteristicLength;
	}

	virtual double getCharacteristicVelocity() const {
		return 0.0;
	}

};
/* class ProblemDescription */

template<size_t dim>
inline ProblemDescription<dim>::ProblemDescription(
		boost::shared_ptr<Mesh<dim> > triangulation, double viscosity,
		double characteristicLength) :
		m_triangulation(triangulation), m_viscosity(viscosity), m_characteristicLength(
				characteristicLength) {
	// make default initial conditions (rho = 1, u = v = 0)
	m_initialRho = boost::make_shared<dealii::ConstantFunction<dim> >(1.0, 1);
	m_initialT = boost::make_shared<dealii::ConstantFunction<dim> >(1.0, 1);
	m_initialU = boost::make_shared<dealii::ConstantFunction<dim> >(0.0, dim);
	/// Create MPI (if not done yet);
	MPIGuard::getInstance();
}

} /* namespace natrium */

#endif /* PROBLEMDESCRIPTION_H_ */
