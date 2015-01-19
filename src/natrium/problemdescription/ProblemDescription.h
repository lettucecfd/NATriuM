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

#include "../utilities/Logging.h"
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

	/// kinematic viscosity
	double m_viscosity;

	/// characteristic length
	double m_characteristicLength;

public:

	/////////////////////////////////
	// CONSTRUCTION // DESTRUCTION //
	/////////////////////////////////

	/// constructor
	ProblemDescription(shared_ptr<dealii::Triangulation<dim> > triangulation,
			double viscosity, double characteristicLength);

	///  destructor
	virtual ~ProblemDescription() {
	}

	/////////////////////////////////
	// GETTER     // SETTER        //
	/////////////////////////////////

	const shared_ptr<dealii::Triangulation<dim> >& getTriangulation() const {
		return m_triangulation;
	}

	const shared_ptr<BoundaryCollection<dim> >& getBoundaries() const {
		return m_boundaries;
	}

	/**
	 * @short set initial densities
	 * @param[out] initialDensities vector of densities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialDensities(distributed_vector& initialDensities,
			const vector<dealii::Point<dim> >& supportPoints) const = 0;

	/**
	 * @short set initial velocities
	 * @param[out] initialVelocities vector of velocities; to be filled
	 * @param[in] supportPoints the coordinates associated with each degree of freedom
	 */
	virtual void applyInitialVelocities(
			vector<distributed_vector>& initialVelocities,
			const vector<dealii::Point<dim> >& supportPoints) const = 0;

	/**
	 * @short check if boundary conditions are uniquely assigned to boundary indicator
	 * @return true, if boundaries OK
	 */
	bool checkBoundaryConditions(){
		// read boundary ids from triangulation
		bool result = true;
		std::vector<dealii::types::boundary_id> tria_boundary_ids(m_triangulation->get_boundary_indicators());

		// read boundary ids from boundary collection
		std::vector<dealii::types::boundary_id> collection_boundary_ids;
		typename BoundaryCollection<dim>::ConstIterator boundary;
		typename BoundaryCollection<dim>::ConstIterator end;
		end = m_boundaries->getBoundaries().end();
		for (boundary = m_boundaries->getBoundaries().begin(); boundary != end; ++ boundary){
			collection_boundary_ids.push_back(boundary->first);
		}
		// check uniqueness in one direction
		std::vector<dealii::types::boundary_id>::iterator it;
		for (size_t i = 0; i < tria_boundary_ids.size(); i++){
			it = std::find(collection_boundary_ids.begin(), collection_boundary_ids.end(), tria_boundary_ids.at(i));
			if (it == collection_boundary_ids.end()){
				LOG(ERROR) << "Found boundary ID " << size_t(tria_boundary_ids.at(i)) << " in mesh, but not in boundaries." << endl;
				result = false;
			} else {
				collection_boundary_ids.erase(it);
			}
		}
		// check uniqueness in other directions
		for (size_t i = 0; i < collection_boundary_ids.size(); i++){
			LOG(ERROR) << "Found boundary ID " << size_t(collection_boundary_ids.at(i)) << " in boundaries, but not in mesh." << endl;
			result = false;
		}
		return result;
	}

	void setTriangulation(
			const shared_ptr<dealii::Triangulation<dim> >& triangulation) {
		m_triangulation = triangulation;
	}

	void setBoundaries(const shared_ptr<BoundaryCollection<dim> >& boundaries) {
		m_boundaries = boundaries;
	}

	double getViscosity() const {
		return m_viscosity;
	}

	void setViscosity(double viscosity) {
		assert (viscosity > 0.0);
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
		shared_ptr<dealii::Triangulation<dim> > triangulation,
		double viscosity, double characteristicLength) :
		m_triangulation(triangulation), m_viscosity(
				viscosity), m_characteristicLength(characteristicLength) {
}

} /* namespace natrium */

#endif /* PROBLEMDESCRIPTION_H_ */
