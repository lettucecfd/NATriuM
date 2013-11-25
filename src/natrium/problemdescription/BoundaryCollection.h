/*
 * BoundaryCollection.h
 *
 *  Created on: 29.10.2013
 *      Author: kraemer
 */

#ifndef BOUNDARYCOLLECTION_H_
#define BOUNDARYCOLLECTION_H_

#include "Boundary.h"

#include "PeriodicBoundary.h"

#include "deal.II/lac/constraint_matrix.h"

#include "../utilities/BasicNames.h"

namespace natrium {

template<size_t dim> class BoundaryCollection {
private:

	/// vector to store boundaries in
	vector<shared_ptr<Boundary<dim> > > m_boundaries;

	/// vector to store periodic boundaries in
	vector<shared_ptr<PeriodicBoundary<dim> > > m_periodicBoundaries;

	/// deal.II constraint matrix
	shared_ptr<dealii::ConstraintMatrix> m_constraintMatrix;

public:
	BoundaryCollection() :
			m_boundaries(0),
			m_periodicBoundaries(0){
		m_constraintMatrix = make_shared<dealii::ConstraintMatrix>();
		m_constraintMatrix->clear();

	}
	virtual ~BoundaryCollection(){

	}
	void addBoundary(shared_ptr<PeriodicBoundary<dim> > boundary) {
		m_boundaries.push_back(boundary);
		m_periodicBoundaries.push_back(boundary);
	}
	void applyBoundaries(
			const shared_ptr<dealii::DoFHandler<dim> > doFHandler) {
		for (size_t i = 0; i < m_boundaries.size(); i++) {
			m_boundaries.at(i)->applyBoundaryConditions(doFHandler,
					m_constraintMatrix);
		}
		m_constraintMatrix->close();
	}

	const vector<shared_ptr<Boundary<dim> > >& getBoundaries() const {
		return m_boundaries;
	}

	const vector<shared_ptr<PeriodicBoundary<dim> > >& getPeriodicBoundaries() const {
		return m_periodicBoundaries;
	}

	const size_t numberOfBoundaries() const {
		return m_boundaries.size();
	}

	const size_t numberOfPeriodicBoundaries() const {
		return m_periodicBoundaries.size();
	}
};

} /* namespace natrium */

#endif /* BOUNDARYCOLLECTION_H_ */
