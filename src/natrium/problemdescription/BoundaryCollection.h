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
#include "MinLeeBoundary.h"

#include "deal.II/lac/constraint_matrix.h"

#include "../utilities/BasicNames.h"

namespace natrium {

/**
 * @short Boundary errors (e.g. duplicate boundary indicators)
 */
class BoundaryCollectionError: public std::exception {
private:
	std::string m_message;
public:
	BoundaryCollectionError() :
			m_message("Error in Boundary collection") {
	}
	BoundaryCollectionError(const std::string& message) :
			m_message("Error in Boundary collection: " + message) {
	}
	virtual const char* what() const throw () {
		return m_message.c_str();
	}
	virtual ~BoundaryCollectionError() throw () {
	}
};

/**
 * @short The BoundaryCollection class defines all boundaries of a flow domain.
 *        Internally, the boundaries are stored in two different std::map, one
 *        for the periodic boundaries and one for the non-periodic ones.
 *        Its keys are the boundary indicators (for periodic boundaries: the first boundary indicator).
 */
template<size_t dim> class BoundaryCollection {
private:
	/// vector to store boundaries in
	std::map<size_t, shared_ptr<Boundary<dim> > > m_boundaries;

	/// vector to store boundaries in
	std::map<size_t, shared_ptr<MinLeeBoundary<dim> > > m_minLeeBoundaries;

	/// vector to store periodic boundaries in
	std::map<size_t, shared_ptr<PeriodicBoundary<dim> > > m_periodicBoundaries;

public:

	typedef typename std::map<size_t, shared_ptr<Boundary<dim> > >::iterator Iterator;
	typedef typename std::map<size_t, shared_ptr<MinLeeBoundary<dim> > >::iterator MinLeeIterator;
	typedef typename std::map<size_t, shared_ptr<MinLeeBoundary<dim> > >::const_iterator ConstMinLeeIterator;
	typedef typename std::map<size_t, shared_ptr<PeriodicBoundary<dim> > >::iterator PeriodicIterator;
	typedef typename std::map<size_t, shared_ptr<PeriodicBoundary<dim> > >::const_iterator ConstPeriodicIterator;

	BoundaryCollection() {
	}

	virtual ~BoundaryCollection() {
	}

	/**
	 * @short Add a boundary to the flow definition. This definition of addBoundary applies to periodic boundaries.
	 * 		  As periodic boundaries have two boundary indicators, they are added twice to the boundary map.
	 * 		  Internally, they are additionally added to the vector periodicBoundaries, which can be accessed separate from non-periodic boundaries.
	 * @param boundary a periodic boundary
	 * @throws BoundaryCollectionError, e.g. if boundary indicators are not unique
	 */
	void addBoundary(shared_ptr<PeriodicBoundary<dim> > boundary) {
		bool success1 = m_boundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator1(), boundary)).second;
		bool success2 = m_boundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator2(), boundary)).second;
		if ((not success1) or (not success2)) {
			throw BoundaryCollectionError(
					"Boundary could not be inserted. Boundary indicators must be unique.");
		}
		m_periodicBoundaries.insert(std::make_pair(boundary->getBoundaryIndicator1(), boundary));
		m_periodicBoundaries.insert(std::make_pair(boundary->getBoundaryIndicator2(), boundary));
	}

	/**
	 * @short Add a boundary to the flow definition.
	 * @param boundary a periodic boundary
	 * @throws BoundaryCollectionError, e.g. if boundary indicators are not unique
	 */
	void addBoundary(shared_ptr<MinLeeBoundary<dim> > boundary) {
		bool success = m_boundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator(), boundary)).second;
		if (not success) {
			throw BoundaryCollectionError(
					"Boundary could not be inserted. Boundary indicators must be unique.");
		}
		m_minLeeBoundaries.insert(std::make_pair(boundary->getBoundaryIndicator1(), boundary));
	}

	/**
	 * @short get a specific boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const shared_ptr<Boundary<dim> >& getBoundary(size_t boundaryIndicator) const{
		if (m_boundaries.count(boundaryIndicator) == 0){
			throw BoundaryCollectionError("in getBoundary: This boundary collection does not contain a boundary with the specified boundary indicator.");
		}
		return m_boundaries.at(boundaryIndicator);
	}

	/**
	 * @short get a specific periodic boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const shared_ptr<PeriodicBoundary<dim> >& getPeriodicBoundary(size_t boundaryIndicator) const {
		assert (isPeriodic(boundaryIndicator));
		if (m_periodicBoundaries.count(boundaryIndicator) == 0){
			throw BoundaryCollectionError("in getPeriodicBoundary: This boundary collection does not contain a periodic boundary with the specified boundary indicator.");
		}
		return m_periodicBoundaries.at(boundaryIndicator);
	}

	/**
	 * @short get a specific MinLee boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const shared_ptr<MinLeeBoundary<dim> >& getMinLeeBoundary(size_t boundaryIndicator) const{
		assert (not isPeriodic(boundaryIndicator));
		if (m_minLeeBoundaries.count(boundaryIndicator) == 0){
			throw BoundaryCollectionError("in minLeeBoundary: This boundary collection does not contain a periodic boundary with the specified boundary indicator.");
		}
		return m_minLeeBoundaries.at(boundaryIndicator);
	}

	/**
	 * @short test if the boundary with the given boundary indicator is periodic
	 */
	bool isPeriodic(size_t boundaryIndicator) const {
		return getBoundary(boundaryIndicator)->isPeriodic();
	}

	const size_t numberOfBoundaries() const {
		return m_boundaries.size();
	}

	const size_t numberOfPeriodicBoundaries() const {
		return m_periodicBoundaries.size();
	}

	const std::map<size_t, shared_ptr<Boundary<dim> > >& getBoundaries() const {
		return m_boundaries;
	}

	const std::map<size_t, shared_ptr<PeriodicBoundary<dim> > >& getPeriodicBoundaries() const {
		return m_periodicBoundaries;
	}

	const std::map<size_t, shared_ptr<MinLeeBoundary<dim> > >& getMinLeeBoundaries() const {
		return m_minLeeBoundaries;
	}
};

} /* namespace natrium */

#endif /* BOUNDARYCOLLECTION_H_ */
