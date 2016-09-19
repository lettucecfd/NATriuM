/*
 * BoundaryCollection.h
 *
 *  Created on: 29.10.2013
 *      Author: kraemer
 */

#ifndef BOUNDARYCOLLECTION_H_
#define BOUNDARYCOLLECTION_H_

#include "deal.II/lac/constraint_matrix.h"

#include "../boundaries/Boundary.h"
#include "../boundaries/PeriodicBoundary.h"
#include "../boundaries/LinearFluxBoundary.h"
#include "../boundaries/DoFBoundary.h"
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"
#include "../solver/DistributionFunctions.h"

namespace natrium {

/**
 * @short Exception class for Boundaries
 */
class BoundaryCollectionException: public NATriuMException {
private:
	std::string message;
public:
	BoundaryCollectionException(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	BoundaryCollectionException(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~BoundaryCollectionException() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

/**
 * @short The BoundaryCollection class is a container for all boundaries of a flow domain.
 *        Internally, the boundaries are stored in two different std::map, one
 *        for the periodic boundaries and one for the non-periodic ones.
 *        Its keys are the boundary indicators (for periodic boundaries: the first boundary indicator).
 */
template<size_t dim> class BoundaryCollection {
private:
	/// vector to store boundaries in
	std::map<size_t, boost::shared_ptr<Boundary<dim> > > m_boundaries;

	/// vector to store linear boundaries in
	std::map<size_t, boost::shared_ptr<LinearFluxBoundary<dim> > > m_linearBoundaries;

	/// vector to store dof boundaries in
	std::map<size_t, boost::shared_ptr<DoFBoundary<dim> > > m_dofBoundaries;

	/// vector to store periodic boundaries in
	std::map<size_t, boost::shared_ptr<PeriodicBoundary<dim> > > m_periodicBoundaries;

public:

	typedef typename std::map<size_t, boost::shared_ptr<Boundary<dim> > >::iterator Iterator;
	typedef typename std::map<size_t, boost::shared_ptr<Boundary<dim> > >::const_iterator ConstIterator;
	typedef typename std::map<size_t, boost::shared_ptr<LinearFluxBoundary<dim> > >::iterator LinearIterator;
	typedef typename std::map<size_t, boost::shared_ptr<LinearFluxBoundary<dim> > >::const_iterator ConstLinearIterator;
	typedef typename std::map<size_t, boost::shared_ptr<DoFBoundary<dim> > >::iterator DoFBoundaryIterator;
	typedef typename std::map<size_t, boost::shared_ptr<DoFBoundary<dim> > >::const_iterator ConstDoFBoundaryIterator;
	typedef typename std::map<size_t, boost::shared_ptr<PeriodicBoundary<dim> > >::iterator PeriodicIterator;
	typedef typename std::map<size_t, boost::shared_ptr<PeriodicBoundary<dim> > >::const_iterator ConstPeriodicIterator;

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
	void addBoundary(boost::shared_ptr<PeriodicBoundary<dim> > boundary) {
		bool success1 =
				m_boundaries.insert(
						std::make_pair(boundary->getBoundaryIndicator1(),
								boundary)).second;
		bool success2 =
				m_boundaries.insert(
						std::make_pair(boundary->getBoundaryIndicator2(),
								boundary)).second;
		if ((not success1) or (not success2)) {
			throw BoundaryCollectionException(
					"Boundary could not be inserted. Boundary indicators must be unique.");
		}
		m_periodicBoundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator1(), boundary));
		m_periodicBoundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator2(), boundary));
	}

	/**
	 * @short Add a boundary to the flow definition.
	 * @param boundary a periodic boundary
	 * @throws BoundaryCollectionError, e.g. if boundary indicators are not unique
	 */
	void addBoundary(boost::shared_ptr<LinearFluxBoundary<dim> > boundary) {
		bool success =
				m_boundaries.insert(
						std::make_pair(boundary->getBoundaryIndicator(),
								boundary)).second;
		if (not success) {
			throw BoundaryCollectionException(
					"Boundary could not be inserted. Boundary indicators must be unique.");
		}
		m_linearBoundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator(), boundary));
	}

	/**
	 * @short Add a boundary to the flow definition.
	 * @param boundary a periodic boundary
	 * @throws BoundaryCollectionError, e.g. if boundary indicators are not unique
	 */
	void addBoundary(boost::shared_ptr<DoFBoundary<dim> > boundary) {
		bool success =
				m_boundaries.insert(
						std::make_pair(boundary->getBoundaryIndicator(),
								boundary)).second;
		if (not success) {
			throw BoundaryCollectionException(
					"Boundary could not be inserted. Boundary indicators must be unique.");
		}
		m_dofBoundaries.insert(
				std::make_pair(boundary->getBoundaryIndicator(), boundary));
	}

	/**
	 * @short get a specific boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const boost::shared_ptr<Boundary<dim> >& getBoundary(
			size_t boundaryIndicator) const {
		if (m_boundaries.count(boundaryIndicator) == 0) {
			throw BoundaryCollectionException(
					"in getBoundary: This boundary collection does not contain a boundary with the specified boundary indicator.");
		}
		return m_boundaries.at(boundaryIndicator);
	}

	/**
	 * @short get a specific periodic boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const boost::shared_ptr<PeriodicBoundary<dim> >& getPeriodicBoundary(
			size_t boundaryIndicator) const {
		assert(isPeriodic(boundaryIndicator));
		if (m_periodicBoundaries.count(boundaryIndicator) == 0) {
			throw BoundaryCollectionException(
					"in getPeriodicBoundary: This boundary collection does not contain a periodic boundary with the specified boundary indicator.");
		}
		return m_periodicBoundaries.at(boundaryIndicator);
	}

	/**
	 * @short get a specific MinLee boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const boost::shared_ptr<LinearFluxBoundary<dim> >& getLinearFluxBoundary(
			size_t boundaryIndicator) const {
		assert(not isPeriodic(boundaryIndicator));
		if (m_linearBoundaries.count(boundaryIndicator) == 0) {
			throw BoundaryCollectionException(
					"in LinearBoundary: This boundary collection does not contain a linear boundary with the specified boundary indicator.");
		}
		return m_linearBoundaries.at(boundaryIndicator);
	}

	/**
	 * @short get a specific dof boundary
	 * @throws BoundaryCollectionError, if the specified boundary indicator does not exist
	 */
	const boost::shared_ptr<DoFBoundary<dim> >& getDoFBoundary(
			size_t boundaryIndicator) const {
		assert(not isPeriodic(boundaryIndicator));
		if (m_dofBoundaries.count(boundaryIndicator) == 0) {
			throw BoundaryCollectionException(
					"in getDoFBoundary: This boundary collection does not contain a dof boundary with the specified boundary indicator.");
		}
		return m_dofBoundaries.at(boundaryIndicator);
	}

	/**
	 * @short test if the boundary with the given boundary indicator is periodic
	 */
	bool isPeriodic(size_t boundaryIndicator) const {
		return getBoundary(boundaryIndicator)->isPeriodic();
	}

	size_t numberOfBoundaries() const {
		return m_boundaries.size();
	}

	size_t numberOfPeriodicBoundaries() const {
		return m_periodicBoundaries.size();
	}

	const std::map<size_t, boost::shared_ptr<Boundary<dim> > >& getBoundaries() const {
		return m_boundaries;
	}

	const std::map<size_t, boost::shared_ptr<PeriodicBoundary<dim> > >& getPeriodicBoundaries() const {
		return m_periodicBoundaries;
	}

	const std::map<size_t, boost::shared_ptr<LinearFluxBoundary<dim> > >& getLinearFluxBoundaries() const {
		return m_linearBoundaries;
	}

	const std::map<size_t, boost::shared_ptr<DoFBoundary<dim> > >& getDoFBoundaries() const {
		return m_dofBoundaries;
	}

	bool hasDoFBoundaries(){
		return m_dofBoundaries.size() != 0;
	}

};

} /* namespace natrium */

#endif /* BOUNDARYCOLLECTION_H_ */
