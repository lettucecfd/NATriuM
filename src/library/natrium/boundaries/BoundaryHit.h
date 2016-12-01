/*
 * BoundaryHit.h
 *
 *  Created on: 22.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYHIT_H_
#define LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYHIT_H_

namespace natrium {

/**
 * @short Defines the point in space and time where a Lagrangian path hit the boundary.
 * 		  Boundary conditions in the Semi-Lagrangian LBM reconstruct the solution at the boundary hits.
 */
template<size_t dim>
struct BoundaryHit {
private:
	LagrangianPathDestination destination; // global degree of freedom at destination (includes direction)
	//size_t beta; // current directions (later: direction of departure degree of freedom (across boundary))
	double dtHit; // the point in time of the boundary hit is (t - dt)
	//dealii::Point<dim> departurePoint; // Lagrangian departure point x^(t-dt)
	dealii::Point<dim> currentPoint; // Current point x^(t-(dt-timeLeft))
	//typename dealii::DoFHandler<dim>::active_cell_iterator currentCell; // cell of currentPoint
	size_t boundaryID;
	typename dealii::DoFHandler<dim>::active_cell_iterator currentCell;

public:
	BoundaryHit(const BoundaryHit<dim>& hit) :
			destination(hit.destination), dtHit(hit.getDtHit()), currentPoint(
					hit.getCurrentPoint()), boundaryID(hit.getBoundaryId()), currentCell(
					hit.getCurrentCell()) {
	}
	BoundaryHit(const LagrangianPathTracker<dim>& tracker,
			const Stencil& stencil, size_t boundary_id) :
			destination(tracker.destination) {

		dtHit = (tracker.currentPoint.distance(tracker.departurePoint)
				/ stencil.getDirection(tracker.beta).l2_norm());
		currentPoint = tracker.currentPoint;
		boundaryID = boundary_id;
		currentCell = tracker.currentCell;
	}
	// TODO Distinguish between BoundaryHits that coincide with support points and those that do not
	//		while the latter require a separate quadrature to obtain the gradients and values at those points
	// 		the former can use the existing quadrature
	// TODO There is a possibility to define the BCs purely at the cell of the destination point
	//		(through Taylor-series expansion in space) -> does that make sense?

	size_t getBoundaryId() const {
		return boundaryID;
	}

	const dealii::Point<dim>& getCurrentPoint() const {
		return currentPoint;
	}

	const LagrangianPathDestination& getDestination() const {
		return destination;
	}

	double getDtHit() const {
		return dtHit;
	}

	const typename dealii::DoFHandler<dim>::active_cell_iterator& getCurrentCell() const {
		return currentCell;
	}
};

/**
 * @short container for boundary hits at a certain point and a certain point in time
 */
template<size_t dim>
struct HitListAtPoint: public std::vector<BoundaryHit<dim> > {
private:
	typedef std::vector<BoundaryHit<dim> > Base;
public:
	void addHit(const BoundaryHit<dim>& hit) {
		Base::push_back(hit);
	}
	size_t n_hits() const {
		return Base::size();
	}
};

template<size_t dim>
struct PointComp {
	bool operator()(const dealii::Point<dim>& lhs,
			const dealii::Point<dim>& rhs) const {
		return lhs.distance(rhs) < 1e-10;
	}
};

/**
 * @short Container for all boundary hits in a cell
 */
template<size_t dim>
struct HitListAtCell: public std::map<dealii::Point<dim>, HitListAtPoint<dim>,
		PointComp<dim> > {
private:
	typedef typename std::map<dealii::Point<dim>, HitListAtPoint<dim>,
			PointComp<dim> > Base;
public:
	typedef typename Base::iterator iterator;
	typedef typename Base::const_iterator const_iterator;
	HitListAtCell() {

	}
	HitListAtCell(const HitListAtCell<dim> & other) :
			Base(other) {
		//hitListSupportPoints = other.hitListSupportPoints;
	}
	virtual ~HitListAtCell() {

	}

	void addHit(const BoundaryHit<dim>& hit) {
		iterator it;
		it = Base::find(hit.getCurrentPoint());
		if (Base::end() == it) {
			HitListAtPoint<dim> hl_p;
			it =
					Base::insert(std::make_pair(hit.getCurrentPoint(), hl_p)).first;
		}
		it->second.addHit(hit);
	}

	size_t n_hits() const {
		size_t result = 0;
		const_iterator it = Base::begin();
		const_iterator end = Base::end();
		for (; it != end; ++it) {
			result += it->second.n_hits();
		}
		return result;
	}

};

/**
 * @short stores all boundary hits that are required for the locally owned dofs
 * Structure of HitList:
 * 	- HitList
 * 		>> HitListAtCell A
 * 			>> HitListAtPoint 1
 * 				>> ...
 * 			>> HitListAtPoint 2
 * 				>> ...
 * 			>> ...
 * 		>> HitListAtCell B
 * 			>> ...
 * 		>> ...
 */
template<size_t dim>
struct HitList: public std::map<
		typename dealii::DoFHandler<dim>::active_cell_iterator,
		HitListAtCell<dim> > {
private:
	typedef std::map<
			typename dealii::DoFHandler<dim>::active_cell_iterator,
			HitListAtCell<dim> > Base;
public:
	typedef typename Base::iterator iterator;
	typedef typename Base::const_iterator const_iterator;
	void addHit(const BoundaryHit<dim>& hit) {
		iterator it;
		it = Base::find(hit.getCurrentCell());
		if (Base::end() == it) {
			HitListAtCell<dim> hl_c;
			typename dealii::DoFHandler<dim>::active_cell_iterator c(
					hit.getCurrentCell());
			it = Base::insert(std::make_pair(c, hl_c)).first;
		}
		it->second.addHit(hit);
	}
	size_t n_hits() const {
		size_t result = 0;
		const_iterator it = Base::begin();
		const_iterator end = Base::end();
		for (; it != end; ++it) {
			result += it->second.n_hits();
		}
		return result;
	}

	size_t n_cells() const {
		return Base::size();
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYHIT_H_ */
