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
	size_t dtHit; // the point in time of the boundary hit is (t - dt)
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

	size_t getDtHit() const {
		return dtHit;
	}

	const typename dealii::DoFHandler<dim>::active_cell_iterator& getCurrentCell() const {
		return currentCell;
	}
};

/**
 * @short Special case of BoundaryHit, where the boundary hit coincides with a support point.
 */
/*template<size_t dim>
 struct BoundaryHitAtSupportPoint {
 private:
 LagrangianPathDestination destination; // global degree of freedom at destination (includes direction)
 //size_t beta; // current directions (later: direction of departure degree of freedom (across boundary))
 size_t dtHit; // the point in time of the boundary hit is (t - dt)
 //dealii::Point<dim> departurePoint; // Lagrangian departure point x^(t-dt)
 //dealii::Point<dim> currentPoint; // Current point x^(t-(dt-timeLeft))
 size_t supportQPoint; // id of the quadrature point
 size_t boundaryID;
 //typename dealii::DoFHandler<dim>::active_cell_iterator currentCell; // cell of currentPoint

 public:
 BoundaryHitAtSupportPoint(const BoundaryHitAtSupportPoint<dim>& hit) :
 destination(hit.destination) {
 dtHit = hit.dtHit;
 supportQPoint = hit.supportQPoint;
 boundaryID = hit.boundaryID;
 }
 BoundaryHitAtSupportPoint(const LagrangianPathTracker<dim>& tracker,
 const Stencil&, size_t boundary_id, size_t q_point) :
 destination(tracker.destination) {
 //assert(tracker.currentPoint.distance(tracker.departurePoint) < 1e-12);
 dtHit = 0;
 supportQPoint = q_point;
 boundaryID = boundary_id;
 }
 // TODO Distinguish between BoundaryHits that coincide with support points and those that do not
 //		while the latter require a separate quadrature to obtain the gradients and values at those points
 // 		the former can use the existing quadrature
 // TODO There is a possibility to define the BCs purely at the cell of the destination point
 //		(through Taylor-series expansion in space) -> does that make sense?

 size_t getBoundaryId() const {
 return boundaryID;
 }

 const LagrangianPathDestination& getDestination() const {
 return destination;
 }

 size_t getDtHit() const {
 return dtHit;
 }

 size_t getSupportQPoint() const {
 return supportQPoint;
 }
 };
 */

/**
 * @short container for boundary hits at a certain point and a certain point in time
 */
template<size_t dim>
struct HitListAtPointAndTime {
	std::vector<BoundaryHit<dim> > m_hitlist;
	const std::vector<BoundaryHit<dim> >& getHitlist() const {
		return m_hitlist;
	}
	void addHit(const BoundaryHit<dim>& hit) {
		m_hitlist.push_back(hit);
	}
	size_t n_hits() const {
		return m_hitlist.size();
	}
};

struct DoubleComp {
	bool operator()(const double& lhs, const double& rhs) const {
		return fabs(lhs - rhs) < 1e-10;
	}
};

/**
 * @short container for boundary hits at a certain point
 */
template<size_t dim>
struct HitListAtPoint {
private:
	std::map<double, HitListAtPointAndTime<dim>, DoubleComp> m_hitlist;
public:
	typedef typename std::map<double, HitListAtPointAndTime<dim>, DoubleComp>::iterator iterator;
	typedef typename std::map<double, HitListAtPointAndTime<dim>, DoubleComp>::const_iterator const_iterator;
	const std::map<double, HitListAtPointAndTime<dim>, DoubleComp>& getHitlist() const {
		return m_hitlist;
	}
	void addHit(const BoundaryHit<dim>& hit) {
		iterator it;
		it = m_hitlist.find(hit.getDtHit());
		if (m_hitlist.end() == it) {
			HitListAtPointAndTime<dim> hl_p;
			it = m_hitlist.insert(std::make_pair(hit.getDtHit(), hl_p)).first;
		}
		it->second.addHit(hit);
	}
	size_t n_hits() const {
		size_t result = 0;
		const_iterator it = m_hitlist.begin();
		const_iterator end = m_hitlist.end();
		for (; it != end; ++it) {
			result += it->second.n_hits();
		}
		return result;
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
struct HitListAtCell {
private:
	std::map<dealii::Point<dim>, HitListAtPoint<dim>, PointComp<dim> > m_hitlist;
public:
	typedef typename std::map<dealii::Point<dim>, HitListAtPoint<dim>,
			PointComp<dim> >::iterator iterator;
	typedef typename std::map<dealii::Point<dim>, HitListAtPoint<dim>,
			PointComp<dim> >::const_iterator const_iterator;
	HitListAtCell() {

	}
	HitListAtCell(const HitListAtCell<dim> & other) :
			m_hitlist(other.getHitlist()) {
		//hitListSupportPoints = other.hitListSupportPoints;
	}
	virtual ~HitListAtCell() {

	}

	const std::map<dealii::Point<dim>, HitListAtPoint<dim> >& getHitList() {
		return m_hitlist;
	}

	void addHit(const BoundaryHit<dim>& hit) {
		iterator it;
		it = m_hitlist.find(hit.getCurrentPoint());
		if (m_hitlist.end() == it) {
			HitListAtPoint<dim> hl_p;
			it =
					m_hitlist.insert(
							std::make_pair(hit.getCurrentPoint(), hl_p)).first;
		}
		it->second.addHit(hit);
	}

	size_t n_hits() const {
		size_t result = 0;
		const_iterator it = m_hitlist.begin();
		const_iterator end = m_hitlist.end();
		for (; it != end; ++it) {
			result += it->second.n_hits();
		}
		return result;
	}

	const std::map<dealii::Point<dim>, HitListAtPoint<dim>, PointComp<dim> >& getHitlist() const {
		return m_hitlist;
	}
};

/**
 * @short stores all boundary hits that are required for the locally owned dofs
 * Structure of HitList:
 * 	- HitList
 * 		>> HitListAtCell A
 * 			>> HitListAtPoint 1
 * 				>> HitListAtPointAndTime 1,t
 * 				>> HitListAtPointAndTime 1,t
 * 				>> ...
 * 			>> HitListAtPoint 2
 * 				>> ...
 * 			>> ...
 * 		>> HitListAtCell B
 * 			>> ...
 * 		>> ...
 */
template<size_t dim>
struct HitList {
private:
	std::map<typename dealii::DoFHandler<dim>::active_cell_iterator,
			HitListAtCell<dim> > m_hitlist;
public:
	typedef typename std::map<
			typename dealii::DoFHandler<dim>::active_cell_iterator,
			HitListAtCell<dim> >::iterator iterator;
	typedef typename std::map<
			typename dealii::DoFHandler<dim>::active_cell_iterator,
			HitListAtCell<dim> >::const_iterator const_iterator;
	const std::map<typename dealii::DoFHandler<dim>::active_cell_iterator,
			HitListAtCell<dim> >& getHitlist() const {
		return m_hitlist;
	}
	void addHit(const BoundaryHit<dim>& hit) {
		iterator it;
		it = m_hitlist.find(hit.getCurrentCell());
		if (m_hitlist.end() == it) {
			HitListAtCell<dim> hl_c;
			typename dealii::DoFHandler<dim>::active_cell_iterator c(
					hit.getCurrentCell());
			it = m_hitlist.insert(std::make_pair(c, hl_c)).first;
		}
		it->second.addHit(hit);
	}
	size_t n_hits() const {
		size_t result = 0;
		const_iterator it = m_hitlist.begin();
		const_iterator end = m_hitlist.end();
		for (; it != end; ++it) {
			result += it->second.n_hits();
		}
		return result;
	}

	size_t n_cells() const {
		return m_hitlist.size();
	}
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_BOUNDARYHIT_H_ */
