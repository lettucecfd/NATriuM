/*
 * SemiLagrangianBoundaryHandler.h
 *
 *  Created on: 15.11.2016
 *      Author: akraem3m
 */

#ifndef LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_
#define LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_

#include <map>
#include <utility>

#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/base/point.h"

#include "../advection/SemiLagrangianTools.h"
#include "../problemdescription/BoundaryCollection.h"
#include "../utilities/BasicNames.h"

namespace natrium {

// forward declaration
template<size_t dim>
class AdvectionOperator;

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

public:
	BoundaryHit(const BoundaryHit<dim>& hit) :
			destination(hit.destination) {
		dtHit = hit.dtHit;
		currentPoint = hit.currentPoint;
		boundaryID = hit.boundaryID;
	}
	BoundaryHit(const LagrangianPathTracker<dim>& tracker,
			const Stencil& stencil, size_t boundary_id) :
			destination(tracker.destination) {

		dtHit = (tracker.currentPoint.distance(tracker.departurePoint)
				/ stencil.getDirection(tracker.beta).l2_norm());
		currentPoint = tracker.currentPoint;
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

	const dealii::Point<dim>& getCurrentPoint() const {
		return currentPoint;
	}

	const LagrangianPathDestination& getDestination() const {
		return destination;
	}

	size_t getDtHit() const {
		return dtHit;
	}
};

/**
 * @short Special case of BoundaryHit, where the boundary hit coincides with a support point.
 */
template<size_t dim>
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
		assert(tracker.currentPoint.distance(tracker.departurePoint) < 1e-12);
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

/**
 * @short Container for all boundary hits in a cell
 */
template<size_t dim>
struct HitListAtCell {
public:
	HitListAtCell() {

	}
	HitListAtCell(const HitListAtCell<dim> & other) {
		hitListArbitrary = other.hitListArbitrary;
		hitListSupportPoints = other.hitListSupportPoints;
	}
	virtual ~HitListAtCell(){

	}
/*	HitListAtCell<dim> & operator= (const HitListAtCell<dim> & other) {
		hitListArbitrary = other.hitListArbitrary;
		hitListSupportPoints = other.hitListSupportPoints;
		return *this;
	}*/

	std::vector<BoundaryHit<dim> > hitListArbitrary;
	std::vector<BoundaryHitAtSupportPoint<dim> > hitListSupportPoints;
	void addHit(const BoundaryHitAtSupportPoint<dim>& hit) {
		hitListSupportPoints.push_back(hit);
	}
	void addHit(const BoundaryHit<dim>& hit) {
		hitListArbitrary.push_back(hit);
	}
};

/**
 * @short stores all boundary hits that are required for the locally owned dofs
 */
template<size_t dim>
using HitList = std::map<typename dealii::DoFHandler<dim>::active_cell_iterator, HitListAtCell<dim> >;


/**
 * @short Handler for semi-Lagrangian boundary conditions.
 */
template<size_t dim>
class SemiLagrangianBoundaryHandler {
private:
	HitList<dim> m_hitList;
	double m_timeStep;
	const Stencil& m_stencil;
	const BoundaryCollection<dim>& m_boundaries;
public:
	SemiLagrangianBoundaryHandler(double dt, const Stencil& stencil, const BoundaryCollection<dim>& boundaries) :
			m_timeStep(dt), m_stencil(stencil), m_boundaries(boundaries) {

	}
	virtual ~SemiLagrangianBoundaryHandler() {

	}
	void addHit(const LagrangianPathTracker<dim>& tracker, size_t boundary_id, const AdvectionOperator<dim>& sl) ;

	void apply(DistributionFunctions& f_new, const DistributionFunctions& f_old,
			const dealii::DoFHandler<dim>& dof,
			boost::shared_ptr<dealii::FEValues<dim> > fe_values);

	/*void print_out(){
		pout << "Semi Lagrangian Boundary Handler Object:" << endl;
		m_hitList.print_out();
		pout << "----------------------------------------" << endl;
		pout << "#Hits at support points: " << endl;
	}*/
};

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_BOUNDARIES_SEMILAGRANGIANBOUNDARYHANDLER_H_ */
