#ifndef SEMILAGRANGIANTOOLS_H_
#define SEMILAGRANGIANTOOLS_H_

#include <map>

#include "deal.II/base/point.h"
#include "deal.II/base/tensor.h"
#include "deal.II/grid/tria.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/grid/tria_accessor.h"
#include "../utilities/BasicNames.h"
#include "../utilities/NATriuMException.h"

namespace natrium {

/**
 * @short Exception class for SemiLagrangian advection operator
 */
class FaceCrossedFirstFailed: public NATriuMException {
private:
	std::string message;
public:
	FaceCrossedFirstFailed(const char *msg) :
			NATriuMException(msg), message(msg) {
	}
	FaceCrossedFirstFailed(const string& msg) :
			NATriuMException(msg), message(msg) {
	}
	~FaceCrossedFirstFailed() throw () {
	}
	const char *what() const throw () {
		return this->message.c_str();
	}
};

struct LagrangianPathDestination {
	size_t index;
	size_t direction; // in case of a boundary: the outgoing direction
	LagrangianPathDestination(size_t i, size_t alpha) :
			index(i), direction(alpha) {
	}
	LagrangianPathDestination(const LagrangianPathDestination& other) :
			index(other.index), direction(other.direction) {
	}
};

template<size_t dim>
struct LagrangianPathTracker {

	LagrangianPathDestination destination; // global degree of freedom at destination (includes direction)
	size_t currentDirection; // current directions (later: direction of departure degree of freedom (across boundary))
	size_t lifeTimeCounter = 0; // counts the number of cell face crossings, needed for corner detection
	dealii::Point<dim> departurePoint; // Lagrangian departure point x^(t-dt)
	dealii::Point<dim> currentPoint; // Current point x^(t-(dt-timeLeft))
	typename dealii::DoFHandler<dim>::active_cell_iterator currentCell; // cell of currentPoint
    bool hit_do_nothing = false;


    LagrangianPathTracker(size_t dof, size_t a, size_t b,
			const dealii::Point<dim>& x,
			const dealii::Point<dim>& current_point,
			typename dealii::DoFHandler<dim>::active_cell_iterator current_cell) :
			destination(dof, a), currentDirection(b), departurePoint(x), currentPoint(
					current_point), currentCell(current_cell) {

	}
	LagrangianPathTracker(LagrangianPathDestination& dest, size_t b,
			const dealii::Point<dim>& x,
			const dealii::Point<dim>& current_point,
			typename dealii::DoFHandler<dim>::active_cell_iterator current_cell) :
			destination(dest), currentDirection(b), departurePoint(x), currentPoint(
					current_point), currentCell(current_cell) {

	}
	LagrangianPathTracker(const LagrangianPathTracker& other) :
			destination(other.destination), currentDirection(
					other.currentDirection), departurePoint(
					other.departurePoint), currentPoint(other.currentPoint), currentCell(
					other.currentCell) {
	}
	LagrangianPathTracker& operator=(const LagrangianPathTracker& other) {
		destination = other.destination;
		currentDirection = other.currentDirection;
		departurePoint = other.departurePoint;
		currentPoint = other.currentPoint;
		currentCell = other.currentCell;
		return *this;
	}
};

/**
 * @short A list that stores cell-specific information for assembly
 */
template<size_t dim>
using DeparturePointList = std::vector<LagrangianPathTracker<dim> >;

/**
 * @short List of neighbors
 */
template<size_t dim>
using Neighborhood = std::vector<typename dealii::DoFHandler<dim>::cell_iterator>;

/**
 * @short Calculates the shape values for arbitrary points.
 * @param[in] cell the cell which contains the points
 * @param[in] a vector of points
 * @param[in] a vector of vector<double> that will contain the shape function values at the points.
 * 				Has to have the size n_points x dofs_per_cell
 * @note This function is similar to FEFieldFunction<dim, DH, VECTOR>::vector_value
 */
template<size_t dim>
void shapeFunctionValue(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const std::vector<dealii::Point<dim> >& points,
		std::vector<std::vector<double> >&values,
		const dealii::Mapping<dim>& mapping);

template<size_t dim>
int supportPointNr(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const dealii::Point<dim>& point, const dealii::Quadrature<dim>& quad,
		const dealii::MappingQ<dim,dim>& mapping);

/**
 * @short Create an FEValues object to calculate shape functions, etc. on other points than the support points
 * @note the returned object is already reinitialized for the cell, i.e. you don't need need to call reinit(cell) for the returned object.
 */
template<size_t dim>
dealii::Quadrature<dim> makeQuadratureAtPoints(
		const typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		const std::vector<dealii::Point<dim> >& points,
		const dealii::Mapping<dim>& mapping);

template<size_t dim>
void getNeighborhood(
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell,
		Neighborhood<dim>& neighborhood, size_t n_shells = 1);

template<size_t dim>
typename dealii::DoFHandler<dim>::active_cell_iterator recursivelySearchInNeighborhood(
		const dealii::Point<dim>& p,
		typename dealii::DoFHandler<dim>::active_cell_iterator& cell);

template<size_t dim>
dealii::Tensor<1, dim, double> normal_vector(
		const typename dealii::TriaIterator<
				dealii::TriaAccessor<dim - 1, dim, dim> >& face);
// forward declare explicit specializations
template<>
dealii::Tensor<1, 2, double> normal_vector<2>(
		const typename dealii::TriaIterator<dealii::TriaAccessor<1, 2, 2> >& face) ;
template<>
dealii::Tensor<1, 3, double> normal_vector<3>(
		const typename dealii::TriaIterator<dealii::TriaAccessor<2, 3, 3> >& face) ;

} /* namespace natrium */



#endif /* includeguard*/
