/*
 * ContinuousBoundaryGradient.h
 *
 *  Created on: 28.09.2016
 *      Author: Andreas Kraemer
 */

#ifndef LIBRARY_NATRIUM_UTILITIES_CONTINUOUSBOUNDARYGRADIENT_H_
#define LIBRARY_NATRIUM_UTILITIES_CONTINUOUSBOUNDARYGRADIENT_H_

#include <vector>
#include <set>

#include "boost/shared_ptr.hpp"
#include "boost/make_shared.hpp"


#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_system.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/fe/mapping.h"
#include "deal.II/lac/trilinos_vector.h"
#include "deal.II/fe/fe_q.h"
#include "deal.II/fe/fe_dgq.h"

// using hp::DoFHandler would be much more appropriate
// but it is not implemented to work with parallel::distributed in deal.II 8.4.2.
// so we switch to the normal DoFHandler
// If one day the other solution is available, just comment this line out ;-)
#define NO_SUPPORT_HP_PARALLEL

#ifdef NO_SUPPORT_HP_PARALLEL
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/fe/fe_values.h"
using dealii::DoFHandler;
using dealii::FEFaceValues;

#else
#include "deal.II/hp/dof_handler.h"
#include "deal.II/hp/fe_collection.h"
#include "deal.II/hp/fe_values.h"
using dealii::hp::DoFHandler;
using dealii::hp::FEFaceValues;
#endif


namespace natrium {

/**
 * @short A class to obtain gradients at a boundary that are continuous over the cell interfaces.
 */
template<size_t dim>
class ContinuousBoundaryGradient {
private:
	const dealii::DoFHandler<dim>& m_dof;
	const dealii::Mapping<dim>& m_mapping;
	const dealii::Quadrature<1> m_quadrature1D;
	const dealii::Quadrature<dim - 1> m_faceQuadrature;
	const dealii::Quadrature<dim> m_cellQuadrature;

	std::set<size_t> m_boundaryIds;
#ifdef NO_SUPPORT_HP_PARALLEL
	dealii::FE_DGQArbitraryNodes<dim> m_FeDg;
	dealii::FE_Q<dim> m_FeConti;
#else
	dealii::hp::FECollection<dim> m_FeDgCollection;
	dealii::hp::FECollection<dim> m_FeContiCollection;
#endif
	boost::shared_ptr<DoFHandler<dim> > m_discontinuousBoundaryDoF;
	boost::shared_ptr<DoFHandler<dim> > m_continuousBoundaryDoF;
	boost::shared_ptr<FEFaceValues<dim> > m_feFaceValues;
	std::vector<dealii::types::global_dof_index> m_dofIndizes;
	dealii::TrilinosWrappers::MPI::Vector m_discontinuousGradient;
	dealii::TrilinosWrappers::MPI::Vector m_continuousGradient;


public:
	/**
	 * constructor
	 */
	ContinuousBoundaryGradient(const dealii::DoFHandler<dim>& dof,
			const dealii::Mapping<dim>& mapping,
			const dealii::Quadrature<1>& quadrature,
			const dealii::Quadrature<dim-1>& face_quadrature,
			const dealii::Quadrature<dim>& cell_quadrature,
			std::set<size_t> boundary_ids);
	virtual ~ContinuousBoundaryGradient();

	typedef typename dealii::DoFHandler<dim>::active_cell_iterator active_cell_iterator;
	typedef typename dealii::DoFHandler<dim>::cell_iterator cell_iterator;

	/**
	 * @short reinitialize the data structures of this boundary gradient
	 */
	void reinit();

	/**
	 * @short Calculate the gradients of the solution
	 */
	void calculateGradients(
			const dealii::TrilinosWrappers::MPI::Vector& solution);

	/**
	 * @short first active cell of this gradient
	 */
	active_cell_iterator begin_active() const {
		return m_continuousBoundaryDoF->begin_active();
	}

	/**
	 * @short last active cell of this gradient
	 */
	cell_iterator end() const {
		return m_continuousBoundaryDoF->end();
	}

	/**
	 * @short Reinit the face values
	 * @note The cell iterator given to this functions has to be a different one than that of your regular dof handler.
	 * 		 In order to hand this function a valid cell iterator, you have to initialize one with begin_active() and
	 * 		 increment it in line with your other cell iterator. The cells will come in the same order.
	 */
	void face_reinit(
			const typename dealii::hp::DoFHandler<dim>::active_cell_iterator& cell,
			size_t face_id) {
		assert(cell->is_locally_owned());
		assert(face_id < dealii::GeometryInfo<dim>::faces_per_cell);
		assert(cell->active_fe_index() == 1);
		m_feFaceValues->reinit(cell, face_id);
		cell->get_dof_indices(m_dofIndizes);
	}

	/**
	 * @short get component i of the (continuous) gradient at a quadrature point
	 */
	double get_gradient_component(size_t q_point, size_t component);

	const dealii::TrilinosWrappers::MPI::Vector& getContinuousGradient() const {
		return m_continuousGradient;
	}

	const dealii::TrilinosWrappers::MPI::Vector& getDiscontinuousGradient() const {
		return m_discontinuousGradient;
	}
} /* class */
;

} /* namespace natrium */

#endif /* LIBRARY_NATRIUM_UTILITIES_CONTINUOUSBOUNDARYGRADIENT_H_ */
