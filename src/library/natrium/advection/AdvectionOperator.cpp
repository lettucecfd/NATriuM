/*
 * AdvectionOperator.cpp
 *
 *  Created on: 26.11.2016
 *      Author: akraem3m
 */

#include "AdvectionOperator.h"
#include "AdvectionTools.h"

#include "deal.II/fe/mapping.h"
#include "deal.II/fe/mapping_q1.h"
#include "deal.II/fe/mapping_cartesian.h"

#include "../problemdescription/ProblemDescription.h"
#include "../stencils/Stencil.h"

namespace natrium {

/*
 template<size_t dim>
 SEDGMinLee<dim>::SEDGMinLee(boost::shared_ptr<Mesh<dim> > triangulation,
 boost::shared_ptr<BoundaryCollection<dim> > boundaries,
 size_t orderOfFiniteElement, boost::shared_ptr<Stencil> Stencil,
 bool useCentralFlux) :
 m_mesh(triangulation), m_boundaries(boundaries), m_mapping(
 orderOfFiniteElement), m_stencil(Stencil), m_orderOfFiniteElement(
 orderOfFiniteElement), m_useCentralFlux(useCentralFlux) {
 // assertions
 assert(orderOfFiniteElement >= 1);
 assert(Stencil->getD() == dim);

 // make dof handler
 m_quadrature = boost::make_shared<QGaussLobatto<dim> >(
 orderOfFiniteElement + 1);
 m_faceQuadrature = boost::make_shared<QGaussLobatto<dim - 1> >(
 orderOfFiniteElement + 1);
 m_fe = boost::make_shared<FE_DGQArbitraryNodes<dim> >(
 QGaussLobatto<1>(orderOfFiniteElement + 1));
 m_doFHandler = boost::make_shared<DoFHandler<dim> >(*triangulation);

 */
template<size_t dim>
AdvectionOperator<dim>::AdvectionOperator(ProblemDescription<dim>& problem,
		size_t fe_order, QuadratureName quad_name,
		SupportPointsName points_name, boost::shared_ptr<Stencil> stencil,
		double delta_t, bool dg) :
		m_problem(problem), m_stencil(stencil), m_orderOfFiniteElement(
				fe_order), m_deltaT(delta_t) {

	// assertions
	assert(fe_order >= 1);
	assert(stencil->getD() == dim);
	assert(delta_t >= 0);
	assert(dim > 1);
	assert(dim < 4);

	// make dof handler
	if (problem.isCartesian() and !dg) {
		m_mapping = boost::make_shared<dealii::MappingCartesian<dim> >();
	} else {
		m_mapping = boost::make_shared<dealii::MappingQ1<dim> >();
		if (!dg)
			LOG(WARNING)
					<< "Semi-Lagrangian advection on non-Cartesian grid is not well tested, yet."
							"We encountered problems with deal.II's multilinear Q1-mapping, when "
							"the length of the computational domain was something like 10.3 or 2pi."
							"You should at least make sure that your domain dimensions are integers"
							"to improve your chances of getting a stable simulation."
							"Check the mass conservation in the results table to make sure that your simulation"
							"is not going badly wrong." << endl;
	}
	m_quadrature = AdvectionTools::make_quadrature_by_name<dim>(quad_name,
			fe_order + 1);
	m_faceQuadrature = AdvectionTools::make_quadrature_by_name<dim - 1>(
			quad_name, fe_order + 1);
	m_fe = AdvectionTools::make_fe_by_name<dim>(points_name, fe_order, dg);
	vector<double> weights(m_fe->dofs_per_cell,
			1. / ((double) m_fe->dofs_per_cell));
	m_supportPointEvaluation = boost::make_shared<dealii::Quadrature<dim> >(
			m_fe->get_unit_support_points(), weights);
	m_doFHandler = boost::make_shared<dealii::DoFHandler<dim> >(*getMesh());

} /* constructor */

template class AdvectionOperator<2> ;
template class AdvectionOperator<3> ;

} /* namespace natrium */

