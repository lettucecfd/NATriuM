/**
 * @file DataMinLee2011.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DataMinLee2011.h"

#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

using namespace dealii;


namespace natrium {

template<size_t dim>
DataMinLee2011<dim>::DataMinLee2011(
		shared_ptr<Triangulation<dim> > triangulation, size_t orderOfFiniteElement):
		m_tria(triangulation)
{
	// make dof handler
	m_fe = make_shared<FE_Q<dim> >(orderOfFiniteElement);
	m_doFHandler = make_shared<DoFHandler<dim> >(*triangulation);

	// distribute degrees of freedom over mesh
	m_doFHandler->distribute_dofs(*m_fe);

	//make sparse matrix
	CompressedSparsityPattern cSparse(m_doFHandler->n_dofs());
	//reorder degrees of freedom
	DoFRenumbering::Cuthill_McKee(*m_doFHandler);
	DoFTools::make_sparsity_pattern(*m_doFHandler, cSparse);
	m_sparsityPattern.copy_from(cSparse);

	//reinitialize system
	//m_systemMatrix.reinit(m_sparsityPattern);
	//m_systemRhs.reinit(m_doFHandler->n_dofs());
	//m_solution.reinit(m_doFHandler.n_dofs());

	// assemble system
	reassemble();

} /* DataMinLee2011<dim>::DataMinLee2011 */
/// The template parameter must be made explicit in order for the code to compile
template DataMinLee2011<2>::DataMinLee2011(
		shared_ptr<Triangulation<2> > triangulation, size_t orderOfFiniteElement);
template DataMinLee2011<3>::DataMinLee2011(
		shared_ptr<Triangulation<3> > triangulation, size_t orderOfFiniteElement);


template<size_t dim>
void DataMinLee2011<dim>::stream() {
}
// The template parameter must be made explicit in order for the code to compile.
template void DataMinLee2011<2>::stream();
template void DataMinLee2011<3>::stream();



template<size_t dim>
void DataMinLee2011<dim>::reassemble() {
	// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way

}
/// The template parameter must be made explicit in order for the code to compile
template void DataMinLee2011<2>::reassemble();
template void DataMinLee2011<3>::reassemble();

} /* namespace natrium */
