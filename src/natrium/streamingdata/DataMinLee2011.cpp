/**
 * @file DataMinLee2011.cpp
 * @short 
 * @date 29.05.2013
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#include "DataMinLee2011.h"

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

using namespace dealii;


namespace natrium {

template<size_t dim>
DataMinLee2011<dim>::DataMinLee2011(
		shared_ptr<Triangulation<dim> > triangulation, size_t orderOfFiniteElement):
		m_tria(triangulation),
		m_fe(orderOfFiniteElement),
		m_doFHandler(*triangulation)
{
	// distribute degrees of freedom over mesh
	m_doFHandler.distribute_dofs(m_fe);

	//make sparse matrix
	CompressedSparsityPattern cSparse(m_doFHandler.n_dofs(),
			m_doFHandler.n_dofs());
	//reorder degrees of freedom
	DoFRenumbering::Cuthill_McKee(m_doFHandler);
	DoFTools::make_sparsity_pattern(m_doFHandler, cSparse);
	m_sparsityPattern.copy_from(cSparse);

	//reinitialize system
	//m_systemMatrix.reinit(m_sparsityPattern);
	m_systemRhs.reinit(m_doFHandler.n_dofs());
	//m_solution.reinit(m_doFHandler.n_dofs());

	// assemble system
	reassemble();

}

template<size_t dim>
void DataMinLee2011<dim>::reassemble() {
	// TODO: if Triangulation changed: reinit dof-handler and sparsity pattern in some way




}

} /* namespace natrium */
