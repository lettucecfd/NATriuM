/**
 * @file MinLeeBoundary.h
 * @short Description of a boundary as described by Min and Lee
 * @date 26.03.2014
 * @author Andreas Kraemer, Bonn-Rhein-Sieg University of Applied Sciences, Sankt Augustin
 */

#ifndef MINLEEBOUNDARY_H_
#define MINLEEBOUNDARY_H_

#include <deal.II/lac/compressed_sparsity_pattern.h>
#include "deal.II/dofs/dof_handler.h"

#include "Boundary.h"
#include "../boltzmannmodels/BoltzmannModel.h"

namespace natrium {

/**
 * @short 	The boundary described by Min and Lee.
 * 			For outgoing particle distribution functions the fluxes are set to 0
 * 		  	For incoming particle distributions fluxes are set to
 * 		  	\f[ f_{\alpha} - f^{+}_{\alpha} = f_{\alpha} - f_{\alpha^{*}} - 2w_{\alpha} \rho_{0} (e_{\alpha}\cdot u_{b})/c^{2}_{s}\f]
 *
 */
template<size_t dim> class MinLeeBoundary: public Boundary<dim> {
private:

	size_t m_boundaryIndicator;
public:

	/// constructor
	MinLeeBoundary(size_t boundaryIndicator);

	/// destructor
	virtual ~MinLeeBoundary(){};

	/**
	 * @short modify sparsity pattern so that the fluxes over periodic boundary can be incorporated
	 * @param cSparse the block-sparsity pattern
	 * @param n_blocks the number of blocks in cSparse
	 * @param n_dofs_per_row number of degrees of freedom per block (normally: overall degrees of freedom on grid)
	 * @param dofs_per_cell number of degrees of freedom per cell
	 */
	void addToSparsityPattern(dealii::BlockCompressedSparsityPattern&  cSparse, const dealii::DoFHandler<dim>& doFHandler, const BoltzmannModel& boltzmannModel) const;


};

} /* namespace natrium */

#endif /* MINLEEBOUNDARY_H_ */
