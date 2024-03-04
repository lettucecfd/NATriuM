
#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/utilities/BasicNames.h"

using namespace natrium;


class EnstrophySubdomain: public DataProcessor<2> {
private:
	double m_result;
public:
	EnstrophySubdomain(const CFDSolver<2> & solver) :
			DataProcessor<2>(solver), m_result(0) {
	}
	virtual ~EnstrophySubdomain() {
	}
	virtual void apply() {
		if (m_solver.getIteration() % m_solver.getConfiguration()->getOutputTableInterval() == 0){
			calculate();
			pout << m_solver.getTime() << " " << m_result << endl;
		}
	}
	void calculate(){
		const vector<distributed_vector>& u = m_solver.getVelocity();
		boost::shared_ptr<AdvectionOperator<2> > advection = m_solver.getAdvectionOperator();
		const size_t n_dofs = u.at(0).size();
		assert(n_dofs == u.at(1).size());

		const distributed_vector& ux = u.at(0);
		const distributed_vector& uy = u.at(1);

		// Integrate ux over whole domain
		const dealii::UpdateFlags cellUpdateFlags = dealii::update_JxW_values
				| dealii::update_gradients;
		const dealii::DoFHandler<2> & dof_handler =
				*(advection->getDoFHandler());
		dealii::FEValues<2> feCellValues(advection->getMapping(),
				*(advection->getFe()), *(advection->getQuadrature()),
				cellUpdateFlags);
		double result = 0.0;
		double vorticity = 0.0;
		size_t local_i;
		size_t q_point;
		size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
		size_t n_q_points = advection->getQuadrature()->size();

		typename dealii::DoFHandler<2>::active_cell_iterator cell =
				dof_handler.begin_active(), endc = dof_handler.end();
		for (; cell != endc; ++cell) {
			if (cell->is_locally_owned()) {
				// restrict to lower right subsquare of domain
				if (cell->center()(0) < 0.5)
					continue;
				if (cell->center()(1) > 0.5)
					continue;


				// get global degrees of freedom
				std::vector<dealii::types::global_dof_index> localDoFIndices(
						dofs_per_cell);
				cell->get_dof_indices(localDoFIndices);
				// calculate the fe values for the cell
				feCellValues.reinit(cell);

				// Enstrophy = int w^2 = int ( dv/dx - du/dy)^2  = w_q ( v_i dphi/dx (x_q) - u_j dphi/dy (x_q) )
				for (size_t q = 0; q < n_q_points; q++) {
					vorticity = 0.0;
					for (size_t i = 0; i < dofs_per_cell; i++) {
						local_i = localDoFIndices.at(i);
						vorticity +=
								(uy(local_i)
										* feCellValues.shape_grad(i, q_point)[0]
										- ux(local_i)
												* feCellValues.shape_grad(i,
														q_point)[1]);
					}
					result += vorticity * vorticity * feCellValues.JxW(q_point);
				}
			} /* if is locally owned */
		} /* for cells */

		// communicate among MPI processes
		dealii::Utilities::MPI::MinMaxAvg global_res =
				dealii::Utilities::MPI::min_max_avg(result, MPI_COMM_WORLD);
		// the enstrophy can be defined with and without factor 1/2; here: without
		m_result = global_res.sum;

	}

	double getResult() {
		calculate();
		return m_result;
	}
};
