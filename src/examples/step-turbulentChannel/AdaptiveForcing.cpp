/*
 * AdaptiveForcing.cpp
 *
 *  Created on: 06.10.2022
 *      Author: dwilde3m
 */

#include "AdaptiveForcing.h"

#include "mpi.h"
#include <utility>
#include "deal.II/base/mpi.h"

namespace natrium {

AdaptiveForcing::AdaptiveForcing(CFDSolver<3> & solver,
		std::string outdir, double target) :
		DataProcessor<3>(solver), m_outDir(outdir), m_u(solver.getVelocity()), m_rho(
				solver.getDensity()), m_targetRhoU(target),
                m_filename(outfile(solver.getConfiguration()->getOutputDirectory())) {

        if (solver.getIterationStart() > 0) {
            if (is_MPI_rank_0()) {
                m_tableFile = boost::make_shared<std::fstream>(m_filename,
                                                               std::fstream::out | std::fstream::app);
            }
        } else {
            if (is_MPI_rank_0()) {
                m_tableFile = boost::make_shared<std::fstream>(m_filename,
                                                               std::fstream::out);
            }
        }


}


void AdaptiveForcing::apply() {
	if (m_solver.getIteration() % 10 == 0) {
        getRhoU();

        write();
	}
}

    void AdaptiveForcing::getRhoU() {
        //////////////////////////
        // Calculate averages ////
        //////////////////////////
        double value=0.0;
        double average=0.0;
        int number_values=0;

        boost::shared_ptr<AdvectionOperator<3> > advection =
                m_solver.getAdvectionOperator();
        const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
                                                 | dealii::update_gradients;
        const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
        dealii::FEValues<3> fe_values(advection->getMapping(),
                                      *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
        size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
        std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);

        // loop
        typename dealii::DoFHandler<3>::active_cell_iterator cell =
                dof_handler.begin_active(), endc = dof_handler.end();
        double y;
        size_t dof_ind;
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {

                cell->get_dof_indices(local_indices);

                // get averages
                fe_values.reinit(cell);
                const std::vector<dealii::Point<3> >& quad_points =
                        fe_values.get_quadrature_points();


                for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                    number_values += 1;
                    // fill value vector
                    value = m_rho(dof_ind) * m_u.at(0)(dof_ind);				// rho


                    // add to averages:
                        average += value;
                } /* for all quadrature points */
            } /* if locally owned */
        } /* for all cells */

        // communicate
        double total_value = dealii::Utilities::MPI::sum(number_values, MPI_COMM_WORLD);;
        int total_number = dealii::Utilities::MPI::sum(average, MPI_COMM_WORLD);


        total_value /= total_number;
        m_currentValue = total_value;
    }

    void AdaptiveForcing::write() {
        if (is_MPI_rank_0()) {


            *m_tableFile << this->m_solver.getIteration() << " ";
            *m_tableFile << this->m_solver.getTime() << " ";

            *m_tableFile << m_currentValue << " " << m_force << "\n";


        }
    }


AdaptiveForcing::~AdaptiveForcing() {
	// TODO Auto-generated destructor stub
}

} /* namespace natrium */
