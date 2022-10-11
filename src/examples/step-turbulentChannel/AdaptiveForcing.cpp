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
		FinalChannelStatistics(solver, outdir), m_outDir(outdir), m_u(solver.getVelocity()), m_rho(
				solver.getDensity()), m_targetRhoU(target), m_lastRhoU(target),
                m_filename(outfile(solver.getConfiguration()->getOutputDirectory())) {

        m_yCoordsUpToDate = false;

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
    if (!isMYCoordsUpToDate())
        updateYValues();

	if (m_solver.getIteration() % 1 == 0) {
        getRhoU();
        calculateForce();
        write();
	}
}

    void AdaptiveForcing::getRhoU() {
        //////////////////////////
        // Calculate averages ////
        //////////////////////////
        double value=0.0;
        vector<double> average;
        average.resize(m_nofCoordinates);
        vector<size_t> number;
        number.resize(m_nofCoordinates);

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
        size_t y_ind;

        size_t dof_ind;
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {



                cell->get_dof_indices(local_indices);

                // get averages
                fe_values.reinit(cell);
                const std::vector<dealii::Point<3> >& quad_points =
                        fe_values.get_quadrature_points();


                for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                    y = quad_points.at(i)(1);
                    assert(
                            m_yCoordinateToIndex.find(y)
                            != m_yCoordinateToIndex.end());
                    y_ind = m_yCoordinateToIndex.at(y);

                    dof_ind = local_indices.at(i);
                    number.at(y_ind) += 1;
                    // fill value vector
                    // add to averages:
                    average.at(y_ind) += m_rho(dof_ind) * m_u.at(0)(dof_ind);				// rho


                } /* for all quadrature points */
            } /* if locally owned */
        } /* for all cells */

        // communicate
        for (size_t i = 0; i < m_nofCoordinates; i++) {
            number.at(i) = dealii::Utilities::MPI::sum(number.at(i), MPI_COMM_WORLD);
            average.at(i) = dealii::Utilities::MPI::sum(average.at(i), MPI_COMM_WORLD);
            average.at(i)/=number.at(i);
        }


        double integral = 0;
        double total_window = 0;
        for (size_t i = 0; i < m_nofCoordinates-1; i++) {
            cout << m_yCoordinates.at(i) << " " << average.at(i) << std::endl;
            double window_size = std::abs( m_yCoordinates.at(i+1) -m_yCoordinates.at(i));
            total_window +=window_size;
            integral += window_size*0.5*(average.at(i)+average.at(i+1));
        }
        cout << total_window << "HIER" << endl;
        m_currentValue = integral/2.0;

    }

    void AdaptiveForcing::write() {
        if (is_MPI_rank_0()) {


            *m_tableFile << this->m_solver.getIteration() << " ";
            *m_tableFile << this->m_solver.getTime() << " ";
            *m_tableFile << m_targetRhoU << " " << m_currentValue << " " << m_force << endl;


        }
    }


AdaptiveForcing::~AdaptiveForcing() {
	// TODO Auto-generated destructor stub
}

    void AdaptiveForcing::calculateForce() {
        const double currentForce = this->m_solver.getProblemDescription()->getExternalForce()->getForce()[0];
        const double timeStepSize = this->m_solver.getTimeStepSize();
        double newForce = currentForce + 1./timeStepSize*(m_targetRhoU-2*m_currentValue+m_lastRhoU);


        if (newForce/currentForce>1.01)
            newForce = 1.01*currentForce;
        if (newForce/currentForce<0.99 or newForce<0)
            newForce = 0.99*currentForce;
        dealii::Tensor<1,3> forceTensor;
        forceTensor[0]=newForce;
        forceTensor[1]=0.0;
        forceTensor[2]=0.0;
        //this->m_solver.getProblemDescription()->setExternalForceTensor(forceTensor);
        m_force = newForce;
    }

} /* namespace natrium */
