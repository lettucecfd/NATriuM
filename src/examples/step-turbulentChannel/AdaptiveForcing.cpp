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

    AdaptiveForcing::AdaptiveForcing(CompressibleCFDSolver<3> &solver,
                                     std::string outdir, double target, bool restart) :
            FinalChannelStatistics(solver, outdir), m_outDir(outdir), m_T(solver.getTemperature()), m_targetRhoU(target),
            m_lastRhoU(target), m_starting_force(this->m_solver.getProblemDescription()->getExternalForce()->getForce()[0]),
            m_restart(restart),
            m_filename(outfile(solver.getConfiguration()->getOutputDirectory())) {


        m_names.push_back("T");
        m_names.push_back("dT/dx");
        m_names.push_back("dT/dy");
        m_names.push_back("dT/dz");

        m_yCoordsUpToDate = false;

        m_nofObservables = m_names.size();
        m_averages.resize(m_nofObservables);
        m_EX3.resize(m_nofObservables);
        m_EX4.resize(m_nofObservables);
        m_correlations.resize(m_nofObservables);
        for (size_t i = 0; i < m_nofObservables; i++) {
            m_correlations.at(i).resize(i + 1);
        }

        m_averages_time.resize(m_nofObservables);
        m_EX3_time.resize(m_nofObservables);
        m_EX4_time.resize(m_nofObservables);
        m_correlations_time.resize(m_nofObservables);
        for (size_t i = 0; i < m_nofObservables; i++) {
            m_correlations_time.at(i).resize(i + 1);
        }
        n_steps = 0;
        m_nofCoordinates = 0;

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

    if(m_restart) {
        if (m_solver.getIteration() % 10 == 0) {
            // add to statistics
            if (!isMYCoordsUpToDate())
                updateYValues();
            updateCompressibleAverages();
            addToTemporalAverages();
        }
        if (m_solver.getIteration() % 5000 == 0) {
            write_to_file();
        }
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
            double window_size = std::abs( m_yCoordinates.at(i+1) -m_yCoordinates.at(i));
            total_window +=window_size;
            integral += window_size*0.5*(average.at(i)+average.at(i+1));
        }
        m_currentValue = integral/2.0;

    }

    void AdaptiveForcing::write() {
        if (is_MPI_rank_0()) {


            *m_tableFile << this->m_solver.getIteration() << " ";
            *m_tableFile << this->m_solver.getTime() << " ";
            *m_tableFile << m_targetRhoU << " " << m_currentValue << " " << m_force << endl;


        }
    }

    void AdaptiveForcing::updateCompressibleAverages() {
        boost::shared_ptr<AdvectionOperator<3> > advection =
                m_solver.getAdvectionOperator();

        // prepare local vectors
        vector<size_t> l_number;
        vector<double> l_values;
        vector<vector<double> > l_averages;
        vector<vector<vector<double> > > l_correlations;
        vector<vector<double> > l_EX3; // for skewness
        vector<vector<double> > l_EX4; // for kurtosis
        l_averages.resize(m_nofObservables);
        l_EX3.resize(m_nofObservables);
        l_EX4.resize(m_nofObservables);
        l_correlations.resize(m_nofObservables);
        for (size_t i = 0; i < m_nofObservables; i++) {
            l_correlations.at(i).resize(i + 1);
        }
        l_values.resize(m_nofObservables);
        for (size_t i = 0; i < m_nofObservables; i++) {
            l_averages.at(i).resize(m_nofCoordinates);
            l_EX3.at(i).resize(m_nofCoordinates);
            l_EX4.at(i).resize(m_nofCoordinates);
            for (size_t j = 0; j <= i; j++) {
                l_correlations.at(i).at(j).resize(m_nofCoordinates);
            }
        }
        l_number.resize(m_nofCoordinates);

        //////////////////////////
        // Calculate averages ////
        //////////////////////////
        const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
                                                 | dealii::update_gradients;
        const dealii::DoFHandler<3> & dof_handler = *(advection->getDoFHandler());
        dealii::FEValues<3> fe_values(advection->getMapping(),
                                      *(advection->getFe()), advection->getSupportPointEvaluation(), update_flags);
        size_t dofs_per_cell = advection->getFe()->dofs_per_cell;
        std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);
        std::vector<dealii::Tensor<1, 3, double> > ux_gradients;
        std::vector<dealii::Tensor<1, 3, double> > uy_gradients;
        std::vector<dealii::Tensor<1, 3, double> > uz_gradients;
        std::vector<dealii::Tensor<1, 3, double> > rho_gradients;
        std::vector<dealii::Tensor<1, 3, double> > T_gradients;

        ux_gradients.resize(dofs_per_cell);
        uy_gradients.resize(dofs_per_cell);
        uz_gradients.resize(dofs_per_cell);
        rho_gradients.resize(dofs_per_cell);
        T_gradients.resize(dofs_per_cell);
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

                // calculate gradients (for w and strain rate)
                fe_values.get_function_gradients(m_u.at(0), ux_gradients);
                fe_values.get_function_gradients(m_u.at(1), uy_gradients);
                fe_values.get_function_gradients(m_u.at(2), uz_gradients);
                fe_values.get_function_gradients(m_rho, rho_gradients);
                fe_values.get_function_gradients(m_T, T_gradients);

                for (size_t i = 0; i < fe_values.n_quadrature_points; i++) {
                    y = quad_points.at(i)(1);
                    assert(
                            m_yCoordinateToIndex.find(y)
                            != m_yCoordinateToIndex.end());
                    y_ind = m_yCoordinateToIndex.at(y);
                    dof_ind = local_indices.at(i);
                    l_number.at(y_ind) += 1;

                    // fill value vector
                    l_values.at(0) = m_rho(dof_ind);				// rho
                    l_values.at(1) = rho_gradients.at(i)[0];		// drho/dx
                    l_values.at(2) = rho_gradients.at(i)[1]; 		// drho/dy
                    l_values.at(3) = rho_gradients.at(i)[2];		// drho/dz

                    l_values.at(4) = m_u.at(0)(dof_ind); 			// ux
                    l_values.at(5) = m_u.at(1)(dof_ind);			// uy
                    l_values.at(6) = m_u.at(2)(dof_ind);			// uz
                    l_values.at(7) = ux_gradients.at(i)[0];			// dux/dx
                    l_values.at(8) = ux_gradients.at(i)[1];			// dux/dy
                    l_values.at(9) = ux_gradients.at(i)[2];			// dux/dz
                    l_values.at(10) = uy_gradients.at(i)[0];		// duy/dx
                    l_values.at(11) = uy_gradients.at(i)[1];		// duy/dy
                    l_values.at(12) = uy_gradients.at(i)[2];		// duy/dz
                    l_values.at(13) = uz_gradients.at(i)[0];		// duz/dx
                    l_values.at(14) = uz_gradients.at(i)[1];		// duz/dy
                    l_values.at(15) = uz_gradients.at(i)[2];		// duz/dz
                    l_values.at(16) = m_T(dof_ind);				// T
                    l_values.at(17) = T_gradients.at(i)[0];		// dT/dx
                    l_values.at(18) = T_gradients.at(i)[1]; 	// dT/dy
                    l_values.at(19) = T_gradients.at(i)[2];		// dT/dz

                    // add to averages:
                    for (size_t j = 0; j < m_nofObservables; j++) {
                        l_averages.at(j).at(y_ind) += l_values.at(j);
                    }
                    // add to correlations
                    for (size_t j = 0; j < m_nofObservables; j++) {
                        for (size_t k = 0; k < j + 1; k++) {
                            l_correlations.at(j).at(k).at(y_ind) += (l_values.at(j)
                                                                     * l_values.at(k));
                        }
                    }
                    // add to third moment
                    for (size_t j = 0; j < m_nofObservables; j++) {
                        l_EX3.at(j).at(y_ind) += pow(l_values.at(j), 3);
                        l_EX4.at(j).at(y_ind) += pow(l_values.at(j), 4);
                    }

                } /* for all quadrature points */
            } /* if locally owned */
        } /* for all cells */

        // communicate
        dealii::Utilities::MPI::sum(l_number, MPI_COMM_WORLD, m_number);
        for (size_t i = 0; i < m_nofObservables; i++) {
            dealii::Utilities::MPI::sum(l_averages.at(i), MPI_COMM_WORLD, m_averages.at(i));
            dealii::Utilities::MPI::sum(l_EX3.at(i), MPI_COMM_WORLD, m_EX3.at(i));
            dealii::Utilities::MPI::sum(l_EX4.at(i), MPI_COMM_WORLD, m_EX4.at(i));
            for (size_t j = 0; j < i + 1; j++) {
                dealii::Utilities::MPI::sum(l_correlations.at(i).at(j),
                                            MPI_COMM_WORLD, m_correlations.at(i).at(j));
            }
        }

        // divide by number of replicates
        for (size_t i = 0; i < m_nofCoordinates; i++) {
            for (size_t j = 0; j < m_nofObservables; j++) {
                m_averages.at(j).at(i) /= m_number.at(i);
                m_EX3.at(j).at(i) /= m_number.at(i);
                m_EX4.at(j).at(i) /= m_number.at(i);
                for (size_t k = 0; k < j + 1; k++) {
                    m_correlations.at(j).at(k).at(i) /= m_number.at(i);
                }
            }
        }
    }


AdaptiveForcing::~AdaptiveForcing() {
	// TODO Auto-generated destructor stub
}

    void AdaptiveForcing::calculateForce() {
        const double currentForce = this->m_solver.getProblemDescription()->getExternalForce()->getForce()[0];
        const double timeStepSize = this->m_solver.getTimeStepSize();
        double newForce = m_force; //currentForce + 1./timeStepSize*(m_targetRhoU-2*m_currentValue+m_lastRhoU);
        bool forceChanged = false;

        if (m_currentValue/m_targetRhoU>1.005)
            {newForce = 0.4*m_starting_force;
            forceChanged = true;}
        if (m_currentValue/m_targetRhoU<0.995)
            {newForce = 0.85*m_starting_force;
            forceChanged = true;}
        if (m_currentValue/m_targetRhoU<0.98)
	    {newForce = 2.0 * m_starting_force;
	    forceChanged = true;}

        dealii::Tensor<1,3> forceTensor;
        forceTensor[0]=newForce;
        forceTensor[1]=0.0;
        forceTensor[2]=0.0;
	if (forceChanged){
        this->m_solver.getProblemDescription()->setExternalForceTensor(forceTensor);
        m_force = newForce;}
    }

    void AdaptiveForcing::changeTemperatureProfile() {

    }

} /* namespace natrium */
