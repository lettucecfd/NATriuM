//
// Created by dwilde3m on 04.02.21.
//

#include "deal.II/grid/tria_iterator.h"
#include "deal.II/fe/fe_update_flags.h"
#include "deal.II/fe/fe_dgq.h"
#include "deal.II/fe/fe_values.h"
#include "deal.II/dofs/dof_handler.h"
#include "deal.II/grid/tria_accessor.h"
#include "deal.II/grid/tria_iterator.h"
#include "CompressibleTurbulenceStats.h"

namespace natrium {
template<size_t dim>
CompressibleTurbulenceStats<dim>::CompressibleTurbulenceStats(CompressibleCFDSolver<dim> &solver) :m_solver(solver),
        m_gamma(solver.getConfiguration()->getHeatCapacityRatioGamma()),
        m_filename(outfile(solver.getConfiguration()->getOutputDirectory())), m_legendFilename(
        legendfile(solver.getConfiguration()->getOutputDirectory())), m_outputOff(
        solver.getConfiguration()->isSwitchOutputOff()) {
    m_names.push_back("timestep");
    m_names.push_back("time");
    m_names.push_back("dilatation");
    m_names.push_back("solenoidal");
    m_names.push_back("maxMach");
    m_names.push_back("totalEnergy");

    // make table file
    if (not m_outputOff) {
        // create file (if necessary)
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
            printHeaderLine();
        }
    }


}


    template<size_t dim>
    void CompressibleTurbulenceStats<dim>::printHeaderLine() {
        assert(not m_outputOff);
        if (is_MPI_rank_0()) {
            (*m_tableFile) << "# A description of this data is found in "
                           << m_legendFilename << endl;
            std::ofstream legend_file(m_legendFilename, std::fstream::out);

            int k = 0;
            legend_file << k << "  <dilatation>" << endl;
            k++;
            legend_file << k << "  <solenoidal>" << endl;
            k++;
            legend_file << k << "  maxMach" << endl;
            legend_file.close();

        } /* is_MPI_rank 0 */
    }

    template<size_t dim>
    void CompressibleTurbulenceStats<dim>::writeToFile() {
        if ((is_MPI_rank_0()) and (not m_outputOff)) {

            *m_tableFile << this->m_solver.getIteration() << " ";
            *m_tableFile << this->m_solver.getTime() << " ";

            *m_tableFile << m_dilatation << " ";
            *m_tableFile << m_solenoidal<< " ";
            *m_tableFile << m_maxMach << " ";
            *m_tableFile << m_totalEnergy << " ";
            *m_tableFile << endl;

        } /* is mpi rank 0 */
    }

    template<size_t dim>
    void CompressibleTurbulenceStats<dim>::calculate() {

        const vector<distributed_vector> & u(this->m_solver.getVelocity());
        const distributed_vector & rho(this->m_solver.getDensity());
        const distributed_vector & T(this->m_solver.getTemperature());

        //////////////////////////
        // Calculate averages ////
        //////////////////////////
        const dealii::UpdateFlags update_flags = dealii::update_values | dealii::update_gradients
                                                 | dealii::update_JxW_values;
        const dealii::DoFHandler<dim> & dof_handler =
                *(this->m_solver.getAdvectionOperator()->getDoFHandler());
        dealii::FEValues<dim> fe_values(
                this->m_solver.getAdvectionOperator()->getMapping(),
                *(this->m_solver.getAdvectionOperator()->getFe()),
                *(this->m_solver.getAdvectionOperator()->getQuadrature()),
                update_flags);
        size_t dofs_per_cell =
                this->m_solver.getAdvectionOperator()->getFe()->dofs_per_cell;
        size_t n_q_points = this->m_solver.getAdvectionOperator()->getQuadrature()->size();
        std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);
        std::vector<double> uxs;
        std::vector<double> uys;
        std::vector<double> uzs;
        std::vector<double> rhos;
        std::vector<double> Ts;

        std::vector<dealii::Tensor<1, dim, double> > ux_gradients;
        std::vector<dealii::Tensor<1, dim, double> > uy_gradients;
        std::vector<dealii::Tensor<1, dim, double> > uz_gradients;
        std::vector<dealii::Tensor<1, dim, double> > rho_gradients;
        std::vector<dealii::Tensor<1, dim, double> > T_gradients;
        double dilatation = 0.0;
        double solenoidal = 0.0;
        double globalMach = 0.0;
        double totalEnergy = 0.0;

        uxs.resize(n_q_points);
        uys.resize(n_q_points);
        uzs.resize(n_q_points);
        rhos.resize(n_q_points);
        ux_gradients.resize(n_q_points);
        uy_gradients.resize(n_q_points);
        uz_gradients.resize(n_q_points);
        rho_gradients.resize(n_q_points);
        Ts.resize(n_q_points);
        T_gradients.resize(n_q_points);


        // loop
            typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                    dof_handler.begin_active(), endc = dof_handler.end();
            for (; cell != endc; ++cell) {
                if (cell->is_locally_owned()) {

                    cell->get_dof_indices(local_indices);

                    // get averages
                    fe_values.reinit(cell);
                    const std::vector<double>& weights = fe_values.get_JxW_values();

                    // calculate gradients (for w and strain rate)
                    fe_values.get_function_gradients(u.at(0), ux_gradients);
                    fe_values.get_function_gradients(u.at(1), uy_gradients);
                    fe_values.get_function_values(u.at(0), uxs);
                    fe_values.get_function_values(u.at(1), uys);
                    if (3 == dim){
                        fe_values.get_function_gradients(u.at(2), uz_gradients);
                        fe_values.get_function_values(u.at(2), uzs);
                    }
                    fe_values.get_function_gradients(rho, rho_gradients);
                    fe_values.get_function_values(rho, rhos);
                    fe_values.get_function_gradients(T, T_gradients);
                    fe_values.get_function_values(T, Ts);

                    for (size_t q = 0; q < n_q_points; q++) {
                        /*	for (size_t i = 0; i < dofs_per_cell; i++) {*/
                        //dof_ind = local_indices.at(i);
                        totalEnergy += Ts.at(q) * weights.at(q);
                        double sunderland_factor = 1.402*pow(Ts.at(q),1.5) / (Ts.at(q) + 0.40417);
                        // add to energies and enstrophies

                        double e1 = ux_gradients.at(q)[0] + uy_gradients.at(q)[1];

                        double e2 = pow(
                                uy_gradients.at(q)[0] - ux_gradients.at(q)[1], 2);
                        if (3 == dim) {
                            e1 += uz_gradients.at(q)[2];
                            e2 += pow(uz_gradients.at(q)[1] - uy_gradients.at(q)[2],
                                      2);
                            e2 += pow(ux_gradients.at(q)[2] - uz_gradients.at(q)[0],
                                      2);
                        }
                        e1 *= e1;
                        dilatation += sunderland_factor * e1 * weights.at(q);
                        solenoidal += sunderland_factor* e2 * weights.at(q);

                        double u_magnitude = uxs.at(q)*uxs.at(q)+uys.at(q)*uys.at(q);
                            if (3 == dim)
                                u_magnitude+=uzs.at(q)*uzs.at(q);
                        u_magnitude=sqrt(u_magnitude);
                        double localMach = u_magnitude/sqrt(m_gamma*Ts.at(q));
                        if (localMach > globalMach)
                            globalMach=localMach;

                        //} /* for all dof indices */
                    } /* for all quadrature nodes */
                } /* if locally owned */
            } /* for all cells */

            // communicate
            //for (size_t i = 0; i < m_nofObservables; i++) {


            m_dilatation = dealii::Utilities::MPI::sum(dilatation, MPI_COMM_WORLD);
            m_solenoidal = dealii::Utilities::MPI::sum(solenoidal, MPI_COMM_WORLD);
            m_maxMach = dealii::Utilities::MPI::max(globalMach,MPI_COMM_WORLD);
            m_totalEnergy = dealii::Utilities::MPI::sum(totalEnergy, MPI_COMM_WORLD);
            //}
        }

        template<size_t dim>
        void CompressibleTurbulenceStats<dim>::apply() {
            {
                if ((this->m_solver.getIteration()
                     % this->m_solver.getConfiguration()->getOutputTableInterval()
                     == 0) and (not m_outputOff)) {
                    calculate();
                    writeToFile();
                }
            }
        }

//    template class CompressibleTurbulenceStats<2> ;
    template class CompressibleTurbulenceStats<3> ;

}