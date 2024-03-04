//
// Created by dominik on 31.10.22.
//

#include "BoundaryTemperature.h"
namespace natrium {

    BoundaryTemperature::BoundaryTemperature(CompressibleCFDSolver<3> & solver) : DataProcessor<3>(solver), m_compressibleSolver(solver)
    {

    }

    BoundaryTemperature::~BoundaryTemperature() {
        // TODO Auto-generated destructor stub
    }

    void BoundaryTemperature::apply() {
        setBoundaryTemperature();
    }

    void BoundaryTemperature::setBoundaryTemperature() {

        std::array<double, 45> w;
        int length;
        for (size_t i = 0; i < 45; i++) {
            w[i] = m_solver.getStencil()->getWeight(i);
        }

        const dealii::UpdateFlags update_flags = dealii::update_quadrature_points
                                                 | dealii::update_gradients;
        const dealii::DoFHandler<3> &dof_handler = *(m_solver.getAdvectionOperator()->getDoFHandler());
        typename dealii::DoFHandler<3>::active_cell_iterator cell =
                dof_handler.begin_active(), endc = dof_handler.end();
        dealii::FEValues<3> fe_values(m_solver.getAdvectionOperator()->getMapping(),
                                      *(m_solver.getAdvectionOperator()->getFe()),
                                      m_solver.getAdvectionOperator()->getSupportPointEvaluation(), update_flags);
        size_t dofs_per_cell = m_solver.getAdvectionOperator()->getFe()->dofs_per_cell;
        std::vector<dealii::types::global_dof_index> local_indices(dofs_per_cell);
        std::set<int> already_set;


        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_indices);
                fe_values.reinit(cell);

                const std::vector<dealii::Point<3> > &quad_points =
                        fe_values.get_quadrature_points();

                for (size_t q = 0; q < fe_values.n_quadrature_points; q++) {
                    size_t dof = local_indices.at(q);
                    if (not m_solver.getAdvectionOperator()->getLocallyOwnedDofs().is_element(
                            local_indices.at(q))) {
                        continue;
                    }
                    if (already_set.count(dof)>0)
                        continue;

                    const double tol_distance = 0.0001;
                    double T_wall = 0.0;

                    if (quad_points.at(q)(2) < tol_distance or
                        std::abs(quad_points.at(q)(2) - 1.0) < tol_distance) {
                        if (quad_points.at(q)(2) < tol_distance)
                            T_wall = 0.8;
                        if (std::abs(quad_points.at(q)(2) - 1.0) < tol_distance)
                            T_wall = 1.2;


                        const double scaling = m_solver.getStencil()->getScaling();
                        const double cs2 = m_solver.getStencil()->getSpeedOfSoundSquare() / (scaling * scaling);
                        const double gamma = 1.4;
                        assert(m_solver.getStencil()->getQ() == 45);
                        std::array<double, 45> f_destination={0.0};
                        std::array<double, 45> g_destination={0.0};;
                        std::array<double, 45> feq={0.0};
                        std::array<double, 45> geq={0.0};

                        for (int i = 0; i < 45; i++) {
                            f_destination[i] = m_solver.getF().at(i)(dof);
                            g_destination[i] = m_compressibleSolver.getG().at(i)(dof);
                        }

                        const double rho = calculateDensity<45>(f_destination);
                        std::array<double, 3> u_local;
                        std::array<std::array<double, 3>, 45> e = getParticleVelocitiesWithoutScaling<3, 45>(
                                *m_solver.getStencil());
                        calculateVelocity<3, 45>(f_destination, u_local, rho, e);

                        const double T_local = calculateTemperature<3, 45>(f_destination, g_destination, u_local, rho,
                                                                           e,
                                                                           cs2,
                                                                           gamma);
                            QuarticEquilibrium<3, 45> eq(cs2, e);
                            eq.polynomial(feq, rho, u_local, T_local, e, w, cs2);
                            calculateGeqFromFeq<3, 45>(feq, geq, T_local, gamma);



                            for (int i = 0; i < 45; i++) {
                                f_destination[i] -= feq[i];
                                g_destination[i] -= geq[i];
                                feq[i] = 0.0;
                                geq[i] = 0.0;
                            }
                            eq.polynomial(feq, rho, u_local, T_wall, e, w, cs2);
                            calculateGeqFromFeq<3, 45>(feq, geq, T_wall, gamma);

                            for (int i = 0; i < 45; i++) {
                                m_solver.getF().at(i)(dof) =
                                        f_destination[i] + feq[i];

                                m_compressibleSolver.getG().at(i)(dof) =
                                        geq[i];

                                already_set.insert(dof);
                            }
                        }
                    }
                }
            }

        m_solver.getF().updateGhosted();
        m_compressibleSolver.getG().updateGhosted();

    } /* if is locally owned */
} // namespace natrium