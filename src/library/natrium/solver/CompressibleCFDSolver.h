/*
 * CompressibleCFDSolver.h
 *
 *  Created on: 25.10.2018
 *      Author: natrium
 */

#ifndef COMPRESSIBLECFDSOLVER_H_
#define COMPRESSIBLECFDSOLVER_H_

#include "CFDSolver.h"
#include "../utilities/BasicNames.h"
#include "../utilities/CFDSolverUtilities.h"
#include "../collision_advanced/CollisionSelection.h"
#include "../collision_advanced/AuxiliaryCollisionFunctions.h"
#include "../collision_advanced/Equilibria.h"
#include "../smoothing/VmultLimiter.h"

namespace natrium {

template<size_t dim>
class CompressibleCFDSolver: public CFDSolver<dim> {
private:
	/// macroscopic density
	distributed_vector m_temperature;
    distributed_vector m_tmpTemperature;


    /// particle distribution functions for internal energy
	DistributionFunctions m_g;

public:
	CompressibleCFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
			boost::shared_ptr<ProblemDescription<dim> > problemDescription) :
			CFDSolver<dim>(configuration, problemDescription) {

		m_temperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
				this->getAdvectionOperator()->getLocallyRelevantDofs(),
				MPI_COMM_WORLD);
		m_tmpTemperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
			MPI_COMM_WORLD);

		this->m_maskShockSensor.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
                this->getAdvectionOperator()->getLocallyRelevantDofs(),
                MPI_COMM_WORLD);

		distributed_vector writeable_temperature;
		// In this case, the density function fulfills the same purpose for the temperature
		CFDSolverUtilities::getWriteableDensity(writeable_temperature, m_temperature,
						this->getAdvectionOperator()->getLocallyOwnedDofs());
		this->applyInitialTemperatures(writeable_temperature, this->getSupportPoints());
		CFDSolverUtilities::applyWriteableDensity(writeable_temperature, m_temperature);
		LOG(BASIC) << "Speed of Sound Square: " << this->m_stencil->getSpeedOfSoundSquare() << endl;

		m_g.reinit(this->m_stencil->getQ(),
				this->getAdvectionOperator()->getLocallyOwnedDofs(),
				this->getAdvectionOperator()->getLocallyRelevantDofs(),
				MPI_COMM_WORLD, (SEDG == configuration->getAdvectionScheme()));

	}

	
void gStream() {



	// no streaming in direction 0; begin with 1
	distributed_block_vector& g = m_g.getFStream();
	const distributed_sparse_block_matrix& systemMatrix =
		this->getAdvectionOperator()->getSystemMatrix();

	if (SEMI_LAGRANGIAN == this->getConfiguration()->getAdvectionScheme()) {

		DistributionFunctions g_tmp(m_g);
		systemMatrix.vmult(m_g.getFStream(), g_tmp.getFStream());
		this->getAdvectionOperator()->applyBoundaryConditions(g_tmp, m_g, this->m_time);
		/*distributed_block_vector f_tmp(f.n_blocks());
		 reinitVector(f_tmp, f);
		 f_tmp = f;
		 systemMatrix.vmult(f, f_tmp);*/

		//m_advectionOperator->applyBoundaryConditions( f_tmp, f,  m_time);
		if (this->getConfiguration()->isVmultLimiter()) {
			TimerOutput::Scope timer_section(Timing::getTimer(), "Limiter");
			VmultLimiter::apply(systemMatrix, m_g.getFStream(),
					g_tmp.getFStream());
		}
	}
}

void compressibleFilter() {

// start timer
	TimerOutput::Scope timer_section(Timing::getTimer(), "Filter");

	if (this->m_configuration->isFiltering()) {
		if (this->m_i % this->m_configuration->getFilterInterval() == 0) {
			for (size_t i = 0; i < this->m_stencil->getQ(); i++) {
				this->m_filter->applyFilter(*(this->m_advectionOperator)->getDoFHandler(),
						this->m_f.at(i));
				this->m_filter->applyFilter(*(this->m_advectionOperator)->getDoFHandler(),
						m_g.at(i));
			}
		}
		this->m_f.updateGhosted();
		m_g.updateGhosted();
	}

}

	void applyInitialTemperatures(
			distributed_vector& initialTemperatures,
			const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints) const {

		// get Function instance
			const boost::shared_ptr<dealii::Function<dim> >& f_T =
					this->getProblemDescription()->getInitialTFunction();
			const unsigned int dofs_per_cell =
					this->getAdvectionOperator()->getFe()->dofs_per_cell;
			vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
			typename dealii::DoFHandler<dim>::active_cell_iterator cell =
					this->getAdvectionOperator()->getDoFHandler()->begin_active(), endc =
							this->getAdvectionOperator()->getDoFHandler()->end();
			for (; cell != endc; ++cell) {
				if (cell->is_locally_owned()) {
					cell->get_dof_indices(local_dof_indices);
					for (size_t i = 0; i < dofs_per_cell; i++) {
						if (not this->getAdvectionOperator()->getLocallyOwnedDofs().is_element(
								local_dof_indices.at(i))) {
							continue;
						}
						assert(
								supportPoints.find(local_dof_indices.at(i))
										!= supportPoints.end());
						initialTemperatures(local_dof_indices.at(i)) = f_T->value(
								supportPoints.at(local_dof_indices.at(i)));
					}
				} /* if is locally owned */
			} /* for all cells */


		}

	const distributed_vector& getTemperature() const {
		return m_temperature;
	}
		void calculateTemperature()
		{
			std::vector<distributed_vector> writeable_u;
			distributed_vector writeable_rho;
			distributed_vector writeable_T;
			distributed_vector temporary;
			CFDSolverUtilities::getWriteableVelocity(writeable_u, this->getVelocity(),
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			CFDSolverUtilities::getWriteableDensity(writeable_rho, this->getDensity(),
					this->getAdvectionOperator()->getLocallyOwnedDofs());

			CFDSolverUtilities::getWriteableDensity(writeable_T, m_temperature,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			size_t Q = this->getStencil()->getQ();
			writeable_T = 0;


			for (size_t i = 0; i < Q; i++) {
					//writeable_T.add(m_f.at(i));
					for (size_t j = 0; j < dim; j++) {
						writeable_T.add(this->getStencil()->getDirection(i)(j), this->getStencil()->getDirection(i)(j));
					}
				}

		}


	inline distributed_vector& getWriteableTemperature(distributed_vector& writeable, const distributed_vector& member, const dealii::IndexSet& locally_owned){
		TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
		writeable.reinit(locally_owned, member.get_mpi_communicator(), true);
		assert (not writeable.has_ghost_elements());
		writeable = member;
		assert (not writeable.has_ghost_elements());
		return writeable;
	}

	/**
	 * @short copy the changes in the writeable copy to the global density, see getWritableVelocity for a detailed explanation
	 */
	inline void applyWriteableTemperature(const distributed_vector& writeable, distributed_vector& member){
		TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
		member = writeable;
	}

	void secondOutput(size_t iteration, bool is_final)
		{
			// no output if solution interval > 10^8
			if (((iteration % this->m_configuration->getOutputSolutionInterval() == 0)
					and this->m_configuration->getOutputSolutionInterval() <= 1e8)
					or (is_final)) {
				// save local part of the solution
				std::stringstream str;
				str << this->m_configuration->getOutputDirectory().c_str() << "/t_"
						<< this->m_problemDescription->getMesh()->locally_owned_subdomain()
						<< "." << iteration << ".vtu";
				std::string filename = str.str();
				std::ofstream vtu_output(filename.c_str());
				dealii::DataOut<dim> data_out;
				data_out.attach_dof_handler(*(this->m_advectionOperator)->getDoFHandler());
				data_out.add_data_vector(this->m_density, "rho");
				data_out.add_data_vector(m_temperature, "T");
                data_out.add_data_vector(this->m_maskShockSensor, "shockSensor");
				if (dim == 2) {
					data_out.add_data_vector(this->m_velocity.at(0), "ux");
					data_out.add_data_vector(this->m_velocity.at(1), "uy");
				} else { //dim == 3
					data_out.add_data_vector(this->m_velocity.at(0), "ux");
					data_out.add_data_vector(this->m_velocity.at(1), "uy");
					data_out.add_data_vector(this->m_velocity.at(2), "uz");
				}

				/// For Benchmarks: add analytic solution
				this->addAnalyticSolutionToOutput(data_out);
				/// For turbulent flows: add turbulent statistics
				if (this->m_configuration->isOutputTurbulenceStatistics()) {
					this->m_turbulenceStats->addReynoldsStatisticsToOutput(data_out);
				}

				// tell the data processor the locally owned cells
				dealii::Vector<float> subdomain(
						this->m_problemDescription->getMesh()->n_active_cells());
				for (unsigned int i = 0; i < subdomain.size(); ++i)
					subdomain(i) =
							this->m_problemDescription->getMesh()->locally_owned_subdomain();
				data_out.add_data_vector(subdomain, "subdomain");

				// Write vtu file

                data_out.build_patches(
                        (this->m_configuration->getSedgOrderOfFiniteElement() == 1) ?
                                this->m_configuration->getSedgOrderOfFiniteElement() :
                                this->m_configuration->getSedgOrderOfFiniteElement() + 1);
				data_out.write_vtu(vtu_output);

				// Write pvtu file (which is a master file for all the single vtu files)
				if (is_MPI_rank_0()) {
					// generate .pvtu filename
					std::stringstream pvtu_filename;
					pvtu_filename << this->m_configuration->getOutputDirectory().c_str()
							<< "/t_"
							<< this->m_problemDescription->getMesh()->locally_owned_subdomain()
							<< "." << iteration << ".pvtu";
					std::ofstream pvtu_output(pvtu_filename.str().c_str());

					// generate all other filenames
					std::vector<std::string> filenames;
					for (unsigned int i = 0;
							i < dealii::Utilities::MPI::n_mpi_processes(
							MPI_COMM_WORLD); ++i) {
						std::stringstream vtu_filename_i;
						vtu_filename_i
						//<< m_configuration->getOutputDirectory().c_str() << "/"
						<< "t_" << i << "." << iteration << ".vtu";
						filenames.push_back(vtu_filename_i.str());
					}
					data_out.write_pvtu_record(pvtu_output, filenames);
				}
			}
		}

    void applyShockSensor()  {
        const unsigned int dofs_per_cell =
                this->getAdvectionOperator()->getFe()->dofs_per_cell;
        vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                this->getAdvectionOperator()->getDoFHandler()->begin_active(), endc =
                this->getAdvectionOperator()->getDoFHandler()->end();
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_dof_indices);
                double cell_average = 0.0;

                distributed_vector writeable_Mask;
                CFDSolverUtilities::getWriteableDensity(writeable_Mask, this->m_maskShockSensor,
                					this->getAdvectionOperator()->getLocallyOwnedDofs());

                for (size_t i = 0; i < dofs_per_cell; i++) {

                	writeable_Mask(local_dof_indices.at(i)) = 0.0;
                    cell_average += this->m_density(local_dof_indices.at(i))*m_temperature(local_dof_indices.at(i));

                }
                cell_average /= dofs_per_cell;
                double sum_mse = 0.0;

                for (size_t i = 0; i < dofs_per_cell; i++) {
                    if ((this->m_density(local_dof_indices.at(i))*m_temperature(local_dof_indices.at(i))-cell_average)>sum_mse)
                    sum_mse = ((this->m_density(local_dof_indices.at(i))*m_temperature(local_dof_indices.at(i)) - cell_average));
                }
                //sum_mse /= dofs_per_cell;



                CFDSolverUtilities::getWriteableDensity(writeable_Mask, this->m_maskShockSensor,
                					this->getAdvectionOperator()->getLocallyOwnedDofs());

                for (size_t i = 0; i < dofs_per_cell; i++)
                {
                    writeable_Mask(local_dof_indices.at(i)) = sum_mse;
                }

                CFDSolverUtilities::applyWriteableDensity(writeable_Mask, this->m_maskShockSensor);




            } /* if is locally owned */
        } /* for all cells */


    }


	void collide()  {
	// start timer
		TimerOutput::Scope timer_section(Timing::getTimer(), "Collision");

		try {
			// get writeable copies of density and velocity and temperature
			std::vector<distributed_vector> writeable_u;
			distributed_vector writeable_rho;
			distributed_vector writeable_T;
            distributed_vector writeable_mSS;
			CFDSolverUtilities::getWriteableVelocity(writeable_u, this->m_velocity,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			CFDSolverUtilities::getWriteableDensity(writeable_rho, this->m_density,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			CFDSolverUtilities::getWriteableDensity(writeable_T, this->m_temperature,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
            CFDSolverUtilities::getWriteableDensity(writeable_mSS, this->m_maskShockSensor,
                    this->getAdvectionOperator()->getLocallyOwnedDofs());


			double delta_t = CFDSolverUtilities::calculateTimestep<dim>(
						*(this->getProblemDescription()->getMesh()),
						this->getConfiguration()->getSedgOrderOfFiniteElement(), *(this->m_stencil),
						this->getConfiguration()->getCFL());





	// TODO member function collisionModel
            selectCollision(*(this->m_configuration), *(this->m_problemDescription), this->m_f, m_g, writeable_rho, writeable_u, writeable_T, writeable_mSS,
				this->m_advectionOperator->getLocallyOwnedDofs(), this->m_problemDescription->getViscosity(), delta_t, *(this->getStencil()), false);

			 //perform collision
			//m_collisionModel->collideAll(m_f, writeable_rho, writeable_u,
			//		m_advectionOperator->getLocallyOwnedDofs(), false);

			// copy back to ghosted vectors and communicate across MPI processors
			CFDSolverUtilities::applyWriteableDensity(writeable_rho, this->m_density);
			CFDSolverUtilities::applyWriteableVelocity(writeable_u, this->m_velocity);
			applyWriteableTemperature(writeable_T, m_temperature);
			this->m_f.updateGhosted();
            this->m_g.updateGhosted();

		} catch (CollisionException& e) {
			natrium_errorexit(e.what());
		}
	}

	void initializeDistributions()  {
	// PRECONDITION: vectors already created with the right sizes

		LOG(BASIC) << "Initialize distribution functions: ";
		array<double,25> feq;
		array<double,2> u;

	// save starting time
		double t0 = this->m_time;

	// Initialize f with the equilibrium distribution functions
	//for all degrees of freedom on current processor
		const dealii::IndexSet& locally_owned_dofs =
				this->m_advectionOperator->getLocallyOwnedDofs();
		dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
		for (; it != end; it++) {
			size_t i = *it;
			for (size_t j = 0; j < dim; j++) {
				u[j] = this->m_velocity.at(j)(i);
			}
        GeneralCollisionData<2,25> data(*(this->m_configuration), *(this->m_problemDescription), this->m_stencil->getScaling(), this->m_problemDescription->getViscosity(), *(this->m_stencil), this->m_stencil->getSpeedOfSoundSquare(), 0.0);
        data.density = this->m_density[i];
        data.temperature = this->m_temperature[i];
        data.velocity = u;
        BGKEquilibrium<2,25> eq;
              eq.calc(feq, data);
              double gamma = 1.4;
              double C_v = 1. / (gamma - 1.0);
			//m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
			for (size_t j = 0; j < this->m_stencil->getQ(); j++) {
				this->m_f.at(j)(i) = feq[j];
                // ONLY FOR 2D! (due to -2.0)
                this->m_g.at(j)(i) = feq[j]*(data.temperature)*(2.0*C_v-2.0);
			}
		}

		LOG(BASIC) << "Equilibrium distribution functions" << endl;
			// do nothing else




		this->m_time = t0;

		LOG(BASIC) << "Initialize distribution functions: done." << endl;

	}

    void run()
    {
        this->m_i = this->m_iterationStart;
        initializeDistributions();
        collide();
        while (true) {
            if (this->stopConditionMet()) {
                break;
            }
            //applyShockSensor();
            this->output(this->m_i);
            this->secondOutput(this->m_i,false);
            this->m_i++;
            this->stream();
            gStream();
            compressibleFilter();
            this->collide();
            for (size_t i = 0; i < this->m_dataProcessors.size(); i++) {
                this->m_dataProcessors.at(i)->apply();
            }
        }
        this->output(this->m_i, true);
        this->secondOutput(this->m_i,true);

        // Finalize
        if (is_MPI_rank_0()) {
            Timing::getTimer().print_summary();
        }
        LOG(BASIC) << "NATriuM run complete." << endl;
        LOG(BASIC) << "Summary: " << endl;
        LOG(BASIC) << Timing::getOutStream().str() << endl;
    }



};



} //namespace natrium
#endif /* COMPRESSIBLECFDSOLVER_H_ */
