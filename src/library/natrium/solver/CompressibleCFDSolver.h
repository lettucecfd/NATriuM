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

namespace natrium {

template<size_t dim>
class CompressibleCFDSolver: public CFDSolver<dim> {
private:
	/// macroscopic density
	distributed_vector m_temperature;
	distributed_vector m_tmpTemperature;

public:
	CompressibleCFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
			boost::shared_ptr<ProblemDescription<dim> > problemDescription) :
			CFDSolver<dim>(configuration, problemDescription) {

		m_temperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
				this->getAdvectionOperator()->getLocallyRelevantDofs(),
				MPI_COMM_WORLD);
		m_tmpTemperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
			MPI_COMM_WORLD);
		distributed_vector writeable_temperature;
		// In this case, the density function fulfills the same purpose for the temperature
		CFDSolverUtilities::getWriteableDensity(writeable_temperature, m_temperature,
						this->getAdvectionOperator()->getLocallyOwnedDofs());
		this->applyInitialDensities(writeable_temperature, this->getSupportPoints());
		CFDSolverUtilities::applyWriteableDensity(writeable_temperature, m_temperature);

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

	void run()
	{
		this->m_i = this->m_iterationStart;
		collide();
		while (true) {
			if (this->stopConditionMet()) {
				break;
			}
			this->output(this->m_i);
			this->m_i++;
			this->stream();
			this->filter();
			this->collide();
			for (size_t i = 0; i < this->m_dataProcessors.size(); i++) {
				this->m_dataProcessors.at(i)->apply();
			}
		}
		this->output(this->m_i, true);

	// Finalize
		if (is_MPI_rank_0()) {
			Timing::getTimer().print_summary();
		}
		LOG(BASIC) << "NATriuM run complete." << endl;
		LOG(BASIC) << "Summary: " << endl;
		LOG(BASIC) << Timing::getOutStream().str() << endl;
	}


	void collide()  {
cout << "COLLISIONTEST" << endl;
	// start timer
		TimerOutput::Scope timer_section(Timing::getTimer(), "Collision");

		try {
			// get writeable copies of density and velocity and temperature
			std::vector<distributed_vector> writeable_u;
			distributed_vector writeable_rho;
			distributed_vector writeable_T;
			CFDSolverUtilities::getWriteableVelocity(writeable_u, this->m_velocity,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			CFDSolverUtilities::getWriteableDensity(writeable_rho, this->m_density,
					this->getAdvectionOperator()->getLocallyOwnedDofs());
			CFDSolverUtilities::getWriteableDensity(writeable_T, this->m_temperature,
					this->getAdvectionOperator()->getLocallyOwnedDofs());

			double delta_t = CFDSolverUtilities::calculateTimestep<dim>(
						*(this->getProblemDescription()->getMesh()),
						this->getConfiguration()->getSedgOrderOfFiniteElement(), *(this->m_stencil),
						this->getConfiguration()->getCFL());





	// TODO member function collisionModel
			selectCollision(*(this->m_configuration), *(this->m_problemDescription), this->m_f, writeable_rho, writeable_u, writeable_T,
				this->m_advectionOperator->getLocallyOwnedDofs(), this->m_problemDescription->getViscosity(), delta_t, *(this->getStencil()), false);

			 //perform collision
	//		m_collisionModel->collideAll(m_f, writeable_rho, writeable_u,
	//				m_advectionOperator->getLocallyOwnedDofs(), false);

			// copy back to ghosted vectors and communicate across MPI processors
			CFDSolverUtilities::applyWriteableDensity(writeable_rho, this->m_density);
			CFDSolverUtilities::applyWriteableVelocity(writeable_u, this->m_velocity);
			applyWriteableTemperature(writeable_T, m_temperature);
			this->m_f.updateGhosted();

		} catch (CollisionException& e) {
			natrium_errorexit(e.what());
		}
	}

};



} //namespace natrium
#endif /* COMPRESSIBLECFDSOLVER_H_ */
