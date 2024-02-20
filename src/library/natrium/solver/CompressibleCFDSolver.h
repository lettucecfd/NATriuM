/*
 * CompressibleCFDSolver.h
 *
 *  Created on: 25.10.2018
 *      Author: natrium
 */

#ifndef COMPRESSIBLECFDSOLVER_H_
#define COMPRESSIBLECFDSOLVER_H_

#include <filesystem>
#include "deal.II/grid/grid_out.h"
#include "CFDSolver.h"
#include "../utilities/BasicNames.h"
#include "../utilities/CFDSolverUtilities.h"
#include "../collision_advanced/CollisionSelection.h"
#include "../collision_advanced/AuxiliaryCollisionFunctions.h"
#include "../collision_advanced/Equilibria.h"
#include "../smoothing/VmultLimiter.h"
#include "../dataprocessors/CompressibleTurbulenceStats.h"
#include "Checkpoint.h"

namespace natrium {

template<size_t dim>
class CompressibleCFDSolver: public CFDSolver<dim> {
    template<size_t dim6> friend class CompressibleTurbulenceStats;
private:
	/// macroscopic density
	distributed_vector m_temperature;
    distributed_vector m_tmpTemperature;
    boost::shared_ptr<CompressibleTurbulenceStats<dim> > m_compressibleTurbulenceStats;
    // starting time
    time_t m_tstart;
    time_t m_tstart3;
    string m_tstart2;

    /// particle distribution functions for internal energy
	DistributionFunctions m_g;

public:
	CompressibleCFDSolver(boost::shared_ptr<SolverConfiguration> configuration,
			boost::shared_ptr<ProblemDescription<dim> > problemDescription) :
            CFDSolver<dim>(configuration, problemDescription) {

        m_tstart = clock();
        m_tstart3 = time(nullptr);
        struct tm* ltm = localtime(&m_tstart3);
        m_tstart2 = string(asctime(ltm));

        bool checkpointExists = applyCheckpointToG();

        m_temperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
                             this->getAdvectionOperator()->getLocallyRelevantDofs(),
                             MPI_COMM_WORLD);
        m_tmpTemperature.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
                                MPI_COMM_WORLD);
        if (configuration->isOutputCompressibleTurbulenceStatistics()) {
            m_compressibleTurbulenceStats = boost::make_shared<CompressibleTurbulenceStats<dim> >(*this);
        }
        this->m_maskShockSensor.reinit(this->getAdvectionOperator()->getLocallyOwnedDofs(),
                                       this->getAdvectionOperator()->getLocallyRelevantDofs(),
                                       MPI_COMM_WORLD);

        if(!checkpointExists) {
            distributed_vector writeable_temperature;
            // In this case, the density function fulfills the same purpose for the temperature
            CFDSolverUtilities::getWriteableDensity(writeable_temperature, m_temperature,
                                                    this->getAdvectionOperator()->getLocallyOwnedDofs());
            this->applyInitialTemperatures(writeable_temperature, this->getSupportPoints());
            CFDSolverUtilities::applyWriteableDensity(writeable_temperature, m_temperature);
            LOG(BASIC) << "Speed of Sound Square: " << this->m_stencil->getSpeedOfSoundSquare() << endl;
            LOG(BASIC) << "Mach compressible:     " << this->m_configuration->getMachNumber() << endl << endl;

            m_g.reinit(this->m_stencil->getQ(),
                       this->getAdvectionOperator()->getLocallyOwnedDofs(),
                       this->getAdvectionOperator()->getLocallyRelevantDofs(),
                       MPI_COMM_WORLD, (SEDG == configuration->getAdvectionScheme()));
        }
        if (checkpointExists)
        {
            this->calculateDensitiesAndVelocities();
            distributed_vector writeable_T;
            CFDSolverUtilities::getWriteableDensity(writeable_T, m_temperature,
                                                    this->m_advectionOperator->getLocallyOwnedDofs());
            //for all degrees of freedom on current processor
            dealii::IndexSet::ElementIterator it(
                    this->m_advectionOperator->getLocallyOwnedDofs().begin());
            dealii::IndexSet::ElementIterator end(
                    this->m_advectionOperator->getLocallyOwnedDofs().end());
            for (; it != end; it++) {
                double temperature = 0.0;
                for (size_t i = 0; i < this->m_stencil->getQ(); i++) {
                    double sum = 0.0;
                    for (size_t a = 0; a < dim; a++) {
                        sum += (this->getStencil()->getDirection(i)[a] - this->m_velocity.at(a)(i)) * (this->getStencil()->getDirection(i)[a] - this->m_velocity.at(a)(i));
                    }
                    temperature += sum * this->m_f.at(i)(*it) / this->getStencil()->getSpeedOfSoundSquare() +
                            m_g.at(i)(*it);
                }
                const double C_v = 1./(this->getConfiguration()->getHeatCapacityRatioGamma()-1.0);
                temperature = temperature * 0.5 / (this->m_density(*it)*C_v);
                writeable_T(*it)=temperature;
            }
            m_g.updateGhosted();

            CFDSolverUtilities::applyWriteableDensity(writeable_T, m_temperature);
        }
	}

    bool applyCheckpointToG() {
        boost::filesystem::path out_dir(this->m_configuration->getOutputDirectory());
        boost::filesystem::path log_file = out_dir / "natrium.log";
        boost::filesystem::path checkpoint_dir = out_dir / "checkpoint";

        // Restart from checkpoint?
        boost::shared_ptr<Checkpoint<dim> > checkpoint;
//        CheckpointStatus checkpoint_status;
        size_t restart_i = this->m_configuration->getRestartAtIteration();
        if (0 != restart_i) {
            checkpoint = boost::make_shared<Checkpoint<dim> >(restart_i,
                    checkpoint_dir, true);
            if (not checkpoint->exists()) {
                std::stringstream msg;
                msg << "You want to restart from iteration " << restart_i
                    << ", but I could not find the required checkpoint files "
                    << checkpoint->getStatusFile().string() << " and "
                    << checkpoint->getDataFile().string() << ".";
                natrium_errorexit(msg.str().c_str());
            } else {
                LOG(BASIC) << "Restart at iteration (also for g)" << restart_i << endl;
            }
        }

        if (checkpoint) {
            m_g.reinit(this->m_stencil->getQ(),
                       this->getAdvectionOperator()->getLocallyOwnedDofs(),
                       this->getAdvectionOperator()->getLocallyRelevantDofs(),
                       MPI_COMM_WORLD, (SEDG == this->m_configuration->getAdvectionScheme()));
            dealii::IndexSet::ElementIterator it(
                    this->m_advectionOperator->getLocallyOwnedDofs().begin());
            dealii::IndexSet::ElementIterator end(
                    this->m_advectionOperator->getLocallyOwnedDofs().end());

            const double scaling = this->m_stencil->getScaling();
                for (; it != end; it++) {
                    double density = 0.0;
                    std::array<double,dim> vel;
                    for (size_t i = 0; i < this->m_stencil->getQ(); i++) {
                        density +=  this->m_f.at(i)(*it);
                        for (size_t a = 0; a < dim; a++) {
                            vel[a] += this->m_stencil->getDirection(i)(a)/scaling *this->m_f.at(i)(*it);
                        }
                    }
                        double temperature = 0.0;
                    for (size_t i = 0; i < this->m_stencil->getQ(); i++) {
                        double sum = 0.0;
                        for (size_t a = 0; a < dim; a++) {
                            sum += (this->m_stencil->getDirection(i)(a)/scaling - vel[a]) * (this->m_stencil->getDirection(i)(a)/scaling - vel[a]);
                        }
                        temperature += sum * this->m_f.at(i)(*it) / (this->m_stencil->getSpeedOfSoundSquare()/scaling/scaling);
                    }
                    // Only one distribution --> gamma = d
                    const double gamma_old =  (static_cast<double>(dim)+2) / static_cast<double>(dim);
                    const double C_v_old = 1. / (gamma_old - 1.0);
                    const double gamma_new = this->getConfiguration()->getHeatCapacityRatioGamma();
                    const double C_v_new = 1. / (gamma_new - 1.0);
                    temperature = temperature * 0.5 / (density*C_v_old);
                for (size_t i = 0; i < this->m_stencil->getQ(); i++) {
                   m_g.at(i)(*it) = this->m_f.at(i)(*it)*temperature*(2.0*C_v_new-dim);
                }
            }
           // checkpoint->load(m_g, *(this->m_problemDescription), *(this->m_advectionOperator),
            //                 checkpoint_status);
        }
        m_g.updateGhosted();
        return 0 != restart_i;
	}

    void stream() {
        // start timer
        TimerOutput::Scope timer_section(Timing::getTimer(), "Stream");

        // no streaming in direction 0; begin with 1
        distributed_block_vector& f = this->m_f.getFStream();
        const distributed_sparse_block_matrix& systemMatrix =
                this->m_advectionOperator->getSystemMatrix();

        if (SEMI_LAGRANGIAN == this->m_configuration->getAdvectionScheme()) {

            DistributionFunctions f_tmp(this->m_f);
            systemMatrix.vmult(this->m_f.getFStream(), f_tmp.getFStream());
            this->m_advectionOperator->applyBoundaryConditions(f_tmp, this->m_f, m_g, this->m_time);
            /*distributed_block_vector f_tmp(f.n_blocks());
             reinitVector(f_tmp, f);
             f_tmp = f;
             systemMatrix.vmult(f, f_tmp);*/

            //m_advectionOperator->applyBoundaryConditions( f_tmp, f,  m_time);
            if (this->m_configuration->isVmultLimiter()) {
                TimerOutput::Scope timer_section(Timing::getTimer(), "Limiter");
                VmultLimiter::apply(systemMatrix, this->m_f.getFStream(),
                                    f_tmp.getFStream());
            }

            if ((BGK_MULTI_AM4 == this->m_configuration->getCollisionScheme()
                 || (BGK_MULTI_BDF2 == this->m_configuration->getCollisionScheme()))
                && (this->m_i - this->m_iterationStart) > 1) {
                //distributed_block_vector& formerF =
                //		m_multistepData->getFormerF().getFStream();
                f_tmp = this->m_multistepData->getFormerF();
                assert(this->m_multistepData != NULL);
                systemMatrix.vmult(this->m_multistepData->getFormerF().getFStream(),
                                   f_tmp.getFStream());
                this->m_advectionOperator->applyBoundaryConditions(f_tmp,
                                                                   this->m_multistepData->getFormerF(), this->m_time);
                //m_advectionOperator->applyBoundaryConditions( f_tmp, formerF,  m_time);
                if (this->m_configuration->isVmultLimiter()) {
                    TimerOutput::Scope timer_section(Timing::getTimer(), "Limiter");
                    VmultLimiter::apply(systemMatrix,
                                        this->m_multistepData->getFormerF().getFStream(),
                                        f_tmp.getFStream());
                }
                //distributed_block_vector& formerFEq =
                //		m_multistepData->getFormerFEq().getFStream();
                f_tmp = this->m_multistepData->getFormerFEq();
                systemMatrix.vmult(this->m_multistepData->getFormerFEq().getFStream(),
                                   f_tmp.getFStream());
                this->m_advectionOperator->applyBoundaryConditions(f_tmp,
                                                                   this->m_multistepData->getFormerFEq(), this->m_time);
                //m_advectionOperator->applyBoundaryConditions( f_tmp, formerFEq,  m_time);
                if (this->m_configuration->isVmultLimiter()) {
                    TimerOutput::Scope timer_section(Timing::getTimer(), "Limiter");
                    VmultLimiter::apply(systemMatrix,
                                        this->m_multistepData->getFormerFEq().getFStream(),
                                        f_tmp.getFStream());
                }
            }
            this->m_time += this->getTimeStepSize();

        } else {
            this->m_boundaryVector = this->m_advectionOperator->getSystemVector();
            double new_dt = this->m_timeIntegrator->step(f, systemMatrix,
                                                         this->m_boundaryVector, 0.0, this->m_timeIntegrator->getTimeStepSize());

            // For multistep methods, the former PDF for f and feq have also to be streamed
            if ((BGK_MULTI_AM4 == this->m_configuration->getCollisionScheme()
                 || (BGK_MULTI_BDF2 == this->m_configuration->getCollisionScheme()))
                && (this->m_i - this->m_iterationStart) > 1) {

                distributed_block_vector& formerF =
                        this->m_multistepData->getFormerF().getFStream();
                distributed_block_vector& formerFEq =
                        this->m_multistepData->getFormerFEq().getFStream();

                this->m_timeIntegrator->step(formerF, systemMatrix, this->m_boundaryVector, 0.0,
                                             this->m_timeIntegrator->getTimeStepSize());

                this->m_timeIntegrator->step(formerFEq, systemMatrix, this->m_boundaryVector,
                                       0.0, this->m_timeIntegrator->getTimeStepSize());
            }
            this->m_timeIntegrator->setTimeStepSize(new_dt);
            this->m_time += new_dt;
            this->m_collisionModel->setTimeStep(this->m_timeIntegrator->getTimeStepSize());
        }

        // communicate
        this->m_f.updateGhosted();
    }

void gStream() {

	// no streaming in direction 0; begin with 1
//	distributed_block_vector& g = m_g.getFStream();
	const distributed_sparse_block_matrix& systemMatrix =
		this->getAdvectionOperator()->getSystemMatrix();

	if (SEMI_LAGRANGIAN == this->getConfiguration()->getAdvectionScheme()) {
		DistributionFunctions g_tmp(m_g);
		systemMatrix.vmult(m_g.getFStream(), g_tmp.getFStream());
		const double gamma = this->getConfiguration()->getHeatCapacityRatioGamma();
        this->getAdvectionOperator()->applyBoundaryConditionsToG(this->m_f,this->m_g,this->m_time,gamma);
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
	else
    {
        this->m_boundaryVector = this->m_advectionOperator->getSystemVector();
//        double new_dt = this->m_timeIntegrator->step(g, systemMatrix,
//                                               this->m_boundaryVector, 0.0, this->m_timeIntegrator->getTimeStepSize());
    }
}

std::string secs_to_stream(int secs) {
    int h = int(secs/3600);
    int m = int((secs - h*3600)/60);
    int s = secs - h*3600 - m*60;
    std::stringstream result;
    result << std::setfill('0') << std::setw(3) << h << ":"
           << std::setfill('0') << std::setw(2) << m << ":"
           << std::setfill('0') << std::setw(2) << s; // << " / " << secs << " seconds";
    return result.str();
};

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
    void calculateTemperature() {
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

/*	inline distributed_vector& getWriteableTemperature(distributed_vector& writeable, const distributed_vector& member, const dealii::IndexSet& locally_owned){
		TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
		writeable.reinit(locally_owned, member.get_mpi_communicator(), true);
		assert (not writeable.has_ghost_elements());
		writeable = member;
		assert (not writeable.has_ghost_elements());
		return writeable;
	}*/

	/**
	 * @short copy the changes in the writeable copy to the global density, see getWritableVelocity for a detailed explanation
	 */
	inline void applyWriteableTemperature(const distributed_vector& writeable, distributed_vector& member){
		TimerOutput::Scope timer_section(Timing::getTimer(), "Copy vectors");
		member = writeable;
	}

	void compressibleOutput(size_t iteration, bool is_final) {
        if (is_final and is_MPI_rank_0()) LOG(DETAILED) << "Printing last outputs." << endl;

// start timer
        TimerOutput::Scope timer_section(Timing::getTimer(), "Output");

// sync MPI processes
        MPI_sync();

// output time estimates

        if ((iteration == 1) or (iteration == 10) or (iteration == 50) or (iteration == 100) or (iteration == 500)
            or ((iteration % 1000 == 0) and iteration > 0)) {
            int secs2 = int(clock() / CLOCKS_PER_SEC);
            time_t t_now = time(nullptr);
            struct tm *ltm = localtime(&t_now);
            LOG(DETAILED) << "Iteration " << iteration << ", t = " << this->m_time << ", server-time = "
                          << secs_to_stream(secs2) << "." << endl;
//                LOG(DETAILED) << ":Started simulation after " << m_tstart << " seconds." << endl;
            LOG(DETAILED) << ":Started at " << m_tstart2;
            LOG(DETAILED) << ":::::::Now, it is " << string(asctime(ltm));

            double secs = 1e-10 + (clock() - this->m_tstart) / CLOCKS_PER_SEC;
            LOG(DETAILED) << ":::::::Average Performance: "
                          << 1.0 * this->m_advectionOperator->getDoFHandler()->n_dofs()
                             * (iteration - this->m_iterationStart) / secs / 1000000.0
                          << " million DoF updates per second" << endl;
            Timing::getTimer().print_summary();

            if ((this->m_configuration->getNumberOfTimeSteps() < 1e8) or (this->m_configuration->getSimulationEndTime() < 1e8)) {
                double factor;
                string base;
                if (this->m_configuration->getNumberOfTimeSteps() < 1e8) {
                    base = "max iterations";
                    // Calculating done iterations
                    int done_iterations = iteration - this->m_iterationStart;
                    // Calculating to-be-done iterations
                    double tobedone_iterations = this->m_configuration->getNumberOfTimeSteps() - this->m_iterationStart;
                    // Calculating iterations factor
                    factor = tobedone_iterations / done_iterations;
                }
                else if (this->m_configuration->getSimulationEndTime() < 1e8) {
                    base = "physical time";
                    // Calculating done time
                    double done_time_ph = this->getTime() - this->m_tstart_ph;
                    // Calculating to-be-done iterations
                    double tobedone_time_ph = this->m_configuration->getSimulationEndTime() - this->m_tstart_ph;
                    // Calculating time factor
                    factor = tobedone_time_ph / done_time_ph;
                }
                else {
                    factor = 1;
                }
                time_t start = m_tstart3;
                time_t done_time = clock()/CLOCKS_PER_SEC;
                time_t tobedone_time = done_time * factor;
                time_t estimated_end = start + tobedone_time;
                struct tm * ltm2 = localtime(&estimated_end);
                //                struct tm * ltm1 = localtime(&start);
                LOG(DETAILED) << ":Finished " << int(100.0/factor * 10000) / 10000 << " % based on " << base << ". Estimated end: " << string(asctime(ltm2));
                time_t server_max = this->m_configuration->getServerEndTime();
                time_t estimated_server_end = start + server_max;
                struct tm * ltm3 = localtime(&estimated_server_end);
                LOG(DETAILED) << ":::::::Server-time left: " << secs_to_stream(server_max - done_time) << ".   Estimated server-end: " << string(asctime(ltm3));
//                    LOG(DETAILED) << "Server time: " << clock()/CLOCKS_PER_SEC << endl;
//                    LOG(DETAILED) << "Calculated done_time: " << done_time << endl;
//                    LOG(DETAILED) << "Calculated tobedone_t " << tobedone_time << endl;
                //                LOG(DETAILED) << "Started at: " << string(asctime(ltm1));
                //                LOG(DETAILED) << "Already did: " << done_time << "[time_t]";
                //                LOG(DETAILED) << "Overall needs: " << tobedone_time << "[time_t]";
            }
        }

// output: vector fields as .vtu files
        if ((iteration == 0) and (not this->m_configuration->isSwitchOutputOff())) {
            std::filesystem::path out_dir(this->m_configuration->getOutputDirectory() + "/vtk");
            std::filesystem::create_directory(out_dir);
        }
        int no_out = this->m_configuration->getNoOutputInterval();
        if ((not this->m_configuration->isSwitchOutputOff())
            and (no_out < int(iteration))) {
            // output elapsed time, server-time, and estimated runtime after iterations 1, 10, 100, 1000, ... and after every 1000
            // add turbulence statistics to output
            if (this->m_configuration->isOutputTurbulenceStatistics())
                this->m_turbulenceStats->addToReynoldsStatistics(this->m_velocity);

            // create vtk-subdirectory
            // no output if solution interval > 10^8
//            double maxP = this->m_solverStats->getMaxP();
            if (((iteration % this->m_configuration->getOutputSolutionInterval() == 0)
                 and (this->m_configuration->getOutputSolutionInterval() <= 1e8))
                or (is_final)
//                or (maxP > 0.6)) {
//                if (maxP > 0.6) {
//                    LOG(DETAILED) << "Doing VTK output at " << iteration << " because MaxP = " << maxP
//                                  << "; MinP = " << this->m_solverStats->getMinP() << endl;
//                }
                    ){
                // save local part of the solution
                std::stringstream str;
                str << this->m_configuration->getOutputDirectory().c_str() << "/vtk/t_"
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
                                  << "/vtk/t_"
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
            // output: table
            // calculate information + physical properties
            if (iteration % this->m_configuration->getOutputTableInterval() == 0) {
                this->m_solverStats->printNewLine();
                if (this->m_configuration->isOutputTurbulenceStatistics()) {
                    assert(this->m_turbulenceStats);
                    this->m_turbulenceStats->printNewLine();
                }
            }
        }
        // output: checkpoint
        // no output if checkpoint interval > 10^8
        if (((iteration % this->m_configuration->getOutputCheckpointInterval() == 0) or is_final)
                and (this->m_configuration->getOutputCheckpointInterval() <= 1e8)
                and (no_out < int(iteration))
                and (this->m_iterationStart != this->m_i)) {

            if (int(iteration) < no_out) {
                if (is_MPI_rank_0()) LOG(DETAILED) << "No output at " << iteration << ". Only after " << no_out << endl;
            } else {
                boost::filesystem::path checkpoint_dir(
                        this->m_configuration->getOutputDirectory());
                checkpoint_dir /= "checkpoint";
                Checkpoint<dim> checkpoint(this->m_i, checkpoint_dir, false);
                Checkpoint<dim> checkpointG(this->m_i, checkpoint_dir, true);
                CheckpointStatus checkpoint_status;
                checkpoint_status.iterationNumber = this->m_i;
                checkpoint_status.stencilScaling = this->m_stencil->getScaling();
                checkpoint_status.time = this->m_time;
                checkpoint_status.feOrder =
                        this->m_configuration->getSedgOrderOfFiniteElement();
                checkpoint.write(*(this->m_problemDescription->getMesh()), this->m_f,
                                 *(this->m_advectionOperator->getDoFHandler()), checkpoint_status);
                checkpointG.write(*(this->m_problemDescription->getMesh()), m_g,
                                  *(this->m_advectionOperator->getDoFHandler()), checkpoint_status);
            }
        } /*if checkpoint interval*/
        if (is_final) {
            LOG(DETAILED) << "Total runtime: " << secs_to_stream(int(clock()/CLOCKS_PER_SEC)) << ". Server-end had been set to " << secs_to_stream(this->m_configuration->getServerEndTime()) << endl;
        }
    }

    void applyShockSensor() {

        const vector<distributed_vector> &u(this->getVelocity());
        const distributed_vector &rho(this->getDensity());
        const distributed_vector &T(this->getTemperature());

        const unsigned int dofs_per_cell =
                this->getAdvectionOperator()->getFe()->dofs_per_cell;
        vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);
        typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                this->getAdvectionOperator()->getDoFHandler()->begin_active(), endc =
                this->getAdvectionOperator()->getDoFHandler()->end();
        const dealii::UpdateFlags update_flags = dealii::update_values | dealii::update_gradients
                                                 | dealii::update_JxW_values;
        const dealii::DoFHandler<dim> &dof_handler =
                *(this->getAdvectionOperator()->getDoFHandler());
        dealii::FEValues<dim> fe_values(
                this->getAdvectionOperator()->getMapping(),
                *(this->getAdvectionOperator()->getFe()),
                *(this->getAdvectionOperator()->getQuadrature()),
                update_flags);
        size_t n_q_points = this->getAdvectionOperator()->getQuadrature()->size();
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

        distributed_vector writeable_Mask;

        CFDSolverUtilities::getWriteableDensity(writeable_Mask, this->m_maskShockSensor,
                                                this->getAdvectionOperator()->getLocallyOwnedDofs());

        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_dof_indices);

                cell->get_dof_indices(local_dof_indices);

                // get averages
                fe_values.reinit(cell);
                const std::vector<double> &weights = fe_values.get_JxW_values();

                // calculate gradients (for w and strain rate)
                fe_values.get_function_gradients(u.at(0), ux_gradients);
                fe_values.get_function_gradients(u.at(1), uy_gradients);
                fe_values.get_function_values(u.at(0), uxs);
                fe_values.get_function_values(u.at(1), uys);
                if (3 == dim) {
                    fe_values.get_function_gradients(u.at(2), uz_gradients);
                    fe_values.get_function_values(u.at(2), uzs);
                }
                fe_values.get_function_gradients(rho, rho_gradients);
                fe_values.get_function_values(rho, rhos);
                fe_values.get_function_gradients(T, T_gradients);
                fe_values.get_function_values(T, Ts);

                for (size_t q = 0; q < n_q_points; q++) {

                    double dilatation = ux_gradients.at(q)[0] + uy_gradients.at(q)[1];

                    if (3 == dim) {
                        dilatation += uz_gradients.at(q)[2];
                    }
                    //double u_magnitude = uxs.at(q)*uxs.at(q)+uys.at(q)*uys.at(q);
                    //if (3 == dim)
                    //    u_magnitude+=uzs.at(q)*uzs.at(q);
                    //    u_magnitude=sqrt(u_magnitude);

                    writeable_Mask(local_dof_indices.at(q)) = dilatation;

                }
            } /* if is locally owned */
        }/* for all cells */
        CFDSolverUtilities::applyWriteableDensity(writeable_Mask, this->m_maskShockSensor);
    }

        void smoothDensities(
            const map<dealii::types::global_dof_index, dealii::Point<dim> >& supportPoints)  {
        array<double,25> feq;
        array<double,2> u;
// get Function instance
        const dealii::IndexSet& locally_owned_dofs =
                this->m_advectionOperator->getLocallyOwnedDofs();
        dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
        dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());
        for (; it != end; it++) {
            size_t i = *it;
 //       const boost::shared_ptr<dealii::Function<dim> >& f_rho =
 //               this->m_problemDescription->getInitialRhoFunction();
/*        const unsigned int dofs_per_cell =
                this->m_advectionOperator->getFe()->dofs_per_cell;
        vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

        typename dealii::DoFHandler<dim>::active_cell_iterator cell =
                this->m_advectionOperator->getDoFHandler()->begin_active(), endc =
                this->m_advectionOperator->getDoFHandler()->end();
        for (; cell != endc; ++cell) {
            if (cell->is_locally_owned()) {
                cell->get_dof_indices(local_dof_indices);
                for (size_t i = 0; i < dofs_per_cell; i++) {
                    if (not this->m_advectionOperator->getLocallyOwnedDofs().is_element(
                            local_dof_indices.at(i))) {
                        continue;
                    }
                    assert(
                            supportPoints.find(local_dof_indices.at(i))
                            != supportPoints.end());*/
                    if(supportPoints.at(i)(0)<27.90) {
                            u[0] = -0.611022;
                            u[1] = -0.0;

                        GeneralCollisionData<2,25> data(*(this->m_configuration), *(this->m_problemDescription), this->m_stencil->getScaling(), this->m_problemDescription->getViscosity(), *(this->m_stencil), this->m_stencil->getSpeedOfSoundSquare(), 0.0);
                        //this->m_density(i) = 1.34206;
                        //m_temperature(i) = 1.129;
                        data.density = 1.34161;
                        data.temperature = 1.12799;
                        data.velocity = u;
                        BGKEquilibrium<2,25> eq;
                        eq.calc(feq, data);
                        double gamma = this->getConfiguration()->getHeatCapacityRatioGamma();
                        double C_v = 1. / (gamma - 1.0);
                        //m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
                        for (size_t j = 0; j < this->m_stencil->getQ(); j++) {
                            this->m_f.at(j)(i) = feq[j];
                            m_g.at(j)(i) = feq[j]*(data.temperature)*(2.0*C_v-2.0);
                        }
                    }
                }
            } /* if is locally owned */
 //       } /* for all cells */
 //   }

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
            CFDSolverUtilities::applyWriteableDensity(writeable_mSS, this->m_maskShockSensor);

            applyWriteableTemperature(writeable_T, m_temperature);
			this->m_f.updateGhosted();
            this->m_g.updateGhosted();

		} catch (CollisionException& e) {
			natrium_errorexit(e.what());
		}
	}

    inline void calcQuarticEquilibrium(vector<double>& feq, size_t T_Q, double density, array<double,dim> velocity,
                                       double temperature, double cs2, vector<array<double,dim>> e,
                                       vector<double> weight) {
        const array<array<size_t,dim>, dim> eye = unity_matrix<dim>();
        double uu_term = 0.0;
        for (size_t j = 0; j < dim; j++) {
            uu_term += -(velocity[j] * velocity[j])
                       / (2.0 * cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {

            double T1 = cs2*(temperature-1);

            double ue_term = 0.0;
            for (size_t j = 0; j < dim; j++) {
                ue_term += (velocity[j] * e[i][j]) / cs2;
            }
            feq[i] = weight[i] * density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
            for (size_t alp = 0; alp < dim; alp++){
                for (size_t bet = 0; bet < dim; bet++){
                    feq[i] += density * weight[i] / (2.0*cs2) * (temperature - 1) *
                            (eye[alp][bet] * e[i][alp] * e[i][bet] - cs2 * eye[alp][bet]);
                    for (size_t gam = 0; gam < dim; gam++){
                        feq[i] += weight[i] * density / (6. * cs2 * cs2 * cs2) *
                                   (velocity[alp] * velocity[bet] * velocity[gam]
                                    + T1 * (eye[alp][bet] * velocity[gam]
                                          + eye[bet][gam] * velocity[alp]
                                          + eye[alp][gam] * velocity[bet])) *
                                          // Hermite-3
                                    (e[i][alp] * e[i][bet] * e[i][gam]
                                        - cs2 * (e[i][gam] * eye[alp][bet]
                                               + e[i][bet] * eye[alp][gam]
                                               + e[i][alp] * eye[bet][gam]));

                        for (size_t det = 0; det < dim; det++)
                        {
                            double power4 = e[i][alp]*e[i][bet]*e[i][gam]*e[i][det];
                            double power2 = e[i][alp]*e[i][bet]*eye[gam][det]
                                            +e[i][alp]*e[i][gam]*eye[bet][det]
                                            +e[i][alp]*e[i][det]*eye[bet][gam]
                                            +e[i][bet]*e[i][gam]*eye[alp][det]
                                            +e[i][bet]*e[i][det]*eye[alp][gam]
                                            +e[i][gam]*e[i][det]*eye[alp][bet];
                            double power0 = eye[alp][bet]*eye[gam][det]
                                           +eye[alp][gam]*eye[bet][det]
                                           +eye[alp][det]*eye[bet][gam];
                            double u4    = velocity[alp]*velocity[bet]*velocity[gam]*velocity[det];
                            double u2 = velocity[alp]*velocity[bet]*eye[gam][det]
                                       +velocity[alp]*velocity[gam]*eye[bet][det]
                                       +velocity[alp]*velocity[det]*eye[bet][gam]
                                       +velocity[bet]*velocity[gam]*eye[alp][det]
                                       +velocity[bet]*velocity[det]*eye[alp][gam]
                                       +velocity[gam]*velocity[det]*eye[alp][bet];
                            double multieye= eye[alp][bet]*eye[gam][det]
                                            +eye[alp][gam]*eye[bet][det]
                                            +eye[alp][det]*eye[bet][gam];

                            feq[i]+= weight[i] * density / (24.*cs2*cs2*cs2*cs2) *
                                    (power4 - cs2 * power2 + cs2 * cs2 * power0) * (u4 + T1 * (u2 + T1 * multieye));
                        }
                    }
                }
            }
        }
    }

	void initializeDistributions()  {
	// PRECONDITION: vectors already created with the right sizes

		LOG(BASIC) << "Initialize distribution functions: ";
        size_t T_Q = this->getStencil()->getQ();
        std::vector<double> feq(T_Q);
		std::array<double,dim> u;
        std::vector<std::array<double,dim>> e(T_Q);
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < this->getStencil()->getQ(); ++j) {
                e[j][i] = this->getStencil()->getDirections().at(j)(i) / this->getStencil()->getScaling();
            }
        }
        std::vector<double> weight(T_Q);
        for (size_t j = 0; j < this->getStencil()->getQ(); ++j) {
            weight[j] = this->getStencil()->getWeight(j);
        }
        double descaled_cs2 = this->m_stencil->getSpeedOfSoundSquare() / (this->getStencil()->getScaling()*this->getStencil()->getScaling());

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
                u[j] = this->m_velocity.at(j)(i)/this->m_stencil->getScaling();
            }
            //GeneralCollisionData<2,25> data(*(this->m_configuration), *(this->m_problemDescription), this->m_stencil->getScaling(), this->m_problemDescription->getViscosity(), *(this->m_stencil), this->m_stencil->getSpeedOfSoundSquare(), 0.0);

            double density = this->m_density[i];
            double temperature = this->m_temperature[i];
            //BGKEquilibrium<2,25> eq;
            //      eq.calc(feq, data);
            calcQuarticEquilibrium(feq,T_Q,density,u,temperature,descaled_cs2,e,weight);


            double gamma = this->getConfiguration()->getHeatCapacityRatioGamma();
                  double C_v = 1. / (gamma - 1.0);
                //m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
                for (size_t j = 0; j < this->m_stencil->getQ(); j++) {
                    this->m_f.at(j)(i) = feq[j];
                    this->m_g.at(j)(i) = feq[j]*(temperature)*(2.0*C_v-dim);
                }
            }

		LOG(BASIC) << "Equilibrium distribution functions" << endl;
			// do nothing else

        switch (this->m_configuration->getInitializationScheme()) {
            case EQUILIBRIUM: {
                LOG(BASIC) << "Equilibrium distribution functions" << endl;
                // do nothing else
                break;
            }
            case COMPRESSIBLE_ITERATIVE: {
                LOG(BASIC) << "Compressible iterative procedure" << endl;
                LOG(DETAILED) << "residual = "
                              << this->m_configuration->getIterativeInitializationResidual();
                LOG(DETAILED) << ", max iterations = "
                              << this->m_configuration->getIterativeInitializationNumberOfIterations()
                              << endl;
                // Iterative procedure; leading to consistent initial values
                size_t loopCount = 0;
                double residual = 1000000000;
//                const bool inInitializationProcedure = true;
                distributed_vector oldDensities;
                while (loopCount
                       < this->m_configuration->getIterativeInitializationNumberOfIterations()) {

                    distributed_vector rho;
                    vector<distributed_vector> u;
                    distributed_vector T;
                    distributed_vector mask;
                    // get writeable copies of rho and u
                    CFDSolverUtilities::getWriteableDensity(oldDensities, this->m_density,
                                                            locally_owned_dofs);
                    CFDSolverUtilities::getWriteableDensity(rho, this->m_density,
                                                            locally_owned_dofs);
                    CFDSolverUtilities::getWriteableDensity(mask, this->m_maskShockSensor,
                                                            locally_owned_dofs);
                    CFDSolverUtilities::getWriteableVelocity(u, this->m_velocity,
                                                             locally_owned_dofs);
                    CFDSolverUtilities::getWriteableDensity(T, m_temperature,
                                                            locally_owned_dofs);
                    try {
                        this->stream();
                    } catch (std::exception& e) {
                        natrium_errorexit(e.what());
                    }
                    // collide without recalculating velocities
                    try {
                        // collide
                        double delta_t = CFDSolverUtilities::calculateTimestep<dim>(
                                *(this->getProblemDescription()->getMesh()),
                                this->getConfiguration()->getSedgOrderOfFiniteElement(), *(this->m_stencil),
                                this->getConfiguration()->getCFL());

                        selectCollision(*(this->m_configuration), *(this->m_problemDescription), this->m_f, m_g, rho, u, T, mask,
                                        this->m_advectionOperator->getLocallyOwnedDofs(), this->m_problemDescription->getViscosity(), delta_t, *(this->getStencil()), true);
                        //m_collisionModel->collideAll(m_f, rho, m_velocity, locally_owned_dofs,
                        //                             inInitializationProcedure);
                        // copy back
                    } catch (CollisionException& e) {
                        natrium_errorexit(e.what());
                    }
                    //oldDensities -= rho;
                    residual = oldDensities.norm_sqr();
                    //CFDSolverUtilities::applyWriteableDensity(rho, this->m_density);
                    this->m_f.updateGhosted();
                    this->m_g.updateGhosted();
                    //CFDSolverUtilities::applyWriteableVelocity(u, m_velocity);
                    loopCount++;
                }
                LOG(DETAILED) << "Residual " << residual << " reached after "
                              << loopCount << " iterations." << endl;

                //for all degrees of freedom on current processor
                /*for (it = locally_owned_dofs.begin(); it != end; it++) {
                 size_t i = *it;
                 for (size_t j = 0; j < dim; j++) {
                 u(j) = m_velocity.at(j)(i);
                 }
                 m_collisionModel->getEquilibriumDistributions(feq, u, m_density(i));
                 for (size_t j = 0; j < m_stencil->getQ(); j++) {
                 m_f.at(j)(i) = feq.at(j);
                 }
                 }*/
                break;
            }
            default: {
                throw CFDSolverException(
                        "Error in CFDSolver::InitializeDistributions. A part of the code was reached, which should never be reached. With ITERATIVE you enforce it for compressible flows, however");
                break;
            }
        }

        this->m_time = t0;

		LOG(BASIC) << "Initialize distribution functions: done." << endl;
	}

    void run()
    {
        this->m_i = this->m_iterationStart;
        if(this->m_i<=1) {
            initializeDistributions();
        }
        collide();
            while (true) {
            if (this->stopConditionMet()) {
                break;
            }
            //applyShockSensor();
            // Deactivated this->output(this->m_i);
            this->compressibleOutput(this->m_i, false);
            this->m_i++;
            this->stream();
            gStream();
            compressibleFilter();
            if (this->m_configuration->isOutputCompressibleTurbulenceStatistics()) {
                this->m_compressibleTurbulenceStats->apply();
            }
            if(this->m_i==200) {
                //smoothDensities(this->m_supportPoints);
            }

            this->collide();

            for (size_t i = 0; i < this->m_dataProcessors.size(); i++) {
                this->m_dataProcessors.at(i)->apply();
            }
        }
        //this->output(this->m_i, true);
        this->compressibleOutput(this->m_i, true);

        // Finalize
        if (is_MPI_rank_0()) {
            Timing::getTimer().print_summary();
        }
        LOG(BASIC) << "NATriuM run complete." << endl;
        LOG(BASIC) << "Summary: " << endl;
        LOG(BASIC) << Timing::getOutStream().str() << endl;
    }

    const DistributionFunctions& getG() const {
        return m_g;
    }

    DistributionFunctions& getG() {
        return m_g;
    }

//protected:
//    /// the physical time passed (normally initialized with 0.0, except for restart at a checkpoint)
//    double m_time;

};



} //namespace natrium
#endif /* COMPRESSIBLECFDSOLVER_H_ */
