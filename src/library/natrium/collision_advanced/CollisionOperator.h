/*
 * CollisionOperator.h
 *
 *  Created on: 06.01.2017
 *      Author: dominik
 */

#ifndef COLLISIONOPERATOR_H_
#define COLLISIONOPERATOR_H_

#include <vector>
#include "../stencils/Stencil.h"
#include "AuxiliaryCollisionFunctions.h"
#include "Equilibria.h"
#include "CollisionSchemes.h"

namespace natrium {


template<int T_D, int T_Q, template<int, int> class T_equilibrium,
		template<int, int,
				template<int, int> class > class T_collision>
class CollisionOperator {
public:

	void collideAll(DistributionFunctions& f, distributed_vector& densities,
			vector<distributed_vector>& velocities,
			const dealii::IndexSet& locally_owned_dofs,
			const bool inInitializationProcedure,
			GeneralCollisionData<T_D, T_Q>& genData,
			typename T_collision<T_D, T_Q, T_equilibrium>::SpecificCollisionData& specData) const {
		//for all degrees of freedom on current processor
		//dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		//dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

		T_collision<T_D, T_Q, T_equilibrium> collisionScheme;

		std::array< double*, T_Q> f_raw;
		int length;
		for (size_t i = 0; i < T_Q; i++){
			f.at(i).trilinos_vector().ExtractView(&f_raw[i], &length);
		}
		double* rho_raw;
		densities.trilinos_vector().ExtractView(&rho_raw, &length);
		std::array< double*, T_D> u_raw;
		for (size_t i = 0; i < T_D; i++){
			velocities.at(i).trilinos_vector().ExtractView(&u_raw[i], &length);
		}


#pragma omp simd
		for (int ii = 0; ii < length; ii++) {

			// Variable that stores the local distribution function values of every node

			// Copy the needed global distribution function values into the local variable
			for (int p = 0; p < T_Q; ++p) {
				genData.fLocal[p] = f_raw[p][ii];
			}
			//copyGlobalToLocalF<T_Q>(genData.fLocal, f, i); // done

			//Calculate the local density and store it into the Parameter Handling System
			genData.density = calculateDensity<T_Q>(genData.fLocal); // done

			//Set temperature so that Equilibrium terms cancel out
			genData.temperature = 1.0;

            genData.H3 = calculateH3<T_D,T_Q>(genData.cs2,genData.e);
            genData.H4 = calculateH4<T_D,T_Q>(genData.cs2,genData.e);

			//Write the local density to the global density vector
			rho_raw[ii] = genData.density; // write local density to global density vector

			//Calculate the local velocity and store it into the Parameter Handling System
			calculateVelocity<T_D, T_Q>(genData.fLocal, genData.velocity,
					genData.density, genData); // TODO velocities for other stencils

			//Write the local density to the global velocity matrix
			if (not inInitializationProcedure) {
				for (size_t j = 0; j < T_D; ++j) {
					u_raw[j][ii] = genData.velocity[j] * genData.scaling;
				}
				if (genData.problemDescription.hasExternalForce()) {
					applyMacroscopicForces<T_D, T_Q>(u_raw, ii, genData);
					applyForces<T_D, T_Q>(genData);
				}
			} else {
				for (size_t j = 0; j < T_D; ++j) {
					genData.velocity[j] = u_raw[j][ii] / genData.scaling;
				}
			}

			//Initialize an object of the desired collision scheme and run the relaxation process

			collisionScheme.relax(genData.fLocal, genData, specData);
            if (genData.problemDescription.hasExternalForce()) {
                postCollisionApplyForces<T_D, T_Q, T_equilibrium>(u_raw, ii, genData);
            }
			//Finally copy the updated distribution function back to the global distribution function
			for (int p = 0; p < T_Q; ++p) {
				 f_raw[p][ii] = genData.fLocal[p];
			}

		}

        if (1e-10 >= genData.density) {
            throw DensityZeroException(
                    "Density too small in collision. Decrease time step size.");
        }

	}

	void collideAll(DistributionFunctions& f, DistributionFunctions& g, distributed_vector& densities,
            vector<distributed_vector>& velocities, distributed_vector& temperature, distributed_vector& maskShockSensor,
			const dealii::IndexSet& locally_owned_dofs,
			const bool inInitializationProcedure,
			GeneralCollisionData<T_D, T_Q>& genData,
			typename T_collision<T_D, T_Q, T_equilibrium>::SpecificCollisionData& specData) const {
		//for all degrees of freedom on current processor
		//dealii::IndexSet::ElementIterator it(locally_owned_dofs.begin());
		//dealii::IndexSet::ElementIterator end(locally_owned_dofs.end());

		T_collision<T_D, T_Q, T_equilibrium> collisionScheme;

		std::array< double*, T_Q> f_raw;
		std::array< double*, T_Q> g_raw;
		int length;
		for (size_t i = 0; i < T_Q; i++){
			f.at(i).trilinos_vector().ExtractView(&f_raw[i], &length);
		}

		for (size_t i = 0; i < T_Q; i++){
			g.at(i).trilinos_vector().ExtractView(&g_raw[i], &length);
		}
		double* rho_raw;
		densities.trilinos_vector().ExtractView(&rho_raw, &length);
		std::array< double*, T_D> u_raw;
		for (size_t i = 0; i < T_D; i++){
			velocities.at(i).trilinos_vector().ExtractView(&u_raw[i], &length);
		}

		double* T_raw;
				temperature.trilinos_vector().ExtractView(&T_raw, &length);

        double* mSS_raw;
                maskShockSensor.trilinos_vector().ExtractView(&mSS_raw, &length);

		genData.H3 = calculateH3<T_D,T_Q>(genData.cs2,genData.e);
		genData.H4 = calculateH4<T_D,T_Q>(genData.cs2,genData.e);


		for (int ii = 0; ii < length; ii++) {

			// Variable that stores the local distribution function values of every node


			// Copy the needed global distribution function values into the local variable
            for (int p = 0; p < T_Q; ++p) {
				genData.fLocal[p] = f_raw[p][ii];
			}
			for (int p = 0; p < T_Q; ++p) {
				genData.gLocal[p] = g_raw[p][ii];
			}

			//copyGlobalToLocalF<T_Q>(genData.fLocal, f, i); // done

			//Calculate the local density and store it into the Parameter Handling System
			genData.density = calculateDensity<T_Q>(genData.fLocal); // done

			//Write the local density to the global density vector
			rho_raw[ii] = genData.density; // write local density to global density vector

			//Calculate the local velocity and store it into the Parameter Handling System
			calculateVelocity<T_D, T_Q>(genData.fLocal, genData.velocity,
					genData.density, genData); // TODO velocities for other stencils

			genData.temperature = calculateTemperature<T_D, T_Q>(genData.fLocal, genData.gLocal, genData.velocity,
                                                                 genData.density, genData.temperature, genData, 0);
			T_raw[ii] = genData.temperature;

			//Write the local density to the global velocity matrix
			if (not inInitializationProcedure) {
				for (size_t j = 0; j < T_D; ++j) {
					u_raw[j][ii] = genData.velocity[j] * genData.scaling;
				}
				if (genData.problemDescription.hasExternalForce()) {
					applyMacroscopicForces<T_D, T_Q>(u_raw, ii, genData);
					applyForces<T_D, T_Q>(genData);
				}
			} else {
				for (size_t j = 0; j < T_D; ++j) {
					genData.velocity[j] = u_raw[j][ii] / genData.scaling;
				}
			}

            //Store the ShockSensorValue
            //genData.maskShockSensor = mSS_raw[ii];

			//Initialize an object of the desired collision scheme and run the relaxation process

           // collisionScheme.relaxWithG(genData.fLocal, genData.gLocal, genData, specData);
            collisionScheme.relaxWithG(genData.fLocal, genData.gLocal, genData, specData);

            mSS_raw[ii] = genData.maskShockSensor;
 			//reApplyForces<T_Q>(fLocal); // TODO

			//Finally copy the updated distribution function back to the global distribution function
			for (int p = 0; p < T_Q; ++p) {
				 f_raw[p][ii] = genData.fLocal[p];
				 g_raw[p][ii] = genData.gLocal[p];
			}

		}

        if (1e-10 >= genData.density) {
            throw DensityZeroException(
                    "Density too small in collision. Decrease time step size.");
        }


    }

};

} /* namespace natrium */

#endif /* COLLISIONOPERATOR_H_ */
