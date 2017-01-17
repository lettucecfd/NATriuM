/*
 * AuxiliaryCollisionFunctions.h
 *
 *  Created on: 16.01.2017
 *      Author: natrium
 */
#ifndef AUXILIARYCOLLISIONFUNCTIONS_H_
#define AUXILIARYCOLLISIONFUNCTIONS_H_

namespace natrium {

template <int T_Q>
inline void calculateDensity(double& density, double fLocal[])
{
	density=0;
	for (int p = 0; p < T_Q; ++p)
	{
		density+=fLocal[p];
	}
}

template <int T_Q>
inline void copyGlobalToLocalF(double fLocal[], DistributionFunctions& f, size_t i)
{
	for (int p = 0; p < T_Q; ++p)
		{
			fLocal[p] = f[p][i];
		}
}


} /* namespace natrium */

#endif /* AUXILIARYCOLLISIONFUNCTIONS_H_ */
