/*
 * Equilibria.h
 *
 *  Created on: 19.01.2017
 *      Author: natrium
 */

#ifndef LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_
#define LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_
//#include "CollisionOperator.h"
#include "AuxiliaryCollisionFunctions.h"
#include "Equilibria.h"
#include "CollisionSchemes.h"
#include <array>


namespace natrium {
    template<int T_D, int T_Q>
    class BGKEquilibrium {

    public:
        //BGKEquilibrium();
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };


    template<>
    inline void BGKEquilibrium<2, 9>::calc(std::array<double, 9> &feq, const GeneralCollisionData<2, 9> &params) {
        double scalar_product;
        double weighting;
        double uSquareTerm;
        double mixedTerm;
        double prefactor = 1. / params.cs2;
        // calculate equilibrium distribution
        scalar_product = params.velocity[0] * params.velocity[0]
                         + params.velocity[1] * params.velocity[1];
        uSquareTerm = -scalar_product / (2 * params.cs2);
        // direction 0
        weighting = 4. / 9. * params.density;
        feq[0] = weighting * (1 + uSquareTerm);
        // directions 1-4
        weighting = 1. / 9. * params.density;
        mixedTerm = prefactor * (params.velocity[0]);
        feq[1] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
        feq[3] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
        mixedTerm = prefactor * (params.velocity[1]);
        feq[2] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
        feq[4] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
        // directions 5-8
        weighting = 1. / 36. * params.density;
        mixedTerm = prefactor * (params.velocity[0] + params.velocity[1]);
        feq[5] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
        feq[7] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
        mixedTerm = prefactor * (-params.velocity[0] + params.velocity[1]);
        feq[6] = weighting * (1 + mixedTerm * (1 + 0.5 * mixedTerm) + uSquareTerm);
        feq[8] = weighting * (1 - mixedTerm * (1 - 0.5 * mixedTerm) + uSquareTerm);
    }


    template<int T_D, int T_Q>
    inline void BGKEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                               const GeneralCollisionData<T_D, T_Q> &params) {

        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(params.velocity[j] * params.velocity[j])
                       / (2.0 * params.cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {
            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
            }
            feq[i] = params.weight[i] * params.density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
        }
    } /* BGKEquilibrium<T_D, T_Q>::calc */

} /* class BGKEquilibrium */

namespace natrium{
template <int T_D,int T_Q>
class QuarticEquilibrium
{

public:
    //QuarticEquilibrium();
    void calc(std::array<double, T_Q>& feq, const GeneralCollisionData<T_D,T_Q> & params);
};



template<int T_D, int T_Q>
inline void QuarticEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q>& feq,
                                           const GeneralCollisionData<T_D, T_Q> & p) {
    double eye[2][2]={{1,0},{0,1}};
    double uu_term = 0.0;
    double T1 = p.cs2*(p.temperature-1);


    std::array<std::array<std::array<double, T_D>, T_D>, T_D> vel3;
    for (int alp = 0; alp < T_D; alp++) {
        for (int bet = 0; bet < T_D; bet++) {
            for (int gam = 0; gam < T_D; gam++) {
                vel3[alp][bet][gam]   = (p.velocity[alp] * p.velocity[bet] * p.velocity[gam]
                                         + T1 *
                                           (eye[alp][bet] * p.velocity[gam] + eye[bet][gam] * p.velocity[alp] +
                                            eye[alp][gam] * p.velocity[bet]));
            }
        }
    }

    std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_D> vel4;
    for (int alp = 0; alp < T_D; alp++) {
        for (int bet = 0; bet < T_D; bet++) {
            for (int gam = 0; gam < T_D; gam++) {
                for (int det = 0; det < T_D; det++) {

                    double u4    = p.velocity[alp]*p.velocity[bet]*p.velocity[gam]*p.velocity[det];
                    double u2 = p.velocity[alp]*p.velocity[bet]*eye[gam][det]+p.velocity[alp]*p.velocity[gam]*eye[bet][det]+p.velocity[alp]*p.velocity[det]*eye[bet][gam]+p.velocity[bet]*p.velocity[gam]*eye[alp][det]+p.velocity[bet]*p.velocity[det]*eye[alp][gam]+p.velocity[gam]*p.velocity[det]*eye[alp][bet];
                    double multieye= eye[alp][bet]*eye[gam][det]+eye[alp][gam]*eye[bet][det]+eye[alp][det]*eye[bet][gam];
                    vel4[alp][bet][gam][det] = (u4+T1*(u2+T1*multieye));

                }
            }
        }
    }




    for (size_t j = 0; j < T_D; j++) {
        uu_term += -(p.velocity[j] * p.velocity[j])
                   / (2.0 * p.cs2);
    }

    for (size_t i = 0; i < T_Q; i++) {

        double ue_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            ue_term += (p.velocity[j] * p.e[i][j]) / p.cs2;
        }
        feq[i] = p.weight[i] * p.density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
        for (int alp = 0; alp < T_D; alp++){
            for (int bet = 0; bet < T_D; bet++){
                feq[i]+=p.density*p.weight[i]/(2.0*p.cs2)*((p.temperature-1)*eye[alp][bet]*p.e[i][alp]*p.e[i][bet]-p.cs2*eye[alp][bet]*(p.temperature-1));
                for (int gam = 0; gam < T_D; gam++){

                    feq[i] += p.weight[i] * p.density / (6. * p.cs2 * p.cs2 * p.cs2) *
                              vel3[alp][bet][gam] * (p.e[i][alp] * p.e[i][bet] * p.e[i][gam] - p.cs2 * (p.e[i][gam] *
                                                                                                                  eye[alp][bet] +
                                                                                                                  p.e[i][bet] *
                                                                                                                  eye[alp][gam] +
                                                                                                                  p.e[i][alp] *
                                                                                                                  eye[bet][gam]));

                    for (int det = 0; det < T_D; det++)
                    {
                        double power4 = p.e[i][alp]*p.e[i][bet]*p.e[i][gam]*p.e[i][det];
                        double power2 = p.e[i][alp]*p.e[i][bet]*eye[gam][det]
                                        +p.e[i][alp]*p.e[i][gam]*eye[bet][det]
                                        +p.e[i][alp]*p.e[i][det]*eye[bet][gam]
                                        +p.e[i][bet]*p.e[i][gam]*eye[alp][det]
                                        +p.e[i][bet]*p.e[i][det]*eye[alp][gam]
                                        +p.e[i][gam]*p.e[i][det]*eye[alp][bet];
                        double power0 = eye[alp][bet]*eye[gam][det]+eye[alp][gam]*eye[bet][det]+eye[alp][det]*eye[bet][gam];

                        feq[i]+= p.weight[i] * p.density /(24.*p.cs2*p.cs2*p.cs2*p.cs2)*(power4-p.cs2*power2+p.cs2*p.cs2*power0)*vel4[alp][bet][gam][det];


                    }

                }
            }

        }
    }
} /* QuarticEquilibrium<T_D, T_Q>::calc */

} /* class QuarticEquilibrium */




// ==========================================================================================

namespace natrium {
    template<int T_D, int T_Q>
    class SteadyStateEquilibrium {

    public:
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };

    template<int T_D, int T_Q>
    inline void SteadyStateEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                       const GeneralCollisionData<T_D, T_Q> &params) {

        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(params.velocity[j] * params.velocity[j])
                       / (2.0 * params.cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {
            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
            }
            feq[i] = params.weight[i] * params.density * (1 + ue_term * (1 + 0.5 * (ue_term) / params.gamma_steadystate)
                                                          + uu_term / params.gamma_steadystate);
        }
    } /* SteadyStateEquilibrium<T_D, T_Q>::calc */

} /* class SteadyStateEquilibrium */

// ==========================================================================================

namespace natrium {
    template<int T_D, int T_Q>
    class IncompressibleEquilibrium {

    public:
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };

    template<int T_D, int T_Q>
    inline void IncompressibleEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                          const GeneralCollisionData<T_D, T_Q> &params) {

        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(params.velocity[j] * params.velocity[j])
                       / (2.0 * params.cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {
            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
            }
            feq[i] = params.weight[i] * (params.density + ue_term * (1 + 0.5 * (ue_term))
                                         + uu_term);
        }
    } /* IncompressibleEquilibrium<T_D, T_Q>::calc */

} /* class IncompressibleEquilibrium */


// ==========================================================================================

namespace natrium {
    template<int T_D, int T_Q>
    class EntropicEquilibrium {

    public:
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };

    template<int T_D, int T_Q>
    inline void EntropicEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                    const GeneralCollisionData<T_D, T_Q> &params) {

        throw NotImplementedException("Entropic equilibrium is not implemented, yet.");
        /*double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(params.velocity[j] * params.velocity[j])
                    / (2.0 * params.cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {
            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (params.velocity[j] * params.e[i][j]) / params.cs2;
            }
            feq[i] = params.weight[i] * (params.density + ue_term * (1 + 0.5 * (ue_term) )
                    + uu_term );
        }*/
    } /* EntropicEquilibrium<T_D, T_Q>::calc */

} /* class EntropicEquilibrium */


#endif /* LIBRARY_NATRIUM_COLLISION_EQUILIBRIA_H_ */
