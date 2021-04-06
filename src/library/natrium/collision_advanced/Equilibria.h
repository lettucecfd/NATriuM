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
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };


    template<int T_D, int T_Q>
    inline void QuarticEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                   const GeneralCollisionData<T_D, T_Q> &p) {
        const std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();

        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(p.velocity[j] * p.velocity[j])
                       / (2.0 * p.cs2);
        }

        const double T1 = p.cs2 * (p.temperature - 1);

        const double a_xxx = p.velocity[0] * p.velocity[0] * p.velocity[0]
                        + T1 * (p.velocity[0] + p.velocity[0] + p.velocity[0]);
        const double a_xxy = p.velocity[0] * p.velocity[0] * p.velocity[1]
                        + T1 * (p.velocity[1]);
        const double a_xyy = p.velocity[0] * p.velocity[1] * p.velocity[1]
                        + T1 * (p.velocity[0]);
        const double a_yyy = p.velocity[1] * p.velocity[1] * p.velocity[1]
                        + T1 * (p.velocity[1] + p.velocity[1] + p.velocity[1]);

        const double a_xxxx = p.velocity[0] * p.velocity[0] * p.velocity[0] * p.velocity[0] +
                         T1 * p.velocity[0] * p.velocity[0] * 6.0 + T1 * T1 * 3.0;
        const double a_yyyy = p.velocity[1] * p.velocity[1] * p.velocity[1] * p.velocity[1] +
                         T1 * p.velocity[1] * p.velocity[1] * 6.0 + T1 * T1 * 3.0;
        const double a_xxxy = p.velocity[0] * p.velocity[0] * p.velocity[0] * p.velocity[1] +
                         T1 * (p.velocity[0] * p.velocity[1] * 3.0 );
        const double a_xyyy = p.velocity[0] * p.velocity[1] * p.velocity[1] * p.velocity[1] +
                         T1 * (p.velocity[0] * p.velocity[1] * 3.0 );
        const double a_xxyy = p.velocity[0] * p.velocity[0] * p.velocity[1] * p.velocity[1] +
                        T1 * (p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1]) + T1 * T1;

        double a_zzz= 0.0, a_xxz= 0.0,a_xzz= 0.0,a_yzz= 0.0,a_yyz= 0.0,a_xyz  = 0.0;
        double a_zzzz= 0.0, a_xzzz= 0.0, a_xxzz= 0.0, a_xxxz= 0.0, a_yzzz= 0.0, a_yyzz= 0.0, a_yyyz= 0.0, a_xxyz= 0.0, a_xyyz= 0.0, a_xyzz = 0.0;

        if(T_D==3) {
            a_zzz = p.velocity[2] * p.velocity[2] * p.velocity[2]
                           + T1 * (p.velocity[2] + p.velocity[2] + p.velocity[2]);
            a_xxz = p.velocity[0] * p.velocity[0] * p.velocity[2]
                           + T1 * (p.velocity[2]);
            a_xzz = p.velocity[0] * p.velocity[2] * p.velocity[2]
                           + T1 * (p.velocity[0]);
            a_yzz = p.velocity[1] * p.velocity[2] * p.velocity[2]
                           + T1 * (p.velocity[1]);
            a_yyz = p.velocity[1] * p.velocity[1] * p.velocity[2]
                           + T1 * (p.velocity[2]);
            a_xyz = p.velocity[0] * p.velocity[1] * p.velocity[2];

            a_zzzz = p.velocity[2] * p.velocity[2] * p.velocity[2] * p.velocity[2] +
                            T1 * p.velocity[2] * p.velocity[2] * 6.0 + T1 * T1 * 3.0;
            a_xxxz = p.velocity[0] * p.velocity[0] * p.velocity[0] * p.velocity[2] +
                            T1 * (p.velocity[0] * p.velocity[2] * 3.0 );
            a_yyyz = p.velocity[1] * p.velocity[1] * p.velocity[1] * p.velocity[2] +
                            T1 * (p.velocity[1] * p.velocity[2] * 3.0 );
            a_xzzz = p.velocity[0] * p.velocity[2] * p.velocity[2] * p.velocity[2] +
                            T1 * (p.velocity[0] * p.velocity[2] * 3.0 );
            a_yzzz = p.velocity[1] * p.velocity[2] * p.velocity[2] * p.velocity[2] +
                            T1 * (p.velocity[1] * p.velocity[2] * 3.0 );
            a_xxzz = p.velocity[0] * p.velocity[0] * p.velocity[2] * p.velocity[2] +
                            T1 * (p.velocity[0] * p.velocity[0] + p.velocity[2] * p.velocity[2]) + T1 * T1;
            a_yyzz = p.velocity[1] * p.velocity[1] * p.velocity[2] * p.velocity[2] +
                            T1 * (p.velocity[1] * p.velocity[1] + p.velocity[2] * p.velocity[2]) + T1 * T1;
            a_xyzz = p.velocity[0] * p.velocity[1] * p.velocity[2] * p.velocity[2] +
                            T1 * (p.velocity[0] * p.velocity[1]);
            a_xyyz = p.velocity[0] * p.velocity[1] * p.velocity[1] * p.velocity[2] +
                            T1 * (p.velocity[0] * p.velocity[2]);
            a_xxyz = p.velocity[0] * p.velocity[0] * p.velocity[1] * p.velocity[2] +
                            T1 * (p.velocity[1] * p.velocity[2]);

        }
        for (size_t i = 0; i < T_Q; i++) {

            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (p.velocity[j] * p.e[i][j]) / p.cs2;
            }
            feq[i] = p.weight[i] * p.density * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
            for (int alp = 0; alp < T_D; alp++) {
                for (int bet = 0; bet < T_D; bet++) {
                    feq[i] += p.density * p.weight[i] / (2.0 * p.cs2) *
                              ((p.temperature - 1) * eye[alp][bet] * p.e[i][alp] * p.e[i][bet] -
                               p.cs2 * eye[alp][bet] * (p.temperature - 1));

                }
            }

            const double H_xxx = p.H3[i][0][0][0];
            const double H_xxy = p.H3[i][0][0][1];
            const double H_xyy = p.H3[i][0][1][1];
            const double H_yyy = p.H3[i][1][1][1];

            feq[i] += p.weight[i] * p.density / (6. * p.cs2 * p.cs2 * p.cs2) *
                      (a_xxx * H_xxx + 3 * (a_xxy * H_xxy + a_xyy * H_xyy) + a_yyy * H_yyy);

            if (T_D == 3) {
                const double H_zzz = p.H3[i][2][2][2];
                const double H_xxz = p.H3[i][0][0][2];
                const double H_xzz = p.H3[i][0][2][2];
                const double H_yzz = p.H3[i][1][2][2];
                const double H_yyz = p.H3[i][1][1][2];
                const double H_xyz = p.H3[i][0][1][2];




                feq[i] += p.weight[i] * p.density / (6. * p.cs2 * p.cs2 * p.cs2) *
                          (a_zzz * H_zzz + 3 * (a_xxz * H_xxz + a_xzz * H_xzz + a_yzz * H_yzz + a_yyz * H_yyz) +
                           6.0*a_xyz * H_xyz);
            }

            const double H_xxxx = p.H4[i][0][0][0][0];
            const double H_yyyy = p.H4[i][1][1][1][1];
            const double H_xxxy = p.H4[i][0][0][0][1];
            const double H_xyyy = p.H4[i][0][1][1][1];
            const double H_xxyy = p.H4[i][0][0][1][1];


            feq[i] += p.weight[i] * p.density / (24. * p.cs2 * p.cs2 * p.cs2 * p.cs2) *
                     (H_xxxx * a_xxxx + H_yyyy * a_yyyy + 6.0 * H_xxyy * a_xxyy + 4.0 * H_xyyy * a_xyyy +
                       4.0 * H_xxxy * a_xxxy);

            if(T_D==3){
                const double H_zzzz = p.H4[i][2][2][2][2];
                const double H_xzzz = p.H4[i][0][2][2][2];
                const double H_xxzz = p.H4[i][0][0][2][2];
                const double H_xxxz = p.H4[i][0][0][0][2];
                const double H_yzzz = p.H4[i][1][2][2][2];
                const double H_yyzz = p.H4[i][1][1][2][2];
                const double H_yyyz = p.H4[i][1][1][1][2];
                const double H_xxyz = p.H4[i][0][0][1][2];
                const double H_xyyz = p.H4[i][0][1][1][2];
                const double H_xyzz = p.H4[i][0][1][2][2];

                feq[i] += p.weight[i] * p.density / (24. * p.cs2 * p.cs2 * p.cs2 * p.cs2) *
                          (H_zzzz * a_zzzz + 4.0 * (H_xzzz * a_xzzz + H_yzzz * a_yzzz+ H_xxxz * a_xxxz + H_yyyz * a_yyyz) + 6.0 * (H_xxzz * a_xxzz + H_yyzz * a_yyzz)
                          + 12.0 * (H_xxyz*a_xxyz + H_xyyz*a_xyyz+H_xyzz*a_xyzz));

            }
        }
    }

}







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
