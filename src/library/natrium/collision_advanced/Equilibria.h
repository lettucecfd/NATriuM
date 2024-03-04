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
        BGKEquilibrium(double cs2, std::array<std::array<double, T_D>, T_Q> e)
        {
            (void) cs2;
            (void) e;
        }

        BGKEquilibrium(const GeneralCollisionData<T_D, T_Q> &params)
        {
            (void) params;
        }
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
        void polynomial(std::array<double, T_Q> &feq,
                        double rho, std::array<double, T_D> u, double T,
                        std::array<std::array<double, T_D>, T_Q> e, std::array<double, T_Q> w,
                        double cs2);
    const std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_Q> H3;
    const std::array<std::array<std::array<std::array<std::array<double, T_D>, T_D>, T_D>,T_D>,T_Q> H4;

    QuarticEquilibrium(double cs2, std::array<std::array<double, T_D>, T_Q> e) : H3(calculateH3(cs2,e)), H4(calculateH4(cs2,e))
    { }

    QuarticEquilibrium(const GeneralCollisionData<T_D, T_Q> &params) : H3(params.H3), H4(params.H4)
    { }
private:
    };

    template<int T_D, int T_Q>
    inline void QuarticEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                   const GeneralCollisionData<T_D, T_Q> &p) {
        polynomial(feq,p.density,p.velocity,p.temperature,p.e,p.weight,p.cs2);
    }

    template<int T_D, int T_Q>
    inline void QuarticEquilibrium<T_D, T_Q>::polynomial(std::array<double, T_Q> &feq,
                                               double rho, std::array<double, T_D> u, double T,
                                               std::array<std::array<double, T_D>, T_Q> e, std::array<double, T_Q> w,
                                               double cs2) {

        const std::array<std::array<size_t,T_D>, T_D> eye = unity_matrix<T_D>();

        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(u[j] * u[j])
                       / (2.0 * cs2);
        }

        const double T1 = cs2 * (T - 1);

        const double a_xxx = u[0] * u[0] * u[0]
                        + T1 * (u[0] + u[0] + u[0]);
        const double a_xxy = u[0] * u[0] * u[1]
                        + T1 * (u[1]);
        const double a_xyy = u[0] * u[1] * u[1]
                        + T1 * (u[0]);
        const double a_yyy = u[1] * u[1] * u[1]
                        + T1 * (u[1] + u[1] + u[1]);

        const double a_xxxx = u[0] * u[0] * u[0] * u[0] +
                         T1 * u[0] * u[0] * 6.0 + T1 * T1 * 3.0;
        const double a_yyyy = u[1] * u[1] * u[1] * u[1] +
                         T1 * u[1] * u[1] * 6.0 + T1 * T1 * 3.0;
        const double a_xxxy = u[0] * u[0] * u[0] * u[1] +
                         T1 * (u[0] * u[1] * 3.0 );
        const double a_xyyy = u[0] * u[1] * u[1] * u[1] +
                         T1 * (u[0] * u[1] * 3.0 );
        const double a_xxyy = u[0] * u[0] * u[1] * u[1] +
                        T1 * (u[0] * u[0] + u[1] * u[1]) + T1 * T1;

        double a_zzz= 0.0, a_xxz= 0.0,a_xzz= 0.0,a_yzz= 0.0,a_yyz= 0.0,a_xyz  = 0.0;
        double a_zzzz= 0.0, a_xzzz= 0.0, a_xxzz= 0.0, a_xxxz= 0.0, a_yzzz= 0.0, a_yyzz= 0.0, a_yyyz= 0.0, a_xxyz= 0.0, a_xyyz= 0.0, a_xyzz = 0.0;

        if(T_D==3) {
            a_zzz = u[2] * u[2] * u[2]
                           + T1 * (u[2] + u[2] + u[2]);
            a_xxz = u[0] * u[0] * u[2]
                           + T1 * (u[2]);
            a_xzz = u[0] * u[2] * u[2]
                           + T1 * (u[0]);
            a_yzz = u[1] * u[2] * u[2]
                           + T1 * (u[1]);
            a_yyz = u[1] * u[1] * u[2]
                           + T1 * (u[2]);
            a_xyz = u[0] * u[1] * u[2];

            a_zzzz = u[2] * u[2] * u[2] * u[2] +
                            T1 * u[2] * u[2] * 6.0 + T1 * T1 * 3.0;
            a_xxxz = u[0] * u[0] * u[0] * u[2] +
                            T1 * (u[0] * u[2] * 3.0 );
            a_yyyz = u[1] * u[1] * u[1] * u[2] +
                            T1 * (u[1] * u[2] * 3.0 );
            a_xzzz = u[0] * u[2] * u[2] * u[2] +
                            T1 * (u[0] * u[2] * 3.0 );
            a_yzzz = u[1] * u[2] * u[2] * u[2] +
                            T1 * (u[1] * u[2] * 3.0 );
            a_xxzz = u[0] * u[0] * u[2] * u[2] +
                            T1 * (u[0] * u[0] + u[2] * u[2]) + T1 * T1;
            a_yyzz = u[1] * u[1] * u[2] * u[2] +
                            T1 * (u[1] * u[1] + u[2] * u[2]) + T1 * T1;
            a_xyzz = u[0] * u[1] * u[2] * u[2] +
                            T1 * (u[0] * u[1]);
            a_xyyz = u[0] * u[1] * u[1] * u[2] +
                            T1 * (u[0] * u[2]);
            a_xxyz = u[0] * u[0] * u[1] * u[2] +
                            T1 * (u[1] * u[2]);
        }
        for (size_t i = 0; i < T_Q; i++) {

            double ue_term = 0.0;
            for (size_t j = 0; j < T_D; j++) {
                ue_term += (u[j] * e[i][j]) / cs2;
            }
            feq[i] = w[i] * rho * (1 + ue_term * (1 + 0.5 * (ue_term)) + uu_term);
            for (int alp = 0; alp < T_D; alp++) {
                for (int bet = 0; bet < T_D; bet++) {
                    feq[i] += rho * w[i] / (2.0 * cs2) *
                              ((T - 1) * eye[alp][bet] * e[i][alp] * e[i][bet] -
                               cs2 * eye[alp][bet] * (T - 1));

                }
            }

            const double H_xxx = H3[i][0][0][0];
            const double H_xxy = H3[i][0][0][1];
            const double H_xyy = H3[i][0][1][1];
            const double H_yyy = H3[i][1][1][1];

            feq[i] += w[i] * rho / (6. * cs2 * cs2 * cs2) *
                      (a_xxx * H_xxx + 3 * (a_xxy * H_xxy + a_xyy * H_xyy) + a_yyy * H_yyy);

            if (T_D == 3) {
                const double H_zzz = H3[i][2][2][2];
                const double H_xxz = H3[i][0][0][2];
                const double H_xzz = H3[i][0][2][2];
                const double H_yzz = H3[i][1][2][2];
                const double H_yyz = H3[i][1][1][2];
                const double H_xyz = H3[i][0][1][2];

                feq[i] += w[i] * rho / (6. * cs2 * cs2 * cs2) *
                          (a_zzz * H_zzz + 3 * (a_xxz * H_xxz + a_xzz * H_xzz + a_yzz * H_yzz + a_yyz * H_yyz) +
                           6.0*a_xyz * H_xyz);
            }

            const double H_xxxx = H4[i][0][0][0][0];
            const double H_yyyy = H4[i][1][1][1][1];
            const double H_xxxy = H4[i][0][0][0][1];
            const double H_xyyy = H4[i][0][1][1][1];
            const double H_xxyy = H4[i][0][0][1][1];


            feq[i] += w[i] * rho / (24. * cs2 * cs2 * cs2 * cs2) *
                            (H_xxxx * a_xxxx
                           + H_yyyy * a_yyyy
                     + 6.0 * H_xxyy * a_xxyy
                     + 4.0 * H_xyyy * a_xyyy
                     + 4.0 * H_xxxy * a_xxxy);

            if(T_D==3){
                const double H_zzzz = H4[i][2][2][2][2];
                const double H_xzzz = H4[i][0][2][2][2];
                const double H_xxzz = H4[i][0][0][2][2];
                const double H_xxxz = H4[i][0][0][0][2];
                const double H_yzzz = H4[i][1][2][2][2];
                const double H_yyzz = H4[i][1][1][2][2];
                const double H_yyyz = H4[i][1][1][1][2];
                const double H_xxyz = H4[i][0][0][1][2];
                const double H_xyyz = H4[i][0][1][1][2];
                const double H_xyzz = H4[i][0][1][2][2];

                feq[i] += w[i] * rho / (24. * cs2 * cs2 * cs2 * cs2) *
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
        SteadyStateEquilibrium(double cs2, std::array<std::array<double, T_D>, T_Q> e)
        {
        }
        SteadyStateEquilibrium(const GeneralCollisionData<T_D, T_Q> &params)
        {
        }
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
        IncompressibleEquilibrium(double cs2, std::array<std::array<double, T_D>, T_Q> e)
        {
        }
        IncompressibleEquilibrium(const GeneralCollisionData<T_D, T_Q> &params)
        {
        }
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
        EntropicEquilibrium(double cs2, std::array<std::array<double, T_D>, T_Q> e)
        {
        }
        EntropicEquilibrium(const GeneralCollisionData<T_D, T_Q> &params)
        {
        }
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
