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


namespace natrium {
    template<int T_D, int T_Q>
    class QuarticEquilibrium {
        std::array<std::array<std::array<std::array<double, T_Q>, T_D>, T_D>, T_D> m_H3 = {{{{0.0}}}};

    public:
        void calc(std::array<double, T_Q> &feq, const GeneralCollisionData<T_D, T_Q> &params);
    };


    template<int T_D, int T_Q>
    inline void QuarticEquilibrium<T_D, T_Q>::calc(std::array<double, T_Q> &feq,
                                                   const GeneralCollisionData<T_D, T_Q> &p) {
        double eye[2][2] = {{1, 0},
                            {0, 1}};
        double uu_term = 0.0;
        for (size_t j = 0; j < T_D; j++) {
            uu_term += -(p.velocity[j] * p.velocity[j])
                       / (2.0 * p.cs2);
        }

        for (size_t i = 0; i < T_Q; i++) {

            double T1 = p.cs2 * (p.temperature - 1);

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
                    /*            for (int gam = 0; gam < T_D; gam++) {
                                    for (int det = 0; det < T_D; det++) {
                                        double power4 = p.e[i][alp] * p.e[i][bet] * p.e[i][gam] * p.e[i][det];
                                        double power2 = p.e[i][alp] * p.e[i][bet] * eye[gam][det]
                                                        + p.e[i][alp] * p.e[i][gam] * eye[bet][det]
                                                        + p.e[i][alp] * p.e[i][det] * eye[bet][gam]
                                                        + p.e[i][bet] * p.e[i][gam] * eye[alp][det]
                                                        + p.e[i][bet] * p.e[i][det] * eye[alp][gam]
                                                        + p.e[i][gam] * p.e[i][det] * eye[alp][bet];
                                        double power0 = eye[alp][bet] * eye[gam][det] + eye[alp][gam] * eye[bet][det] +
                                                        eye[alp][det] * eye[bet][gam];
                                        double u4 = p.velocity[alp] * p.velocity[bet] * p.velocity[gam] * p.velocity[det];
                                        double u2 = p.velocity[alp] * p.velocity[bet] * eye[gam][det] +
                                                    p.velocity[alp] * p.velocity[gam] * eye[bet][det] +
                                                    p.velocity[alp] * p.velocity[det] * eye[bet][gam] +
                                                    p.velocity[bet] * p.velocity[gam] * eye[alp][det] +
                                                    p.velocity[bet] * p.velocity[det] * eye[alp][gam] +
                                                    p.velocity[gam] * p.velocity[det] * eye[alp][bet];
                                        double multieye = eye[alp][bet] * eye[gam][det] + eye[alp][gam] * eye[bet][det] +
                                                          eye[alp][det] * eye[bet][gam];

                                        feq[i] += p.weight[i] * p.density / (24. * p.cs2 * p.cs2 * p.cs2 * p.cs2) *
                                                  (power4 - p.cs2 * power2 + p.cs2 * p.cs2 * power0) *
                                                  (u4 + T1 * (u2 + T1 * multieye));


                                    }

                                }*/
                }
            }

            double H_xxx = p.e[i][0] * p.e[i][0] * p.e[i][0] - p.cs2 * (p.e[i][0] * 3.0);
            double H_xxy = p.e[i][0] * p.e[i][0] * p.e[i][1] - p.cs2 * p.e[i][1];
            double H_xyy = p.e[i][0] * p.e[i][1] * p.e[i][1] - p.cs2 * p.e[i][0];
            double H_yyy = p.e[i][1] * p.e[i][1] * p.e[i][1] - p.cs2 * (p.e[i][1] * 3.0);
            double a_xxx = p.velocity[0] * p.velocity[0] * p.velocity[0]
                           + T1 * (p.velocity[0] + p.velocity[0] + p.velocity[0]);
            double a_xxy = p.velocity[0] * p.velocity[0] * p.velocity[1]
                           + T1 * (p.velocity[1]);
            double a_xyy = p.velocity[0] * p.velocity[1] * p.velocity[1]
                           + T1 * (p.velocity[0]);
            double a_yyy = p.velocity[1] * p.velocity[1] * p.velocity[1]
                           + T1 * (p.velocity[1] + p.velocity[1] + p.velocity[1]);
            feq[i] += p.weight[i] * p.density / (6. * p.cs2 * p.cs2 * p.cs2) *
                      (a_xxx * H_xxx + 3 * (a_xxy * H_xxy + a_xyy * H_xyy) + a_yyy * H_yyy);

            if (T_D == 3) {
                double H_zzz = p.e[i][2] * p.e[i][2] * p.e[i][2] - p.cs2 * (p.e[i][2] * 3.0);
                double H_xxz = p.e[i][0] * p.e[i][0] * p.e[i][2] - p.cs2 * p.e[i][2];
                double H_xzz = p.e[i][0] * p.e[i][2] * p.e[i][2] - p.cs2 * p.e[i][2];
                double H_yzz = p.e[i][1] * p.e[i][2] * p.e[i][2] - p.cs2 * p.e[i][1];
                double H_yyz = p.e[i][1] * p.e[i][1] * p.e[i][2] - p.cs2 * p.e[i][2];
                double H_xyz = p.e[i][0] * p.e[i][1] * p.e[i][2];

                double a_zzz = p.velocity[2] * p.velocity[2] * p.velocity[2]
                               + T1 * (p.velocity[2] + p.velocity[2] + p.velocity[2]);
                double a_xxz = p.velocity[0] * p.velocity[0] * p.velocity[2]
                               + T1 * (p.velocity[2]);
                double a_xzz = p.velocity[0] * p.velocity[2] * p.velocity[2]
                               + T1 * (p.velocity[0]);
                double a_yzz = p.velocity[1] * p.velocity[2] * p.velocity[2]
                               + T1 * (p.velocity[1]);
                double a_yyz = p.velocity[1] * p.velocity[1] * p.velocity[2]
                               + T1 * (p.velocity[2]);
                double a_xyz = p.velocity[0] * p.velocity[1] * p.velocity[2];


                feq[i] += p.weight[i] * p.density / (6. * p.cs2 * p.cs2 * p.cs2) *
                          (a_zzz * H_zzz + 3 * (a_xxz * H_xxz + a_xzz * H_xzz + a_yzz * H_yzz + a_yyz * H_yyz) +
                           a_xyz * H_xyz);
            }

            double H_xxxx = p.e[i][0] * p.e[i][0] * p.e[i][0] * p.e[i][0] - p.cs2 * p.e[i][0] * p.e[i][0] * 6.0 +
                            p.cs2 * p.cs2 * 3.0;
            double H_yyyy = p.e[i][1] * p.e[i][1] * p.e[i][1] * p.e[i][1] - p.cs2 * p.e[i][1] * p.e[i][1] * 6.0 +
                            p.cs2 * p.cs2 * 3.0;
            double H_xxxy = p.e[i][0] * p.e[i][0] * p.e[i][0] * p.e[i][1] -
                            p.cs2 * (p.e[i][0] * p.e[i][1] * 3.0);
            double H_xyyy = p.e[i][0] * p.e[i][1] * p.e[i][1] * p.e[i][1] -
                            p.cs2 * (p.e[i][0] * p.e[i][1] * 3.0);
            double H_xxyy = p.e[i][0] * p.e[i][0] * p.e[i][1] * p.e[i][1] -
                            p.cs2 * (p.e[i][0] * p.e[i][0] + p.e[i][1] * p.e[i][1]) + p.cs2 * p.cs2;

            /*            for (int gam = 0; gam < T_D; gam++) {
                                    for (int det = 0; det < T_D; det++) {
                                        double power4 = p.e[i][alp] * p.e[i][bet] * p.e[i][gam] * p.e[i][det];
                                        double power2 = p.e[i][alp] * p.e[i][bet] * eye[gam][det]
                                                        + p.e[i][alp] * p.e[i][gam] * eye[bet][det]
                                                        + p.e[i][alp] * p.e[i][det] * eye[bet][gam]
                                                        + p.e[i][bet] * p.e[i][gam] * eye[alp][det]
                                                        + p.e[i][bet] * p.e[i][det] * eye[alp][gam]
                                                        + p.e[i][gam] * p.e[i][det] * eye[alp][bet];
                                        double power0 = eye[alp][bet] * eye[gam][det] + eye[alp][gam] * eye[bet][det] +
                                                        eye[alp][det] * eye[bet][gam];
                                        double u4 = p.velocity[alp] * p.velocity[bet] * p.velocity[gam] * p.velocity[det];
                                        double u2 = p.velocity[alp] * p.velocity[bet] * eye[gam][det] +
                                                    p.velocity[alp] * p.velocity[gam] * eye[bet][det] +
                                                    p.velocity[alp] * p.velocity[det] * eye[bet][gam] +
                                                    p.velocity[bet] * p.velocity[gam] * eye[alp][det] +
                                                    p.velocity[bet] * p.velocity[det] * eye[alp][gam] +
                                                    p.velocity[gam] * p.velocity[det] * eye[alp][bet];
                                        double multieye = eye[alp][bet] * eye[gam][det] + eye[alp][gam] * eye[bet][det] +
                                                          eye[alp][det] * eye[bet][gam];

                                        feq[i] += p.weight[i] * p.density / (24. * p.cs2 * p.cs2 * p.cs2 * p.cs2) *
                                                  (power4 - p.cs2 * power2 + p.cs2 * p.cs2 * power0) *
                                                  (u4 + T1 * (u2 + T1 * multieye));


                                    }

                                }*/

            double a_xxxx = p.velocity[0] * p.velocity[0] * p.velocity[0] * p.velocity[0] +
                            T1 * p.velocity[0] * p.velocity[0] * 6.0 + T1 * T1 * 3.0;
            double a_yyyy = p.velocity[1] * p.velocity[1] * p.velocity[1] * p.velocity[1] +
                            T1 * p.velocity[1] * p.velocity[1] * 6.0 + T1 * T1 * 3.0;
            double a_xxxy = p.velocity[0] * p.velocity[0] * p.velocity[0] * p.velocity[1] +
                            T1 * (p.velocity[0] * p.velocity[1] * 3.0 );
            double a_xyyy = p.velocity[0] * p.velocity[1] * p.velocity[1] * p.velocity[1] +
                            T1 * (p.velocity[0] * p.velocity[1] * 3.0 );
            double a_xxyy = p.velocity[0] * p.velocity[0] * p.velocity[1] * p.velocity[1] +
                            T1 * (p.velocity[0] * p.velocity[0] + p.velocity[1] * p.velocity[1]) + T1 * T1;


            feq[i] += p.weight[i] * p.density / (24. * p.cs2 * p.cs2 * p.cs2 * p.cs2) *
                     (H_xxxx * a_xxxx + H_yyyy * a_yyyy + 6.0 * H_xxyy * a_xxyy + 4.0 * H_xyyy * a_xyyy +
                       4.0 * H_xxxy * a_xxxy);
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
