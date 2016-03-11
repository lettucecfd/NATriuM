/*
 * D3Q15Model.h
 *
 *  Created on: Nov 20, 2015
 *      Author: kraemer
 */

#ifndef D3Q15MODEL_H_
#define D3Q15MODEL_H_


#include "Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {

  /** @short D3Q15 Model. The lattice velocities are: 0. 0.0 0.0 0.0
   *
   * -# scaling, 0.0, 0.0
   * -# -scaling, 0.0, 0.0
   * -# 0.0, scaling, 0.0
   * -# 0.0, -scaling, 0.0
   * -# 0.0, 0.0, scaling
   * -# 0.0, 0.0, -scaling
   * -# scaling, scaling, scaling
   * -# -scaling, -scaling, -scaling
   * -# scaling, scaling, -scaling
   * -# -scaling, -scaling, scaling
   * -# scaling, -scaling, scaling
   * -# -scaling, scaling, -scaling
   * -# scaling, -scaling, -scaling
   * -# -scaling, scaling, scaling
   *
   * The weights are:
   * - 0) 2./9.
   * - 1-6) 1./9.
   * - 7-18) 1./72.
   *
   * Taken from http://arxiv.org/pdf/comp-gas/9611001v1.pdf
   */
  class D3Q15: public Stencil {

  private:

    /** @short function to create the vector of directions
     *  @return the vector of directions
     */
    vector<numeric_vector> makeDirections(double scaling);

    /** @short function to create the vector of weights
     *  @return the vector of weights
     */
    vector<double> makeWeights();

    numeric_matrix makeMomentBasis(vector<numeric_vector> e);

  protected:

    /// speed of sound
    const double m_speedOfSound;

    /// (speed of sound)^2
    const double m_speedOfSoundSquare;

    /// scaling of the stencil
    const double m_scaling;

  public:

    /// D
    static const size_t D;

    /// Q
    static const size_t Q;

    /// constructor
    D3Q15(double scaling = 1.0);

    /// destructor
    virtual ~D3Q15();

    virtual double getSpeedOfSound() const {
      return m_speedOfSound;
    }
    virtual double getSpeedOfSoundSquare() const {
      return m_speedOfSoundSquare;
    }

    virtual size_t getIndexOfOppositeDirection(size_t index) const {
        switch (index){
          case 0:
            return 0;
            break;
          case 1:
            return 2;
            break;
          case 2:
            return 1;
            break;
          case 3:
            return 4;
            break;
          case 4:
            return 3;
            break;
          case 5:
            return 6;
            break;
          case 6:
            return 5;
            break;
          case 7:
            return 8;
            break;
          case 8:
            return 7;
            break;
          case 9:
            return 10;
            break;
          case 10:
            return 9;
            break;
          case 11:
            return 12;
            break;
          case 12:
            return 11;
            break;
          case 13:
            return 14;
            break;
          case 14:
            return 13;
            break;
          default:
            return 0;
        }
      }

      virtual double getMaxParticleVelocityMagnitude() const {
        return m_scaling*sqrt(3);
      }

      virtual double getScaling() const {
        return m_scaling;
      }


  };



}/* namespace natrium */
#endif /* D3Q15MODEL_H_ */
