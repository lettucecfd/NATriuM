/*
 * D3Q19Model.h
 *
 *  Created on: Sep 16, 2014
 *      Author: bajat
 */

#ifndef D3Q19MODEL_H_
#define D3Q19MODEL_H_


#include "Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {

  /** @short D3Q19 Model
   */
  class D3Q19: public Stencil {

  private:

    /** @short function to create the vector of directions
     *  @return the vector of directions
     */
    vector<numeric_vector> makeDirections(double scaling);

    /** @short function to create the vector of weights
     *  @return the vector of weights
     */
    vector<double> makeWeights();

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
    D3Q19(double scaling = 1.0);

    /// destructor
    virtual ~D3Q19();

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
            return 3;
            break;
          case 2:
            return 4;
            break;
          case 3:
            return 1;
            break;
          case 4:
            return 2;
            break;
          case 5:
            return 6;
            break;
          case 6:
            return 5;
            break;
          case 7:
            return 9;
            break;
          case 8:
            return 10;
            break;
          case 9:
            return 7;
            break;
          case 10:
            return 8;
            break;
          case 11:
            return 13;
            break;
          case 12:
            return 14;
            break;
          case 13:
            return 11;
            break;
          case 14:
            return 12;
            break;
          case 15:
            return 17;
            break;
          case 16:
            return 18;
            break;
          case 17:
            return 15;
            break;
          case 18:
            return 16;
            break;
          default:
            return 0;
        }
      }

      virtual double getMaxParticleVelocityMagnitude() const {
        return m_scaling*sqrt(2);
      }

      virtual double getScaling() const {
        return m_scaling;
      }


  };



}/* namespace natrium */
#endif /* D3Q19MODEL_H_ */
