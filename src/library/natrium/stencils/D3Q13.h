/*
 * D3Q13Model.h
 *
 *  Created on: Nov 20, 2015
 *      Author: kraemer
 */

#ifndef D3Q13MODEL_H_
#define D3Q13MODEL_H_


#include "Stencil.h"
#include "../utilities/BasicNames.h"

namespace natrium {

  /** @short D3Q13 Model.
   *
   * This model is not based on a cubic grid, but on a icosahedron.
   *
   * The lattice velocities are: 0. 0.0 0.0 0.0
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
   * - 0) 2./5.
   * - 1-12) 1./20.
   *
   *
   * Taken from Krüger et al.: The Lattice Boltzmann Method, 2017
   */
  class D3Q13: public Stencil {

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
    D3Q13(double scaling = 1.0);

    /// destructor
    virtual ~D3Q13();

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
          default:
            return 0;
        }
      }

      virtual double getMaxParticleVelocityMagnitude() const {
    	  return m_scaling*sqrt((5+sqrt(5.))/2.);
      }

      virtual double getScaling() const {
        return m_scaling;
      }


  };



}/* namespace natrium */
#endif /* D3Q13MODEL_H_ */
