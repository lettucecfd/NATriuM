/*
 * D3Q19.cpp

 *
 *  Created on: Sep 16, 2014
 *      Author: bajat
 */

#include "D3Q19.h"
#include <math.h>

// enable + operator for filling vectors
#include "boost/assign/std/vector.hpp"

// enable + operator for filling vectors
using namespace boost::assign;

namespace natrium {

  /////////////////////////////
  // ASSIGN STATIC VARIABLES //
  /////////////////////////////

  // assign D and Q
  // has to be done outside the class, because function calls are not allowed in initialization of statics
  /// D
  const size_t D3Q19::D = 3;
  /// Q
  const size_t D3Q19::Q = 19;

  /// constructor
  D3Q19::D3Q19(double scaling):
      Stencil(3, 19, makeDirections(scaling), makeWeights(), Stencil_D3Q19),
      m_speedOfSound(scaling*pow(3, -0.5)),
      m_speedOfSoundSquare(scaling*scaling/3.),
      m_scaling(scaling){
  } //constructor


  /// destructor
  D3Q19::~D3Q19() {
  } /// destructor


  // make weights
  vector<double> D3Q19::makeWeights()  {
    vector<double> result;
    result += 1./3., 1./18.,1./18.,1./18.,1./18.,1./18.,1./18.,1./36.,1./36.,1./36.,
        1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.;
    return result;
  }/// make weights


  /// make directions
  /* Following the definition by Mohamad LBM book */
  vector<numeric_vector> D3Q19::makeDirections(double scaling) {
    const double directionsArray[][3] = { { 0.0, 0.0 , 0.0 },{scaling , 0.0, 0.0},{0.0 , 0.0, scaling},
        { -scaling, 0.0 , 0.0 },{ 0.0, 0.0 , -scaling },{ 0.0, -scaling , 0.0 },{ 0.0, scaling , 0.0 },
        { scaling, 0.0 , scaling }, { -scaling , 0.0 , scaling }, { -scaling , 0.0 , -scaling }, {scaling , 0.0 , -scaling },
        {scaling, -scaling, 0.0 },{ scaling, scaling , 0.0 },{ -scaling, scaling , 0.0 },{ -scaling, -scaling , 0.0 },
        { 0.0, -scaling , scaling },{ 0.0, scaling , scaling },{ 0.0, scaling , -scaling },{ 0.0, -scaling , -scaling }};
    vector<numeric_vector> result;
    for (size_t i = 0; i < Q; i++) {
      numeric_vector direction(D);
      direction(0) = directionsArray[i][0];
      direction(1) = directionsArray[i][1];
      direction(2) = directionsArray[i][2];
      result += direction;
    }
    return result;
  }/// make directions







}/* namespace natrium */

