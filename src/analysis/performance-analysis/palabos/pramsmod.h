
/*
 * prams.h
 *
 *  Created on: 14.05.2014
 *      Author: dominik
 */

#ifndef PRAMSMOD_H_
#define PRAMSMOD_H_

#include "core/globalDefs.h"
#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"
#include <string>
#include <fstream>

namespace plb {
template<typename T>
class ModIncomprFlowParam {
public:
    /// Constructor
    /** \param latticeU_  Characteristic velocity in lattice units (proportional to Mach number).
     *  \param Re_ Reynolds number.
     *  \param N_  Resolution (a lattice of size 1 has N_+1 cells).
     *  \param lx_ x-length in dimensionless units (e.g. 1).
     *  \param ly_ y-length in dimensionless units (e.g. 1).
     *  \param lz_ z-length in dimensionless units (e.g. 1).
     */


    ModIncomprFlowParam(T latticeU_, T Re_, plint resolution_, T lx_, T ly_, T lz_=T())
    : latticeU(latticeU_), Re(Re_), resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    {
        physicalU      = (T)1;
        physicalLength = (T)1;
    }
    /// velocity in lattice units (proportional to Mach number)
    T getLatticeU() const { return latticeU; }
    /// velocity in physical units
    T getPhysicalU() const { return physicalU; }
    /// Reynolds number
    T getRe() const      { return Re; }
    /// physical resolution
    T getPhysicalLength() const { return physicalLength; }
    /// resolution
    plint getResolution() const { return resolution; }
    /// x-length in dimensionless units
    T getLx() const      { return getPhysicalLength()*lx; }
    /// y-length in dimensionless units
    T getLy() const      { return getPhysicalLength()*ly; }
    /// z-length in dimensionless units
    T getLz() const      { return getPhysicalLength()*lz; }
    /// lattice spacing in physical units
    T getDeltaX() const  { return (T)lx/(T)getResolution(); }
    /// time step in physical units
    T getDeltaT() const  { return getDeltaX()*getLatticeU()/getPhysicalU(); }
    /// conversion from dimensionless to lattice units for space coordinate
    T getLatticeNu() const { return getLatticeU()*(T)getResolution()/Re; }
    /// viscosity in dimensionless units
    T getPhysicalNu() const { return lx*getPhysicalU()/Re; }
    /// relaxation time
    T getTau() const       { return (T)3*getLatticeNu()+(T)0.5; }
    /// relaxation frequency
    T getOmega() const     { return (T)1 / getTau(); }


public:
    T physicalU, latticeU, physicalLength, Re;
    plint resolution;
    T lx, ly, lz;
};

template<typename T>
void ModwriteLogFile(ModIncomprFlowParam<T> const& parameters, plb_ofstream &parameters_log)
{

	parameters_log << "Mach Number                Ma=" << parameters.latticeU * sqrt(3) << "\n";
    parameters_log << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
    parameters_log << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    parameters_log << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
    parameters_log << "Relaxation frequency:      omega=" << parameters.getOmega() << "\n";
    parameters_log << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    parameters_log << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    parameters_log << "Extent of the system:      lz=" << parameters.getLz() << "\n";
    parameters_log << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    parameters_log << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
    parameters_log << "Physical Nu:               PhyNu=" << parameters.getPhysicalNu() << "\n";
    parameters_log << "Lattice Nu:                LatNu=" << parameters.getLatticeNu() << "\n\n\n";


}

}

#endif /* PRAMS_H_ */
