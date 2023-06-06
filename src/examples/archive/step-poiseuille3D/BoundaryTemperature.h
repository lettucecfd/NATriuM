//
// Created by dominik on 31.10.22.
//

#ifndef NATRIUM_BOUNDARYTEMPERATURE_H
#define NATRIUM_BOUNDARYTEMPERATURE_H


#include "natrium/solver/TurbulenceStats.h"
#include "natrium/dataprocessors/DataProcessor.h"
#include "natrium/solver/CompressibleCFDSolver.h"
#include "natrium/utilities/BasicNames.h"
#include "natrium/utilities/Math.h"
#include "boost/filesystem.hpp"



namespace natrium {


    class BoundaryTemperature: public DataProcessor<3>{
    private:
        CompressibleCFDSolver<3> m_compressibleSolver;

    public:
        BoundaryTemperature(CompressibleCFDSolver<3> & solver);
        void setBoundaryTemperature();
        virtual void apply();
        virtual ~BoundaryTemperature();

    };

} /* namespace natrium */



#endif //NATRIUM_BOUNDARYTEMPERATURE_H
