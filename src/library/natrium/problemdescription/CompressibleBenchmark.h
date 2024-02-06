//
// Created by dwilde3m on 24.02.20.
//

#ifndef NATRIUM_COMPRESSIBLEBENCHMARK_H
#define NATRIUM_COMPRESSIBLEBENCHMARK_H



#include "deal.II/grid/tria.h"
#include "deal.II/base/function.h"

#include "BoundaryCollection.h"
#include "Benchmark.h"

#include "../problemdescription/ProblemDescription.h"

#include "../utilities/BasicNames.h"

namespace natrium{

template<size_t dim> class CompressibleBenchmark: public Benchmark<dim> {
private:
    /// function to define initial densities
    boost::shared_ptr<dealii::Function<dim> > m_analyticT;

protected:
    void setAnalyticT(boost::shared_ptr<dealii::Function<dim> > ana_T) {
        m_analyticT = ana_T;
    }

public:

    CompressibleBenchmark(boost::shared_ptr<Mesh < dim>

    > triangulation,
    double viscosity,
    double characteristicLength
    ) :
    Benchmark<dim>(triangulation, viscosity, characteristicLength
    ) {

        m_analyticT = boost::make_shared<DealIIExtensions::Functions::ConstantFunction<dim> >(1.0, 1);

    }

    const boost::shared_ptr<dealii::Function<dim> > &getAnalyticTFunction(double time) const {
        m_analyticT->set_time(time);
        return m_analyticT;
    }

    virtual const boost::shared_ptr<dealii::Function<dim> > &getInitialTFunction() const {
        m_analyticT->set_time(0);
        return m_analyticT;
    }

};
} // namespace natrium

#endif //NATRIUM_COMPRESSIBLEBENCHMARK_H