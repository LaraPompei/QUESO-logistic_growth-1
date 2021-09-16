#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <cmath>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace QUESO;
using namespace boost::numeric::odeint;

// Defining a shorthand for the type of the mathematical state
typedef vector<double> state_type;
typedef runge_kutta4<state_type> rk4;

template<class V = GslVector, class M = GslMatrix>
class Likelihood : public BaseScalarFunction<V,M>{
    private:
        vector<double> m_data_mean; //exponencial model results
        vector<double> m_t; //values of t 
        vector<double> m_stdDevs; //account for uncertainties
        double* model;
        double k;
        double m;
        int poi;
        int position;
        const QUESO::BaseEnvironment* m_env;
   
   public:
        Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* data_mean, const double* t, const double* stdDevs, const int& poi, const double& s);
        void my_observer(const state_type& x, const double& t);
        void my_system(const state_type& x , state_type& dxdt, const double& t);
        virtual ~Likelihood();
        double norm(const double* vector) const;
        virtual double lnValue(const V& domainVector) const;
        virtual double actualValue(const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
};

#endif
