#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/Environment.h>
#include <queso/Defines.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <cmath>

#include <boost/numeric/odeint.hpp>

#include "compute.h"
#include "func.h"

using namespace std;
using namespace QUESO;
using namespace boost::numeric::odeint;

template<class V = GslVector, class M = GslMatrix>
class Likelihood : public BaseScalarFunction<V,M>{
    private:
        vector<double> m_data_mean; //exponencial model results
        vector<double> m_t; //values of t 
        vector<double> m_stdDevs; //account for uncertainties
        vector<double> m_poi;
       	const QUESO::BaseEnvironment* m_env;
  	
	int numPoi;
	double Y0;
	int dim;

	//Defining the time interval of the model to be used by the function integrate
	double t0;
	double t1;
	double dt;


    public:
	Likelihood(const char* prefix, const VectorSet<V, M>& domain, int& poi, double& t0, double& t1, double& dt, double* points, double* data_mean, double* std_data, double& Y0, int& dim);
	//Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* data_mean, const double* t, const double* stdDevs, const int& poi, const double& s);
        //void my_observer(const state_type& x, state_type& dxdt, const double& t);
        
	//int my_system(double t, const double* y, double* f, void *info);
	
	//void my_system(const state_type& x , state_type& dxdt, const double& t);
        virtual ~Likelihood();
        double norm(const double* vector) const;
        virtual double lnValue(const V& domainVector) const;
	virtual double actualValue(const V&  paramValues, const V* paramDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
};

#endif
