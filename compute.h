#ifndef COMPUTE_H
#define COMPUTE_H

#include <queso/Environment.h>
#include <cmath>
#include <queso/GslMatrix.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GaussianVectorRV.h>
#include <queso/LogNormalVectorRV.h>
#include <queso/GenericVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <queso/VectorSet.h>

#include <boost/numeric/odeint.hpp>

#include <sys/time.h>
#include <random>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cmath>

using namespace QUESO;
using namespace std;
using namespace boost::numeric::odeint;

// Defining a shorthand for the type of the mathematical state
typedef vector<double> state_type;
typedef runge_kutta4<state_type> rk4;

void my_system_boost(const state_type& x , state_type& dxdt, const double& t);
void my_observer(const state_type& x, const double& t);
void compute(const FullEnvironment& env);
void filling_matrix(double* t, double** data, double* values_m, double* values_k);
void save_data(double* model, double** baseModel, double* data, double* values_of_a);
#endif
