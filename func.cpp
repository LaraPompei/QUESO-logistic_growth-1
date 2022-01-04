#include "func.h"

int my_system(double t, const double* y, double* f, void *info){
    	double* params = (double *)info;
    	double K = params[0];
    	double m = params[1];
    	f[0] = K*y[0]*(1-y[0]/m);
    	return GSL_SUCCESS;
    }
