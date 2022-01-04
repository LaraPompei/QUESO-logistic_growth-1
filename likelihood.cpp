#include "likelihood.h"
//Construtor
template<class V, class M>
Likelihood<V,M>::Likelihood(const char* prefix, const VectorSet<V, M>& domain, int& poi, double& t0, double& t1, double& dt, double* points, double* data_mean, double* std_data, double& Y0, int& dim):BaseScalarFunction<V,M>(prefix,domain), m_data_mean(poi), m_stdDevs(poi){
    for(int i = 0; i<poi; i++){
    	this->m_data_mean.push_back(data_mean[i]); 
	this->m_stdDevs.push_back(std_data[i]); 
	this->m_poi.push_back(points[i]); 
    }
    this->numPoi = poi; 
    this->t0 = t0; 
    this->t1 = t1;
    this->dt = dt;
    this->Y0 = Y0; 
    this->dim = dim; 
}

/*####################################################################################
 *    		Function used to calculate the norm of a vector			     *
 *####################################################################################*/
template<class V, class M>
double Likelihood<V,M>::norm(const double* vector) const{
    double norm = 0;
    for(int i=0; i<this->numPoi; i++){
        norm += pow(vector[i],2);
    }
    return sqrt(norm);
}

/************************************************************************************
 * 		Function used to define the system of equations			    *
 ************************************************************************************
 */
/*
template<class V, class M>
int Likelihood<V,M>::my_system(double , const double* y, double* f, void *info){
    double* params = (double *)info;
    double K = params[0];
    double m = params[1];
    f[0] = K*y[0]*(1-y[0]/m);
    return GSL_SUCCESS;
}
*/
/************************************************************************************
 * 		Function that is called by the compute				    *
 ************************************************************************************/
template<class V, class M>
double Likelihood<V, M>::lnValue(const V& domainVector) const{
    //definindo parametro para calculos
    double k = domainVector[0];
    double m = domainVector[1];
    
    double params[]={k, m};

    //definindo parametros para a integral
    const gsl_odeiv_step_type *T   = gsl_odeiv_step_rkf45; //rkf45; //gear1;
    gsl_odeiv_step *s   	   = gsl_odeiv_step_alloc(T,dim);
    gsl_odeiv_control *c           = gsl_odeiv_control_y_new(1e-6,0.0); //(erro_abs, erro_rel)
    gsl_odeiv_evolve *e            = gsl_odeiv_evolve_alloc(dim);
    gsl_odeiv_system sys           = {my_system, NULL, 1, (void *)params};
    
    //criando um vetor para armazenar o erro em cada ponto
    double* erro = new double[this->numPoi]();
    
    double t = this->t0;
    double t1 = this->t1;
    double dt = this->dt;
    double* Y0 = new double[1];
    double *std_Devs = new double[numPoi];

    //calcular o erro em cada ponto
    for(int i=0; i<numPoi; i++){
	Y0[0] = this->Y0;
        int status = gsl_odeiv_evolve_apply(e, c, s, &sys, &t, t1, &dt, Y0);
	std_Devs[i] = m_stdDevs[i];
	
	//calculando o erro
        erro[i] = Y0[0] - m_data_mean[i];
    }

    //norma do erro
    double misfitValue = pow(norm(erro)/norm(std_Devs),2);
    
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free   (s);
    return -1.5*misfitValue;
}


template<class V, class M>
Likelihood<V,M>::~Likelihood(){}

template<class V, class M>
double Likelihood<V,M>::actualValue(const V& paramValues, const V* paramDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const{
   return std::exp(this->lnValue(paramValues));
}

template class Likelihood<GslVector,GslMatrix>;
