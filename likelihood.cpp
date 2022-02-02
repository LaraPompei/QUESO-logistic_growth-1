#include "likelihood.h"
//Construtor
template<class V, class M>
Likelihood<V,M>::Likelihood(const char* prefix, const VectorSet<V, M>& domain, int& poi, double& t0, double& t1, double& dt, double* points, double* data_mean, double* std_data, double& Y0, int& dim):BaseScalarFunction<V,M>(prefix,domain), m_data_mean(poi), m_stdDevs(poi){
    for(int i = 0; i<poi; i++){
    	cout<<"Desvio_construtor recebido: "<<std_data[i]<<endl;
	cout<<"Media construtor recebida: "<<data_mean[i]<<endl;
	cout<<"poi construtor recebido: "<<points[i]<<endl;
	this->m_data_mean.insert(m_data_mean.begin()+i,data_mean[i]); 
	this->m_stdDevs.insert(m_stdDevs.begin()+i,std_data[i]); 
	this->m_poi.insert(m_poi.begin()+i,points[i]); 
	cout<<"Desvio_constutor armazenado: "<<m_stdDevs.at(i)<<endl;
	cout<<"poi construtor armazenado: "<<m_poi.at(i)<<endl;
    	cout<<"media construtor armazenado: "<<m_data_mean.at(i)<<endl;
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
    cout<<"lnValue"<<endl;
    double k = domainVector[0];
    double m = domainVector[1];
    cout<<"Parametro likelihood: "<<endl;
    cout<<"k\t"<<k<<"\tm\t"<<m<<endl;    
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
	cout<<"gsl_status_integrate: "<<status<<endl;
	cout<<"Y0: "<<Y0[0]<<endl;
	std_Devs[i] = m_stdDevs.at(i);
	cout<<"Desvio: "<<std_Devs[i]<<endl;
	cout<<"Desvio main: "<<m_stdDevs.at(i)<<endl;	
	//calculando o erro
        erro[i] = Y0[0] - m_data_mean[i];
	cout<<"erro "<<i<<" "<<erro[i]<<endl;
    }

    //norma do erro
    double misfitValue=0;
    misfitValue = norm(erro)/norm(std_Devs);
    misfitValue = pow(misfitValue,2);

    cout<<"gsl_odeiv_evolve_free"<<endl;
    gsl_odeiv_evolve_free (e);
    cout<<"gsl_odeiv_control_free"<<endl;
    gsl_odeiv_control_free(c);
    cout<<"gsl_odeiv_step_free"<<endl;
    gsl_odeiv_step_free   (s);
    return -.5*misfitValue;
}


template<class V, class M>
Likelihood<V,M>::~Likelihood(){}

template<class V, class M>
double Likelihood<V,M>::actualValue(const V& paramValues, const V* paramDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const{
   cout<<"Erro no actualValue?"<<endl;
   return std::exp(this->lnValue(paramValues));
}

template class Likelihood<GslVector,GslMatrix>;
