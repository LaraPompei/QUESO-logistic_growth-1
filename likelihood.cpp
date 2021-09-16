#include "likelihood.h"

//Construtor
template<class V, class M>
Likelihood<V,M>::Likelihood(const char* prefix, const VectorSet<V, M>& domain, const double* data_mean, const double* t, const double* stdDevs, const int& poi, const double& s):BaseScalarFunction<V,M>(prefix,domain), m_data_mean(poi), m_t(s), m_stdDevs(poi){
    m_data_mean.assign(data_mean,data_mean+poi); 
    size_t const sizet = s;
    m_t.assign(t,t+sizet);
    m_stdDevs.assign(stdDevs,stdDevs+poi);
    model = new double[poi]{0};
    this->poi = poi;
}

template<classV, class M>
void Likelihood::my_observer(const state_type& x, const double& t){
    model[position][0] = x[0];
}

template<classV, class M>
void Likelihood::my_system(const state_type& x, state_type& dxdt, const double& t){
    dxdt[0] = this->k*x[0]*(1-x[0]/this->m);
}

template<classV, class M>
double Likelihood<V,M>::norm(const double* vector) const{
    double norm = 0;
    for(int i=0; i<poi; i++){
        norm += pow(vector[i],2);
    }
    return sqrt(norm);
}

template<class V, class M>
double Likelihood<V, M>::lnValue(const V& domainVector) const{
    this->k = domainVector[0];
    this->m = domainVector[1];
    //resolver pelo integrate
    double* erro = new double[poi];
    for(int i=0; i<poi; i++){
        position = i;
        integrate_const(rk4(), my_system, x, t0, t1, dt, my_observer);
        m_stdDevs[i] = pow(m_data_mean[i]);
        erro[i] = model[i] - m_data_mean[i];
    }
    double misfitValue = pow(norm(erro)/norm(m_stdDevs),2);
    
    return -1.5*misfitValue;
}


template<class V, class M>
Likelihood<V,M>::~Likelihood(){}

template<class V, class M>
double Likelihood<V,M>::actualValue(const V & domainVector, const V * domainDirection, V * gradVector, M * hessianMatrix, V * hessianEffect) const{
    return exp(this->lnValue(domainVector));
}

 template class Likelihood<GslVector,GslMatrix>;
