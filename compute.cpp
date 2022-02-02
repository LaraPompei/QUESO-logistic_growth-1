#include "compute.h"
#include "likelihood.h"
 
#define tam 100
#define num_variaveis 2
#define K 0
#define M 1
#define number_samples 10000
#define spacing_points 0.5

//Defining the time interval of the model to be used by the function integrate
double t0 = 0.0;
double t1 = 50.0;
double dt = 0.1;


double value_k = 0;
double value_M = 0;
int position = 0;
double** data = new double*[number_samples];

// Defining a shorthand for the type of the mathematical state
typedef std::vector< double > state_type;
typedef runge_kutta4<state_type> rk4; 


// Defining the ODE and the params

// System to be solved: dx/dt = k*x*(1-x/M) ~ Logistic Growth
void my_system_boost(const state_type& x , state_type& dxdt, const double& t){
    dxdt[0] = value_k*x[0]*(1-x[0]/value_M);
}

void my_observer(const state_type& x, const double& t){
    
    if(t == 1*spacing_points){
        data[position][0] = x[0]; 
    }
    else if(t == 3*spacing_points){
        data[position][1] = x[0]; 
    }
    else if(t == 6*spacing_points){
        data[position][2] = x[0]; 
    }
    else if(t == 30*spacing_points){
        data[position][3] = x[0]; 
    }
	
}


void filling_matrix(double* t, double** data, double* values_m, double* values_k){
    mt19937 engine; // uniform random bit engine
    
    // seed the URBG
    random_device dev{};
    engine.seed(dev());

    //setting up a lognormal distribution:
    double mu    = 0.2; //mean value of the lognormal distribution
    double sigma = 0.3; //standard deviation of the lognormal distribution
    lognormal_distribution<double> lognormal_dist(mu, sigma);
    //setting up a normal distribution:
    mu    = 15; //mean value of the normal distribution
    sigma = 1.0; //standard deviation of the normal distribution
    normal_distribution<double> normal_dist(mu,sigma);
    cout<< "Filling Matrix"<<endl;
    //filling the matrix
    cout<<"k\tm\tdata[i][0]\tdata[i][1]\tdata[i][2]\tdata[i][3]"<<endl;
    for (int i =0; i<number_samples; i++) {
	state_type x = {1.0};
        values_k[i] = lognormal_dist(engine);
        values_m[i] = normal_dist(engine);
	value_k = values_k[i];
	value_M = values_m[i];
	position = i;
        //chamar o odeint para preencher a matriz
	boost::numeric::odeint::integrate_const(rk4(), my_system_boost, x, t0, t1, dt, my_observer);
	cout  << value_k << "\t" << value_M << "\t" << data[i][0] << "\t" << data[i][1] << "\t" << data[i][2] << "\t" << data[i][3] << "\n";
    }
}

void compute(const FullEnvironment& env){
    struct timeval timevalNow;

    gettimeofday(&timevalNow, NULL);
    if (env.fullRank() == 0){
        cerr<<"\nBeginning run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "<<ctime(&timevalNow.tv_sec)<<endl;
    }
    env.fullComm().Barrier();
    env.subComm().Barrier();

    //instantiating the parameter space
    cerr<<"instantiating the parameter space.."<<endl;
    VectorSpace<> paramSpace(env, "param_", 2, NULL);

    //defining the domain of k e m 
    cerr<<"instantiating the parameter domain"<<endl;
    GslVector paramMinValues(paramSpace.zeroVector());
    GslVector paramMaxValues(paramSpace.zeroVector());
    paramMinValues[K] = 0.0;
    paramMaxValues[K] = 5.0;
    paramMinValues[M] = 10.0;
    paramMaxValues[M] = 20.0;
    BoxSubset<> paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

    //instantiating the likelihood function object
    int num_poi = 4;
    cerr<< "Instantiating the likelihood function object and generating the samples"<<endl;
    double *t             = new double[tam];
    double* values_m      = new double[number_samples];
    double* values_k      = new double[number_samples];
    double* data_mean     = new double[num_poi];
    double* data_std      = new double[num_poi];
    double* poi           = new double[num_poi]{1,3,6,30}; //defining the points in the time vector that will be used to fit the parameters (k and m) of the model
    double Y0 = 1.0;
    int dim = 1;
    //allocating memory for the data matrix[number_samples,tam]
    cout<<"allocating memory for the data matrix"<<endl;
    for(int i=0; i<number_samples; i++){
        data[i] = new double[tam];
    }

    //vector t is filled with [0.0,50] interval catching each number after 0.5
    cout<<"Filling vector t"<<endl<<"t[ ";
    for (int i = 0; i<tam; i++){
        t[i] = i*spacing_points;
        cout<<t[i]<<" ";
    }
    cout<<"]"<<endl;
    
    //generating and filling the matrix
    cout<<"Generating and filling the matrix"<<endl;
    filling_matrix(t,data,values_m,values_k);

    //instantiating the data_mean and data_std arrays
    for (int i=0; i<num_poi; i++){
        data_mean[i] = 0.0;
        data_std[i]  = 0.0;
    }

    //data mean
    cout<<"Calculating the mean of the data"<<endl;
    for(int i = 0; i < num_poi;i++){
        for(int j=0; j<number_samples; j++){
            data_mean[i] += data[j][i]/number_samples;
        }
    }
    //standard deviation
    cout<<"Calculating the standard deviation"<<endl;
    for(int i = 0; i < num_poi;i++){
        for(int j=0; j<number_samples; j++){
            data_std[i] += pow((data[j][i] - data_mean[i]),2);
        }
        data_std[i] = sqrt(data_std[i]/number_samples);
        cout<<"Para t = "<<i<<endl<<"Data mean: "<<data_mean[i]<<endl<<"Data standard deviation: "<<data_std[i]<<endl;
    }

    cerr<<"Creating the likelihood object"<<endl;   
    Likelihood<> lhood("like_", paramDomain, num_poi, t0, t1, dt, poi, data_mean, data_std, Y0, dim); 
     
    //defining the prior RV
    cerr<<"Defining the prior RV"<<endl;
    UniformVectorRV<> priorRV("prior_",paramDomain);

    //instantiating the inverse problem
    cerr<<"Instantiating the inverse problem"<<endl;
    GenericVectorRV<> postRv("post_", paramSpace);
    StatisticalInverseProblem<> ip("", NULL, priorRV, lhood, postRv);

    //Solving the inverse problem
    cerr<<"Solving the inverse problem"<<endl;
    GslVector paramInitials(paramSpace.zeroVector());
    priorRV.realizer().realization(paramInitials);

    GslMatrix proposalCovMatrix(paramSpace.zeroVector());

    proposalCovMatrix(0,0) = pow(abs(paramInitials[0])/20.0, 2.0);
    proposalCovMatrix(1,1) = pow(abs(paramInitials[1])/20.0, 2.0);

    cerr<<"Inicio do mcmc"<<endl;
    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
    cerr<<"Fim do mcmc"<<endl;

    if (env.fullRank() == 0) {
        cerr << "Ending run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
    }
}


// ************************* COMENTANDO PARA TESTAR O CÓDIGO BASE ********************

/*

void save_data(double* model, double** baseModel, double* data, double* values_of_a){
    fstream a_data, a_mcmc, model_data, model_mcmc;
    char fileName[50];

    //writing original a on file
    sprintf(fileName, "a_data.m");
    a_data.open(fileName, ios_base::out);
    a_data<<"a_data = [";
    for(int i=0 ; i<number_samples; i++){
        a_data<<data[i]<<endl;
    }
    a_data<<"];"<<endl;
    a_data.close();

    //writing original model on file
    sprintf(fileName, "model_data.m");
    model_data.open(fileName, ios_base::out);
    model_data<<"model_data = [";
    for(int i=0 ; i<number_samples; i++){
        model_data<<baseModel[i][poi]<<endl;
    }
    model_data<<"];"<<endl;
    model_data.close();

    //writing mcmc a on file

    //writing mcmc model on file
}


*/


