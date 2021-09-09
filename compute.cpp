#include "compute.h"
#include "likelihood.h"

#define tam 100
#define num_variaveis 2
#define num_poi 4
#define K 0
#define M 1
#define Y0 1.0
#define number_samples 10000

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

    //filling the matrix
    for (int i =0; i<number_samples; i++) {
        values_k[i] = lognormal_dist(engine);
        values_m[i] = normal_dist(engine);
        //chamar o odeint para preencher a matriz
    }
}

void compute(const FullEnvironment& env){
    struct timeval timevalNow;

    gettimeofday(&timevalNow, NULL);
    if (env.fullRank() == 0){
        cout<<"\nBeginning run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "<<ctime(&timevalNow.tv_sec)<<endl;
    }
    env.fullComm().Barrier();
    env.subComm().Barrier();

    //instantiating the parameter space
    cerr<<"instantiating the parameter space.."<<endl;
    VectorSpace<> paramSpace(env, "param_", 2, NULL);

    //defining the domain of a
    cerr<<"instantiating the parameter domain"<<endl;
    GslVector paramMinValues(paramSpace.zeroVector());
    GslVector paramMaxValues(paramSpace.zeroVector());
    paramMinValues[K] = 0.0;
    paramMaxValues[K] = 5.0;
    paramMinValues[M] = 10.0;
    paramMaxValues[M] = 20.0;
    BoxSubset<> paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

    //instantiating the likelihood function object
    cerr<< "Instantiating the likelihood function object and generating the samples"<<endl;
    double spacing_points = 0.5;
    double *t             = new double[tam];
    double* values_m      = new double[number_samples];
    double* values_k      = new double[number_samples];
    double** data         = new double*[number_samples];
    double* data_mean     = new double[num_poi];
    double* data_std      = new double[num_poi];
    double* poi           = new double[num_poi]{1,3,6,30}; //defining the points in the time vector that will be used to fit the parameters (k and m) of the model

    //allocating memory for the data matrix[number_samples,tam]
    cerr<<"allocating memory for the data matriz"<<endl;
    for(int i=0; i<number_samples; i++){
        data[i] = new double[tam];
    }

    //vector t is filled with [0.0,50] interval catching each number after 0.5
    cerr<<"Filling vector t"<<endl<<"t[ ";
    for (int i = 0; i<tam; i++){
        t[i] = i*spacing_points;
        cerr<<t[i]<<" ";
    }
    cerr<<"]"<<endl;
    
    //generating and filling the matrix
    cerr<<"Generating and filling the matrix"<<endl;
    filling_matrix(t,data,values_m,values_k);

    //instantiating the data_mean and data_std arrays
    for (int i=0; i<num_poi; i++){
        data_mean[i] = 0.0;
        data_std[i]  = 0.0;
    }

    //data mean
    cerr<<"Calculating the mean of the data"<<endl;
    for(int i = 0; i < num_poi;i++){
        for(int j=0; j<number_samples; j++){
            data_mean[i] += data[j][i];
        }
    }
    //standard deviation
    cerr<<"Calculating the standard deviation"<<endl;
    for(int i = 0; i < num_poi;i++){
        for(int j=0; j<number_samples; j++){
            data_std[i] += pow((data[j][i] - data_mean[i]),2);
        }
        data_std[i] = sqrt(data_std[i]/number_samples);
        cerr<<"Para t = "<<i<<endl<<"Data mean: "<<data_mean[i]<<endl<<"Data standard deviation: "<<data_std[i]<<endl;
    }

    cerr<<"Creating the likelihood object"<<endl;
    /*************************************************************************************************************************************
     *                                              MEXER NA LIKELIHOO                                                                   *
     *************************************************************************************************************************************/
    Likelihood<> lhood("like_", paramDomain, data_mean, t, data_std, poi, tam);
    
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

    ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

    if (env.fullRank() == 0) {
        cout << "Ending run of 'Example 1: Log-Normal Distribution Function (0.5,.2)' example at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
    }
}

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
