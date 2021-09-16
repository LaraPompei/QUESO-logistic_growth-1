#include "compute.h"

int main(int argc, char* argv[]) {
    // Initialize QUESO environment

#ifndef QUESO_HAS_MPI
    cerr<<"Iniciando MPI"<<endl;
    MPI_Init(&argc,&argv);

    cerr<<"Criando o environment"<<endl;
    FullEnvironment* env = new FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
    cerr<<"Criando o environment sem MPI"<<endl;
    FullEnvironment* env = new FullEnvironment(argv[1],"",NULL);
#endif
    cerr<<"compute"<<endl;
    //call apllication
    compute(*env);
    cerr<<"acabou o compute"<<endl;

    //Finalize QUESO environment
    cerr<<"Apagando o environment"<<endl;
    delete env;

#ifndef QUESO_HAS_MPI
    cerr<<"Finalizando o MPI"<<endl;
    MPI_Finalize();
#endif
    
    return 0;
} 
