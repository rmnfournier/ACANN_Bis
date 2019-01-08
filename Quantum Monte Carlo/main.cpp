#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <ctime>       /* time */
#include <fstream>
#include <omp.h>
#include <mpi/mpi.h>

using namespace std;
#define N (8*8)
#define beta 0.25
#define omega 20.
#define m 1.
#define NBLOCS 300 //75 worked, but with high variance
#define bloc_size 20000
#define CONFIGURATIONS NBLOCS*bloc_size
#define warmup 100000
#define NB_UPDATES 1000
#define INITIAL_SPREAD 0.0025
#define POINTS_TO_UPDATE 1
#define XI_RESOLUTION 0.001
#define XI_T_MAX 5
#define J_OMEGA_MAX 75
#define J_OMEGA_RES 0.001

double STEP_SIZE (0.25);
// *************** Physics ***************** //
double V(const vector<double>& x){
    double sum(0);
    #pragma simd
    for(size_t i(0);i<x.size();i++){
        sum+=x[i]*x[i];
    }
    return 0.5*omega*omega*m*sum;
}
double K(const vector<double>& x,const double dtau){
    double sum(0);
    #pragma simd
    for(size_t i(0);i<x.size()-1;i++){
        sum+=(x[i+1]-x[i])*(x[i+1]-x[i]);
    }
    sum+=(x[x.size()-1]-x[0])*(x[x.size()-1]-x[0]);
    return sum*0.5/(m+0.0)/(dtau*dtau);
}
double H(const vector<double>& x,const double dtau){
    return V(x)+K(x,dtau);
}

// ************** Monte Carlo *****************//
double p(const vector<double>& x,const double dtau){
    return exp(-dtau*H(x,dtau));
}

void init (vector<double>& x){
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> distribution(0.0,INITIAL_SPREAD);

    for (auto& el : x){
        el = distribution(generator);
    }
}

int metropolis(vector<double>& x,double dtau){
    int accepted(0);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 generator(rd()); //Standard mersenne_twister_engine seeded with rd()

    std::normal_distribution<double> distribution(0.0,STEP_SIZE);
    std::uniform_real_distribution<double> uniform(0,1);
    std::uniform_int_distribution<int> uniint(0,N-1);
    vector<double>x_new(x);
    for (size_t i(0);i<POINTS_TO_UPDATE;i++){
        //select slice to move
        int slice = uniint(generator);
        // generate a move
        double delta_xi = distribution(generator);

        x_new[slice]+=delta_xi;
    }
    if(p(x_new,dtau)/p(x,dtau)>=uniform(generator)){
        //cout<<"change accepted "<<endl;
        x=std::move(x_new);
        accepted=1;
    }
    return accepted;
}


/* Sampling with a Bath */
vector<double> build_xi(){
    double xi_0(225),a1(1.486*10000),a2(285),alpha1(903),alpha2(75),f(0.2);
    vector<double>xi(XI_T_MAX/XI_RESOLUTION);
    double t(0);
#pragma simd
    for (size_t c(0);c<xi.size();c++){
        xi[c]= xi_0*( exp(-alpha1*pow(f*t,2))*(1+a1*pow(f*t,4)) + a2*pow(f*t,4)*exp(-alpha2*pow(f*t,2)));
        t+=XI_RESOLUTION;
    }
    return xi;
}

vector<double> build_J(){
        vector<double> xi (build_xi());
        vector<double> J ((J_OMEGA_MAX/J_OMEGA_RES));
        double freq(0);
        for (size_t c(0);c<J_OMEGA_MAX/J_OMEGA_RES;c++){
            double t(0),integral(0);
            //Integrate cos(omega t)*xi(t) for each omega
            //first count all intermediate points twice
            for(size_t ii(2);ii<XI_T_MAX/XI_RESOLUTION-1;ii++){
                t+=XI_RESOLUTION;
                integral+=2*cos(freq*t)*xi[ii];
            }

            //Then add the boundary terms once
            integral+=xi[0];
            integral+=cos(XI_T_MAX*freq)*xi[XI_T_MAX/XI_RESOLUTION-1];
            //finally multiply by dt
            integral*=XI_RESOLUTION/2.;
            //save the results
            J[c]=freq*integral;
            //update frequency
            freq+=J_OMEGA_RES;
        }
    return J;
}

vector<double> build_L(){
    vector<double> J(build_J());
    vector<double> L(N);
    //For each component of L
    for (size_t ii(0);ii<N;ii++){
        //Compute the integral defining L (see matlab live script for details)
        // Start by adding twice the inside terms
        double freq(0),sum(0);
        for (size_t w(1);w<J.size()-1;w++){
            freq+=J_OMEGA_RES;
            sum+=2*J[w]*cosh(beta/2.*freq-freq*ii*beta/(N+0.0))/sinh(beta/2.*freq);
        }
        sum+=J[J.size()-1]*cosh(beta/2.*J_OMEGA_MAX-J_OMEGA_MAX*ii*beta/(N+0.0))/sinh(beta/2.*J_OMEGA_MAX);
        freq=1e-6;
        sum+=J[0]*cosh(beta/2.*freq+freq*ii*beta/(N+0.0))/sinh(beta/2.*freq);
        sum*=J_OMEGA_RES/2.;
        L[ii]=sum;
    }
    return L;
}

double I(const vector<double>&x,const vector<double>& L,double dtau){
    double sum(0);
    for (size_t ii(0);ii<x.size();++ii){
        for(size_t jj(0);jj<=ii;jj++){
            if(jj==ii)
                sum+=pow(x[ii]-x[jj],2)*L[ii-jj]*dtau*dtau*0.5;
            else
                sum+=pow(x[ii]-x[jj],2)*L[ii-jj]*dtau*dtau;
        }
    }
    return exp(-sum*0.5);
}




// Main function
int main(int argc, char * argv[]) {
    // Parallelization
    MPI_Init(&argc,&argv);

    int prank,psize;


    MPI_Comm_rank(MPI_COMM_WORLD,&prank);
    MPI_Comm_size(MPI_COMM_WORLD,&psize);
    vector<double> x(N);
    vector<double> L (build_L());

    //For each tau going from 0 to beta
    ofstream my_file,corre;
    my_file.open("J.csv");
    vector<double>j(build_J());
    double axis(0);
    for(size_t i(0);i<j.size();i++){
        my_file<<axis<<","<<j[i]<<endl;
        axis+=J_OMEGA_RES;
    }
    my_file.close();
    my_file.open("xi.csv");
    j=(build_xi());
    axis=0;
    for(size_t i(0);i<j.size();i++){
        my_file<<axis<<","<<j[i]<<endl;
        axis+=XI_RESOLUTION;
    }
    my_file.close();

    for(int N1=(prank);N1<N;N1+=psize){
        double correlations(0);
        double variance(0);
        corre.open("correlations"+to_string(N1)+".csv");
        my_file.open("simulation"+to_string(N1)+".csv");

        init(x);


        double dtau ((beta+0.0)/(N+0.0));

        // warm-up phase
        double ratio=0;
        STEP_SIZE/=1.5;
        do{
            ratio=0;
            STEP_SIZE*=1.25;
            cout<<"Warm up for N1 = "<<N1<<" with step size "<<STEP_SIZE<<endl;
            for (size_t w(0);w<warmup;w++){
                ratio+=metropolis(x,dtau);
            }
        }while(ratio/(warmup+0.0)>0.35);
        cout<<"ratio = "<<ratio/(warmup+0.0)<<endl;
        // sampling phase
        cout<<"Sampling for N1 = "<<N1<<endl;
        // For each bloc
        for(size_t sampling(0);sampling<NBLOCS;sampling++){
            // Initialize the mean and the variance for the current bloc
            double bloc_mean(0);
            double bloc_square(0);
            init(x);
            ratio=0;
            int total=0;
            for(size_t bloc(0);bloc<bloc_size;bloc++) {

                //Update the trajectory
                for (size_t k(0); k < NB_UPDATES; k++){
                    ratio+=metropolis(x, dtau);
                    total++;
                }

                // Compute the mean
                double term(x[0]*x[N1]*I(x,L,dtau));
                bloc_mean+=term;
                // Keep this step in memory for later computation of variance
                bloc_square +=  term*term;
            }
            //ratio/=(NB_UPDATES*bloc_size+0.0);
            // cout<<"Acceptance ratio "<<ratio<<endl;

            // estimate the mean and variance for the current bloc
            bloc_mean/=(bloc_size+0.0);

            corre<<bloc_mean<<",";
            corre.flush();
            correlations+=bloc_mean;
            variance+=1./(bloc_size-1.)* (bloc_square+0.0-bloc_mean*bloc_mean*bloc_size);
            if(prank==0)
                cout<<" End of Block nb "<<sampling<<endl;

        }
        correlations/=(NBLOCS+0.0);
        variance/=(NBLOCS*NBLOCS+0.0);
        //Computing mean and var
        my_file<<beta*N1/(N)<<","<<correlations<<","<<sqrt(variance)<<endl;
        my_file.flush();
        corre.close();
        my_file.close();

    }
    MPI_Finalize();
    return 0;
}