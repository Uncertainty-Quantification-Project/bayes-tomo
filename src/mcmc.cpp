#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include <string.h>
#include<chrono>
#include <armadillo>
#include "fmmiolib.hpp"
#include "fmmlib.hpp"
#include "tomolib.hpp"

using namespace std;

int main(int nargin, char *varargin[])
{
  clock_t tim = clock();
  Parameters param;
  string paramfilename, VfileSource, EfileSource, temp, energy, TfileSource, filevoronoi;
  vector<locations::Location> sensors;
  vector<Data> Data_obs;
  Data Tobs;
  int it_aceite;

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
 
  if (nargin!=2){
    cout << "Wrong number of input arguments ! Needs a parameter file as input.\n";
    exit(1);
  }
  
  cout << "Running " << varargin[0] << "...\n";
  cout << "Using " << varargin[1] << " as input parameter file.\n";

  paramfilename.assign(varargin[1]);
  
  param.initialise(paramfilename);

  // read source location
  sensors = locations::readlocations(param.geosourcefile.c_str());

  for (int k=0; k<sensors.size(); k++){
    temp = std::to_string(k);
    TfileSource = param.tsurveyfile + "__source_" + temp + ".bin";
    Tobs.load(TfileSource.c_str(),param.Nx,param.Ny,param.Nz,param.h);
    Data_obs.push_back(Tobs);
    Tobs.dat.erase(Tobs.dat.begin(),Tobs.dat.end());
  }

  Model modelo(param);
  State state(param);
  state.InitRandomState(param, generator);

  modelo.CreateFildVelocity(state.GetCelulas());
  modelo.CreateFildAnisotropy(state.GetCelulas());
  modelo.SuavizaModel();

  Grid grid(modelo.GetHeight(),modelo.GetWidth(),modelo.GetDeep(),modelo.GetH());  
  grid.initialise(modelo.GetVelocity(),modelo.GetAnisotropy());

  state.Energy(modelo, sensors, grid, Data_obs);  

  Chain chain(state, Data_obs);
  it_aceite = 0;
  for (int n=0; n < param.Nit; n++) { 
    std::cout << "Interation:  "<<n<<std::endl;

    chain.RjMCMC(modelo, sensors, grid, n, generator);
  
    double ratio = ((double)(chain.GetNumAt() - 1)/(double)n);
    modelo.CreateFildVelocity(chain.GetState().GetCelulas());
    modelo.CreateFildAnisotropy(chain.GetState().GetCelulas());
    modelo.SuavizaModel();
    if(chain.GetSaveAceite() == 1){
      it_aceite += 1; 
      temp = to_string(it_aceite);
      VfileSource = param.vinvertedfile + "_iteration_" + temp + ".bin";
      modelo.SaveVelocity(VfileSource.c_str());     
    }
    std::cout << "Taxa de aceite:  "<<ratio<<std::endl;
    chain.DelBeginState();
  }

  cout <<"Energia Saving..."<<endl;
  energy = "Chain_Energy_" + temp + ".bin";
  chain.SaveEnergy(energy.c_str());

  tim = clock() - tim;

  cout << "time elapsed: " << double(tim)/CLOCKS_PER_SEC << " sec"<<endl;

  return 0;
}
