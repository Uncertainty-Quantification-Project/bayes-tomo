#include<iostream>
#include<cstdlib>
#include<cmath>
#include <cstdlib>
#include <random>
#include <algorithm>
# include "tomolib.hpp"
# include "fmmlib.hpp"
# include "raylib.hpp"


Parameters::Parameters(){}
/// class constructor with call to a parameter file.
Parameters::Parameters(const std::string& filename)
{
  initialise(filename);
}

/// initialisation of parameters from an ascii parameter file.
void Parameters::initialise(const std::string& filename)
{
  std::ifstream fis;
  std::string line;
  std::vector<std::string> list;
  std::size_t found;
  
  paramfile = filename;

  fis.open(filename.c_str());
  while(!fis.eof())
    {
      std::getline(fis, line);
      found = line.find_first_of("\t\n");
      if (found<std::string::npos)
	{
	  list.push_back(line.substr(0,found));
	}
    }
  
  // at this point we have a vector of strings containing all parameter values
  // now extract and assign to members of PArameters class:
  if (list.size()!=17)
    {
      std::cout << "Invalid parameter file ! Needs to be 17 lines." << "\n";
      std::cout << "The file " << filename << " contains " <<list.size() << " lines.\n";
      std::cout << "Last line contains: "<< list.back() << " \n";
      exit(EXIT_FAILURE);
    }
  else
    {
      folder = list[0];
      geosourcefile = folder + list[1];
      tsurveyfile = folder + list[2];
      vinvertedfile = list[3];
      Nx = std::stoi(list[4]);
      Ny = std::stoi(list[5]);
      Nz = std::stoi(list[6]);
      h = std::stod(list[7]);
      Nit = std::stoi(list[8]); 
      s_surveys = std::stod(list[9]);
      s_l = std::stoi(list[10]);
      d_V = std::stod(list[11]);
      minVel = std::stod(list[12]);
      maxVel = std::stod(list[13]);
      celulas = std::stoi(list[14]);
      invert_parametrization = std::stoi(list[15]);
      suavization = std::stoi(list[16]);
    }
}


