//safeguard library
#ifndef TOMOLIB_INCLUDED
#define TOMOLIB_INCLUDED

# include<armadillo>
# include<vector>
# include<string>

# include "fmmiolib.hpp"
# include "fmmlib.hpp"

/// a global constant to use in place of -inf values in ln(v) (i.e., a very large negative number but not infinite...)
const double V_SHADOW = -100.0; // remember that the wavespeed will then be exp(-100) = 3.72e-44 ... close enough from zero !

/// the class Parameters is a structure containing inversion parameters, file names, etc.
class Parameters
{
public:
  std::string folder; ///< name fo folder to store all input/saved files.
  std::string paramfile; ///< name of parameter file (ascii)
  std::string geosourcefile; ///< name of event position file (ascii)
  std::string tsurveyfile; ///< name of survey arrival times file (.dat)
  std::string vinvertedfile; ///< name of inverted V structure file (.dat)
  int Nx; // dim x
  int Ny; // dim y
  int Nz; // dim z
  double h; // grid spacing
  int Nit; // Interations number of markov chain monte carlo
  double s_surveys; ///< variance on arrival times from surveys
  int s_l; ///< variance on model lacate voronoi cel
  double d_V; //atualization of model velocidade
  double minVel; // min velocidade a prior
  double maxVel; // max velocidade a prior
  int celulas; //numero iniciais de celulas na parametrização
  int invert_parametrization; // 1 inverte parametrização, 0 mantém fixa a parametrização
  int suavization; // numeros de pontos para suavização do modelos V e E

  Parameters();
  Parameters(const std::string& );
  void initialise(const std::string& );
};

#endif
