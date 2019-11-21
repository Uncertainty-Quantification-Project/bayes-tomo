// safeguard library
#ifndef FMMLIB_INCLUDED
#define FMMLIB_INCLUDED
//

# include <array>
# include <vector>
# include <random>
# include <boost/heap/fibonacci_heap.hpp>
# include "fmmiolib.hpp"

//forward declare class Data:
class Data;

class Parameters;

/// The class Node corresponds to a node in a three dimensional orthonormal grid.
class Node
{
    
public:
  int index; ///< the index of the node in the grid, a unique identifier. It corresponds to the linear index in a 3D matrix, similar to Matlab linear indices.
  int i; ///< the subscript index along the x direction.
  int j; ///< the subscript index along the y direction.
  int k; ///< the subscript index along the z direction.

  /// an array of the 6 indices corresponding to the 6 neighbours of the node in the grid.
  /*! If the index is -1, it means that the neighbour does not exist (point outside the grid).
   * nhb[0] corresponds to the neighbour in +x
   * nhb[1] corresponds to the neighbour in -x
   * nhb[2] corresponds to the neighbour in +y
   * nhb[3] corresponds to the neighbour in -y
   * nhb[4] corresponds to the neighbour in +z
   * nhb[5] corresponds to the neighbour in -z
   */
  int nhb[6]; 
  double T; ///< the arrival time of the wavefront at the node.
  double V; ///< the horizontal wavespeed at this node.
  double E; ///< the anisotropy parameter (ratio of vertical to horizontal wavespeed).
  bool trial; ///< boolean tag for "trial" nodes.
  bool known; ///< boolean tag for "known" nodes.
  bool unknown; ///< boolean tag for "unknown nodes.

  Node();

  /// define comparison operator > for class "node" (compares the arrival time T)
  bool operator > (const Node& pt) const
  {
    return ( T > pt.T);
  }

  /// define comparison operator < for class "node" (compares the arrival time T)
  bool operator < (const Node& pt) const
  {
    return ( T < pt.T);
  }

};


/// The class Tetra embeds the information (velocity, anisotropy, arrival time of neighbours, etc) and methods required to compute the arrival time at a given node. It is named "Tetra" in reference to the tetrahedron that is formed by a given node of interest and three of its neighbours (each in x, y and z direction, resp.).
class Tetra
{
    
public:
  double vel; ///< the horizontal wavespeed at the node of interest (i,j,k).
  double eps; ///< the aniotropy parameter at the node of interest (i,j,k).
  double h; ///< grid spacing.
  bool test_a; ///< true is the x-neighbour (i+/-1,j,k) exists and is known.
  bool test_b; ///< true if the y-neighbour (i,j+/-1,k) exists and is known.
  bool test_c; ///< true if the z-neighbour (i,j,k+/-1) exists and is known. 
  double a; ///< arrival time at the x-neighbour.
  double ap; ///< arrival time at the x-neighbour of the x-neighbour (use for second order FD).
  double sw_a; ///< switch: equal to 1 if the x-neighbour of the x-neighbour exists and is knwon, 0 other wise. 
  double b; ///< same as for a, corresponding to the y-neighbour.
  double bp;  ///< same as for ap, corresponding to the y-neighbour.
  double sw_b; ///< same as for sw_a, corresponding to the y-neighbour.
  double c;  ///< same as for a, corresponding to the z-neighbour.
  double cp;  ///< same as for ap, corresponding to the z-neighbour.
  double sw_c;  ///< same as for sw_a, corresponding to the z-neighbour.
    
  Tetra(const double&, const double&, const double&);
  Tetra();
  double compute_T();
    
private:

  double threeD_aniso();
  double twoD_aniso_a();
  double twoD_aniso_b();

};


//type declaration:
typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<std::greater<Node> > >::handle_type handle_t;


/// The class Grid corresponds to the set of nodes at which arrival time is to be computed, and contains the methods to do that.
class Grid
{

public:
  int Nx; ///< the number of nodes in the x direction
  int Ny; ///< the number of nodes in the y direction
  int Nz; ///< the number of nodes in the z direction
  int Ntot; ///< the total number of nodes
  double h; ///< the node spacing
  std::vector<Node> node; ///< the vector of nodes
  std::vector<handle_t> tab_handle_trial; ///< vector of handles in the trial heap
  int ind_source; ///< linear index of source node
  

  Grid(const int&, const int&, const int&, const double&);
  ~Grid();

  void initialise(const std::vector<double>& , const std::vector<double>& );
  void initialise_log(const std::vector<double>& , const std::vector<double>& );
  void reset_T();
  void reset(const std::vector<double>& , const std::vector<double>&);
  int march(const int& , const bool& );
  std::vector<double> export_Tvector();
  void import_T(const Data& );
  int sub2ind(const int& ,const int& , const int& );
  void ind2sub(const int& , int& , int& , int& );    
  
private:

  boost::heap::fibonacci_heap<Node, boost::heap::compare<std::greater<Node> > > trial_heap; ///< the heap structure in which "trial" nodes are pushed.

  Tetra buildtetra(const Node& , const int& ,const int&, const int&);
  std::array<Tetra,4> findalltetra(const int& ,const int& );
  double update(const int& , const int& );
  void init_box();
    
};


class Point{
  public:
    Point();
    Point(int&, int&, int&);
    void SetX(int&);
    void SetY(int&);
    void SetZ(int&);
    int GetX();
    int GetY();
    int GetZ();
    int DistanceSqrd(int&, int&, int&);
  private:
  int x_, y_, z_;
};

class Celula{
  public:
    Celula();
    Celula(Point&, double&, double&, double&, double&);
    void SetCoord(Point&);
    void SetVel(double&);
    void SetAnis(double&);
    void SetminVel(double&);
    void SetmaxVel(double&);
    void SetminAnis(double&);
    void SetmaxAnis(double&);
    Point GetCoord();
    double GetVel();
    double GetAnis();
    double GetMaxVel();
    double GetMinVel();
    double GetMaxAnis();
    double GetMinAnis();

    bool operator == (Celula& cl){
      bool eq;
      if(coord_.GetX() == cl.GetCoord().GetX() && coord_.GetY() == cl.GetCoord().GetY() && coord_.GetZ() == cl.GetCoord().GetZ())
        eq = true;
      else
        eq = false;
      return eq;
    }
  private:
  Point coord_;
  double vel_, anis_, minVel_, maxVel_, minAnis_, maxAnis_;
};

class Model {
 public:
  Model();  
  Model(Parameters&); 
  void SetVelocity(double&, int&, int&, int&);
  void SetAnisotropy(double&, int&, int&, int&);
  void CreateFildVelocity(std::vector<Celula>&);
  void CreateFildAnisotropy(std::vector<Celula>&);
  double GetVelocitySample(int&);
  double GetAnisotropySample(int&);
  double GetVelocitySampleCoord(int&, int&, int&);
  double GetAnisotropySampleCoord(int&, int&, int&);
  void SuavizaModel();
  int GetInvertParametrization();
  std::vector<double>& GetVelocity();
  std::vector<double>& GetAnisotropy();
  void SaveVelocity(const char*);
  void SaveVelocityVoronoi(const char*);
  void SaveAnisotropy(const char*) ;
  int GetWidth();
  int GetHeight(); 
  int GetDeep();
  double GetH();
  double GetMaxvel(); 
  double GetMinvel(); 
  double GetMaxanis();
  double GetMinanis();
  double GetVarSurv();
  int GetAniso(); 
  std::string GetNameDataTime(); 
 private:
  std::vector<double> datV, datE, datVoronoi;
  int width_, height_, deep_;
  double maxvel_,minvel_, maxanis_, minanis_, h_, varSurv_,sigSuavization_; 
  std::string name_data_time_;
  int aniso_, invert_parametrization_, suavization_;
};

class State{ 
  public:
    State();                                   // To store the energy E of the state
    State(Parameters&);
    void SetCelulas(std::vector<Celula>&);
    void SetCelula(Celula&, int&);
    void AddCelula(Celula&);
    void RemoveCelula(int&);
    std::vector<Celula>& GetCelulas();
    Celula GetCelula(int&);
    void InitRandomState(Parameters&, std::default_random_engine&);
    void Energy(Model&, std::vector<locations::Location>&, Grid&, std::vector<Data>&);
    double GetEnergy();
    int GetVarl();
    double GetVarE();
    double GetPassV(); 
    double GetPassE(); 
  private:
    std::vector<Celula> celulas_;                // To store the parameters that describe the state
    double energy_, varV_, varE_, passV_, passE_;
};


class Chain{       // To store a Markov chain
  public:
    void SetNumAt(int&);
    void SetNumRt(int&);
    int GetSaveAceite();
    int GetNumAt();
    int GetNumRt();
    Chain(State&, std::vector<Data>&);  
    Chain();  
    void SetState(State&);
    State GetState();
    void DelBeginState();
    void SaveEnergy(const char *);
    void RjMCMC(Model&, std::vector<locations::Location>&, Grid&, int&, std::default_random_engine&);
  private:
    std::vector<State> states_;
    std::vector<double>aceiteEnergy_;
    std::vector<Data>Data_obs_;      
    int at_, rt_, save_aceite_;                                   // Number of accepted MH transitions e reject  
};


#endif
