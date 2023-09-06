
#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <array>
#include "vec.hpp"
#include "stencil.hpp"
#include "serialization.hpp"

/** Type used to represent values on the grid */
using field = std::vector<double>;
/** Grid coordinate */
using coord = vec<unsigned, 3>;

/** Model class
 *
 * This class contains the whole program and is mainly used to be able to
 * scatter the implementation to different files without having to declare
 * every single variable extern. The architecture has not been really thought
 * through and one should not create two of these objects.
 * */
struct Model
{
  /** Simulation variables
   * @{ */

  /** List of neighbours
   *
   * The list of neighbours is used for computing the derivatives and is pre-
   * computed during initialization. The neighbours_patch variable has the same
   * role but for a region of the size of the cell patches (all patches have the
   * same size). These variables are computed in Initialize() and do not change
   * at runtime.
   * */
  std::vector<stencil> neighbors, neighbors_patch;

  /** Phase fields and derivatives */
  std::vector<field> phi, phi_dx, phi_dy, phi_dz;
  /** Predicted phi in a PC^n step */
  std::vector<field> phi_old;
  /** Phi difference */
  std::vector<field> dphi;
  /** Predicted phi difference in a P(C)^n step */
  std::vector<field> dphi_old;
  /** V = delta F / delta phi */
  std::vector<field> V;
  /** Sum of phi at each node */
  field sum;
  /** Sum_i S_i phi_i */
  field sumS00, sumS01, sumS02, sumS12, sumS11, sumS22; 
  /** Sum_i Q_i phi_i */
  field sumQ00, sumQ01;
  /** Sum_i vol_i phi_i */
  /** Sum of square, third, and fourth powers of phi at each node */
  field square, thirdp, fourthp, cIds;
  /** Phase-field for the walls and their derivatives */
  field walls, walls_dx, walls_dy, walls_dz, walls_laplace;
  /** Total polarization of the tissue */
  field P0, P1, P2;
  /** Total velocity of the tissue */
  field U0, U1, U2;
  /** vol associated with a phase field */
  std::vector<double> vol; 
  /** Forces */
  std::vector<vec<double, 3>> Fpol, Fnem, Fshape, Fpressure;
  /** Velocity */
  std::vector<vec<double, 3>> velocity;
  /** Structure tensor */
  std::vector<double> S00, S01, S02, S12, S11, S22;

  std::vector<double> gams;
  std::vector<double> omega_ccs;
  std::vector<double> omega_cws;
  std::vector<double> xis;
  std::vector<double> alphas;
  std::vector<double> kappas;
  std::vector<double> mus;
  std::vector<double> Rs;
  std::vector<double> V0;
  std::vector<double> zetaS_field;
  std::vector<double> zetaQ_field;
  std::vector<double> cellTypes;
   
  /** Q-tensor */
  std::vector<double> Q00, Q01;
  /** Polarisation */
  std::vector<vec<double, 3>> polarization;
  /** Direction of the polarisation */
  std::vector<double> theta_pol, theta_pol_old;
  /** Direction of the nematics */
  std::vector<double> theta_nem, theta_nem_old;
  /** Polarisation total torque */
  std::vector<double> delta_theta_pol;
  /** Stress tensor */
  field stress_xx, stress_xy , stress_yy, stress_xz, stress_zz, pressure;
  /** Elastic torque for the nematic */
  std::vector<double> tau;
  /** Vorticity around each cell */
  std::vector<vec<double, 3>> vorticity;
  /** Alignement options (see options.cpp) */
  int align_nematic_to = 0, align_polarization_to = 1;

  /** @} */

  /** Domain managment
   * @{ */

  /** Min of the boundaries of the patches and center of mass */
  std::vector<coord> patch_min;
  /** Max of the boundaries of the patches and center of mass */
  std::vector<coord> patch_max;
  /** Size of the patch (set by margin) */
  coord patch_size, patch_margin;
  /** Total number of nodes in a patch */
  unsigned patch_N;
  /** Memory offset for each patch */
  std::vector<coord> offset;
  /** Center-of-mass */
  std::vector<vec<double, 3>> com, com_prev;
  /** Counter to compute com in Fourier space */
  std::vector<std::complex<double>> com_x, com_y, com_z;
  /** Precomputed tables for sin and cos (as complex numbers) used in the
   * computation of the com.
   * */
  std::vector<std::complex<double>> com_x_table, com_y_table, com_z_table;

  /** @} */

  /** Program options
   * @{ */

  /** verbosity level
   *
   * 0: no output
   * 1: some output
   * 2: extended output (default)
   * 3: debug
   * */
  unsigned verbose = 2;
  /** compress output? (we use zip) */
  bool compress, compress_full;
  /** name of the run */
  std::string runname;
  /** Output dir (or tmp dir before moving files to the archive) */
  std::string output_dir;
  /** write any output? */
  bool no_write = false;
  /** skip runtime warnings? */
  bool no_warning = false;
  /** are the runtime warnings fatal? (i.e. they do stop the simulation) */
  bool stop_at_warning = false;
  /** shall we perform runtime checks ? */
  bool runtime_check = false;
  /** shall we print some runtime stats ? */
  bool runtime_stats = false;
  /** padding for onscreen output */
  unsigned pad;
  /** name of the inpute file */
  std::string inputname = "";
  /** Delete output? */
  bool force_delete;
  /** The random number seed */
  unsigned long seed; // picked using a fair die
  /** Flag if seed was set in the arguments, see options.hpp */
  bool set_seed;
  /** Number of predictor-corrector steps */
  unsigned npc = 1;
  /** Relaxation time at initialization */
  unsigned relax_time = 0;
  /** Value of nsubstep to use for initialization */
  unsigned relax_nsubsteps = 0;
  /** Total time spent writing output */
  std::chrono::duration<double> write_duration;
  /** @} */

  /** Simulation parameters
   * @{ */

  /** size of the system */
  vec<unsigned, 3> Size;
  /** size of the system */
  // vec<unsigned, 3> LDims;  
  /** total number of nodes */
  unsigned N;
  /** Total number of time steps */
  unsigned nsteps;
  /** Time interval between data outputs */
  unsigned ninfo = 10;
  /** Time at which to start the output */
  unsigned nstart = 0;
  /** number of subdivisions for a time step */
  unsigned nsubsteps = 10;
  /** effective time step */
  double time_step, sqrt_time_step;
  /** Number of phases */
  unsigned nphases;
  /** angle in degrees (input variable only) */
  double angle_deg;
  /** Margin for the definition of patches */
  unsigned margin = 25;
  /** Boudaries for cell generation
   *
   * These are the boundaries (min and max x and y components) of the domain in
   * which the cells are created when the initial config 'random' is choosen.
   * */
  std::vector<unsigned> birth_bdries;
  /** @} */

  /** Cell properties
   * @{ */

  /** Elasticity */
  double gam = 0.1;
  /** Energy penalty for vol */
  double mu = 3;
  /** Interface thickness */
  double lambda = 3;
  /**  Interaction stength */
  double kappa = 0.2;
  /** Adhesion */
  double omega = 0;
  /** Activity from shape */
  double zetaS = 0, sign_zetaS = 0;
  /** Activity from internal Q tensor */
  double zetaQ = 0, sign_zetaQ = 0;
  /** Propulsion strength */
  double alpha = 0;
  /** Substrate friction parameter */
  double xi = 1;
  /** Prefered radii (vol = pi*R*R) and radius growth */
  double R = 8;
  /** kij between dissimilar cell types */
  double kij = 5.;
  
  /** Base vol: V0 = Pi*R*R */
  // double V0;
  /** Repuslion by the wall */
  double wall_kappa = 0.5;
  /** Adhesion on the wall */
  double wall_omega = 0;
  /** Elasitc parameters */
  double Knem = 0, Kpol = 0;
  /** Strength of polarity / nematic tensor */
  double Spol = 0, Snem = 0;
  /** Flow alignment strenght */
  double Jpol = 0, Jnem = 0;
  /** Vorticity coupling */
  double Wnem = 0;
  /** Noise strength */
  double Dpol = 0, Dnem = 0;

  /** @} */
  /** Multi-threading parameters
   * @{ */
  #ifdef _OPENMP

  /** number of threads */
  unsigned nthreads;

  #endif
  /** @} */

  // ===========================================================================
  // Options. Implemented in options.cpp

  /** the variables map used to collect options */
  opt::variables_map vm;

  /** Set parameters from input */
  void ParseProgramOptions(int ac, char **av);

  /** Output all program options */
  void PrintProgramOptions();

  // =========================================================================
  // Program managment. Implemented in main.cpp

  /** The main loop */
  void Algorithm();

  /** Setup computation */
  void Setup(int, char**);

  /** Do the computation */
  void Run();

  /** Clean after you */
  void Cleanup();

  // ===========================================================================
  // Configuration. Implemented in configure.cpp

  /** Initial configuration parameters
   * @{ */

  /** Initial configuration name */
  std::string init_config;
  std::string sub_adh;
  /** Boundary conditions flag */
  unsigned BC = 0;
  /** Noise level for the nematic tensor initial configuration */
  double noise = 1;
  /** Wall thickness */
  double wall_thickness = 1.;
  /** Ratio of the cross vs size of the domain (BC=4) */
  double cross_ratio = .25;
  /** Ratio of the wound vs size of the domain (BC=5) */
  double wound_ratio = .50;
  /** Ratio of the tumor vs size of the domain (BC=6) */
  double tumor_ratio = .80;

  /** @} */

  /** Add cell with number n at a certain position */
  void AddCell(unsigned n, const coord& center);
  void AddCellMix(unsigned n, const coord& center, double zetaSi, double zetaQi, double gami, double omegai, double omega_walli, double kappai, double mui, double alphai, double xi, double Ri, double cellType);

  /** Subfunction for AddCell() */
  void AddCellAtNode(unsigned n, unsigned q, const coord& center);

  /** Set initial condition for the fields */
  void Configure();

  /** Set initial configuration for the walls */
  void ConfigureWalls(int BC);

  // ==========================================================================
  // Writing to file. Implemented in write.cpp

  /** Write current state of the system */
  void WriteFrame(unsigned);
  
  /** (added by Siavash) Write phase-field for cell n */
  void Write_phi(unsigned);
  void Write_dphi(unsigned);
  void Write_COM(unsigned);
  void Write_velocities(unsigned);    
  void Write_forces(unsigned);    
  void Write_visData(unsigned);
  void Write_contArea(unsigned);
  void Write_Density(unsigned);
  
  /** Write run parameters */
  void WriteParams();

  /** Remove old files */
  void ClearOutput();

  /** Create output directory */
  void CreateOutputDir();

  // ==========================================================================
  // Initialization. Implemented in init.cpp

  /** Initialize memory for field */
  void Initialize();

  /** Allocate memory for individual cells */
  void SetCellNumber(unsigned new_nphases);

  /** Initialize neighbors list (stencils) */
  void InitializeNeighbors();

  /** Swap two cells in the internal arrays */
  void SwapCells(unsigned n, unsigned m);

  // ===========================================================================
  // Random numbers generation. Implemented in random.cpp

  /** Pseudo random generator */
  std::mt19937 gen;
  //ranlux24 gen;

  /** Return random real, uniform distribution */
  double random_real(double min=0., double max=1.);

  /** Return random real, gaussian distributed */
  double random_normal(double sigma=1.);

  /** Return geometric dist numbers, prob is p */
  unsigned random_geometric(double p);

  /** Return poisson distributed unsigned integers */
  unsigned random_poisson(double lambda);

  /** Return exp distributed unsigned integers */
  int random_exponential(double lambda);

  /** Return random unsigned uniformly distributed */
  unsigned random_unsigned();

  /** Initialize random numbers
   *
   * If CUDA is enabled, alos intialize CUDA random numbers
   * */
  void InitializeRandomNumbers();

  // ==========================================================================
  // Support for cuda. Implemented in cuda.cu

  #ifdef _CUDA_ENABLED

  /** Device(s) propeties
   * @{ */

  /** Obtain (and print) device(s) properties
   *
   * Also checks that the device properties are compatible with the values given
   * in src/cuda.h.
   * */
  void QueryDeviceProperties();

  /** @} */

  /** Pointer to device global memory
   *
   * These pointers reflects the program data strcuture and represents the cor-
   * responding data on the device global memory. All names should be identical
   * to their host counterparts apart from the d_ prefix.
   *
   * @{ */
    
  double *d_phi, *d_phi_old, *d_V, *d_potential, *d_potential_old, *d_sum,
         *d_square, *d_Q00, *d_Q01, *d_walls,
         *d_walls_laplace, *d_walls_dx, *d_walls_dy, *d_walls_dz, *d_sum_cnt, *d_square_cnt,
         *d_Q00_cnt, *d_Q01_cnt, *d_vol, *d_vol_cnt, *d_c, *d_S00, *d_S01, *d_S02, *d_S12, *d_S11, *d_S22, *d_S_order,
         *d_S_angle, *d_theta, *d_alpha, *d_gam, *d_mu, *d_gams, *d_zetaS_field, *d_zetaQ_field, *d_cellTypes, *d_alphas,*d_omega_ccs, *d_omega_cws, *d_xis ,*d_kappas, *d_mus, *d_Rs, *d_V0;
  vec<double, 3>  *d_vel, *d_force_p, *d_force_c, *d_force_f, *d_com, *d_com_prev;
  stencil         *d_neighbors, *d_neighbors_patch;
  coord           *d_patch_min, *d_patch_max, *d_offset;
  cuDoubleComplex *d_com_x, *d_com_y, *d_com_z, *d_com_x_table, *d_com_y_table, *d_com_z_table;

  /** @} */

  /** Random number generation
   * @{ */

  /** Random states on the device */
  curandState *d_rand_states;

  /** Initialization function */
  void InitializeCuda();

  /** Initialization function for random numbers */
  void InitializeCUDARandomNumbers();

  /** @} */
  /** CUDA device memory managment
    * @{ */

  /** In which direction do we copy data? */
  enum class CopyMemory {
    HostToDevice,
    DeviceToHost
  };

  /** Allocate or free memory? */
  enum class ManageMemory {
    Allocate,
    Free
  };

  /** Implementation for AllocDeviceMemory() and FreeDeviceMemory() */
  void _manage_device_memory(ManageMemory);
  /** Implementation for PutToDevice() and GetFromDevice() */
  void _copy_device_memory(CopyMemory);

  /** Copy data to the device global memory
   *
   * This function is called at the begining of the program just before the main
   * loop but after the system has been initialized.
   * */
  void PutToDevice();

  /** Copy data from the device global memory
   *
   * This function is called every time the results need to be dumped on the disk.
   * */
  void GetFromDevice();

  /** Allocate memory for all device arrays */
  void AllocDeviceMemory();

  /** Allocate memory for all device arrays */
  void FreeDeviceMemory();

  /** @} */

  /** Runtime properties
   * @{ */

  int n_total, n_blocks, n_threads;

  /** @} */

  #endif

  // ===========================================================================
  // Run. Implemented in run.cpp

  /** Time step
   *
   * This is the time-stepping function and performs the main computation.
   * */
  void Step();

  /** Prepare before run */
  void Pre();

  /** Prints some stats before running */
  void PreRunStats();

  /** Prints some stats in between ninfo steps */
  void RuntimeStats();

  /** Performs punctual check at runtime */
  void RuntimeChecks();

  /** Post run function */
  void Post();

  double SubAdhesion_field(double,double);
  /** Subfunction for update */    
  void UpdatePotAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdatePhaseFieldAtNode(unsigned, unsigned, bool);

  /** Subfunction for update */
  void UpdateForcesAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdateStructureTensorAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void UpdateSumsAtNode(unsigned, unsigned);

  /** Subfunction for update */
  void ReinitSumsAtNode(unsigned);

  /** Compute center of mass of a given phase field */
  void ComputeCoM(unsigned);

  /** Update polarisation of a given field */
  void UpdatePolarization(unsigned, bool);

  /** Update nematic tensor of a given field */
  void UpdateNematic(unsigned, bool);

  /** Compute shape parameters
   *
   * This function effectively computes the second moment of vol, which ca n be used to
   * fit the shape of a cell to an ellipse.
   * */
  void ComputeShape(unsigned);

  /** Update the moving patch following each cell */
  void UpdatePatch(unsigned);

  /** Update fields
   *
   * The boolean argument is used to differentiate between the predictor step
   * (true) and subsequent corrector steps.
   * */
  void Update(bool, unsigned=0);

  // ===========================================================================
  // Serialization

  /** Serialization of parameters (in and out) */
  template<class Archive>
  void SerializeParameters(Archive& ar)
  {
    ar & auto_name(gam)
       & auto_name(mu)
       & auto_name(lambda)
       & auto_name(nphases)
       & auto_name(init_config)
	   & auto_name(sub_adh)
       & auto_name(kappa)
       & auto_name(xi)
       & auto_name(kij)
       & auto_name(R)
       & auto_name(alpha)
       & auto_name(zetaS)
       & auto_name(zetaQ)
       & auto_name(omega)
       & auto_name(wall_thickness)
       & auto_name(wall_kappa)
       & auto_name(wall_omega)
       & auto_name(walls)
       & auto_name(patch_margin)
       & auto_name(relax_time)
       & auto_name(relax_nsubsteps)
       & auto_name(npc)
       & auto_name(seed)
       & auto_name(Knem)
       & auto_name(Kpol)
       & auto_name(Snem)
       & auto_name(Spol)
       & auto_name(Jnem)
       & auto_name(Jpol)
       & auto_name(Wnem)
       & auto_name(Dpol)
       & auto_name(Dnem)
       & auto_name(margin)
       & auto_name(patch_size);
  }

  /** Serialization of parameters (in and out) */
  template<class Archive>
  void SerializeFrame(Archive& ar)
  {
    ar & auto_name(nphases)
       & auto_name(phi)
       & auto_name(offset)
       & auto_name(com)
       & auto_name(vol)
       & auto_name(S00)
       & auto_name(S01)
       & auto_name(S02)
       & auto_name(S12)
       & auto_name(S11)
       & auto_name(S22)     
       & auto_name(Q00)
       & auto_name(Q01)
       & auto_name(gams)
       & auto_name(zetaS_field)
       & auto_name(zetaQ_field)
       & auto_name(cellTypes)
       & auto_name(alphas)
       & auto_name(omega_ccs)
       & auto_name(omega_cws)
       & auto_name(xis)
       & auto_name(kappas)
       & auto_name(mus)
       & auto_name(Rs)
       & auto_name(V0)
       & auto_name(stress_xx)
       & auto_name(stress_xy)
       & auto_name(stress_yy)
       & auto_name(stress_xz)
       & auto_name(stress_zz)
       & auto_name(velocity)
       & auto_name(vorticity)
       & auto_name(Fpol)
       & auto_name(Fnem)
       & auto_name(Fshape)
       & auto_name(Fpressure)
       & auto_name(theta_pol)
       & auto_name(theta_nem)
       & auto_name(patch_min)
       & auto_name(patch_max);
  }
  
  // ===========================================================================
  // Tools
  
  /** Gives domain coordinate corresponding to a domain index */
  coord GetPosition(unsigned k) const
  { return { GetXPosition(k), GetYPosition(k), GetZPosition(k) }; }

  /** Gives domain x-coordinate corresponding to a domain index */
  unsigned GetXPosition(unsigned k) const
  { return (k/Size[1])%Size[0];}

  /** Gives domain y-coordinate corresponding to a domain index */
  unsigned GetYPosition(unsigned k) const
  {return k%Size[1];}
  
  unsigned GetZPosition(unsigned k) const 
  {return k/(Size[0]*Size[1]);}
    
  /** Get domain index from domain coordinates */
  unsigned GetIndex(const coord& p) const
  {return p[1] + Size[1]*p[0] + Size[0]*Size[1]*p[2];}

  /** Get patch index from domain coordinates */
  unsigned GetPatchIndex(unsigned n, coord p) const
  {
    // get difference to the patch min
    p = (p + Size - patch_min[n])%Size;
    // correct for offset
    p = (p + patch_size - offset[n])%patch_size;
    // remap linearly
    return p[1] + patch_size[1]*p[0] + patch_size[0]*patch_size[1]*p[2];
  }
  
  /** Get patch index from domain index */
  unsigned GetPatchIndex(unsigned n, unsigned k) const
  {
    return GetPatchIndex(n, GetPosition(k));
  }

  /** Get domain index from patch index */
  unsigned GetIndexFromPatch(unsigned n, unsigned q) const
  {
    // position on the patch
    const coord qpos = { (q/patch_size[1])%patch_size[0] , q%patch_size[1]  , q/( patch_size[0]*patch_size[1] ) };
    // position on the domain
    const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
    // return domain index
    return GetIndex(dpos);
  }
  
  coord GetNodePosOnPatch(unsigned n, unsigned q) const
  {
    // position on the patch
    const coord qpos = { (q/patch_size[1])%patch_size[0] , q%patch_size[1]  , q/( patch_size[0]*patch_size[1] ) };
    return qpos;
  }  
  
  coord GetNodePosOnDomain(unsigned n, unsigned q) const
  {
    // position on the patch
    const coord qpos = { (q/patch_size[1])%patch_size[0] , q%patch_size[1]  , q/( patch_size[0]*patch_size[1] ) };
    // position on the domain
    const coord dpos = ( (qpos + offset[n])%patch_size + patch_min[n] )%Size;
    return dpos;
  }    
  
};

#endif//MODEL_HPP_
