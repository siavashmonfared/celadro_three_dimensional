
#include "header.hpp"
#include "model.hpp"

using namespace std;

void Model::ParseProgramOptions(int ac, char **av)
{
  // ===========================================================================
  // Declaration

  // options allowed only in the command line
  opt::options_description generic("Generic options");
  generic.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v", opt::value<unsigned>(&verbose)->implicit_value(2),
     "verbosity level (0=none, 1=little, 2=normal, 3=debug)")
    ("input,i", opt::value<string>(&inputname),
     "input file")
    ("force-delete,f", opt::bool_switch(&force_delete),
     "force deletion of existing output file")
#ifdef _OPENMP
    ("threads,t",
     opt::value<unsigned>(&nthreads)->default_value(0)->implicit_value(1),
     "number of threads (0=no multithreading, 1=OpenMP default, "
     ">1=your favorite number)")
#endif
    ("compress,c", opt::bool_switch(&compress),
     "compress individual files using zip")
    ("compress-full", opt::bool_switch(&compress_full),
     "compress full output using zip (might be slow)")
    ("no-write", opt::bool_switch(&no_write),
     "disable file output (for testing purposes)");

  // options allowed both in the command line and config file
  opt::options_description config("Program options");
  config.add_options()
    ("output,o", opt::value<string>(&runname),
     "output name (if compression is on, then .zip is added automatically)")
    ("seed", opt::value<unsigned long>(&seed),
     "set seed for random number generation (random if unset)")
    ("no-warning", opt::bool_switch(&no_warning),
     "disable model specific runtime warnings")
    ("stop-at-warning", opt::bool_switch(&stop_at_warning),
     "runtime warnings interrupt the algorithm")
    ("check", opt::bool_switch(&runtime_check),
     "perform runtime checks")
    ("stat", opt::bool_switch(&runtime_stats),
     "print runtime stats")
    ("nstart", opt::value<unsigned>(&nstart)->default_value(0u),
     "time at which to start the output")
    ("bc", opt::value<unsigned>(&BC)->default_value(0u),
     "boundary conditions flag (0=pbc, 1=box, 2=channel, 3=ellipse)");

  // model specific options
  opt::options_description simulation("Simulation options");
  simulation.add_options()
    ("LX", opt::value<unsigned>(&Size[0]),
     "# of nodes in the x direction")
    ("LY", opt::value<unsigned>(&Size[1]),
     "# of nodes in the y direction")
    ("LZ", opt::value<unsigned>(&Size[2]),
     "# of nodes in the y direction")
    ("nsteps", opt::value<unsigned>(&nsteps),
     "iterate this many steps in total")
    ("nsubsteps", opt::value<unsigned>(&nsubsteps)->default_value(1u),
     "subdivision of a single time step "
     "(effective time step is 1/nsubsteps)")
    ("ninfo", opt::value<unsigned>(&ninfo),
     "save frame every so many steps")
    ("nphases", opt::value<unsigned>(&nphases),
      "Number of phases")
    ("gamma", opt::value<double>(&gam),
      "Elastic constant of each phase (array)")
    ("mu", opt::value<double>(&mu),
      "Energy penalty for vol of each phase (array)")
    ("lambda", opt::value<double>(&lambda),
      "Interface thickness parameter")
    ("kappa_cc", opt::value<double>(&kappa_cc),
      "Interaction strength")
    ("npc", opt::value<unsigned>(&npc)->default_value(1u),
      "Number of predictor-corrector steps")
    ("margin", opt::value<unsigned>(&margin)->default_value(0u),
      "Margin for the definition of restricted domains (if 0: update full box)")
    //("friction", opt::value<double>(&f),
    //  "Cell-cell friction parameter")
    //("friction-walls", opt::value<double>(&f_walls),
    //  "Cell-wall friction parameter")
    ("xi", opt::value<double>(&xi),
      "Substrate friction parameter")
    ("K-pol", opt::value<double>(&Kpol),
     "elastic constant for the polarisation")
    ("K-nem", opt::value<double>(&Knem),
     "elastic constant for the nematic")
    //("C-pol", opt::value<double>(&Cpol),
    // "Strength of LdG potential for the polarisation")
    ("J-pol", opt::value<double>(&Jpol),
     "Nematic flow alignment strength")
    ("J-nem", opt::value<double>(&Jnem),
     "Nematic flow alignment strength")
    ("W-nem", opt::value<double>(&Wnem),
     "Strength of vorticity torque")
    ("D-nem", opt::value<double>(&Dnem),
      "Nematic noise strength")
    ("D-pol", opt::value<double>(&Dpol),
      "Polarisation noise strength")
    ("S-nem", opt::value<double>(&Snem),
      "Order of the nematic tensors")
    ("S-pol", opt::value<double>(&Spol),
      "Norm of the polarisation vector")
    ("alpha", opt::value<double>(&alpha),
     "Strength of propulsion")
    ("zetaQ", opt::value<double>(&zetaQ),
     "Activity from internal nematic tensor")
    ("zetaS", opt::value<double>(&zetaS),
     "Activity from shape")
    ("omega_cc", opt::value<double>(&omega_cc),
      "Adhesion parameter")
    ("wall-thickness", opt::value<double>(&wall_thickness),
      "Wall thickness (typical decay length)")
    ("kappa_cs", opt::value<double>(&kappa_cs)->default_value(kappa_cc),
      "Wall repulsion")
    ("omega_cs", opt::value<double>(&omega_cs)->default_value(0.),
      "Wall adhesion")
    ("R", opt::value<double>(&R),
      "Preferred radius (defines volume (4/3)*Pi*R*R*R)")
    ("align-polarization-to", opt::value<int>(&align_polarization_to),
     "Align polarization to velocity (=0) or pressure force (=1)")
    ("align-nematic-to", opt::value<int>(&align_nematic_to),
     "Align nematic tensor to velocity (=0), pressure force (=1), or shape (=2)");
     
  // init config options
  opt::options_description init("Initial configuration options");
  init.add_options()
    ("config", opt::value<string>(&init_config),
      "Initial configuration")
    ("relax-time", opt::value<unsigned>(&relax_time)->default_value(0u),
      "Relaxation time steps at initialization.")
    ("noise", opt::value<double>(&noise),
      "Noise level for initial nematic angle, in (0,1).")
    ("cross-ratio", opt::value<double>(&cross_ratio),
      "Ratio of the size of the cross compared to the domain size (for BC=4)")
    ("wound-ratio", opt::value<double>(&wound_ratio),
      "Ratio of the size of the wound open space compared to the domain size (for BC=5)")
    ("tumor-ratio", opt::value<double>(&tumor_ratio),
      "Ratio of the size of the tumor compared to the domain size (for BC=6)")
    ("birth-boundaries", opt::value<vector<unsigned>>(&birth_bdries)->multitoken(),
     "Boundaries in which the cells are created "
     "when the initial configuration 'random' is choosed. "
     "Format: {min x, max, x, min y, max y}")
    ("relax-nsubsteps", opt::value<unsigned>(&relax_nsubsteps)->default_value(0u),
      "Value of nsubsteps to use at initial relaxation (0 means use nsubsteps).");

  // ===========================================================================
  // Parsing

  // command line options
  opt::options_description cmdline_options;
  cmdline_options.add(generic).add(config).add(simulation).add(init);
  // config file options
  opt::options_description config_file_options;
  config_file_options.add(config).add(simulation).add(init);

  // first unnamed argument is the input file
  opt::positional_options_description p;
  p.add("input", 1);

  // reintialize vm in case we run this function twice
  vm = opt::variables_map();

  // parse first the cmd line
  opt::store(
    opt::command_line_parser(ac, av)
    .options(cmdline_options)
    .positional(p)
    .run(), vm);
  opt::notify(vm);

  // print help msg and exit
  if(vm.count("help"))
  {
    cout << config_file_options << endl;
    exit(0);
  }

  // parse input file (values are not erased, such that cmd line args
  // are 'stronger')
  if(inputname.empty())
    throw error_msg("please provide an input file / type -h for help.");
  else
  {
    std::fstream file(inputname.c_str(), std::fstream::in);
    if(!file.good()) throw error_msg("can not open runcard file ", inputname);
    opt::store(opt::parse_config_file(file, config_file_options), vm);
    opt::notify(vm);
  }

  // ===========================================================================
  // Fixing some secondary values

  // fix compression mode: if we compress the full archive we do not compress
  // individual files.
  if(compress_full) compress=false;

  // Set default value for runname (depends on compression)
  if(vm.count("output")==0)
  {
    if(compress_full) runname = "output";
    else runname = "./";
  }

  // init random numbers?
  set_seed = vm.count("seed");

  // compute the correct padding
  pad = inline_str(nsteps).length();

  // compute effective time step
  time_step = 1./nsubsteps;

  // set nstart to the next correct frame (round above)
  if(nstart%ninfo) nstart = (1u+nstart/ninfo)*ninfo;
}

/** Print variables from variables_map
  *
  * from: https://gist.github.com/gesquive/8673796
  */
void print_vm(const opt::variables_map& vm, unsigned padding)
{
  for (opt::variables_map::const_iterator it = vm.begin(); it != vm.end(); ++it)
  {
    // pass if defaulted or empty
    if (vm[it->first].defaulted() || it->second.defaulted()) continue;
    if (((boost::any)it->second.value()).empty()) continue;

    std::cout << std::left << std::setw(floor(padding/2)) << it->first;


    std::cout << std::right << std::setw(ceil(padding/2));

    bool is_char;
    try {
      boost::any_cast<const char*>(it->second.value());
      is_char = true;
    } catch (const boost::bad_any_cast &) {
      is_char = false;
    }
    bool is_str;
    try {
      boost::any_cast<std::string>(it->second.value());
      is_str = true;
    } catch (const boost::bad_any_cast &) {
      is_str = false;
    }

    if (((boost::any)it->second.value()).type() == typeid(int)) {
      std::cout << vm[it->first].as<int>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(unsigned)) {
      std::cout << vm[it->first].as<unsigned>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(size_t)) {
      std::cout << vm[it->first].as<size_t>() << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
      std::cout << (vm[it->first].as<bool>() ? "true" : "false") << std::endl;
    } else if (((boost::any)it->second.value()).type() == typeid(double)) {
      std::cout << vm[it->first].as<double>() << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<double>)) {
      std::cout << vec2str(vm[it->first].as<vector<double>>()) << std::endl;
    } else if (((boost::any)it->second.value()).type()
               == typeid(vector<unsigned>)) {
      std::cout << vec2str(vm[it->first].as<vector<unsigned>>()) << std::endl;
    } else if (is_char) {
      std::cout << vm[it->first].as<const char *>() << std::endl;
    } else if (is_str) {
      std::string temp = vm[it->first].as<std::string>();
      if (temp.size()) {
        std::cout << temp << std::endl;
      } else {
        std::cout << "true" << std::endl;
      }
    } else { // Assumes that the only remainder is vector<string>
      try {
        auto vect = vm[it->first].as<std::vector<std::string> >();
        uint i = 0;
        for (auto oit=vect.begin();
            oit != vect.end(); oit++, ++i) {
          std::cout << "\r> " << it->first
                    << "[" << i << "]=" << (*oit) << std::endl;
        }
      } catch (const boost::bad_any_cast &) {
        std::cout << "UnknownType("
                  << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
      }
    }
  }
}

void Model::PrintProgramOptions()
{
  // print the simulation parameters
  print_vm(vm, width);
}
