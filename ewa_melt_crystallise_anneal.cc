#include<iostream>
#include<vector>
#include<map>
#include<cstdlib>
#include<cmath>
#include<fstream>



//=============================================================
/// Namespace for parameters
//=============================================================
namespace GlobalParameters
{

 // Primary crystallisation
 //------------------------

 /// Magnitude of crystallation rate fct
 double G_0_cryst=1.0;

 /// Temperature at which crystallisation rate is max.
 double T_max_cryst=5.0;

 /// Width of Gaussian for crystallisation rate
 double S_cryst=0.1;

 /// Gaussian crystallation rate as a function of current temperature
 double G_cryst_rate(const double& temp)
 {
  return G_0_cryst*exp(-pow((temp-T_max_cryst)/S_cryst,2));
 }


 // Crystal size
 //-------------

 /// Volume factor (1 for cubic crystal)
 double K_crystal=1.0;

 /// Crystal size (volume fraction) as function of scalar size parameter
 double crystal_volume(const double& l)
 {
  return K_crystal*l*l*l;
 }


 // Primary crystal formation/melting characteristics
 //--------------------------------------------------

 /// Crystal size formed at zero temperature
 double L_c_0=0.1;

 /// Temperature at which crystals grow to infinite size
 /// (or the temperature at which even infinitely large crystals melt)
 double T_infinite_cryst=12.0;
 
 /// Scalar size of crystal that forms at given temperature
 double size_of_new_crystal(const double& temperature)
 {
  return L_c_0/(1.0-temperature/T_infinite_cryst);
 }


 /// Melt temperature of crystal of given size
 double melt_temperature(const double& l)
 {
  return T_infinite_cryst*(1.0-L_c_0/l);
 }



 // Secondary crystallisation
 //--------------------------


 /// Scaling parameter for growth rate
 double K_growth_1=1.0;

 /// Exponential parameter for growth rate
 double K_growth_2=1.0;

 /// Boltzmann constant
 double K_boltzmann=1.0;
 
 /// Secondary growth rate (dL/dt) as a function of the current crystal
 /// length, l, the temperature, temp, and the current value of the crystalline
 /// volume fraction, assumed to be related to melt via c = 1-m
 double secondary_growth_rate(const double& l,
                              const double& temp,
                              const double& c)
 {
  double m=1.0-c;
  return m*K_growth_1*exp(-(K_growth_2*l/(K_boltzmann*temp)));
 }
 

 // Temperature variation
 //----------------------

 // Initial temperature at t=0
 double T_init=1.2*T_max_cryst;

 // End temperature at t=1
 double T_end=0.8*T_max_cryst;
 
 /// current temperature
 double current_temperature(const double& time)
 {
  // Cut off
  return std::max(T_end,T_init+(T_end-T_init)*time);
 }


 

 // Info about crystal phase
 //-------------------------
 
 /// Current size (scalar) of crystals initially created at timestep i
 std::vector<double> Current_size_of_crystals_initially_created_at_timestep;
 
 /// Number of crystals created at timestep i (represent as double)
 std::vector<double> Number_of_crystals_initially_created_at_timestep;


 /// Volume currently occupied by crystalline phase
 double volume_of_crystalline_phase()
 {
  double vol=0.0;
  unsigned n=Current_size_of_crystals_initially_created_at_timestep.size();
  for (unsigned i=0;i<n;i++)
   {
    vol+=crystal_volume(
     Current_size_of_crystals_initially_created_at_timestep[i])*
     Number_of_crystals_initially_created_at_timestep[i];
   }
  return vol;
 }

 
}


//------------------------------------------------------------
// Driver
//------------------------------------------------------------
int main()
{



 // Number of timesteps
 unsigned nstep=1000;

 // Max. time
 double t_max=10.0;
 
 // Timestep
 double dt=t_max/double(nstep);

 // Initial doc of stuff
 std::ofstream some_file;
 some_file.open("crystallsation_properties.dat");
 unsigned nplot=100;
 double T_max=GlobalParameters::T_infinite_cryst;
 double dT=T_max/double(nplot);
 for (unsigned i=0;i<nplot;i++)
  {
   double temp=double(i)*dT;
   double l=GlobalParameters::size_of_new_crystal(temp);
   // This should return the temperature again
   double temp_check=GlobalParameters::melt_temperature(l);
   some_file << temp << " "
             << l << " "
             << temp_check << " "
             << std::endl;
  }
 some_file.close();



 // Create storage:
 //----------------

 // Current size (scalar) of crystals initially created at timestep i
 GlobalParameters::Current_size_of_crystals_initially_created_at_timestep.
  reserve(nstep);
 
 // Number of crystals created at timestep i
 GlobalParameters::Number_of_crystals_initially_created_at_timestep.
  reserve(nstep);


 std::ofstream trace_file;
 trace_file.open("trace.dat");
 
 std::ofstream size_distribution_file;
 size_distribution_file.open("size_distribution.dat");

 // Timestepping loop
 double time=0.0;
 for (unsigned i=0;i<nstep;i++)
  {
   time+=dt;

   // Current volume of crystalline phase (from previous timestep)
   double c=GlobalParameters::volume_of_crystalline_phase();
  
   // Current temperature
   double temp=GlobalParameters::current_temperature(time);


   // STAGE 1: CRYSTALLISE STUFF FROM MELT.
   //--------------------------------------
   
   // Get the size (scalar) of the crystals that are currently being formed
   double l_new=GlobalParameters::size_of_new_crystal(temp);

   // Push back
   GlobalParameters::Current_size_of_crystals_initially_created_at_timestep.
    push_back(l_new);
   
   // Incremental volume fraction that is currently being crystallised
   // into crystals of (scalar) size l_new over current timestep;
   // First order explicit Euler scheme.
   double dc_new_cryst = GlobalParameters::G_cryst_rate(temp)*(1.0-c)*dt;

   // How many crystals fit into this incremental volume? (double prec)
   double n_new_cryst=dc_new_cryst/GlobalParameters::crystal_volume(l_new);
   GlobalParameters::Number_of_crystals_initially_created_at_timestep.
    push_back(
    n_new_cryst);

   std::cout << "during primary crystallisation phase, "
             << "total crystalline volume changed from "
             << c << " to " << GlobalParameters::volume_of_crystalline_phase()
             << " corresp to change of " 
             << GlobalParameters::volume_of_crystalline_phase()-c  << " a.k.a. "
             << dc_new_cryst << std::endl;
   

   // STATE 2: GROW THE CRYSTALS THAT WERE CREATED AT PREVIOUS TIMES (FIRST
   //----------------------------------------------------------------------
   // ORDER EULER SCHEME.
   //--------------------
   double c_before_growth=GlobalParameters::volume_of_crystalline_phase();

   // Deal with all the timesteps excluding the current one
   unsigned n=
    GlobalParameters::Number_of_crystals_initially_created_at_timestep.size();
   for (unsigned j=0;j<n-1;j++)
    {
     double num_of_cryst=
      GlobalParameters::Number_of_crystals_initially_created_at_timestep[j];
     if (num_of_cryst>0.0)
      {
       double l=
        GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j];
       
       // Get (linear) growth rate
       double dL_dt=
        GlobalParameters::secondary_growth_rate(l,temp,c_before_growth);
       std::cout << "dL/dt = " << dL_dt << std::endl;
       
       // Bump (first-order Euler scheme)
       GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]
        +=dL_dt*dt;
      }
    }

   double dc_growth=
    GlobalParameters::volume_of_crystalline_phase()-c_before_growth;
   std::cout << "during secondary growth phase, "
             << "total crystalline volume changed from "
             << c_before_growth << " to "
             << GlobalParameters::volume_of_crystalline_phase()
             << " corresp to change of " 
             << dc_growth
             << std::endl;
   

   
   // STAGE 3: MELT STUFF
   //--------------------
   double c_before_melt=GlobalParameters::volume_of_crystalline_phase();
   
   // Loop over all current individual crystal sizes and work 
   // out their melt temperature based on their linear size.
   // If the temperature is above the melt temperature reset the
   // number of these crystals to zero.
   // hierher later just pop them from a collection of pairs.
   double dc_new_melt=0.0;
    for (unsigned j=0;j<n;j++)
    {
     double num_of_cryst=
      GlobalParameters::Number_of_crystals_initially_created_at_timestep[j];
     if (num_of_cryst>0.0)
      {
       double l=
        GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j];
       double T_melt=GlobalParameters::melt_temperature(l);
       if (T_melt<temp)
        {
         dc_new_melt+=num_of_cryst*GlobalParameters::crystal_volume(l);
         std::cout << "melting " << num_of_cryst << " crystals of size " << l
                   << " i.e. of volume " << GlobalParameters::crystal_volume(l)
                   << " so total vol "
                   << num_of_cryst*GlobalParameters::crystal_volume(l)
                   << std::endl;

         // Kill 'em
         GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
          =0.0;
        }
      }
    }
   
   std::cout << "during melting phase, total crystalline volume changed from "
             << c_before_melt << " to "
             << GlobalParameters::volume_of_crystalline_phase()
             << " corresp to change of " 
             << GlobalParameters::volume_of_crystalline_phase()-c_before_melt
             << " a.k.a. "
             << dc_new_melt << std::endl;
   
   // Doc stuff
   trace_file << time << " "
              << temp << " "
              << l_new << " "
              << GlobalParameters::G_cryst_rate(temp) << " " 
              << dc_new_cryst << " "
              << dc_growth << " "
              << dc_new_melt << " " 
              << GlobalParameters::volume_of_crystalline_phase()<< " "
              << std::endl;

   
   // Size distribution
   {
    unsigned n=
     GlobalParameters::Number_of_crystals_initially_created_at_timestep.size();
    for (unsigned j=0;j<n;j++)
     {
      size_distribution_file
       <<  GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j] << " "
       << GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
       << " "
       << std::endl;
     }
    // Double blank line for gnuplot
    size_distribution_file  << std::endl << std::endl;
   }
   
  }
 
  }
