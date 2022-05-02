#include<iostream>
#include<vector>
#include<map>
#include<cstdlib>
#include<cmath>
#include<fstream>


#include<fenv.h>

//#define PARANOID

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
  double t_switch_over=3.0;
  if (time<t_switch_over)
   {
    return std::max(T_end,T_init+(T_end-T_init)*time);
   }
  else
   {
    return std::min(0.9*T_infinite_cryst,
                    T_end-(T_end-T_init)*(time-t_switch_over)); 
   }
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


 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Number of timesteps
 unsigned nstep=10000; // 10000;

 // Max. time
 double t_max=10.0; //10.0;
 
 // Timestep
 double dt=t_max/double(nstep);

 // Initial doc of stuff
 std::ofstream some_file;
 some_file.open("RESLT/crystallsation_properties.dat");
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
 trace_file.open("RESLT/trace.dat");

 
 #ifdef PARANOID
 double max_error_stage1=0.0;
 double max_error_stage2=0.0;
 double max_error_stage3=0.0;
 #endif
 
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
    push_back(n_new_cryst);

   // std::cout << "time, temp, l_new, n_new, dc_new "
   //           << time << " "
   //           << temp << " "
   //           << l_new << " "
   //           << n_new_cryst << " "
   //           << dc_new_cryst << " "
   //           << std::endl;
    
#ifdef PARANOID
   // Check error
   {
    double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c)-dc_new_cryst);
    max_error_stage1=std::max(max_error_stage1,error);
    //std::cout << "Error stage 1: " << error << std::endl;
   }
#endif
  
   

   // STATE 2: GROW THE CRYSTALS THAT WERE CREATED AT PREVIOUS TIMES (FIRST
   //----------------------------------------------------------------------
   // ORDER EULER SCHEME.
   //--------------------

#ifdef PARANOID
   double c_before_growth=GlobalParameters::volume_of_crystalline_phase();
#endif
   
   // Deal with all the timesteps excluding the current one
   double dc_growth=0.0;
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
       
       // Get (linear) growth rate: Forward Euler so base c on crystalline phase
       // at previous timestep
       double dL_dt=
        GlobalParameters::secondary_growth_rate(l,temp,c);
       
       // Bump (first-order Euler scheme)
       double dL=dL_dt*dt;

       double v_old=num_of_cryst*GlobalParameters::crystal_volume(l);
       double v_new=num_of_cryst*GlobalParameters::crystal_volume(l+dL);
       dc_growth+=(v_new-v_old);
       
       GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]
        +=dL_dt*dt;
      }
    }

#ifdef PARANOID
   // Check error
   {
    double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c_before_growth)-dc_growth);
    max_error_stage2=std::max(max_error_stage2,error);
    //std::cout << "Error stage 2: " << error << std::endl;
   }
#endif


   
   // STAGE 3: MELT STUFF
   //--------------------
   double dc_new_melt=0.0;
   bool allow_melting=true;
   if (allow_melting)
    {
#ifdef PARANOID
     double c_before_melt=GlobalParameters::volume_of_crystalline_phase();
#endif
     
     // Loop over all current individual crystal sizes and work 
     // out their melt temperature based on their linear size.
     // If the temperature is above the melt temperature reset the
     // number of these crystals to zero.
     // hierher later just pop them from a collection of pairs.
     for (unsigned j=0;j<n-1;j++) // hierher don't melt newly created crystals straight away
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
           //std::cout << "temp : " << temp << " melting size " << l << std::endl;
           dc_new_melt-=num_of_cryst*GlobalParameters::crystal_volume(l);
           
           // Kill 'em (sizes are set to zero too so they don't clunk up the animations
           // and 0 x 0 is still zero!
           GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
            =0.0;
           GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]
            =0.0;
          }
         else
          {
           //std::cout << "temp : " << temp << " not melting size " << l << std::endl;
          }
         
        }
      }
     
#ifdef PARANOID
     // Check error
     {
      double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c_before_melt)-dc_new_melt);
      max_error_stage3=std::max(max_error_stage3,error);
      //std::cout << "Error stage 3: " << error << std::endl;
     }
#endif
    }
   
   
   // Doc stuff
   trace_file << time << " " //1
              << temp << " " //2
              << l_new << " " //3
              << GlobalParameters::G_cryst_rate(temp) << " "  //4
              << dc_new_cryst << " " //5
              << dc_growth << " " //6
              << dc_new_melt << " " //7
              << GlobalParameters::volume_of_crystalline_phase()<< " " //8
              << std::endl;

   
   // Size distribution and temperature
   {

    // current temperature
    std::ofstream temperature_file;
    std::string filename="RESLT/temperature"+std::__cxx11::to_string(i)+".dat";
    temperature_file.open(filename.c_str());
    temperature_file << temp << std::endl;
    temperature_file.close();

    // Vertical line at current melt/cryst size
    std::ofstream melt_temperature_file;
    filename="RESLT/melt_temperature"+std::__cxx11::to_string(i)+".dat";
    melt_temperature_file.open(filename.c_str());
    double l_melt=GlobalParameters::size_of_new_crystal(temp);
    melt_temperature_file << l_melt << " " << 0.0 << std::endl;
    melt_temperature_file << l_melt << " "
                          << GlobalParameters::T_infinite_cryst << std::endl; 
    melt_temperature_file.close();

    // Plot the size distribution
    std::ofstream size_distribution_file;
    filename="RESLT/size_distribution"+std::__cxx11::to_string(i)+".dat";
    size_distribution_file.open(filename.c_str());
    
    unsigned n=
     GlobalParameters::Number_of_crystals_initially_created_at_timestep.size();
    for (unsigned j=0;j<n;j++)
     {
      //if (GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]>0.0)
       {
        size_distribution_file
         <<  GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j] << " "
         << GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
         << " "
         << std::endl;
       }
     }
    size_distribution_file.close();
    
   }
   
  }

 #ifdef PARANOID
 std::cout << "Max. errors: "
           << max_error_stage1 << " "
           << max_error_stage2 << " "
           << max_error_stage3 << " "
           << std::endl;
#endif

 
}
