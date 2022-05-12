#include<iostream>
#include<vector>
#include<map>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include <cfloat>

#include<fenv.h>


#define PARANOID

//=============================================================
/// Namespace for parameters
//=============================================================
namespace GlobalParameters
{

 // Primary crystallisation
 //------------------------

 /// Magnitude of crystallation rate fct
 double G_0_cryst=0.6055; // min^-1 // 1.0; 

 /// Temperature at which crystallisation rate is max.
 double T_max_cryst=383.902; // K // 5.0;

 /// Width of Gaussian for crystallisation rate
 double S_cryst=4.39142; // K // 0.1;

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
 double L_c_0=4.081e-10; // m // 0.1;

 /// Temperature at which crystals grow to infinite size
 /// (or the temperature at which even infinitely large crystals melt)
 double T_infinite_cryst=415.873; // K // 12.0;
 
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



 /// Melt instantaneously or by modifying the growth rate of crystals
 bool Melt_instantaneously=false;


 // hierher park this for now
 // /// Magnitude of melt rate
 // double M_hat=100.0;

 // /// Sharpness of continous melt rate
 // double Alpha_melt=0.01;
 
 // /// Continuous melting: Melt rate
 // double melt_rate(const double& l, const double& temp)
 // {
 //  double t_melt=melt_temperature(l);
 //  return M_hat*0.5*(1.0+tanh((temp-t_melt)/Alpha_melt));
 // }

 
 // Secondary crystallisation
 //--------------------------

 /// Scaling parameter for growth rate
 double K_growth_1=6.49e-2; // min^{-1} // used to be constant: 1.0; 

 /// Activation energy for secondary crystallisation
 double E_activation=48.7353; // kJ mol^{-1}

 /// Gas constant
 double R_gas_constant=8.31e-3; // kH mol^{-1} K^{-1}
 
 /// Exponential parameter for growth rate
 double K_growth_2=9.11e-14; // m J^{-1} 

 /// Boltzmann constant
 double K_boltzmann=1.38e-23; // J K^{-1} // 1.0;
 
 /// Secondary growth rate (dL/dt) as a function of the current crystal
 /// length, l, the temperature, temp, and the current value of the crystalline
 /// volume fraction, assumed to be related to melt via c = 1-m
 double secondary_growth_rate(const double& l,
                              const double& temp,
                              const double& c)
 {
  double m=1.0-c;
  double k_growth_1=K_growth_1*exp(-E_activation/(R_gas_constant*temp));

  // hierher overflows were probably an artefact of underflow/negative lengths
  // double numer=K_growth_2*l;
  // double denom=K_boltzmann*temp;
  // double fract=std::max(1.0e1,numer/denom);
  // double exp_thing=exp(-fract);
  // return m*k_growth_1*exp_thing;
  return m*k_growth_1*exp(-(K_growth_2*l/(K_boltzmann*temp)));
 }
 

 // Temperature variation
 //----------------------

 // Initial temperature at t=0
 double T_init=2.73e2; //  K // 1.2*T_max_cryst;

 // End temperature at t=1
 double T_end=4.63e2; // K // 0.8*T_max_cryst;

 // Rate of temperature change
 double Rate_of_temperature_change=10.0; // K min^{-1}

 /// current temperature
 double current_temperature(const double& time)
  {
   double t_switch_over=(T_end-T_init)/Rate_of_temperature_change;
   double temperature_threshold=T_infinite_cryst-20.0; // hierher this is what phil suggested -0.1;
   if (time< t_switch_over)
    {
     return std::min(temperature_threshold,
                     T_end-Rate_of_temperature_change*time);
    }
   else
    {
     return std::min(temperature_threshold,
                     T_init+Rate_of_temperature_change*(time-t_switch_over));
    }

  // double t_switch_over=3.0;
  // if (time<t_switch_over)
  //  {
  //   return std::max(T_end,T_init+(T_end-T_init)*time);
  //  }
  // else
  //  {
  //   return std::min(0.9*T_infinite_cryst,
  //                   T_end-(T_end-T_init)*(time-t_switch_over)); 
  //  }

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


 feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Number of timesteps
 unsigned nstep=1000; // 10000;

 // Max. time: Once down once up
 double t_max=2.0*(GlobalParameters::T_end-GlobalParameters::T_init)/
  GlobalParameters::Rate_of_temperature_change; //10.0;
 
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
             << GlobalParameters::T_max_cryst << " " 
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


 // Limits for size distribution
 double l_min=DBL_MAX;
 double l_max=-DBL_MAX;
 double n_max=-DBL_MAX;

 
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

   
   // Incremental volume fraction that is currently being crystallised
   // into crystals of (scalar) size l_new over current timestep;
   // First order explicit Euler scheme.
   double dc_new_cryst = GlobalParameters::G_cryst_rate(temp)*(1.0-c)*dt;

   // How many crystals fit into this incremental volume? (double prec)
   double threshold_for_number=1.0e-20; // hierher
   double n_new_cryst=dc_new_cryst/GlobalParameters::crystal_volume(l_new);
   if (n_new_cryst<threshold_for_number)
    {
     n_new_cryst=0.0;
     l_new=0.0;
    }
   
   // Push back
   GlobalParameters::Current_size_of_crystals_initially_created_at_timestep.
    push_back(l_new);
   GlobalParameters::Number_of_crystals_initially_created_at_timestep.
    push_back(n_new_cryst);
  
    
#ifdef PARANOID
   // Check error
   {
    double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c)-dc_new_cryst);
    max_error_stage1=std::max(max_error_stage1,error);
   }
#endif
  
   

   // STATE 2: GROW THE CRYSTALS THAT WERE CREATED AT PREVIOUS TIMES (FIRST
   //----------------------------------------------------------------------
   // ORDER EULER SCHEME; POSSIBLY COMBINED WITH MELTING
   //---------------------------------------------------

   // Change in crystalline volume fraction due to secondary growth and melting
   double dc_melt=0.0;
   double dc_growth=0.0;
   
#ifdef PARANOID
   double c_before_growth=GlobalParameters::volume_of_crystalline_phase();
#endif
   
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

       // Did we melt during this timestep?
       bool have_melted=false;
       
       // Get (linear) growth rate: Forward Euler so base c on crystalline phase
       // at previous timestep
       double dL_dt=
        GlobalParameters::secondary_growth_rate(l,temp,c);

       // Combine/override with melt?
       if (!GlobalParameters::Melt_instantaneously)
        {
         double t_melt=GlobalParameters::melt_temperature(l);
         // hierher pretty ad hoc
         double Alpha_exponential_melt_rate=1.0;
         if (temp>t_melt)
          {
           have_melted=true;
           dL_dt=-Alpha_exponential_melt_rate*l;
          }
        }
       
       // Bump (first-order Euler scheme)
       double dL=dL_dt*dt;

       // Stop lengths from becoming negative!
       double l_threshold=1.0e-30;
       if (l+dL<l_threshold)
        {
         //std::cout << "rescuing " << l << " " << dL << std::endl;
         dL=-l;
         GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]=0.0;
         GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
            =0.0;
        }
       else
        {
         GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]
          +=dL;
        }
       
       double v_old=num_of_cryst*GlobalParameters::crystal_volume(l);
       double v_new=num_of_cryst*GlobalParameters::crystal_volume(l+dL);
       double dc=(v_new-v_old);

       if (have_melted)
        {
         dc_melt+=dc;
        }
       else
        {
         dc_growth+=dc;
        }
      }
    }

#ifdef PARANOID
   // Check error
   {
    double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c_before_growth)
                      -(dc_growth+dc_melt));
    max_error_stage2=std::max(max_error_stage2,error);
   }
#endif


   
   // STAGE 3: MELT STUFF INSTANTANEOUSLY
   //------------------------------------
   if (GlobalParameters::Melt_instantaneously)
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
           dc_melt-=num_of_cryst*GlobalParameters::crystal_volume(l);
           
           // Kill 'em (sizes are set to zero too so they don't clunk up the animations
           // and 0 x 0 is still zero!
           GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]
            =0.0;
           GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]
            =0.0;
          }         
        }
      }
     
#ifdef PARANOID
     // Check error
     {
      double error=fabs((GlobalParameters::volume_of_crystalline_phase()-c_before_melt)-dc_melt);
      max_error_stage3=std::max(max_error_stage3,error);
     }
#endif
    }
   
   
   // Doc stuff
   trace_file << time << " " //1
              << temp << " " //2
              << l_new << " " //3
              << GlobalParameters::G_cryst_rate(temp) << " "  //4
              << dc_new_cryst/dt << " " //5
              << dc_growth/dt << " " //6
              << dc_melt/dt << " " //7
              << GlobalParameters::volume_of_crystalline_phase()<< " " //8
              << (dc_new_cryst+dc_growth+dc_melt)/dt << " " // 9
              << std::endl;

   
   // Size distribution and temperature
   {

    // current temperature
    std::ofstream temperature_file;
    std::string filename="RESLT/temperature"+std::__cxx11::to_string(i)+".dat";
    temperature_file.open(filename.c_str());
    temperature_file << temp << std::endl;
    temperature_file.close();

    // Vertical line at current time
    std::ofstream time_file;
    filename="RESLT/time"+std::__cxx11::to_string(i)+".dat";
    time_file.open(filename.c_str());
    time_file << time << " " << 0.0 << std::endl;
    time_file << time << " " << 1.0 << std::endl;
    time_file.close();

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

        if (GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]<l_min)
         {
          l_min=GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j];
         }
        if (GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j]>l_max)
         {
          l_max=GlobalParameters::Current_size_of_crystals_initially_created_at_timestep[j];
         }
        if (GlobalParameters::Number_of_crystals_initially_created_at_timestep[j]>n_max)
         {
          n_max=GlobalParameters::Number_of_crystals_initially_created_at_timestep[j];
         }

         
       }
     }
    size_distribution_file.close();
    
   }
   
  }

 
 // Limits for size distribution
 std::ofstream lim_file;
 lim_file.open("RESLT/size_distribution_limits.dat");
 lim_file << l_min << " "
          << l_max << " "
          << n_max << " "
          << std::endl;
 lim_file.close();
 
 #ifdef PARANOID
 std::cout << "Max. errors: "
           << max_error_stage1 << " "
           << max_error_stage2 << " "
           << max_error_stage3 << " "
           << std::endl;
#endif

 
}
