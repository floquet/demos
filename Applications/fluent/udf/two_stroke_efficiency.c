/*********************************************************************************************
                   UDF for two stroke engine scavenging efficiency calculation

		   For 2-stroke engines, scaverning is critical since there is no separate stroke
	   to drive the burnt gas out of the combustion chamber.  This udf provides a tool to
	   calcuate different scaverging efficiency.  This udf calculates the following:

	   - mass trapping efficiency, TEm = fresh charge mass trapped in the cylinder at EPC
	     divided by mass of fresh charge delivered to the cylinder.
	   - volume scavenging ratio, SRv = volume of fresh charge supplied to the cylinder
	     at EPC divided by volume of the cylinder at BDC.
	   - volume scavenging efficiency, SEv = volume of fresh charge trapped in the cylinder
	     at EPC divided by volume of the cylinder at BDC.

	   Note that the above volume is calculated by the corresponding mass divided by the
	     average density inside the cylinder at EPC.

	   Mesh and setup requirement:
	   - Chamber, tranfer port including boost port, and exhaust port will each need
	     a separate cell zone
	   - Fresh charge has to be the first species in the mixture.
	   - Can not restart in Windows.  The calculation involves time integral.  In order to
	     restart, the value has to be saved in the data file.  Windows does not support
	     that functionality.

	   How to use the udf:
	   - Set up your IC case
	   - Modify the user inputs part of the udf.
	   - Build the library
	   - Hook the DEFINE_EXECUTE_AT_END, DEFINE_RW_FILE udf
	   - Initialize your flow field, and patch the correct the values
	   - Solve the flow
	   - DEFINE_ON_DEMAND gives the results

       Note:
	   - UDF works in 2d and 3d.
	   - UDF works in both serial and parallel.

	   Written by : Xiao HU (Fluent Inc)
	   Last updated : 7/24/2006
***********************************************************************************************/

#include "udf.h"

# define RPM RP_Get_Real("dynamesh/in-cyn/crank-rpm")

#if RP_DOUBLE
#  define REAL_FMT "%le"
#else  /* RP_DOUBLE */
#  define REAL_FMT "%e"
#endif /* RP_DOUBLE */

# define ADVANCED 0

/********************************* User input starts *****************************************/

static int Cylinder_Inlet_ID[]={31, -1};  /* from (rpgetvar 'sliding-interfaces) */
static int Chamber_Zone_ID[]= {6, 3, -1};
static real EPC_CA = 260;

/********************************** User input ends ******************************************/

/********************************* Optional/advanced inputs **********************************

     WARNING: using this option only when you are absoluately sure about it.  Otherwise, it may
     have detrimental impact on your results.

     These inputs are optional only for confirmation purpose.  To calcuate the delivered mass
  into cylinder using the interior face zone from interface has a very tiny numerical error.
  This is due to the fact that flux terms are not saved for interior face zones.  So, to calculate
  the flux, an averaged value from adjacent cells are used.  This will cause a tiny numerical error.
  To use the following approach to calcuate the delivered mass, on the other hand, does not
  calulate the delivered mass at the interface.  Instead, it calcuates the delivered mass at
  the inlet and the accumulation at the transport ports and boost ports.  Both numbers are from
  Fluent and thus are true values from Fluent as opposed to interpolated values.  But to use
  this option is more tedious and thus error-prone.  So, use it only if you are very comfortable
  with the inputs.  */

/* To use this option, you will first need to change the ADVANCED from 0 to 1 above !  Note that
  if you already have a data file, switching the flag can cause trouble because one additional
  variable is written in the data file.  You should have no problem if you define this before
  saving the data file.  Again, use this option with caution. */

#if ADVANCED

static int Inlet_Zone_ID[] = {12, -1};
static int Tranfer_port_ID[]= {6, 2, -1};
static real Init_Charge_in_Transfer_Port = 5.030717e-4;

#endif

/********************************** Advanced input ends ******************************************/

static real BDC_CA = 180;
static real charge_delivered_to_cyl=0, cyl_vol_BTC, density_cyl_EPC;
static int counter_EPC=0, counter_BDC=0;

#if ADVANCED

static real charge_delivered_to_inlet=0;

#endif

/* Given a cell zone IC, calculates volume */

static real f_volume(int * ID)
{
 real volume=0;
 Thread *tc;
 cell_t c;
 Domain* domain;
 int i;

 domain=Get_Domain(1);
 i=0;
 while(ID[i]>=0)
  {
   tc=Lookup_Thread(domain, ID[i]);

   begin_c_loop_int(c, tc)
    {
      volume += C_VOLUME(c,tc);
    }
   end_c_loop_int(c, tc)

   i++;
  }

 volume = PRF_GRSUM1(volume);

 return volume;
}

/* Given a cell zone IC, calculates trapped charge mass */

static real trapped_charge_mass(int * ID)
{
 real trapped_charge=0;
 Thread *tc;
 cell_t c;
 Domain* domain;
 int i;

 domain=Get_Domain(1);

 i=0;
 while(ID[i]>=0)
  {
   tc=Lookup_Thread(domain, ID[i]);

   begin_c_loop_int(c, tc)
    {
      trapped_charge += C_R(c,tc) * C_VOLUME(c,tc) * C_YI(c, tc, 0);
    }
   end_c_loop_int(c, tc)

   i++;
  }

 trapped_charge = PRF_GRSUM1(trapped_charge);

 return trapped_charge;
}

/* Given a cell zone IC, calculates trapped total mass */

static real trapped_total_mass(int * ID)
{
 real total_mass=0;
 Thread *tc;
 cell_t c;
 Domain* domain;
 int i;

 domain=Get_Domain(1);

 i=0;
 while(ID[i]>=0)
  {
   tc=Lookup_Thread(domain, ID[i]);

   begin_c_loop_int(c, tc)
    {
      total_mass += C_R(c,tc) * C_VOLUME(c,tc);
    }
   end_c_loop_int(c, tc)

   i++;
  }

 total_mass = PRF_GRSUM1(total_mass);

 return total_mass;
}

DEFINE_ON_DEMAND(print_out_results)
{
 real charge_trapped_in_cyl;

#if ADVANCED

 real build_up_charge;

#endif

 charge_trapped_in_cyl=trapped_charge_mass(Chamber_Zone_ID);

#if ADVANCED

 build_up_charge = trapped_charge_mass(Tranfer_port_ID) - Init_Charge_in_Transfer_Port;

#endif

 Message0("\n\n************************* Results at CA = %6.2f deg *************************\n",
                  CURRENT_TIME*RPM*6.0+RP_Get_Real("dynamesh/in-cyn/crank-start-angle"));
 Message0("\n Fresh charge mass delivered to the cylinder (kg) = %11.4e", fabs(charge_delivered_to_cyl));
 Message0("\n Fresh charge mass trapped in the cylinder   (kg) = %11.4e", charge_trapped_in_cyl);

 if(counter_BDC>0)
    Message0("\n\n Cylinder volume at BDC (m^3)      = %11.4e", cyl_vol_BTC);
 if(counter_EPC>0)
    Message0("\n Charge density at EPC  (kg/m^3)   = %11.4e", density_cyl_EPC);

 Message0("\n\n Mass trapping efficiency, TEm     = %11.4e", charge_trapped_in_cyl/fabs(charge_delivered_to_cyl));

 if(counter_EPC>0)
   {
    Message0("\n Volume scavenging ratio, SRv      = %11.4e", fabs(charge_delivered_to_cyl)/density_cyl_EPC/cyl_vol_BTC);
    Message0("\n Volume scavenging efficiency, SEv = %11.4e", charge_trapped_in_cyl/density_cyl_EPC/cyl_vol_BTC);
   }

#if ADVANCED

 Message0("\n\n The following results are from the advaned option, which should be simlar to the above:\n");

 Message0("\n Fresh charge mass delivered to the inlet    (kg) = %11.4e", charge_delivered_to_inlet);
 Message0("\n Fresh charge mass initially in transfer port(kg) = %11.4e", Init_Charge_in_Transfer_Port);
 Message0("\n Fresh charge mass currently in transfer port(kg) = %11.4e", trapped_charge_mass(Tranfer_port_ID));
 Message0("\n Fresh charge mass built up in transfer port (kg) = %11.4e", build_up_charge);
 Message0("\n Fresh charge mass delivered to the cylinder (kg) = %11.4e", charge_delivered_to_inlet - build_up_charge);
 Message0("\n Fresh charge mass trapped in the cylinder   (kg) = %11.4e", charge_trapped_in_cyl);


 Message0("\n\n Mass trapping efficiency, TEm     = %11.4e", charge_trapped_in_cyl/(charge_delivered_to_inlet - build_up_charge));

 if(counter_EPC>0)
   {
    Message0("\n Volume scavenging ratio, SRv      = %11.4e", (charge_delivered_to_inlet - build_up_charge)/density_cyl_EPC/cyl_vol_BTC);
    Message0("\n Volume scavenging efficiency, SEv = %11.4e", charge_trapped_in_cyl/density_cyl_EPC/cyl_vol_BTC);
   }
#endif

 Message0("\n\n******************************************************************************\n");
}

DEFINE_EXECUTE_AT_END(cal_charge)
{
#if !RP_HOST

 real charge_one_dt, CA, yi;
 Thread *tf;
 face_t f;
 Domain* domain;
 int i;

 domain=Get_Domain(1);

 /* Calcualte fresh charge delivered to the cylinder */

 i=0;
 while(Cylinder_Inlet_ID[i]>=0)
  {
    tf=Lookup_Thread(domain, Cylinder_Inlet_ID[i]);

    charge_one_dt=0;
    begin_f_loop(f, tf)
     if (PRINCIPAL_FACE_P(f,tf))
      {
		if (THREAD_T1(tf) == NULL)
		    yi = F_YI(f, tf, 0);
		else
		    yi = (C_YI(F_C0(f, tf), THREAD_T0(tf), 0) +
		          C_YI(F_C1(f, tf), THREAD_T1(tf), 0))/2.;

	    charge_one_dt += CURRENT_TIMESTEP * (-F_FLUX(f, tf)) * yi;
      }
    end_f_loop(f, tf)

	i++;

    charge_one_dt = PRF_GRSUM1(charge_one_dt);

    charge_delivered_to_cyl += charge_one_dt;
   }

#if ADVANCED
 /* Calcualte fresh charge delivered to the inlet */

 i=0;
 while(Inlet_Zone_ID[i]>=0)
  {
    tf=Lookup_Thread(domain, Inlet_Zone_ID[i]);

    charge_one_dt=0;
    begin_f_loop(f, tf)
     if (PRINCIPAL_FACE_P(f,tf))
      {
		yi = F_YI(f, tf, 0);

	    charge_one_dt += CURRENT_TIMESTEP * (-F_FLUX(f, tf)) * yi;
      }
    end_f_loop(f, tf)

	i++;

    charge_one_dt = PRF_GRSUM1(charge_one_dt);

    charge_delivered_to_inlet += charge_one_dt;
   }
#endif

#endif

#if ADVANCED

 node_to_host_real_1(charge_delivered_to_inlet);

#endif

 node_to_host_real_1(charge_delivered_to_cyl);

#if !RP_HOST

 CA = CURRENT_TIME*RPM*6.0+RP_Get_Real("dynamesh/in-cyn/crank-start-angle");

 if((CA>=EPC_CA)&&(counter_EPC==0))
   {
	   density_cyl_EPC = trapped_total_mass(Chamber_Zone_ID)/f_volume(Chamber_Zone_ID);
	   counter_EPC++;
   }

 if((CA>=BDC_CA)&&(counter_BDC==0))
   {
	   cyl_vol_BTC = f_volume(Chamber_Zone_ID);
	   counter_BDC++;
   }

#endif

 node_to_host_real_2(density_cyl_EPC, cyl_vol_BTC);
 node_to_host_int_2(counter_EPC, counter_BDC);
}

DEFINE_RW_FILE(write_data, fp)
{
  Message0("\nWriting user defined data to the data file...\n");

 #if PARALLEL
   if(I_AM_NODE_HOST_P)
 #endif
    {
     fprintf(fp, "\n" REAL_FMT, charge_delivered_to_cyl);

   #if ADVANCED
    fprintf(fp, "\n" REAL_FMT, charge_delivered_to_inlet);
   #endif

     fprintf(fp, "\n" REAL_FMT, cyl_vol_BTC);
     fprintf(fp, "\n" REAL_FMT, density_cyl_EPC);
     fprintf(fp, "\n%d", counter_EPC);
     fprintf(fp, "\n%d", counter_BDC);
    }
}

DEFINE_RW_FILE(read_data, fp)
{
  Message0("\nReading user defined data from the data file...\n");

 #if PARALLEL
   if(I_AM_NODE_HOST_P)
 #endif
    {
     fscanf(fp, REAL_FMT, &charge_delivered_to_cyl);

   #if ADVANCED
     fscanf(fp, REAL_FMT, &charge_delivered_to_inlet);
   #endif

     fscanf(fp, REAL_FMT, &cyl_vol_BTC);
     fscanf(fp, REAL_FMT, &density_cyl_EPC);
     fscanf(fp, "%d", &counter_EPC);
     fscanf(fp, "%d", &counter_BDC);
    }

#if ADVANCED

 host_to_node_real_4(charge_delivered_to_cyl, charge_delivered_to_inlet, cyl_vol_BTC, density_cyl_EPC);

#else

 host_to_node_real_3(charge_delivered_to_cyl, cyl_vol_BTC, density_cyl_EPC);

#endif

 host_to_node_int_2(counter_EPC, counter_BDC);
}
