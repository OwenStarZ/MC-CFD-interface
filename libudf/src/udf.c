#include "udf.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

extern void ReadRMCPowerTally(real b[], real p[]);
extern void WriteFuelData(real data[]);
extern void WriteModeratorData(real data[]);
extern void WriteCoolantData(real data1[], real data2[]);

/*******************change this dims for reading power*******************/

#define DIM0 8
#define DIM1 8
#define DIM2 100

/*******************change this for output*******************/
/*******************also change this for multi-output meshes!*******************/

#define OUTDIMF0 15
#define OUTDIMF1 15
#define OUTDIMF2 40
//Output for Fuel temp

#define OUTDIMC0 15
#define OUTDIMC1 15
#define OUTDIMC2 40
//Output for Coolant temp and dens

#define OUTDIMM0 15
#define OUTDIMM1 15
#define OUTDIMM2 40
//Output for Moderator temp

real rmc_bndry[6];
real rmc_power[DIM0*DIM1*DIM2];
int bin_x = DIM0;
int bin_y = DIM1;
int bin_z = DIM2;
real fluent_bndry_x[2];
real fluent_bndry_y[2];
real fluent_bndry_z[2];
real hx, hy, hz;

/*******************no need to change functions！*******************/

/*
// Binary search
// Get mesh index by real position, suitable for meshtype2
int GetIndexByPos(int bin_num, real* nums, real pos)
{
	int index_left = 0;
	int index_right = bin_num; //bin_num = sizeof(nums) - 1
	
	if ( pos < nums[0] ) return 0;
	
	if ( pos > nums[bin_num] ) return bin_num;
	
	while( index_left <= index_right ){
		int mid = ( index_left + index_right ) / 2;
		if ( nums[mid] == pos) {
			return mid;
		}else if ( nums[mid] < pos ) {
			index_left = mid + 1;
		}else if ( nums[mid] > pos ) {
			index_right = mid - 1;
		}
	}
	return index_right;
}*/

//suitable for meshtype1
int GetIndexByPos(real h, real xmin, real pos)
{
	int index;
	index = (int)((pos-xmin)/h);
	return index;
}

DEFINE_INIT(read_power, d)
{
	// On all nodes

	real iwork[DIM0*DIM1*DIM2];
    real total = 0.0;
	real move_ori_x = 0.;
	real move_ori_y = 0.;
	real move_ori_z = 0.; //unit in cm
	int ix, iy, iz, i;

	for (i = 0; i < DIM0*DIM1*DIM2; i++){
		rmc_power[i]=0.0;
		iwork[i]=0.0;
	}
	for (i = 0; i < 6; i++){
		rmc_bndry[i] = 0.0;
	}

	#if RP_HOST
		int max_iter = 1;
		real power_all = 0.;
		real frac_direct_heat = 0.;
		
		// Read from coupling.dat
		Message("\n Read from coupling.dat \n");
		FILE *fp = fopen ("coupling.dat", "r");
		if(!fp)exit(0);
		
		fscanf(fp, "%d", &max_iter); // max iteration number
		fscanf(fp, "%lf", &power_all); // total power of the problem
		fscanf(fp, "%lf", &frac_direct_heat); // fraction of direct heat into the coolant, %
		fscanf(fp, "%lf", &move_ori_x);
		fscanf(fp, "%lf", &move_ori_y);
		fscanf(fp, "%lf", &move_ori_z);
		fclose(fp);
		
		Message("\n start reading power... \n");
		ReadRMCPowerTally(rmc_bndry, rmc_power);
		Message("\n power transfered \n");

		for (iz = 0; iz < DIM2; iz++){
            for (iy = 0; iy < DIM1; iy++){
                for(ix = 0; ix < DIM0; ix++){
                    rmc_power[ix+DIM0*iy+DIM0*DIM1*iz] *= 1.e6 * power_all;
                //power_all in MW，not power density but power
                }
            }
		}
		
	#endif//RP_HOST
		
	// Send from host to node
	host_to_node_real(rmc_bndry, 6);
	host_to_node_real(rmc_power, DIM0*DIM1*DIM2);
	host_to_node_real_3(move_ori_x, move_ori_y, move_ori_z);

    #if !RP_HOST

	real scale_factor = 100.0;
	
	// convert from rmc mesh boundary to fluent mesh boundary
	for (i = 0; i < 2; i++){
		fluent_bndry_x[i] = rmc_bndry[3 * i] + move_ori_x;
		fluent_bndry_x[i] /= scale_factor;
	}
	for (i = 0; i < 2; i++){
		fluent_bndry_y[i] = rmc_bndry[3 * i + 1] + move_ori_y;
		fluent_bndry_y[i] /= scale_factor;
	}
	for (i = 0; i < 2; i++){
		fluent_bndry_z[i] = rmc_bndry[3 * i + 2] + move_ori_z;
		fluent_bndry_z[i] /= scale_factor;
	}

    hx = (fluent_bndry_x[1] - fluent_bndry_x[0])/bin_x;
    hy = (fluent_bndry_y[1] - fluent_bndry_y[0])/bin_y;
    hz = (fluent_bndry_z[1] - fluent_bndry_z[0])/bin_z;

	cell_t c;
	real x[ND_ND];
	Thread *t;
	int ID_fuel[1] = {7}; /******************* Change this !*******************/
	real vol = 0.0;
	
	real vol_fuel[DIM0*DIM1*DIM2];
	for (i = 0; i < DIM0*DIM1*DIM2; i++){
		vol_fuel[i]=0.0;
	}

	/******************* Change this !*******************/
	for (i = 0; i < 1; i++){
        t = Lookup_Thread(d, ID_fuel[i]);
		begin_c_loop_int(c, t)
		{
            C_CENTROID(x, c, t);
            ix = GetIndexByPos(hx, fluent_bndry_x[0], x[0]);
            iy = GetIndexByPos(hy, fluent_bndry_y[0], x[1]);
            iz = GetIndexByPos(hz, fluent_bndry_z[0], x[2]);
			vol = C_VOLUME(c, t);
		    vol_fuel[ix+DIM0*iy+DIM0*DIM1*iz] += vol;
		}
		end_c_loop_int(c, t)
    }
	
	PRF_GRSUM(vol_fuel, DIM0*DIM1*DIM2, iwork);

	for (iz=0;iz<DIM2;iz++){
		for(iy=0;iy<DIM1;iy++){
			for(ix=0;ix<DIM0;ix++){
				if (vol_fuel[ix+DIM0*iy+DIM0*DIM1*iz]>1.e-14){
					rmc_power[ix+DIM0*iy+DIM0*DIM1*iz]/=vol_fuel[ix+DIM0*iy+DIM0*DIM1*iz];
				}
			}
		}
	}

	#endif
}

DEFINE_SOURCE(fuel_power, cell, thread, dS, eqn)
{
    real source;
	int ix, iy, iz;
	real x[ND_ND];
    C_CENTROID(x, cell, thread);
    ix = GetIndexByPos(hx, fluent_bndry_x[0], x[0]);
    iy = GetIndexByPos(hy, fluent_bndry_y[0], x[1]);
    iz = GetIndexByPos(hz, fluent_bndry_z[0], x[2]);
    source = rmc_power[ix+DIM0*iy+DIM0*DIM1*iz];
	dS[eqn] = 0.;
	return source;
}

/*DEFINE_SOURCE(direct_heat, cell, thread, dS, eqn)
{
	int max_iter = 1;
	real power_all = 0.;
	real frac_direct_heat = 0.;
	real vol_coolant_all = 32.73; // in[cm3]
	FILE *fp = fopen ("coupling.dat", "r");
	if(!fp)exit(0);
	fscanf(fp, "%d", &max_iter); // max iteration number
	fscanf(fp, "%lf", &power_all); // total power of the problem
	fscanf(fp, "%lf", &frac_direct_heat); // fraction of direct heat into the coolant, %
	fclose(fp);
	real src_direct_heat = power_all * 1.e6  * frac_direct_heat/100.0 / vol_coolant_all;
	real source = src_direct_heat;
	dS[eqn] = 0.;
	return source;
}*/

DEFINE_EXECUTE_AT_END(cal_th)
{
    int ix, iy, iz;

	/*******************Change the dimensions for multi-output meshes!*******************/

	real cnt_vol_fuel[OUTDIMF0*OUTDIMF1*OUTDIMF2];
	real cnt_vol_fluid[OUTDIMC0*OUTDIMC1*OUTDIMC2];
	real cnt_vol_moderator[OUTDIMM0*OUTDIMM1*OUTDIMM2];
	real cnt_temp_fuel[OUTDIMF0*OUTDIMF1*OUTDIMF2];
	real cnt_temp_fluid[OUTDIMC0*OUTDIMC1*OUTDIMC2];
	real cnt_temp_moderator[OUTDIMM0*OUTDIMM1*OUTDIMM2];
	real cnt_dens_fluid[OUTDIMC0*OUTDIMC1*OUTDIMC2];
	real F_iwork[OUTDIMF0*OUTDIMF1*OUTDIMF2];
	real C_iwork[OUTDIMC0*OUTDIMC1*OUTDIMC2];
	real M_iwork[OUTDIMM0*OUTDIMM1*OUTDIMM2];
	
	#if RP_HOST // RP_HOST
		// OUTPUT
		/*******************Change the dimensions for multi-output meshes!*******************/
		real temp_fuel[OUTDIMF0*OUTDIMF1*OUTDIMF2];
		real temp_fluid[OUTDIMC0*OUTDIMC1*OUTDIMC2];
		real temp_moderator[OUTDIMM0*OUTDIMM1*OUTDIMM2];
		real dens_fluid[OUTDIMC0*OUTDIMC1*OUTDIMC2];
		real foo; 
		// unit in K, kg/m^3
		real default_temp_fuel;
		real default_temp_moderator;
		real default_temp_fluid;
		real default_dens_fluid;
		// Read from coupling.dat
		FILE *fp = fopen ("coupling.dat", "r");
		if(!fp)exit(0);
		fscanf(fp, "%lf", &foo); // &max_iter
		fscanf(fp, "%lf", &foo); // &power_all
		fscanf(fp, "%lf", &foo); // &frac_direct_heat
		fscanf(fp, "%lf", &foo); // &move_ori_x
		fscanf(fp, "%lf", &foo); // &move_ori_y
		fscanf(fp, "%lf", &foo); // &move_ori_z
		fscanf(fp, "%lf", &default_temp_fuel);
		fscanf(fp, "%lf", &default_temp_moderator);
		fscanf(fp, "%lf", &default_temp_fluid);
		fscanf(fp, "%lf", &default_dens_fluid);
		fclose(fp);
	#endif // RP_HOST
	
	/*******************Change the dimensions for multi-output meshes!*******************/

	real size_F_mesh = OUTDIMF0*OUTDIMF1*OUTDIMF2;
	real size_C_mesh = OUTDIMC0*OUTDIMC1*OUTDIMC2;
	real size_M_mesh = OUTDIMM0*OUTDIMM1*OUTDIMM2;
	real zero = 1.e-10;
	
	// initialize
	for (ix = 0; ix < size_F_mesh; ix++){
		cnt_vol_fuel[ix] = 0.0;
		cnt_temp_fuel[ix] = 0.0;
		F_iwork[ix] = 0.0;
	}
	for (ix = 0; ix < size_C_mesh; ix++){
		cnt_vol_fluid[ix] = 0.0;
		cnt_temp_fluid[ix] = 0.0;
		cnt_dens_fluid[ix] = 0.0;
		C_iwork[ix] = 0.0;
	}
	for (ix = 0; ix < size_M_mesh; ix++){
		cnt_vol_moderator[ix] = 0.0;
		cnt_temp_moderator[ix] = 0.0;
		M_iwork[ix] = 0.0;
	}
	
	#if !RP_HOST		
		
	/*******************Change the dimensions for multi-output meshes!*******************/

        int i;
		int F_bin_x = OUTDIMF0;
		int F_bin_y = OUTDIMF1;
		int F_bin_z = OUTDIMF2;
		int C_bin_x = OUTDIMC0;
		int C_bin_y = OUTDIMC1;
		int C_bin_z = OUTDIMC2;
		int M_bin_x = OUTDIMM0;
		int M_bin_y = OUTDIMM1;
		int M_bin_z = OUTDIMM2;
		real F_bndry_x[2] = {-0.0075, 0.0075};
		real F_bndry_y[2] = {-0.0075, 0.0075};
		real F_bndry_z[2] = {0., 0.2};
		real C_bndry_x[2] = {-0.0075, 0.0075};
		real C_bndry_y[2] = {-0.0075, 0.0075};
		real C_bndry_z[2] = {0., 0.2};
		real M_bndry_x[2] = {-0.0075, 0.0075};
		real M_bndry_y[2] = {-0.0075, 0.0075};
		real M_bndry_z[2] = {0., 0.2};
        real F_hx = (F_bndry_x[1] - F_bndry_x[0]) / F_bin_x;
        real F_hy = (F_bndry_y[1] - F_bndry_y[0]) / F_bin_y;
        real F_hz = (F_bndry_z[1] - F_bndry_z[0]) / F_bin_z;
        real C_hx = (C_bndry_x[1] - C_bndry_x[0]) / C_bin_x;
        real C_hy = (C_bndry_y[1] - C_bndry_y[0]) / C_bin_y;
        real C_hz = (C_bndry_z[1] - C_bndry_z[0]) / C_bin_z;
        real M_hx = (M_bndry_x[1] - M_bndry_x[0]) / M_bin_x;
        real M_hy = (M_bndry_y[1] - M_bndry_y[0]) / M_bin_y;
        real M_hz = (M_bndry_z[1] - M_bndry_z[0]) / M_bin_z;

		// Fuel cell 
        Domain *d = Get_Domain(1);
        cell_t c;
        Thread *t;
        int ID_fuel[1] = {7}; /******************* Change this !*******************/
        int ID_fluid[1] = {6}; /******************* Change this !*******************/
        int ID_moderator[1] = {8}; /******************* Change this !*******************/
        real x[ND_ND];
        real vol;
        /******************* Change this !*******************/
        for (i = 0; i < 1; i++){
            t = Lookup_Thread(d, ID_fuel[i]);
		    begin_c_loop_int(c, t)
		    {
                C_CENTROID(x, c, t);
                ix = GetIndexByPos(F_hx, F_bndry_x[0], x[0]);
                iy = GetIndexByPos(F_hy, F_bndry_y[0], x[1]);
                iz = GetIndexByPos(F_hz, F_bndry_z[0], x[2]);
			    vol = C_VOLUME(c, t);
			    cnt_vol_fuel[ix+OUTDIMF0*iy+OUTDIMF0*OUTDIMF1*iz] += vol;
			    cnt_temp_fuel[ix+OUTDIMF0*iy+OUTDIMF0*OUTDIMF1*iz] += C_T(c, t) * vol;
		    }
		    end_c_loop_int(c, t)
        }
        
		// Fluid cell 
		/******************* Change this !*******************/
		for (i = 0; i < 1; i++){
			t = Lookup_Thread(d, ID_fluid[i]);
			begin_c_loop_int(c, t)
			{	
            	C_CENTROID(x, c, t);
           		ix = GetIndexByPos(C_hx, C_bndry_x[0], x[0]);
            	iy = GetIndexByPos(C_hy, C_bndry_y[0], x[1]);
            	iz = GetIndexByPos(C_hz, C_bndry_z[0], x[2]);
				vol = C_VOLUME(c, t);
				cnt_vol_fluid[ix+OUTDIMC0*iy+OUTDIMC0*OUTDIMC1*iz] += vol;
				cnt_temp_fluid[ix+OUTDIMC0*iy+OUTDIMC0*OUTDIMC1*iz] += C_T(c, t) * vol;
				cnt_dens_fluid[ix+OUTDIMC0*iy+OUTDIMC0*OUTDIMC1*iz] += C_R(c, t) * vol;
			}
			end_c_loop_int(c, t)
		}

		//moderator cell 
		/******************* Change this !*******************/
		for (i = 0; i < 1 ; i++){
			t = Lookup_Thread(d, ID_moderator[i]);
			begin_c_loop_int(c, t)
			{
        		C_CENTROID(x, c, t);
        		ix = GetIndexByPos(M_hx, M_bndry_x[0], x[0]);
        		iy = GetIndexByPos(M_hy, M_bndry_y[0], x[1]);
        		iz = GetIndexByPos(M_hz, M_bndry_z[0], x[2]);
				vol = C_VOLUME(c, t);
				cnt_vol_moderator[ix+OUTDIMM0*iy+OUTDIMM0*OUTDIMM1*iz] += vol;
				cnt_temp_moderator[ix+OUTDIMM0*iy+OUTDIMM0*OUTDIMM1*iz] += C_T(c, t) * vol;
			}
			end_c_loop_int(c, t_moderator)
		}
				
	#endif // !RP_HOST
	
	// Global summation and send data from node_0 to host 

	PRF_GRSUM(cnt_vol_fuel, size_F_mesh, F_iwork);
	PRF_GRSUM(cnt_vol_fluid, size_C_mesh, C_iwork);
	PRF_GRSUM(cnt_vol_moderator, size_M_mesh, M_iwork);
	PRF_GRSUM(cnt_temp_fuel, size_F_mesh, F_iwork);
	PRF_GRSUM(cnt_temp_fluid, size_C_mesh, C_iwork);
	PRF_GRSUM(cnt_temp_moderator, size_M_mesh, M_iwork);
	PRF_GRSUM(cnt_dens_fluid, size_C_mesh, C_iwork);
	
	node_to_host_real(cnt_vol_fuel, size_F_mesh);
	node_to_host_real(cnt_vol_fluid, size_C_mesh);
	node_to_host_real(cnt_vol_moderator, size_M_mesh);
	node_to_host_real(cnt_temp_fuel, size_F_mesh);
	node_to_host_real(cnt_temp_fluid, size_C_mesh);
	node_to_host_real(cnt_temp_moderator, size_M_mesh);
	node_to_host_real(cnt_dens_fluid, size_C_mesh);
	
	#if RP_HOST // RP_HOST

    for (ix = 0; ix < OUTDIMF0*OUTDIMF1*OUTDIMF2; ix++){
        if (cnt_vol_fuel[ix] < zero) {
            temp_fuel[ix] = default_temp_fuel;
        }else{
			temp_fuel[ix] = cnt_temp_fuel[ix] / cnt_vol_fuel[ix];
		}
	}
		// Fluid temp and dens
	for (ix = 0; ix < OUTDIMC0*OUTDIMC1*OUTDIMC2; ix++){
		if (cnt_vol_fluid[ix] < zero) {
			temp_fluid[ix] = default_temp_fluid;
			dens_fluid[ix] = default_dens_fluid;
		}else{
			temp_fluid[ix] = cnt_temp_fluid[ix] / cnt_vol_fluid[ix];
			dens_fluid[ix] = cnt_dens_fluid[ix] / cnt_vol_fluid[ix];
		}
	}	
		// moderator temp
	for (ix = 0; ix < OUTDIMM0*OUTDIMM1*OUTDIMM2; ix++){
		if (cnt_vol_moderator[ix] < zero) {
			temp_moderator[ix] = default_temp_moderator;
		}else{
			temp_moderator[ix] = cnt_temp_moderator[ix] / cnt_vol_moderator[ix];
		}
    }
		
		// Write data to h5Files 
		WriteFuelData(temp_fuel);
		WriteModeratorData(temp_moderator);
		WriteCoolantData(temp_fluid, dens_fluid);
		
	#endif // RP_HOST
	
	#if RP_HOST // RP_HOST
		Message("MACRO DEFINE_EXECUTE_AT_END DONE! FROM HOST \n");
	#endif // RP_HOST
	
}
