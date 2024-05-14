#include "stdio.h"
#include "string.h"
#include "hdf5.h"
#include "stdlib.h"
#include "math.h"

/******** change the dimensions, note that they should be same as those in UDF ********/
#define OUTDIMF0 8
#define OUTDIMF1 8
#define OUTDIMF2 100

#define OUTDIMC0 8
#define OUTDIMC1 8
#define OUTDIMC2 100


#define OUTDIMM0 1
#define OUTDIMM1 1
#define OUTDIMM2 100

double lambda1(double T, double R, double q)
{
    double knew, kold, k1d, k1p, k2p, k4r, k0, T0, Tmid;
	double tal = 0.001;
	double p = 0.05;
	knew = 5.0;
	do
	{
		kold = knew;
		T0 = T + q * R * R / 6 / kold;
		Tmid = (T0 + T) / 2;
		k1d = (1.09/pow(tal, 3.265)+0.0643/sqrt(tal)*sqrt(Tmid))*atan(1/(1.09/pow(tal, 3.265)+0.0643/sqrt(tal)*sqrt(Tmid)));
		k1p = 1+0.019*tal/(3-0.019*tal)*(1/(1+exp(12-0.01*Tmid)));
		k2p = (1-p)/(1+2*p);
		k4r = 1-(0.2/(1+exp((Tmid-900)/80)))*(1-exp(-1*tal));
		k0 = 1/(0.0375+2.165e-4*Tmid)+4.715e9/pow(Tmid, 2)*exp(-16361/Tmid);
		knew = k1d*k1p*k2p*k4r*k0;
	} while (fabs(knew - kold) > 0.000001);

	return knew;
}

double lambda4(double T, double R0, double R1, double R2, double q)
{
	double knew, kold, T0, Tmid;
	knew = 42.0;
	do
	{
		kold = knew;
		T0 = T + pow(R0, 3) * q / 3 / kold * (1/R1-1/R2);
		Tmid = (T0 + T) / 2;
		knew = 42.58-15564/Tmid+1.2977E7/pow(Tmid, 2)-1.8458E9/pow(Tmid, 3);
	} while (fabs(knew - kold) > 0.000001);
	
	return knew;
}

int main()
{
	FILE* in;

	hid_t file, space, group, attr, dataset;
	herr_t status;
	hsize_t dim_attr[1] = { 1 }, dim_BinNumber[1] = { 3 }, dim_Boundary[2] = { 3,2 };
	hsize_t dim_F_data[3] = { OUTDIMF0,OUTDIMF1,OUTDIMF2 };
	hsize_t dim_C_data[3] = { OUTDIMC0,OUTDIMC1,OUTDIMC2 };
	hsize_t dim_M_data[3] = { OUTDIMM0,OUTDIMM1,OUTDIMM2 };

	int a[1] = { 1 }; // 1 for uniform mesh, 2 for the other

	double F_geo_binnumber[3] = { OUTDIMF0,OUTDIMF1,OUTDIMF2 };
	double F_geo_boundary[3][2] = { {-0.75,0.75},{-0.75,0.75},{0,20} };         //  change this, unit in [cm]!!!
	double F_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double power_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double F1_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double F2_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double F3_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double F4_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double F5_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];
    double r1 = 0.00025, r2 = 0.000345, r3 = 0.000385, r4 = 0.00042, r5 = 0.00046;
    double T1, T2, T3, T4;
    double k1, k4;
    double k2 = 0.5, k3 = 4.0, k5 = 4.0;

	double C_geo_binnumber[3] = { OUTDIMC0,OUTDIMC1,OUTDIMC2 };
	double C_geo_boundary[3][2] = { {-0.75,0.75},{-0.75,0.75},{0,20} };         //  change this, unit in [cm]!!!
	double C_data[OUTDIMC0][OUTDIMC1][OUTDIMC2];
	
	double M_geo_binnumber[3] = { OUTDIMM0,OUTDIMM1,OUTDIMM2 };
	double M_geo_boundary[3][2] = { {-0.75,0.75},{-0.75,0.75},{0,20} };         //  change this, unit in [cm]!!!
	double M_data[OUTDIMM0][OUTDIMM1][OUTDIMM2];
	

	//Write HDF5 file for fuel
	file = H5Fcreate("info_fuel.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//Group
	group = H5Gcreate(file, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Attribute
	space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(group, "MeshType", H5T_STD_I32BE, space, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr, H5T_NATIVE_INT, a);
	//Geometry_Binnumber
	space = H5Screate_simple(1, dim_BinNumber, NULL);
	dataset = H5Dcreate(file, "/Geometry/BinNumber", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F_geo_binnumber);
	//Geometry_Boundary
	space = H5Screate_simple(2, dim_Boundary, NULL);
	dataset = H5Dcreate(file, "/Geometry/Boundary", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F_geo_boundary);

	//Dataset
	if ((in = fopen("powerfuel.dat", "r")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}
	for (int k = 0; k < OUTDIMF2; k++) {
		for (int j = 0; j < OUTDIMF1; j++) {
			for (int i = 0; i < OUTDIMF0; i++) {
				fscanf(in, "%lf", &power_data[i][j][k]);
			}
		}
	}
	fclose(in);
    if ((in = fopen("tempfuel.dat", "r")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}
	for (int k = 0; k < OUTDIMF2; k++) {
		for (int j = 0; j < OUTDIMF1; j++) {
			for (int i = 0; i < OUTDIMF0; i++) {
				fscanf(in, "%lf", &F_data[i][j][k]);
                F5_data[i][j][k] = F_data[i][j][k] + (pow(r1, 3)*power_data[i][j][k]/(6*k5)) * (r5*r5+r5*r4-2*r4*r4) / (r5*(r5*r5+r5*r4+r4*r4));
                T4 = F_data[i][j][k] + (pow(r1, 3)*power_data[i][j][k]/(3*k5))*(1/r4-1/r5);
                k4 = lambda4(T4, r1, r3, r4, power_data[i][j][k]);
                F4_data[i][j][k] = T4 + (pow(r1, 3)*power_data[i][j][k]/(6*k4)) * (r4*r4+r4*r3-2*r3*r3) / (r4*(r4*r4+r4*r3+r3*r3));
                T3 = T4 + (pow(r1, 3)*power_data[i][j][k]/(3*k4))*(1/r3-1/r4);
                F3_data[i][j][k] = T3 + (pow(r1, 3)*power_data[i][j][k]/(6*k3)) * (r3*r3+r3*r2-2*r2*r2) / (r3*(r3*r3+r3*r2+r2*r2));
                T2 = T3 + (pow(r1, 3)*power_data[i][j][k]/(3*k3))*(1/r2-1/r3);
                F2_data[i][j][k] = T2 + (pow(r1, 3)*power_data[i][j][k]/(6*k2)) * (r2*r2+r2*r1-2*r1*r1) / (r2*(r2*r2+r2*r1+r1*r1));
                T1 = T2 + (pow(r1, 3)*power_data[i][j][k]/(3*k2))*(1/r1-1/r2);
                k1 = lambda1(T1,r1, power_data[i][j][k]);
                F1_data[i][j][k] = T1 + power_data[i][j][k]*r1*r1/(15*k1);
			}
		}
	}
	fclose(in);

    space = H5Screate_simple(3, dim_F_data, NULL);

	dataset = H5Dcreate(file, "/temp_fuel", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F_data);
    dataset = H5Dcreate(file, "/temp_fuel1", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F1_data);
	dataset = H5Dcreate(file, "/temp_fuel2", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F2_data);
	dataset = H5Dcreate(file, "/temp_fuel3", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F3_data);
	dataset = H5Dcreate(file, "/temp_fuel4", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F4_data);
	dataset = H5Dcreate(file, "/temp_fuel5", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F5_data);

	status = H5Gclose(group);
	status = H5Sclose(space);
	status = H5Dclose(dataset);
	status = H5Fclose(file);


	//Write HDF5 file for coolant
	file = H5Fcreate("info_coolant.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//Group
	group = H5Gcreate(file, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Attribute
	space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(group, "MeshType", H5T_STD_I32BE, space, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr, H5T_NATIVE_INT, a);
	//Geometry_Binnumber
	space = H5Screate_simple(1, dim_BinNumber, NULL);
	dataset = H5Dcreate(file, "/Geometry/BinNumber", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C_geo_binnumber);
	//Geometry_Boundary
	space = H5Screate_simple(2, dim_Boundary, NULL);
	dataset = H5Dcreate(file, "/Geometry/Boundary", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C_geo_boundary);

	//Dataset for coolant temperature
	space = H5Screate_simple(3, dim_C_data, NULL);
	dataset = H5Dcreate(file, "/temp_coolant", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ((in = fopen("tempcoolant.dat", "r")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}
	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				fscanf(in, "%lf", &C_data[i][j][k]);
			}
		}
	}
	fclose(in);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C_data);

	//Dataset for coolant density
	dataset = H5Dcreate(file, "/r_coolant", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ((in = fopen("rcoolant.dat", "r")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}
	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				fscanf(in, "%lf", &C_data[i][j][k]);
			}
		}
	}
	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				C_data[i][j][k]/=1000;
			}
		}
	}
	fclose(in);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C_data);

	status = H5Gclose(group);
	status = H5Sclose(space);
	status = H5Dclose(dataset);
	status = H5Fclose(file);

	
	// Write HDF5 file for moderator
	file = H5Fcreate("info_moderator.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//Group
	group = H5Gcreate(file, "/Geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//Attribute
	space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(group, "MeshType", H5T_STD_I32BE, space, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite(attr, H5T_NATIVE_INT, a);
	//Geometry_Binnumber
	space = H5Screate_simple(1, dim_BinNumber, NULL);
	dataset = H5Dcreate(file, "/Geometry/BinNumber", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, M_geo_binnumber);
	//Geometry_Boundary
	space = H5Screate_simple(2, dim_Boundary, NULL);
	dataset = H5Dcreate(file, "/Geometry/Boundary", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, M_geo_boundary);

	//Dataset
	space = H5Screate_simple(3, dim_M_data, NULL);
	dataset = H5Dcreate(file, "/temp_moderator", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if ((in = fopen("tempmoderator.dat", "r")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}
	for (int k = 0; k < OUTDIMM2; k++) {
		for (int j = 0; j < OUTDIMM1; j++) {
			for (int i = 0; i < OUTDIMM0; i++) {
				fscanf(in, "%lf", &M_data[i][j][k]);
			}
		}
	}
	fclose(in);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, M_data);

	status = H5Gclose(group);
	status = H5Sclose(space);
	status = H5Dclose(dataset);
	status = H5Fclose(file);
	

	return 0;
}