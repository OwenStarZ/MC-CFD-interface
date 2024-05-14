#include "stdio.h"
#include "string.h"
#include "hdf5.h"
#include "stdlib.h"
#include "math.h"
#include "h5rw.hpp"

/******** change the dimensions, note that they should be same as those in UDF ********/
#define DIM0 8
#define DIM1 8
#define DIM2 100

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

void ReadRMCPowerTally(double b[], double p[])
{
	hid_t file, space, group, attr, dataset;
	herr_t status;
	double DIMS[3];
	double boundary[3][2];
	double data3D[DIM0][DIM1][DIM2];
	file = H5Fopen("MeshTally1.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/Geometry/BinNumber", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, DIMS);
	dataset = H5Dopen(file, "/Geometry/Boundary", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, boundary);

	/*double*** data3D, ** data2D, * data1D;
	data3D = (double***)malloc((size_t)DIMS[0] * sizeof(double**));
	data2D = (double**)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * sizeof(double*));
	data1D = (double*)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * (size_t)DIMS[2] * sizeof(double));
	for (int i = 0; i < DIMS[0] * DIMS[1]; i++) {
		data2D[i] = &data1D[i * (int)DIMS[2]];
	}
	for (int j = 0; j < DIMS[0]; j++) {
		data3D[j] = &data2D[j * (int)DIMS[1]];
	}*/

	dataset = H5Dopen(file, "/Type2", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3D);

	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 3; i++) {
			b[i + 3 * j] = boundary[i][j];
		}
	}
	
	double total = 0.0;

	for (int k = 0; k < DIMS[2]; k++) {
		for (int j = 0; j < DIMS[1]; j++) {
			for (int i = 0; i < DIMS[0]; i++) {
				p[i + (int)DIMS[0] * j + (int)DIMS[0] * (int)DIMS[1] * k] = data3D[i][j][k];
				total += data3D[i][j][k];
			}
		}
	}

	for (int k = 0; k < DIMS[2]; k++) {
		for (int j = 0; j < DIMS[1]; j++) {
			for (int i = 0; i < DIMS[0]; i++) {
				p[i + (int)DIMS[0] * j + (int)DIMS[0] * (int)DIMS[1] * k] /= total;
			}
		}
	}

	//free(data3D);
	//free(data2D);
	//free(data1D);

	status = H5Dclose(dataset);
	status = H5Fclose(file);
}

void WriteFueldata(double data1[], double data2[])
{
	hid_t file, dataset;
	herr_t status;

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

	//Write HDF5 file for fuel
	file = H5Fopen("info_fuel.h5", H5F_ACC_RDWR, H5P_DEFAULT);

	//Dataset
	
	for (int k = 0; k < OUTDIMF2; k++) {
		for (int j = 0; j < OUTDIMF1; j++) {
			for (int i = 0; i < OUTDIMF0; i++) {
				F_data[i][j][k] = data1[i + OUTDIMF0 * j + OUTDIMF0 * OUTDIMF1 * k];
                power_data[i][j][k] = data2[i + OUTDIMF0 * j + OUTDIMF0 * OUTDIMF1 * k];
                F5_data[i][j][k] = F_data[i][j][k] + (pow(r1, 3)*power_data[i][j][k]/(6*k5)) * (r5*r5+r5*r4-2*r4*r4) / (r5*(r5*r5+r5*r4+r4*r4));
                T4 = F_data[i][j][k] + (pow(r1, 3)*power_data[i][j][k]/(3*k5))*(1/r4-1/r5);
                k4 = lambda4(T4, r1, r3, r4, power_data[i][j][k]);
                F4_data[i][j][k] = T4 + (pow(r1, 3)*power_data[i][j][k]/(6*k4)) * (r4*r4+r4*r3-2*r3*r3) / (r4*(r4*r4+r4*r3+r3*r3));
                T3 = T4 + (pow(r1, 3)*power_data[i][j][k]/(3*k4))*(1/r3-1/r4);
                F3_data[i][j][k] = T3 + (pow(r1, 3)*power_data[i][j][k]/(6*k3)) * (r3*r3+r3*r2-2*r2*r2) / (r3*(r3*r3+r3*r2+r2*r2));
                T2 = T3 + (pow(r1, 3)*power_data[i][j][k]/(3*k3))*(1/r2-1/r3);
                F2_data[i][j][k] = T2 + (pow(r1, 3)*power_data[i][j][k]/(6*k2)) * (r2*r2+r2*r1-2*r1*r1) / (r2*(r2*r2+r2*r1+r1*r1));
                T1 = T2 + (pow(r1, 3)*power_data[i][j][k]/(3*k2))*(1/r1-1/r2);
                k1 = lambda1(T1, r1, power_data[i][j][k]);
                F1_data[i][j][k] = T1 + power_data[i][j][k]*r1*r1/(15*k1);
			}
		}
	}

	dataset = H5Dopen(file, "/temp_fuel", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F_data);
    dataset = H5Dopen(file, "/temp_fuel1", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F1_data);
	dataset = H5Dopen(file, "/temp_fuel2", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F2_data);
	dataset = H5Dopen(file, "/temp_fuel3", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F3_data);
	dataset = H5Dopen(file, "/temp_fuel4", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F4_data);
	dataset = H5Dopen(file, "/temp_fuel5", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F5_data);

	status = H5Dclose(dataset);
	status = H5Fclose(file);
}

void WriteModeratorData(double data[])
{
	hid_t file, dataset;
	herr_t status;

	double data3D[OUTDIMM0][OUTDIMM1][OUTDIMM2];

	// Write HDF5 file for moderator
	file = H5Fopen("info_moderator.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/temp_moderator", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMM2; k++) {
		for (int j = 0; j < OUTDIMM1; j++) {
			for (int i = 0; i < OUTDIMM0; i++) {
				data3D[i][j][k] = data[i + OUTDIMM0 * j + OUTDIMM0 * OUTDIMM1 * k];
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3D);

	status = H5Dclose(dataset);
	status = H5Fclose(file);
}

void WriteCoolantData(double data1[], double data2[])
{
	hid_t file, dataset;
	herr_t status;

	double data3D[OUTDIMC0][OUTDIMC1][OUTDIMC2];

	//Write HDF5 file for coolant
	file = H5Fopen("info_coolant.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/temp_coolant", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				data3D[i][j][k] = data1[i + OUTDIMC0 * j + OUTDIMC0 * OUTDIMC1 * k];
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3D);

	//Dataset for coolant density
	dataset = H5Dopen(file, "/r_coolant", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				data3D[i][j][k] = data2[i + OUTDIMC0 * j + OUTDIMC0 * OUTDIMC1 * k] / 1000;
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3D);

	status = H5Dclose(dataset);
	status = H5Fclose(file);
}
