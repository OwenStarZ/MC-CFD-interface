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

#define OUTDIMR0 1
#define OUTDIMR1 1
#define OUTDIMR2 100

double lambda1(double T, double R, double q)
{
    double knew, kold, k1d, k1p, k2p, k4r, k0, Tave;
	double tal = 0.001;
	double p = 0.05;
	knew = 5.0;
	do
	{
		kold = knew;
		Tave = T + q * R * R / (15 * kold);
		k1d = (1.09/pow(tal, 3.265)+0.0643/sqrt(tal)*sqrt(Tave))*atan(1/(1.09/pow(tal, 3.265)+0.0643/sqrt(tal)*sqrt(Tave)));
		k1p = 1+0.019*tal/(3-0.019*tal)*(1/(1+exp(12-0.01*Tave)));
		k2p = (1-p)/(1+2*p);
		k4r = 1-(0.2/(1+exp((Tave-900)/80)))*(1-exp(-1*tal));
		k0 = 1/(0.0375+2.165e-4*Tave)+4.715e9/pow(Tave, 2)*exp(-16361/Tave);
		knew = k1d*k1p*k2p*k4r*k0;
	} while (fabs(knew - kold) > 0.000001);

	return knew;
}

double lambda4(double T, double R0, double R1, double R2, double q)
{
	double knew, kold, Tave;
	knew = 42.0;
	do
	{
		kold = knew;
		Tave = T + (pow(R0, 3)* q/(6*kold)) * (R2*R2+R2*R1-2*R1*R1) / (R2*(R2*R2+R2*R1+R1*R1));
		knew = 0.0391112*exp(0.00224732*973.15)*(42.58-15564/Tave+1.2977E7/pow(Tave, 2)-1.8458E9/pow(Tave, 3));
	} while (fabs(knew - kold) > 0.000001);
	
	return knew;
}

void ReadRMCPowerTally(double b[], double p[])
{
	hid_t file, space, group, attr, dataset;
	herr_t status;
	double DIMS[3];
	double boundary[3][2];
	//double data3D[DIM0][DIM1][DIM2];
	file = H5Fopen("MeshTally1.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/Geometry/BinNumber", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, DIMS);
	dataset = H5Dopen(file, "/Geometry/Boundary", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, boundary);

	double*** data3D = (double***)malloc((size_t)DIMS[0] * sizeof(double**));
	double** data2D = (double**)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * sizeof(double*));
	double* data1D = (double*)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * (size_t)DIMS[2] * sizeof(double));
	for (int i = 0; i < DIMS[0] * DIMS[1]; i++) {
		data2D[i] = &data1D[i * (int)DIMS[2]];
	}
	for (int j = 0; j < DIMS[0]; j++) {
		data3D[j] = &data2D[j * (int)DIMS[1]];
	}

	dataset = H5Dopen(file, "/Type2", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

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

	free(data3D);
	free(data2D);
	free(data1D);

	status = H5Dclose(dataset);
	status = H5Fclose(file);
}

void WriteFuelData(double data[])
{
	hid_t file, dataset;
	herr_t status;

	double*** data3D = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		data2D[i] = &data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		data3D[j] = &data2D[j * (int)OUTDIMF1];
	}
	
	//Write HDF5 file for fuel
	file = H5Fopen("info_fuel.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/temp_fuel", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMF2; k++) {
		for (int j = 0; j < OUTDIMF1; j++) {
			for (int i = 0; i < OUTDIMF0; i++) {
				data3D[i][j][k] = data[i + OUTDIMF0 * j + OUTDIMF0 * OUTDIMF1 * k];
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	free(data3D);
	free(data2D);
	free(data1D);
}

void WriteFuelData_Multilevel(double data1[], double data2[])
{
	hid_t file, dataset;
	herr_t status;

	double*** F_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F_data2D[i] = &F_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F_data[j] = &F_data2D[j * (int)OUTDIMF1];
	}
    double*** power_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** power_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* power_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		power_data2D[i] = &power_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		power_data[j] = &power_data2D[j * (int)OUTDIMF1];
	}
	double*** F1_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F1_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F1_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F1_data2D[i] = &F1_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F1_data[j] = &F1_data2D[j * (int)OUTDIMF1];
	}
	double*** F2_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F2_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F2_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F2_data2D[i] = &F2_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F2_data[j] = &F2_data2D[j * (int)OUTDIMF1];
	}
	double*** F3_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F3_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F3_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F3_data2D[i] = &F3_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F3_data[j] = &F3_data2D[j * (int)OUTDIMF1];
	}
	double*** F4_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F4_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F4_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F4_data2D[i] = &F4_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F4_data[j] = &F4_data2D[j * (int)OUTDIMF1];
	}
	double*** F5_data = (double***)malloc((size_t)OUTDIMF0 * sizeof(double**));
	double** F5_data2D = (double**)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * sizeof(double*));
	double* F5_data1D = (double*)malloc((size_t)OUTDIMF0 * (size_t)OUTDIMF1 * (size_t)OUTDIMF2 * sizeof(double));
	for (int i = 0; i < OUTDIMF0 * OUTDIMF1; i++) {
		F5_data2D[i] = &F5_data1D[i * (int)OUTDIMF2];
	}
	for (int j = 0; j < OUTDIMF0; j++) {
		F5_data[j] = &F5_data2D[j * (int)OUTDIMF1];
	}

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
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F_data[0][0][0]);
    dataset = H5Dopen(file, "/temp_fuel1", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F1_data[0][0][0]);
	dataset = H5Dopen(file, "/temp_fuel2", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F2_data[0][0][0]);
	dataset = H5Dopen(file, "/temp_fuel3", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F3_data[0][0][0]);
	dataset = H5Dopen(file, "/temp_fuel4", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F4_data[0][0][0]);
	dataset = H5Dopen(file, "/temp_fuel5", H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &F5_data[0][0][0]);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	free(F_data);
	free(F_data2D);
	free(F_data1D);
	free(power_data);
	free(power_data2D);
	free(power_data1D);
	free(F1_data);
	free(F1_data2D);
	free(F1_data1D);
	free(F2_data);
	free(F2_data2D);
	free(F2_data1D);
	free(F3_data);
	free(F3_data2D);
	free(F3_data1D);
	free(F4_data);
	free(F4_data2D);
	free(F4_data1D);
	free(F5_data);
	free(F5_data2D);
	free(F5_data1D);
}

void WriteModeratorData(double data[])
{
	hid_t file, dataset;
	herr_t status;

	double*** data3D = (double***)malloc((size_t)OUTDIMM0 * sizeof(double**));
	double** data2D = (double**)malloc((size_t)OUTDIMM0 * (size_t)OUTDIMM1 * sizeof(double*));
	double* data1D = (double*)malloc((size_t)OUTDIMM0 * (size_t)OUTDIMM1 * (size_t)OUTDIMM2 * sizeof(double));
	for (int i = 0; i < OUTDIMM0 * OUTDIMM1; i++) {
		data2D[i] = &data1D[i * (int)OUTDIMM2];
	}
	for (int j = 0; j < OUTDIMM0; j++) {
		data3D[j] = &data2D[j * (int)OUTDIMM1];
	}

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

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	free(data3D);
	free(data2D);
	free(data1D);
}

void WriteCoolantData(double data1[], double data2[])
{
	hid_t file, dataset;
	herr_t status;

	double*** data3D = (double***)malloc((size_t)OUTDIMC0 * sizeof(double**));
	double** data2D = (double**)malloc((size_t)OUTDIMC0 * (size_t)OUTDIMC1 * sizeof(double*));
	double* data1D = (double*)malloc((size_t)OUTDIMC0 * (size_t)OUTDIMC1 * (size_t)OUTDIMC2 * sizeof(double));
	for (int i = 0; i < OUTDIMC0 * OUTDIMC1; i++) {
		data2D[i] = &data1D[i * (int)OUTDIMC2];
	}
	for (int j = 0; j < OUTDIMC0; j++) {
		data3D[j] = &data2D[j * (int)OUTDIMC1];
	}

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

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

	//Dataset for coolant density
	dataset = H5Dopen(file, "/r_coolant", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				data3D[i][j][k] = data2[i + OUTDIMC0 * j + OUTDIMC0 * OUTDIMC1 * k] / 1000;
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	free(data3D);
	free(data2D);
	free(data1D);
}

void WriteReflectorData(double data[])
{
	hid_t file, dataset;
	herr_t status;

	double*** data3D = (double***)malloc((size_t)OUTDIMR0 * sizeof(double**));
	double** data2D = (double**)malloc((size_t)OUTDIMR0 * (size_t)OUTDIMR1 * sizeof(double*));
	double* data1D = (double*)malloc((size_t)OUTDIMR0 * (size_t)OUTDIMR1 * (size_t)OUTDIMR2 * sizeof(double));
	for (int i = 0; i < OUTDIMR0 * OUTDIMR1; i++) {
		data2D[i] = &data1D[i * (int)OUTDIMR2];
	}
	for (int j = 0; j < OUTDIMR0; j++) {
		data3D[j] = &data2D[j * (int)OUTDIMR1];
	}

	// Write HDF5 file for reflector
	file = H5Fopen("info_reflector.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/temp_reflector", H5P_DEFAULT);

	for (int k = 0; k < OUTDIMR2; k++) {
		for (int j = 0; j < OUTDIMR1; j++) {
			for (int i = 0; i < OUTDIMR0; i++) {
				data3D[i][j][k] = data[i + OUTDIMR0 * j + OUTDIMR0 * OUTDIMR1 * k];
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data3D[0][0][0]);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	free(data3D);
	free(data2D);
	free(data1D);
}