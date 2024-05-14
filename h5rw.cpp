#include "stdio.h"
#include "string.h"
#include "hdf5.h"
#include "stdlib.h"
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

void WriteFuelData(double data[])
{
	hid_t file, dataset;
	herr_t status;

	double data3D[OUTDIMF0][OUTDIMF1][OUTDIMF2];

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

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3D);

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
