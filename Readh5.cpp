#include "stdio.h"
#include "string.h"
#include "hdf5.h"
#include "stdlib.h"


int main()
{
	FILE* in;

	hid_t file, space, group, attr, dataset;
	herr_t status;
	double DIMS[3];
	double boundary[3][2];
	file = H5Fopen("MeshTally1.h5", H5F_ACC_RDWR, H5P_DEFAULT);
	dataset = H5Dopen(file, "/Geometry/BinNumber", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, DIMS);
	dataset = H5Dopen(file, "/Geometry/Boundary", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, boundary);

	double*** data3D, ** data2D, * data1D;
	data3D = (double***)malloc((size_t)DIMS[0] * sizeof(double**));
	data2D = (double**)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * sizeof(double*));
	data1D = (double*)malloc((size_t)DIMS[0] * (size_t)DIMS[1] * (size_t)DIMS[2] * sizeof(double));
	for (int i = 0; i < DIMS[0] * DIMS[1]; i++) {
		data2D[i] = &data1D[i * (int)DIMS[2]];
	}
	for (int j = 0; j < DIMS[0]; j++) {
		data3D[j] = &data2D[j * (int)DIMS[1]];
	}


	dataset = H5Dopen(file, "/Type2", H5P_DEFAULT);
	status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(data3D[0][0][0]));
	if ((in = fopen("RMCPower.dat", "w")) == NULL)
	{
		puts("can't be opened"); exit(0);
	}

	for (int j = 0; j < 2; j++) {
		for (int i = 0; i < 3; i++) {
			fprintf(in, "%lf\n", boundary[i][j]);
		}
	}
	for (int k = 0; k < DIMS[2]; k++) {
		for (int j = 0; j < DIMS[1]; j++) {
			for (int i = 0; i < DIMS[0]; i++) {
				fprintf(in, "%lf\n", data3D[i][j][k]);
			}
		}
	}

	fclose(in);

	free(data3D);
	free(data2D);
	free(data1D);

	status = H5Dclose(dataset);
	status = H5Fclose(file);

	return 0;
}