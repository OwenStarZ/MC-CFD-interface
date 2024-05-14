#include "stdio.h"
#include "string.h"
#include "hdf5.h"
#include "stdlib.h"

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

#define F_xmin -0.75
#define F_xmax 0.75
#define F_ymin -0.75
#define F_ymax 0.75
#define F_zmin 0.0
#define F_zmax 20.0

#define C_xmin -0.75
#define C_xmax 0.75
#define C_ymin -0.75
#define C_ymax 0.75
#define C_zmin 0.0
#define C_zmax 20.0

#define M_xmin -0.75
#define M_xmax 0.75
#define M_ymin -0.75
#define M_ymax 0.75
#define M_zmin 0.0
#define M_zmax 20.0

#define INIT_temp_fuel 700
#define INIT_temp_coolant 500
#define INIT_r_coolant 0.8
#define INIT_temp_moderator 600

int main()
{
	hid_t file, space, group, attr, dataset;
	herr_t status;
	hsize_t dim_attr[1] = { 1 }, dim_BinNumber[1] = { 3 }, dim_Boundary[2] = { 3,2 };
	hsize_t dim_F_data[3] = { OUTDIMF0,OUTDIMF1,OUTDIMF2 };
	hsize_t dim_C_data[3] = { OUTDIMC0,OUTDIMC1,OUTDIMC2 };
	hsize_t dim_M_data[3] = { OUTDIMM0,OUTDIMM1,OUTDIMM2 };

	int a[1] = { 1 }; // 1 for uniform mesh, 2 for the other

	double F_geo_binnumber[3] = { OUTDIMF0,OUTDIMF1,OUTDIMF2 };
	double F_geo_boundary[3][2] = { {F_xmin,F_xmax},{F_ymin,F_ymax},{F_zmin,F_zmax} };         //  change this, unit in [cm]!!!
	double F_data[OUTDIMF0][OUTDIMF1][OUTDIMF2];

	double C_geo_binnumber[3] = { OUTDIMC0,OUTDIMC1,OUTDIMC2 };
	double C_geo_boundary[3][2] = { {C_xmin,C_xmax},{C_ymin,C_ymax},{C_zmin,C_zmax} };         //  change this, unit in [cm]!!!
	double C_data[OUTDIMC0][OUTDIMC1][OUTDIMC2];
	
	double M_geo_binnumber[3] = { OUTDIMM0,OUTDIMM1,OUTDIMM2 };
	double M_geo_boundary[3][2] = { {M_xmin,M_xmax},{M_ymin,M_ymax},{M_zmin,M_zmax} };         //  change this, unit in [cm]!!!
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
	space = H5Screate_simple(3, dim_F_data, NULL);

	dataset = H5Dcreate(file, "/temp_fuel", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	for (int k = 0; k < OUTDIMF2; k++) {
		for (int j = 0; j < OUTDIMF1; j++) {
			for (int i = 0; i < OUTDIMF0; i++) {
				F_data[i][j][k]=INIT_temp_fuel;
			}
		}
	}
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, F_data);

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

	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				C_data[i][j][k]=INIT_temp_coolant;
			}
		}
	}

	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, C_data);

	//Dataset for coolant density
	dataset = H5Dcreate(file, "/r_coolant", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	for (int k = 0; k < OUTDIMC2; k++) {
		for (int j = 0; j < OUTDIMC1; j++) {
			for (int i = 0; i < OUTDIMC0; i++) {
				C_data[i][j][k]=INIT_r_coolant;
			}
		}
	}


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

	for (int k = 0; k < OUTDIMM2; k++) {
		for (int j = 0; j < OUTDIMM1; j++) {
			for (int i = 0; i < OUTDIMM0; i++) {
				M_data[i][j][k]=INIT_temp_moderator;
			}
		}
	}
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, M_data);

	status = H5Gclose(group);
	status = H5Sclose(space);
	status = H5Dclose(dataset);
	status = H5Fclose(file);
	

	return 0;
}