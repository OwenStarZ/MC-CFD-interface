#ifndef _SO_H5RW_H_
#define _SO_H5RW_H_

#if __cplusplus
extern "C" {
#endif
	void ReadRMCPowerTally(double b[], double p[]);
	void WriteFuelData(double data[]);
	void WriteModeratorData(double data[]);
	void WriteCoolantData(double data1[], double data2[]);
#if __cplusplus
}
#endif

#endif
