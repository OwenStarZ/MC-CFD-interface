/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_ADJUST(Hello, d);
extern DEFINE_INIT(read_power, d);
extern DEFINE_SOURCE(fuel_power, cell, thread, dS, eqn);
extern DEFINE_EXECUTE_AT_END(cal_TH);
UDF_Data udf_data[] = {
{"Hello", (void (*)(void))Hello, UDF_TYPE_ADJUST},
{"read_power", (void (*)(void))read_power, UDF_TYPE_INIT},
{"fuel_power", (void (*)(void))fuel_power, UDF_TYPE_SOURCE},
{"cal_TH", (void (*)(void))cal_TH, UDF_TYPE_EXECUTE_AT_END},
};
int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
