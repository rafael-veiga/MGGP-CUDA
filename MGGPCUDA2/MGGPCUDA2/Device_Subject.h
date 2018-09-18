#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Subject.h"

class Device_Subject
{
public:
	double* d_tree_exp;
	int vp;
	int fp;
	int fn;
	int vn;
	double erro;
	int d_tree_countExp;
	
	void iniDeviceTree(Subject* sub);
	void destDeviceTree();
	__host__ __device__ Device_Subject();
	__host__ __device__ ~Device_Subject() {

	}
	

};

