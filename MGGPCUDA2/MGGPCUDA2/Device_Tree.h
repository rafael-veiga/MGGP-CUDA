#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Device_Tree
{
public:
	int expCounter;
	double* exp;

	__host__ __device__ Device_Tree();
	~Device_Tree();
	
};

