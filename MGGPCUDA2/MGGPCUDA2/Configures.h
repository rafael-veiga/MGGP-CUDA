/*
* File:   Configures.h
* Author: rafae_000
*
* Created on 14 de Janeiro de 2016, 10:16
*/

#ifndef CONFIGURES_INCLUDED
#define CONFIGURES_INCLUDED
#include <stdlib.h>
#include <cuda_runtime.h>

class Configures {
public:
	int treeHigh;
	int popSize;
	int iterations;
	int leastSquare = 0;
	int mono;
	//    double peso0;
	//    double peso1;
	double elitism;

	int complexity = 0; //0=high 1=terminals

	__host__ __device__  Configures();

};

extern Configures* h_conf;
extern Configures* d_conf;

#endif // CONFIGURES_INCLUDED
