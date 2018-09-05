#include "Configures.h"
#include "Device_Subject.h"

Device_Subject::Device_Subject() {
	
};

void Device_Subject::iniDeviceTree( Subject* sub) {
	
	
	size_t tam = sizeof(double)*sub->tree->expCounter;
	this->d_tree_countExp = sub->tree->expCounter;
	cudaMalloc(&this->d_tree_exp, tam);
	cudaMemcpy(this->d_tree_exp, sub->tree->exp, tam, cudaMemcpyHostToDevice);
	

}

void Device_Subject::destDeviceTree() {
	
	cudaFree(this->d_tree_exp);
}