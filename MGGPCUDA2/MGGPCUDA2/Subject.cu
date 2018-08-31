#include <iostream>
#include <cmath>
#include "Configures.h"
#include "Subject.h"
#include "Tree.h"



Subject::Subject() {
	tree = new Tree();
	//    fitness = new double[objectives];
	//    fitnessTest = new double[objectives];
	fitnessLS = INFINITY;
	fitnessTestLS = INFINITY;
	printing = true;
};

void Subject::iniDeviceTree() {
	Device_Tree host;
	size_t tam = sizeof(double)*host.expCounter;
	host.expCounter = this->tree->expCounter;
	cudaMalloc(&host.exp, tam);
	cudaMemcpy(host.exp, this->tree->exp, tam, cudaMemcpyHostToDevice);
	tam = sizeof(Device_Tree);
	cudaMalloc(&this->d_tree, tam);
	cudaMemcpy(this->d_tree, &host, tam, cudaMemcpyHostToDevice);


}

void Subject::destDeviceTree() {
	cudaFree(this->d_tree_exp);
	cudaFree(this->d_tree);
}

Subject::Subject(Tree* n) {
	tree = n;
	fitnessLS = INFINITY;
	fitnessTestLS = INFINITY;
	printing = true;
};

Subject::~Subject() {
	//    tree->print();
	delete tree;
};

double Subject::complexity() {
	if (h_conf->complexity == 0)
		return tree->high;
	else if (h_conf->complexity == 1)
		return tree->terminals;
}

void Subject::print() {
	//    cout << ".";
	tree->print();
};
