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
