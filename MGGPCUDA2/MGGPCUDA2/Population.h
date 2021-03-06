#ifndef POPULATION_H_INCLUDED
#define POPULATION_H_INCLUDED

#include "Gramatica.h"
#include "Subject.h"
#include "Configures.h"

extern Configures* h_conf;
extern Configures* d_conf;

class Population {
public:
	//        Subject* pop[conf->popSize];
	vector<Subject*> pop;

	Population();
	Population(Population* p);
	~Population();
	void print();
};

#endif // POPULATION_H_INCLUDED
