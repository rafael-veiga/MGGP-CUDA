#ifndef SUBJECT_H_INCLUDED
#define SUBJECT_H_INCLUDED

#include "Tree.h"
#include <vector>

class Subject {
public:
	Tree * tree;
	int complex;
	double fitnessLS;
	double fitnessTestLS;
	double fitnessValidLS;

	int treino_vp;// verdadeiro positivo
	int treino_fp;// falso positivo
	int treino_vn;// verdadeiro negativo
	int treino_fn;// falso negativo

	int teste_vp;// verdadeiro positivo
	int teste_fp;// falso positivo
	int teste_vn;// verdadeiro negativo
	int teste_fn;// falso negativo

	int valida_vp;// verdadeiro positivo
	int valida_fp;// falso positivo
	int valida_vn;// verdadeiro negativo
	int valida_fn;// falso negativo
	int ranking = 0;

	bool printing;

	double crowdingDistance = 0;
	double* fitness;
	double* fitnessTest;

	Subject();
	Subject(Tree* n);
	double complexity();
	~Subject();
	void print();
};

extern Configures* h_conf;
extern Configures* d_conf;

#endif // SUBJECT_H_INCLUDED
