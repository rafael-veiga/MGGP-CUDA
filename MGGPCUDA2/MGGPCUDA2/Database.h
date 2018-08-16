#ifndef DATABASE_H_INCLUDED
#define DATABASE_H_INCLUDED
#include "Gramatica.h"
#include "Configures.h"
#include <cstdio>
#include <string>
//#include <iostream>
#include <string>
#include <fstream>

using namespace std;



class Database {
public:
	int countVar; //quantidade de variaveis no arquivo texto
	string* vars;
	int countTestValues;
	double** values; //totalResultados x totalVariaveis
	int countResults;
	double* results;
	int* training;
	int trainCount;
	int* test;
	int testCount;
	int* validation;
	int validCount;

	Database(string base, string groups);
	void loadBase(string base);
	void loadGroups(string groups);
	void print();
	//void selectData();
	double* getVars(int position);


};
//Database* banco_dados;
//

#endif // DATABASE_H_INCLUDED

