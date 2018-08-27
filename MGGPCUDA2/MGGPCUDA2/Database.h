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
	//char** d_vars;
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
	
	//__host__ Database* clone();
	__host__ __device__ Database() = default;
	__host__ __device__ Database(string base, string groups);
	__host__ Database* copyDevice();
	//__host__  ~Database();
	void loadBase(string base);
	void loadGroups(string groups);
	void print();
	//void selectData();
	double* getVars(int position);
private:
//__host__	char** arrayStringAlocate(string* array);
//__host__	double** alocateValues();


};
//Database* banco_dados;
//

#endif // DATABASE_H_INCLUDED

