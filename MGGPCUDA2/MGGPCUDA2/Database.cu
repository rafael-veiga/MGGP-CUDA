#include "Database.h"

//Database* data;

//Database* Database::clone() {
//	Database* novo = new Database();
//	novo->countResults = this->countResults;
//	novo->countTestValues = this->countTestValues;
//	novo->countVar = this->countVar;
//	novo->testCount = this->testCount;
//	novo->trainCount = this->trainCount;
//	novo->validCount = this->validCount;
//}


Database* Database::copyDevice() {
	Database aux;
	Database* device;
	aux.countResults = this->countResults;
	aux.countTestValues = this->countTestValues;
	aux.countVar = this->countVar;
	aux.testCount = this->testCount;
	aux.trainCount = this->trainCount;
	aux.validCount = this->validCount;
	//alocar valores
	cudaMalloc(&aux.values, sizeof(double*)*aux.countResults);
	
	double** inter = new double*[aux.countResults];
	for (int i = 0; i < aux.countResults; i++) {
		cudaMalloc(&inter[i], sizeof(double)*aux.countVar);
		cudaMemcpy(inter[i], this->values[i], sizeof(double)*aux.countVar, cudaMemcpyHostToDevice);
	}
	
	cudaMemcpy(aux.values, inter, sizeof(double**)*aux.countResults,cudaMemcpyHostToDevice);
	delete inter;
	

	//alocar resultados
	cudaMalloc(&aux.results, sizeof(double)*aux.countResults);
	cudaMemcpy(aux.results, this->results, sizeof(double)*aux.countResults, cudaMemcpyHostToDevice);

	//alocar treino
	cudaMalloc(&aux.training, sizeof(int)*aux.trainCount);
	cudaMemcpy(aux.training, this->training, sizeof(int)*aux.trainCount, cudaMemcpyHostToDevice);

	//alocar teste
	cudaMalloc(&aux.test, sizeof(int)*aux.testCount);
	cudaMemcpy(aux.test, this->test, sizeof(int)*aux.testCount, cudaMemcpyHostToDevice);

	//alocar valida
	cudaMalloc(&aux.validation, sizeof(int)*aux.validCount);
	cudaMemcpy(aux.validation, this->validation, sizeof(int)*aux.validCount, cudaMemcpyHostToDevice);

	//carregando objeto
	cudaMalloc(&device, sizeof(Database));
	cudaMemcpy(device, &aux, sizeof(Database), cudaMemcpyHostToDevice);

	return device;
}

Database::Database(string base, string groups) {
	loadBase(base);
	loadGroups(groups);
};

void Database::loadBase(string base) {
	ifstream arq;
	string line;

	arq.open(base.c_str());

	arq >> countVar >> countResults;
	vars = new string[countVar];
	values = new double*[countResults];
	for (int i = 0; i < countResults; i++)
		values[i] = new double[countVar];
	results = new double[countResults];

	//inserido


	for (int i = 0; i < countVar; i++)
		arq >> vars[i];

	arq >> line;

	for (int i = 0; i < countResults; i++) {
		for (int j = 0; j < countVar; j++) {
			arq >> values[i][j];
		}
		arq >> results[i];
	}
	arq.close();
}

void Database::loadGroups(string groups) {
	//    ifstream arq;
	//    string line;
	//
	//    for(int i = 0; i < countVar; i++){
	//        int pos = gram->checkMapV(vars[i]);
	//        if(pos != -1){
	//            for(int j = 0; j < countResults; j++){
	//                swap(values[j][pos], values[j][i]);
	//            }
	//            swap(vars[i], vars[pos]);
	//        }
	//    }
	ifstream arq;
	string line;
	//int temp[30];

	int *temp = new int[countVar];
	for (int i = 0; i < countVar; i++)
		temp[i] = gram->checkMapV(vars[i]);


	for (int i = 0; i < countVar; i++) {
		int pos = temp[i];
		if (pos != -1) {
			for (int j = 0; j < countResults; j++) {
				swap(values[j][pos], values[j][i]);
			}
			swap(vars[i], vars[pos]);
			swap(temp[i], temp[pos]);
		}
	}

	arq.open(groups.c_str());
	arq >> trainCount >> testCount >> validCount;

	//conf->peso1 = 1.0 - conf->peso0;

	training = new int[trainCount];
	test = new int[testCount];
	validation = new int[validCount];

	arq >> line;
	for (int i = 0; i < trainCount; i++) {
		arq >> training[i];
	}
	arq >> line;
	for (int i = 0; i < testCount; i++) {
		arq >> test[i];
	}
	arq >> line;
	for (int i = 0; i < validCount; i++) {
		arq >> validation[i];
	}
	delete temp;
}

void Database::print() {
	cout << "vars: ";
	for (int i = 0; i < countVar; i++)
		cout << vars[i] << "  ";
	cout << endl;
	for (int i = 0; i < trainCount; i++) {
		cout << i << "   ";
		for (int j = 0; j < countVar; j++)
			cout << values[training[i]][j] << "  ";
		cout << "   " << results[training[i]] << endl;
	}
	cout << endl;
	for (int i = 0; i < testCount; i++) {
		cout << i << "   ";
		for (int j = 0; j < countVar; j++)
			cout << values[test[i]][j] << "  ";
		cout << "   " << results[test[i]] << endl;
	}
	cout << endl;
	for (int i = 0; i < validCount; i++) {
		cout << i << "   " << validation[i] << "  ";
		for (int j = 0; j < countVar; j++)
			cout << values[validation[i]][j] << "  ";
		cout << "   " << results[validation[i]] << endl;
	}
	cout << endl;
}

double* Database::getVars(int position) {
	return values[position];
}

