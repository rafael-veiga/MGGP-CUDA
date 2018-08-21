
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
//#include <string>

//#include "Datamaker.h"
#include "No.h"
#include "Tree.h"
#include "Database.h"
#include "Gramatica.h"
#include "Configures.h"
#include "Operators.h"
#include "Search.h"

using namespace std;
Configures* h_conf;
Configures* d_conf;
/*
*
*/
__global__ void teste(Database* d_dados) {
	d_dados->countVar = 20;
}

int main(int argc, char** args) {

	// imprimindo argumentos
	//    for(int i=0; i<argc;i++){
	//    cout << " " << args[i];    
	//    }
	//    cout << endl;
	//string nome_saida = args[8];    
	freopen(args[8], "w", stdout);

	int seed = atoi(args[1]);
	string gramatica = args[2];
	string dados = args[3];
	string grupo = args[4];
	int geracoes = atoi(args[5]);
	int populacao = atoi(args[6]);
	int altura = atoi(args[7]);
	//int complexidade = atoi(args[8]); //0=high 1=terminals

	srand(seed);
	cout.precision(7);

	//configuracoes
	h_conf = new Configures();
	h_conf->treeHigh = altura;
	h_conf->popSize = populacao;
	h_conf->iterations = geracoes;
	h_conf->leastSquare = 1;
	h_conf->elitism = 0.1;
	h_conf->mono = 1; // 0 = monobjetivo;  1 = multiobjetivo
					  //conf->complexity = complexidade; // 0 = high 1 = terminals
	h_conf->complexity = 1;
	//gramatica
	gram = new Gramatica(gramatica);
	// cout<< "fim da gramatica" << endl;
	//operadores
	// cout << "inicio do operadores" << endl;
	// gram->imprimeGramatica();
	op = new Operators();
	// cout << "fim do operadores" << endl << "inicio dos database" << endl;
	//dados
	//    data = new Database("read/base5.txt", "read/grupo5.txt");
	Database *banco_dados = new Database(dados, grupo);
	Database* d_banco_dados =banco_dados->copyDevice();
	teste<<<1, 1>>>(d_banco_dados);
	//banco_dados->print();
	//cout << "fim do database" << endl << "inicio do search" << endl;
	//busca
	Search* s = new Search(banco_dados);
	delete s;

	fclose(stdout);

	return 0;
}

