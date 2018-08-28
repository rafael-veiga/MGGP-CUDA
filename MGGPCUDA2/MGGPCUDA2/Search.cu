#include <algorithm>
#include <math.h>
#include <iomanip>
#include "Search.h"
//#include "Population.h"
#include "Subject.h"
#include "Database.h"
#include "Gramatica.h"
#include "Configures.h"
#include "Parser.h"

using namespace std;

__global__ void teste(Database* d_dados, Configures* d_conf,Subject** d_pop) {
	d_dados->countVar = 20;
}

void Search::GPUcalcFitnessLS(int ini,int fim) {
	Subject** d_pop;
	Subject** aux;
	int tam = fim - ini;
	aux = new Subject*[tam];
	//carregando na GPU
	for (int i = 0; i < tam; i++) {
		//pop[i + ini]->iniDeviceTree();
		cudaMalloc(&aux[i], sizeof(Subject));
		cudaMemcpy(aux[i],pop[i+ini],sizeof(Subject),cudaMemcpyHostToDevice);
	}

	cudaMalloc(&d_pop, sizeof(Subject*)*tam);
	cudaMemcpy(d_pop, aux, sizeof(Subject*)*tam, cudaMemcpyHostToDevice);

	
	//executando
	teste<<<1, 1>>>(this->d_banco_dados, this->d_conf,d_pop);
	//descaregando da GPU
	for (int i = 0; i < tam; i++) {
		//pop[i + ini]->destDeviceTree();
		cudaFree(&d_pop[i]);
	}
	cudaFree(d_pop);
	delete aux;
}

bool mySort(Subject* a, Subject* b) {
	//    return (a->fitness < b->fitness);
	return (a->fitnessLS < b->fitnessLS);
};

bool sortTest(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->fitnessTestLS < b->fitnessTestLS);
};

bool sortRank(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->ranking < b->ranking);
};

bool Search::sortHigh(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->tree->high < b->tree->high);
};

bool Search::sortComplexity(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->complex < b->complex);
};

bool Search::sortSize(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->tree->terminals < b->tree->terminals);
};

bool Search::sortFit(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->fitnessLS < b->fitnessLS);
};

bool Search::sortFitTest(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->fitnessTestLS < b->fitnessTestLS);
};

bool Search::sortFitValid(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->fitnessValidLS < b->fitnessValidLS);
};

bool sortCrow(Subject* a, Subject* b) {
	//    return (a->fitnessTest < b->fitnessTest);
	return (a->crowdingDistance > b->crowdingDistance);
};

Search::Search(Database* banco_dados, Database* d_banco) {
	testError = INFINITY;
	testTr = INFINITY;
	errors = 0;
	this->banco_dados = banco_dados;
	this->d_banco_dados = d_banco;
	cudaMalloc(&this->d_conf, sizeof(Configures));
	cudaMemcpy(this->d_conf, h_conf, sizeof(Configures), cudaMemcpyHostToDevice);
	aux = new int[h_conf->popSize];
	for (int i = 0; i < h_conf->popSize; i++)
		aux[i] = 0;

	h_conf->mono == 0 ? doMonoSearch() : doMultiSearch();
};

void Search::doMonoSearch() {
	pop = new Subject*[h_conf->popSize * 2];
	convergence = new double[h_conf->iterations];

	for (int i = 0; i < h_conf->popSize; i++) { //Cria população
		pop[i] = new Subject(gram->criaArvExp());
	}

	for (int i = 0; i < h_conf->popSize; i++) //Calcula primeira vez
	{
		calcFitnessLS(pop[i]);
	}

	int size = (h_conf->popSize * h_conf->elitism);

	for (int it = 0; it < h_conf->iterations; it++) {
		Operate();
		sort(pop, pop + h_conf->popSize, sortFit);
		sort(pop + h_conf->popSize, pop + h_conf->popSize * 2, sortFit);

		sort(pop, pop + h_conf->popSize + size, sortFit);
		for (int i = h_conf->popSize; i < h_conf->popSize * 2; i++)
			delete pop[i];
		convergence[it] = pop[0]->fitnessLS;

		//        cout << it << endl;
		//        if(it > (int)(conf->iterations * 0.00005))
		//            cout << convergence[it]/convergence[(int)(it - conf->iterations * 0.00005)] << endl;
		//        if(it > (int)(conf->iterations * 0.00005) && convergence[it] > (convergence[(int)(it - conf->iterations * 0.00005)] - convergence[(int)(it - conf->iterations * 0.00005)] * 0.00000000001)){
		//            cout << "Stopped Converging on " << it << "!" << endl;
		//            break;
		//        }
	}

	for (int i = 0; i < h_conf->popSize; i++) {
		calcFitnessTestLS(pop[i]);
	}

	sort(pop, pop + h_conf->popSize, sortFitTest);

	for (int i = 0; i < h_conf->popSize; i++) {
		calcFitnessValidLS(pop[i]);
	}

	sort(pop, pop + h_conf->popSize, sortFitValid);

	for (int i = 0; i < h_conf->popSize; i++) {
		cout << pop[i]->ranking << ", "
			<< pop[i]->tree->high << ", "
			<< pop[i]->fitnessLS << ", "
			<< pop[i]->fitnessTestLS << ", "
			<< pop[i]->fitnessValidLS << ", "
			<< pop[i]->tree->infixPrint() << ", ";
		pop[i]->tree->printAlphas();
		cout << endl;
	}

	//    double diff = abs(convergence[conf->iterations - 1]) + 1;
	//    double max = convergence[conf->iterations - 1] < 0 ? convergence[0] + convergence[conf->iterations - 1] : convergence[0];
	//    cout << diff << " " << max << endl;
	//    cin.get();
	//
	//    for(int i = 0; i < conf->iterations; i += conf->iterations * .1){
	//        cout << convergence[i] << endl;
	////        for(int j = 0; j < convergence[i] + diff; j += max * .1){
	////            cout << "-";
	////        }
	////        cout << "." << endl;
	//    }
};

void Search::doMultiSearch() {
	pop = new Subject*[h_conf->popSize * 2];
	convergence = new double[h_conf->iterations];

	for (int i = 0; i < h_conf->popSize; i++) { //Cria população
		pop[i] = new Subject(gram->criaArvExp());
	}

	for (int i = 0; i < h_conf->popSize; i++) {
		calcFitnessLS(pop[i]);
	}

	//    for(int i = 0; i < 100000000; i++)
	//        calcRank(conf->popSize);

	//    double old = 0;
	for (int it = 0; it < h_conf->iterations; it++) {
		Operate();
		int size = h_conf->popSize * 2;
		calcRank(size);
		// apagar piores
		for (int i = h_conf->popSize; i < size; i++)
			delete pop[i];

		// teste
		if (it % 10 == 0) {
			cout << it << "\t";
			for (int i = 0; i < h_conf->popSize; i++) {
				if (pop[i]->ranking == 0) {
					cout << pop[i]->tree->high << ","
						<< pop[i]->tree->terminals << ","
						<< pop[i]->fitnessLS;
					cout << ";";
				}

			}
			cout << endl;
		}
	}
	// fim teste
	cout << "RESULTADO" << endl;
	for (int i = 0; i < h_conf->popSize; i++) {
		calcFitnessTestLS(pop[i]);
	}

	calcRankTest(h_conf->popSize);

	for (int i = 0; i < h_conf->popSize; i++) {
		calcFitnessValidLS(pop[i]);
	}

	calcRankValid(h_conf->popSize);

	for (int i = 0; i < h_conf->popSize; i++) {
		cout << pop[i]->ranking << ", "
			<< pop[i]->tree->high << ", "
			<< pop[i]->tree->terminals << ", "
			<< pop[i]->fitnessLS << ", "
			<< pop[i]->fitnessTestLS << ", "
			<< pop[i]->fitnessValidLS << ", "
			<< pop[i]->tree->infixPrint() << ", " << "pop[i]->tree->print()"
			<< ",( " << pop[i]->treino_vp << ";" << pop[i]->treino_vn << ";" << pop[i]->treino_fp << ";" << pop[i]->treino_fn << ")"
			<< ",( " << pop[i]->teste_vp << ";" << pop[i]->teste_vn << ";" << pop[i]->teste_fp << ";" << pop[i]->teste_fn << ")"
			<< ",( " << pop[i]->valida_vp << ";" << pop[i]->valida_vn << ";" << pop[i]->valida_fp << ";" << pop[i]->valida_fn << ")";
		//             pop[i]->tree->printAlphas();
		cout << endl;
	}
	// inserir modificação
	cout << "VETORES" << endl;
	cout << "TREINO" << endl;
	for (int i = 0; i<h_conf->popSize; i++) {
		Subject* s = pop[i];
		for (int l = 0; l<banco_dados->trainCount; l++) {

			cout << s->tree->treeResult(banco_dados->values[banco_dados->training[l]], NULL, 0) << ";";
		}

		cout << endl;
	}
	cout << "TESTE" << endl;
	for (int i = 0; i<h_conf->popSize; i++) {
		Subject* s = pop[i];
		for (int l = 0; l<banco_dados->testCount; l++) {

			cout << s->tree->treeResult(banco_dados->values[banco_dados->test[l]], NULL, 0) << ";";
		}

		cout << endl;
	}
	cout << "VALIDATION" << endl;
	for (int i = 0; i<h_conf->popSize; i++) {
		Subject* s = pop[i];
		for (int l = 0; l<banco_dados->validCount; l++) {

			cout << s->tree->treeResult(banco_dados->values[banco_dados->validation[l]], NULL, 0) << ";";
		}
		cout << endl;
	}

	// Encontrado
	cout << "ENCONTRADO" << endl;
	for (int i = 0; i<banco_dados->trainCount; i++) {
		cout << banco_dados->results[banco_dados->training[i]] << ";";
	}
	cout << endl;
	for (int i = 0; i<banco_dados->testCount; i++) {
		cout << banco_dados->results[banco_dados->test[i]] << ";";
	}
	cout << endl;
	for (int i = 0; i<banco_dados->validCount; i++) {
		cout << banco_dados->results[banco_dados->validation[i]] << ";";
	}
};

int Search::tournamentMono(int a, int b) {
	if (pop[a]->fitnessLS < pop[b]->fitnessLS)
		return a;
	else
		return b;
}

int Search::tournamentMulti(int a, int b) {
	if (pop[a]->ranking < pop[b]->ranking)
		return a;
	else if (pop[a]->ranking > pop[b]->ranking)
		return b;
	else {
		if (pop[a]->crowdingDistance > pop[b]->crowdingDistance)
			return a;
		else
			return b;
	}
}

void Search::Operate() {
	//instanciar
	for (int i = h_conf->popSize; i < h_conf->popSize * 2; i += 2) {
		pop[i] = (new Subject());
		pop[i + 1] = (new Subject());
	}
	//mutação e cross
	for (int i = h_conf->popSize; i < h_conf->popSize * 2; i += 2) {
		int a, b;
		int s1 = rand() % h_conf->popSize;
		int s2 = rand() % h_conf->popSize;
		a = h_conf->mono == 0 ? tournamentMono(s1, s2) : tournamentMulti(s1, s2);

		s1 = rand() % h_conf->popSize;
		s2 = rand() % h_conf->popSize;
		b = h_conf->mono == 0 ? tournamentMono(s1, s2) : tournamentMulti(s1, s2);

		aux[a]++;
		aux[b]++;

		op->Cross(pop[a]->tree, pop[b]->tree, pop[i]->tree, pop[i + 1]->tree);

		op->Mutate(pop[i]->tree);
		op->Mutate(pop[i + 1]->tree);
	}

	//calcFitnessLS 
	//paralelizar

	//GPUcalcFitnessLS(h_conf->popSize, h_conf->popSize * 2);
	GPUcalcFitnessLS(0,h_conf->popSize);

	for (int i = h_conf->popSize; i < h_conf->popSize * 2; i++) {
		calcFitnessLS(pop[i]);
	}

	
};

bool Search::dominate(Subject* a, Subject* b) {
	if (a->fitnessLS <= b->fitnessLS && a->complex <= b->complex) {
		if (a->complex < b->complex) {
			return true;
		}
		else if (a->fitnessLS < b->fitnessLS) {
			return true;
		}
	}
	return false;
};

bool Search::dominateTest(Subject* a, Subject* b) {
	if (a->fitnessTestLS <= b->fitnessTestLS && a->complex <= b->complex) {
		if (a->complex < b->complex) {
			return true;
		}
		else if (a->fitnessTestLS < b->fitnessTestLS) {
			return true;
		}
	}
	return false;
};

bool Search::dominateValid(Subject* a, Subject* b) {
	if (a->fitnessValidLS <= b->fitnessValidLS && a->complex <= b->complex) {
		if (a->complex < b->complex) {
			return true;
		}
		else if (a->fitnessValidLS < b->fitnessValidLS) {
			return true;
		}
	}
	return false;
};

void Search::calcRank(int size) {
	/**
	Juntar tudo
	selecionar 1s não dominados rank = 0
	aumentar rank dos outros

	**/

	int ranking = 0;

	for (int i = 0; i < size; i++)
		pop[i]->ranking = 0;

	bool stop = false;
	while (!stop) {
		stop = true;
		for (int i = 0; i < size; i++) // ele
		{
			if (pop[i]->ranking == ranking) // se tiver na hora dele
			{
				for (int j = 0; j < size; j++) // pra cada um candidato
				{
					if (i != j && (pop[j]->ranking == ranking)) // se não for ele e se tiver no mesmo ranking
					{
						// if (pop[i]->tree->infixPrint().compare(pop[j]->tree->infixPrint()) == 0) {
						if (rand() % 100 < 50 && pop[i]->fitnessLS == pop[j]->fitnessLS && pop[i]->complex == pop[j]->complex) {
							pop[j]->ranking = 10000000;
							pop[j]->complex = 10000000;
							pop[j]->fitnessLS = 10000000;
						}
						else if (dominate(pop[i], pop[j])) // se for dominado por j então ele aumenta o ranking
						{
							pop[j]->ranking = ranking + 1;
							stop = false;
						}



					}
				}
			}
		}
		ranking++;
	}
	sort(pop, pop + size, sortRank);

	int i, f, r;
	i = f = r = 0;
	while (f <= h_conf->popSize) {
		while (pop[f]->ranking == r)
			f++;

		sort(pop + i, pop + f, sortFit);
		crowdingDistanceFitness(i, f - 1);
		sort(pop + i, pop + f, sortComplexity);
		crowdingDistanceComplexity(i, f - 1);

		if (f > h_conf->popSize)
			sort(pop + i, pop + f, sortCrow);

		i = f;
		f++;
		r++;
	}
};

void Search::calcRankTest(int size) {
	/**
	Juntar tudo
	selecionar 1s não dominados rank = 0
	aumentar rank dos outros

	**/

	int ranking = 0;

	for (int i = 0; i < size; i++)
		pop[i]->ranking = 0;

	bool stop = false;
	while (!stop) {
		stop = true;
		for (int i = 0; i < size; i++) // ele
		{
			if (pop[i]->ranking == ranking) // se tiver na hora dele
			{
				for (int j = 0; j < size; j++) // pra cada um candidato
				{
					if (i != j && (pop[j]->ranking == ranking)) // se não for ele e se tiver no mesmo ranking
					{
						if (pop[i]->tree->infixPrint() == pop[j]->tree->infixPrint()) {
							pop[j]->ranking = 10000000;
							pop[j]->complex = 10000000;
							pop[j]->fitnessLS = 10000000;
						}
						else if (dominateTest(pop[i], pop[j])) // se for dominado por j então ele aumenta o ranking
						{
							pop[j]->ranking = ranking + 1;
							stop = false;
						}
					}
				}
			}
		}
		ranking++;
	}
	sort(pop, pop + size, sortRank);

	int i, f, r;
	i = f = r = 0;
	while (f <= h_conf->popSize) {
		while (f < h_conf->popSize && pop[f]->ranking == r)
			f++;

		sort(pop + i, pop + f, sortFitTest);
		crowdingDistanceFitnessTest(i, f - 1);
		sort(pop + i, pop + f, sortComplexity);
		crowdingDistanceComplexity(i, f - 1);

		//desordena o rank se o tamanho dele for >= conf->popSize
		if (f > h_conf->popSize)
			sort(pop + i, pop + f - 1, sortCrow);

		i = f;
		f++;
		r++;
	}
};

void Search::calcRankValid(int size) {
	/**
	Juntar tudo
	selecionar 1s não dominados rank = 0
	aumentar rank dos outros

	**/

	int ranking = 0;

	for (int i = 0; i < size; i++)
		pop[i]->ranking = 0;

	bool stop = false;
	while (!stop) {
		stop = true;
		for (int i = 0; i < size; i++) // ele
		{
			if (pop[i]->ranking == ranking) // se tiver na hora dele
			{
				for (int j = 0; j < size; j++) // pra cada um candidato
				{
					if (i != j && (pop[j]->ranking == ranking)) // se não for ele e se tiver no mesmo ranking
					{
						if (pop[i]->tree->infixPrint() == pop[j]->tree->infixPrint()) {
							pop[j]->ranking = 10000000;
							pop[j]->complex = 10000000;
							pop[j]->fitnessLS = 10000000;
						}
						else if (dominateValid(pop[i], pop[j])) // se for dominado por j então ele aumenta o ranking
						{
							pop[j]->ranking = ranking + 1;
							stop = false;
						}
					}
				}
			}
		}
		ranking++;
	}
	sort(pop, pop + size, sortRank);

	int i, f, r;
	i = f = r = 0;
	while (f <= h_conf->popSize) {
		while (f < h_conf->popSize && pop[f]->ranking == r)
			f++;

		sort(pop + i, pop + f, sortFitValid);
		crowdingDistanceFitnessValid(i, f - 1);
		sort(pop + i, pop + f, sortComplexity);
		crowdingDistanceComplexity(i, f - 1);

		if (f > h_conf->popSize)
			sort(pop + i, pop + f - 1, sortCrow);

		i = f;
		f++;
		r++;
	}
};

void Search::crowdingDistanceFitness(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = MAX;
		for (int j = i + 1; j < f; j++) {
			pop[j]->crowdingDistance = (pop[j + 1]->fitnessLS - pop[j - 1]->fitnessLS) / (pop[f]->fitnessLS - pop[i]->fitnessLS);
		}
	}
}

void Search::crowdingDistanceFitnessTest(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = MAX;
		for (int j = i + 1; j < f; j++) {
			pop[j]->crowdingDistance = (pop[j + 1]->fitnessTestLS - pop[j - 1]->fitnessTestLS) / (pop[f]->fitnessTestLS - pop[i]->fitnessTestLS);
		}
	}
}

void Search::crowdingDistanceFitnessValid(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = INFINITY;
		for (int j = i + 1; j < f; j++) {
			pop[j]->crowdingDistance = (pop[j + 1]->fitnessValidLS - pop[j - 1]->fitnessValidLS) / (pop[f]->fitnessValidLS - pop[i]->fitnessValidLS);
		}
	}
}

void Search::crowdingDistanceHigh(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = MAX;
		for (int j = i + 1; j < f; j++) {
			double num = (pop[j + 1]->tree->high - pop[j - 1]->tree->high);
			double den = (pop[f]->tree->high - pop[i]->tree->high);
			if (den == 0)
				pop[j]->crowdingDistance += 0;
			else
				pop[j]->crowdingDistance += num / den;
		}
	}
}

void Search::crowdingDistanceComplexity(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = MAX;
		for (int j = i + 1; j < f; j++) {
			double num = (pop[j + 1]->complex - pop[j - 1]->complex);
			double den = (pop[f]->complex - pop[i]->complex);
			if (den == 0)
				pop[j]->crowdingDistance += 0;
			else
				pop[j]->crowdingDistance += num / den;
		}
	}
}

void Search::crowdingDistanceSize(int i, int f) {
	if (i == f) {
		pop[i]->crowdingDistance = MAX;
	}
	else {
		pop[i]->crowdingDistance = pop[f]->crowdingDistance = MAX;
		for (int j = i + 1; j < f; j++) {
			double num = (pop[j + 1]->tree->terminals - pop[j - 1]->tree->terminals);
			double den = (pop[f]->tree->terminals - pop[i]->tree->terminals);
			if (den == 0)
				pop[j]->crowdingDistance += 0;
			else
				pop[j]->crowdingDistance += num / den;
		}
	}
}

void Search::calcFitnessLS(Subject* s) {
	//double fitness = 0;
	int lines = banco_dados->trainCount;
	// int columns = s->tree->subCounter;
	int dimension = 0;
	//   double** subMat;
	double* alpha = NULL;

	s->treino_vp = s->treino_fp = s->treino_vn = s->treino_fn = 0;

	//classification

	//double total = 0;

	for (int i = 0; i < lines; i++) {
		double res = (s->tree->treeResult(banco_dados->values[banco_dados->training[i]], alpha, dimension));
		//        if (data->results[data->training[i]] == 0) {
		//            total += conf->peso0;
		//        } else {
		//            total += conf->peso1;
		//        }
		if (res != banco_dados->results[banco_dados->training[i]]) {
			if (banco_dados->results[banco_dados->training[i]] == 0) {
				//fitness += conf->peso0;
				s->treino_fp++;
			}
			else {
				//fitness += conf->peso1;
				s->treino_fn++;
			}
		}
		else {
			if (banco_dados->results[banco_dados->training[i]] == 0) {
				s->treino_vn++;
			}
			else {
				s->treino_vp++;
			}
		}
	}

	s->fitnessLS = ((double)(s->treino_fn + s->treino_fp) / lines) * 100;
	s->complex = s->complexity();

	//    //regression
	//    if(conf->leastSquare == 1){
	//        dimension = columns;
	//        subMat = new double*[lines];
	//        for(int i = 0; i < lines; i++)
	//            subMat[i] = new double[columns];
	//
	//        //Avalia cada subExp
	//        for(int i = 0; i < columns; i++)
	//        {
	//            double* exp = &s->tree->sub.at(i).at(0);
	//            int sizeExp = s->tree->sub.at(i).size();
	//            for(int j = 0; j < lines; j++)
	//            {
	//                double* a = new double[sizeExp];
	//                double* var = data->values[data->training[j]];
	//                for(int k = 0; k < sizeExp; k += 2)
	//                {
	//                    a[k] = exp[k];
	//                    a[k + 1] = exp[k + 1];
	//                    if(a[k] == 3.0)
	//                    {
	//                        a[k] = 0.0;
	//                        a[k + 1] = var[(int)a[k + 1]];
	//                    }
	//                    if(a[k] == 5.0)
	//                    {
	//                        a[k] = 0.0;
	//                    }
	//                }
	//
	//                subMat[j][i] = Avalia(a, sizeExp);
	//                delete [] a;
	//            }
	//        }
	//
	//
	//        QRDecomposition* qrDec = new QRDecomposition(subMat, lines, columns);
	//        double* b = new double[lines];
	//        for(int i = 0; i < lines; i++)
	//        {
	//            b[i] = data->results[data->training[i]];
	//        }
	////        se tem algo lá, apagar
	//        alpha = s->tree->alpha;
	//        if(alpha != NULL)
	//            delete [] alpha;
	////        novos coeficientes
	//        alpha = s->tree->alpha = qrDec->solveLeastSquares(b, lines);
	//
	//        if(alpha == NULL){
	//            dimension = 0;
	//        }
	//
	//        for(int i = 0; i < lines; i++)
	//            delete [] subMat[i];
	//        delete [] subMat;
	//        delete [] b;
	//        delete qrDec;
	//    }
	//
	//    for(int i = 0; i < lines ; i++){
	//        double res = (s->tree->treeResult(data->values[data->training[i]], alpha, dimension));
	//
	////        for(int j = 0; j < data->countVar; j++)
	////            cout << data->values[data->training[i]][j] << " ";
	////        if(s->tree->infixPrint() == "(a + ((c * b) + (d / e)))"){
	////            cout << data->results[data->training[i]] << " " << res << " " << data->results[data->training[i]] - res << " " << pow(data->results[data->training[i]] - res, 2) << endl;
	////            cin.get();
	////        }
	//
	//        res = pow(data->results[data->training[i]] - res, 2);
	//
	//        if(res == INFINITY){
	//            fitness = INFINITY;
	//            break;
	//        }
	//        if(res == NAN){
	//            fitness = INFINITY;
	//            break;
	//        }
	//        fitness += res;
	//    }
	//    s->fitnessLS = fitness/lines;
};

void Search::calcFitnessTestLS(Subject* s) {
	//double fitness = 0;
	int lines = banco_dados->testCount;
	int columns = s->tree->subCounter;
	double* alpha = s->tree->alpha;
	int dimension = columns;
	if (h_conf->leastSquare == 0 || alpha == NULL) {
		dimension = 0;
	}

	s->teste_vp = s->teste_fp = s->teste_vn = s->teste_fn = 0;

	//classification
	//double total = 0;

	for (int i = 0; i < lines; i++) {
		double res = (s->tree->treeResult(banco_dados->values[banco_dados->test[i]], alpha, dimension));
		//        if (data->results[data->training[i]] == 0) {
		//            total += conf->peso0;
		//        } else {
		//            total += conf->peso1;
		//        }
		if (res != banco_dados->results[banco_dados->test[i]]) {
			if (banco_dados->results[banco_dados->test[i]] == 0) {
				//fitness += conf->peso0;
				s->teste_fp++;
			}
			else {
				//fitness += conf->peso1;
				s->teste_fn++;
			}
		}
		else {
			if (banco_dados->results[banco_dados->test[i]] == 0) {
				s->teste_vn++;
			}
			else {
				s->teste_vp++;
			}
		}
	}

	s->fitnessTestLS = ((double)(s->teste_fn + s->teste_fp) / lines) * 100;

	//    //regression
	//    for(int i = 0; i < lines ; i++)
	//    {
	//        double res = (s->tree->treeResult(data->values[data->test[i]], alpha, dimension));
	//        res = pow(data->results[data->test[i]] - res, 2);
	//        if(res == INFINITY)
	//        {
	//            fitness = INFINITY;
	//            break;
	//        }
	//        if(res == NAN)
	//        {
	//            fitness = INFINITY;
	//            break;
	//        }
	//        fitness += res;
	//    }
	//    s->fitnessTestLS = fitness/lines;
};

void Search::calcFitnessValidLS(Subject* s) {
	//double fitness = 0;
	int lines = banco_dados->validCount;
	int columns = s->tree->subCounter;
	double* alpha = s->tree->alpha;
	int dimension = columns;
	if (h_conf->leastSquare == 0 || alpha == NULL) {
		dimension = 0;
	}

	s->valida_vp = s->valida_fp = s->valida_vn = s->valida_fn = 0;

	//classification
	//double total = 0;

	for (int i = 0; i < lines; i++) {
		double res = (s->tree->treeResult(banco_dados->values[banco_dados->validation[i]], alpha, dimension));
		//        if (data->results[data->training[i]] == 0) {
		//            total += conf->peso0;
		//        } else {
		//            total += conf->peso1;
		//        }
		if (res != banco_dados->results[banco_dados->validation[i]]) {
			if (banco_dados->results[banco_dados->validation[i]] == 0) {
				//fitness += conf->peso0;
				s->valida_fp++;
			}
			else {
				//fitness += conf->peso1;
				s->valida_fn++;
			}
		}
		else {
			if (banco_dados->results[banco_dados->validation[i]] == 0) {
				s->valida_vn++;
			}
			else {
				s->valida_vp++;
			}
		}
	}

	s->fitnessValidLS = ((double)(s->valida_fn + s->valida_fp) / lines) * 100;
	//      //regression
	//    for(int i = 0; i < lines ; i++)
	//    {
	//        double res = (s->tree->treeResult(data->values[data->validation[i]], alpha, dimension));
	//        res = pow(data->results[data->validation[i]] - res, 2);
	//        if(res == INFINITY)
	//        {
	//            fitness = INFINITY;
	//            break;
	//        }
	//        if(res == NAN)
	//        {
	//            fitness = INFINITY;
	//            break;
	//        }
	//        fitness += res;
	//    }
	//    s->fitnessValidLS = fitness/lines;
};

void Search::Replace() {
	for (int i = h_conf->popSize * 2 - 1; i >= h_conf->popSize; i--)
		delete pop[i];
};

Search::~Search() {
	delete[] pop;
};
