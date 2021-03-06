#include "stdafx.h"

#include<math.h>
#include<time.h>
#include<stdlib.h> 
#include<chrono>
#include<thread>
#include<algorithm>
#include<numeric>
#include<vector>
#include<iostream>
#include<future>
#include<random>
#include<fstream>
#include<queue>
#include<string>

#define INT int

using namespace std;

#include "generalfunctions.h"
#include "threadpool.h"
#include "sequentialfunction.h"
#include "parallelfunction.h"
#include "poolfunction.h"

//used by function testFunctionForRandomInput
int MAX_TEST_SIZE = 1000000;
int MIN_TEST_SIZE = 1000;

//how many times to call function.
//used by functions runSequential and runParallel. 
int TEST_TRIES = 10;

//which functions to test
bool testInitializeFunction = false;
bool testValueFunction = false;
bool testGradFunction = false;
bool testValgradFunction = true;

//which methods correctness should program test
bool TEST_PARALLEL = true;
bool TEST_POOL = true;

//whether program should print all errors or only one
bool printAllErrors = true;

//use for parameter to double_equals function if default epsilon is not suitable
double EPSILON = 0.00000000001;

//Output format
bool OUTPUT_HIGH_PRECISION = false;
bool SHOW_TIMING = true;
bool SHOW_EXACT_TIME = true;

//used to compare results to the correct solution
bool double_equals(double a, double b, double epsilon = 0.00000000001)
{
	if (abs(a - b) < epsilon) {
		return true;
	}
	return abs((a / b) - 1) < epsilon;
}

/*
This function runs original, sequential version of the function from file "oldfunction.h"
x - parameters of the function
g - gradients
n - number of parameters
sum - sum of function's results
init - if set to false "StartingGuess" method will not be called
value - if set to false "value" method will not be called
gradient - if set to false "grad" method will not be called
valgrad - if set to false "valgrad" method will not be called
*/
double runSequential(double* x, double* g, int n, double& sum, bool init, bool value, bool gradient, bool valgrad) {
	sum = 0;
	if (init == true) {
		sequential::initialize(x, n);
	}
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < TEST_TRIES; i++) {
		if (valgrad == true) {
			sum += sequential::getValGrad(x, g, n);
		}
		/*
		if (gradient == true) {
			sequential::EG2grad(g, x, n);
		}
		if (value == true) {
			sum += sequential::EG2value(x, n);
		}
		*/
	}
	auto finish = std::chrono::high_resolution_clock::now();
	return chrono::duration<double>(finish - start).count();
}

/*
This function runs new, parallel version of the function from file "newfunction.h"
x - parameters of the function
g - gradients
n - number of parameters
sum - sum of function's results
init - if set to false "StartingGuess" method will not be called
value - if set to false "value" method will not be called
gradient - if set to false "grad" method will not be called
valgrad - if set to false "valgrad" method will not be called
*/
double runParallel(double* x, double* g, int n, double& sum, bool init, bool value, bool gradient, bool valgrad) {
	sum = 0;
	if (init == true) {
		parallel::initialize(x, n);
	}
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < TEST_TRIES; i++) {
		if (valgrad == true) {
			sum += parallel::getValGrad(x, g, n);
		}
		/*
		if (gradient == true) {
			parallel::EG2grad(g, x, n);
		}
		if (value == true) {
			sum += parallel::EG2value(x, n);
		}
		*/
	}
	auto finish = std::chrono::high_resolution_clock::now();
	return chrono::duration<double>(finish - start).count();
}

/*
This function runs new, parallel version of the function from file "newfunction.h" using pools
x - parameters of the function
g - gradients
n - number of parameters
sum - sum of function's results
init - if set to false "StartingGuess" method will not be called
value - if set to false "value" method will not be called
gradient - if set to false "grad" method will not be called
valgrad - if set to false "valgrad" method will not be called
*/
double runPool(double* x, double* g, int n, double& sum, bool init, bool value, bool gradient, bool valgrad) {
	sum = 0;
	if (init == true) {
		pool::initialize(x, n);
	}
	threadPool = new ThreadPool(3);
	auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < TEST_TRIES; i++) {
		if (valgrad == true) {
			sum += pool::getValGrad(x, g, n);
		}
	}
	auto finish = std::chrono::high_resolution_clock::now();
	delete threadPool;
	return chrono::duration<double>(finish - start).count();
}

void checkCorrectness(double firstSum, double secondSum, double* firstG, double* secondG, double* firstX, double* secondX, int size, string firstName, string secondName) {
	if (!double_equals(firstSum, secondSum, EPSILON)) {
		cout << "Value ERROR!!! (from valgrad, " + firstName + ") " << endl << firstSum << " " << endl << secondSum << endl;
	}
	for (int i = 0; i < size; i++) {
		if (!double_equals(firstG[i], secondG[i])) {
			cout << "Gradient ERROR!!! (from valgrad, " + firstName + ") " << i << " " << endl << firstG[i] << " " << endl << secondG[i] << endl;
			if (!printAllErrors) {
				break;
			}
		}
		if (!double_equals(firstX[i], secondX[i])) {
			cout << "X ERROR!!! (from valgrad, " + firstName + ") " << i << " " << endl << firstX[i] << " " << endl << secondX[i] << endl;
			if (!printAllErrors) {
				break;
			}
		}
	}
}

//runs sequential and parallel versions of function, compares their results and time efficiency
void testFunctionForRandomInput() {
	random_device rd;
	default_random_engine generator(rd());
	uniform_real_distribution<double> distribution(0, 1000);

	if (OUTPUT_HIGH_PRECISION == true) {
		cout.precision(17);
	}
	for (int size = MIN_TEST_SIZE; size <= MAX_TEST_SIZE; size = size * 10) {
		initializer(size, 100);
		double sequentialSum, parallelSum, poolSum;
		double *sequentialX = (double *)malloc(size * sizeof(double));
		double *sequentialG = (double *)malloc(size * sizeof(double));
		double *parallelX = (double *)malloc(size * sizeof(double));
		double *parallelG = (double *)malloc(size * sizeof(double));
		double *poolX = (double *)malloc(size * sizeof(double));
		double *poolG = (double *)malloc(size * sizeof(double));

		for (int i = 0; i < size; i++) {
			double number = distribution(generator);
			sequentialX[i] = number;
			parallelX[i] = number;
			poolX[i] = number;
		}

		double sequentialsTime = 0, parallelsTime = 0, poolTime = 0;
		sequentialsTime = runSequential(sequentialX, sequentialG, size, sequentialSum, testInitializeFunction, testValueFunction, testGradFunction, testValgradFunction);
		parallelsTime = runParallel(parallelX, parallelG, size, parallelSum, testInitializeFunction, testValueFunction, testGradFunction, testValgradFunction);
		poolTime = runPool(poolX, poolG, size, poolSum, testInitializeFunction, testValueFunction, testGradFunction, testValgradFunction);

		//cout << "parallelSum= " << parallelSum << endl;
		//cout << "sequentialSum= " << sequentialSum << endl;
		if (TEST_PARALLEL == true) {
			checkCorrectness(parallelSum, sequentialSum, parallelG, sequentialG, parallelX, sequentialX, size, "parallel", "sequential");
		}
		if (TEST_POOL == true) {
			checkCorrectness(poolSum, sequentialSum, poolG, sequentialG, poolX, sequentialX, size, "pool", "sequential");
		}
		/*
		runParallel(parallelX, parallelG, size, parallelSum, false, testValueFunction, testGradFunction, false);
		runSequential(sequentialX, sequentialG, size, sequentialSum, false, testValueFunction, testGradFunction, false);

		//cout << "SUM= " << parallelSum << endl;
		if (!double_equals(parallelSum, sequentialSum, EPSILON)) {
			cout << "Value ERROR!!! " << parallelSum << " " << sequentialSum << endl;
		}
		for (int i = 0; i < size; i++) {
			if (!double_equals(parallelG[i], sequentialG[i])) {
				cout << "Gradient ERROR!!! " << i << " " << parallelG[i] << " " << sequentialG[i] << endl;
				break;
			}
			if (!double_equals(parallelX[i], sequentialX[i])) {
				cout << "X ERROR!!! " << i << " " << parallelX[i] << " " << sequentialX[i] << endl;
				break;
			}
		}
		*/

		free(sequentialX);
		free(sequentialG);
		free(parallelX);
		free(parallelG);
		free(poolX);
		free(poolG);

		if (SHOW_TIMING == true) {
			if(SHOW_EXACT_TIME == true) {
			cout << "Time taken by sequential " << sequentialsTime << endl;
			cout << "Time taken by parallel " << parallelsTime << endl;
			cout << "Time taken by pool " << poolTime << endl;
			}
			cout << "Test size = " << size << endl;
			cout << "Parallel is " << sequentialsTime / parallelsTime << " times faster than sequential" << endl;
			cout << "Pooling is " << sequentialsTime / poolTime << " times faster than sequential" << endl;
		}
	}
}

//used to check if extremum found by parallel function is correct
void testValue() {
	int size;
	ifstream input;
	input.open("D:\\Desktop\\Phd\\New Test\\CG_DESCENT-6.8_EXP\\CG_DESCENT-6.8_EXP\\input.txt");
	input >> size;
	cout << size << endl;
	double *sequentialX = (double *)malloc(size * sizeof(double));
	double *sequentialG = (double *)malloc(size * sizeof(double));

	for (int i = 0; i < size; i++) {
		input >> sequentialX[i];
	}
	cout << sequentialX[0] << endl;
	cout << sequentialX[1] << endl;
	cout << sequentialX[2] << endl;

	//runSequential(sequentialX, sequentialG, sequentialSum, false);
	//cout << "SUM= " << sequential::EG2value(sequentialX, size) << endl;

	free(sequentialX);
	free(sequentialG);
}

int main(void)
{
	cout << "Number of hardware concurrency: " << HARDWARE_THREADS << endl;
	testFunctionForRandomInput();
	//testValue();
}