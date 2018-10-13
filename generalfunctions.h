const long HARDWARE_THREADS = thread::hardware_concurrency();
long num_threads;
long block_size;

/*
Accumulation with begin and end indices
begin - index where accumulation starts
end - index where accumulation ends
m - index shift, used by function
n - length of the array
*/
double valgradParallelRangedWithMWithN
(
	double(*functionToCompute)(double *x, double *g, INT begin, INT end, INT m, INT placeholder),
	double *x,
	double *g,
	INT begin,
	INT end,
	INT m,
	INT n
)
{
	vector<future<double>> futures(num_threads - 1);
	INT block_start = begin, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, g, block_start, block_end, m, n);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, block_start, end, m, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

/*
Accumulation with begin and end indices
begin - index where accumulation starts
end - index where accumulation ends
m - index shift, used by function
placeholder - just added temporarly because of wierd bug with async
*/
double valgradParallelRangedWithM
(
	double(*functionToCompute)(double *x, double *g, INT begin, INT end, INT m, INT placeholder),
	double *x,
	double *g,
	INT begin,
	INT end,
	INT m
)
{
	vector<future<double>> futures(num_threads - 1);
	INT block_start = begin, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, g, block_start, block_end, m, 0);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, block_start, end, m, 0);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

/*
In case when all threads must write gradient into one index
index - indicates index of congested gradient
n - number of items to sum, not number of elements in array
*/
double valgradParallelOneCongestion
(
	double(*functionToCompute)(double *x, double *g, double* sum, INT index, INT begin, INT end),
	double *x,
	double *g,
	INT n,
	INT index
)
{
	vector<double> sums(num_threads, 0.0);
	vector<future<double>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, g, &sums[i], index, block_start, block_end);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, &sums[num_threads - 1], index, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	for (i = 0; i < num_threads; i++) {
		g[index] += sums[i];
	}
	return sum;
}


/*
In case when all threads must write gradient into one index
index - indicates index of congested gradient
n - number of items to sum, not number of elements in array
*/
void gradParallelOneCongestion
(
	void(*functionToCompute)(double *g, double *x, double *sum, INT index, INT begin, INT end),
	double *g,
	double *x,
	INT n,
	INT index
)
{
	vector<double> sums(num_threads, 0.0);
	vector<future<void>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, g, x, &sums[i], index, block_start, block_end);
		block_start = block_end;
	}
	functionToCompute(g, x, &sums[i], index, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	for (i = 0; i < num_threads; i++) {
		g[index] += sums[i];
	}
}

/*
In case when ranged functions need to know some general index that is out of range
n - number of items to sum, not number of elements in array
*/
double valueParallelOneCongestion
(
	double(*functionToCompute)(double *x, INT index, INT begin, INT end),
	double *x,
	INT n,
	INT index
)
{
	vector<future<double>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, index, block_start, block_end);
		block_start = block_end;
	}
	double sum = functionToCompute(x, index, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

/*
For simplest case of accumulation functions
n - number of items to sum, not number of elements in array
*/
double valgradParallel
(
	double(*functionToCompute)(double *x, double *g, INT begin, INT end),
	double *x,
	double *g,
	INT n
)
{
	vector<future<double>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, g, block_start, block_end);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

/*
For simplest case of accumulation functions
n - number of items to sum, not number of elements in array
*/
void gradParallel
(
	void(*functionToCompute)(double *g, double *x, INT begin, INT end),
	double *g,
	double *x,
	INT n
)
{
	vector<future<void>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, g, x, block_start, block_end);
		block_start = block_end;
	}
	functionToCompute(g, x, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
}

/*
For simplest case of accumulation functions
n - number of items to sum, not number of elements in array
*/
double valueParallel
(
	double(*functionToCompute)(double *x, INT begin, INT end),
	double *x,
	INT n
)
{
	vector<future<double>> futures(num_threads - 1);
	INT block_start = 0, i;
	for (i = 0; i < (num_threads - 1); ++i)
	{
		INT block_end = block_start + block_size;
		futures[i] = async(functionToCompute, x, block_start, block_end);
		block_start = block_end;
	}
	double sum = functionToCompute(x, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

//used to initialize block_size
void initializer
(
	INT elements,
	INT minPerThread
)
{
	const long max_threads = (elements + minPerThread - 1) / minPerThread;
	num_threads = std::min(HARDWARE_THREADS != 0 ? HARDWARE_THREADS : 2, max_threads);
	//num_threads = 2;
	//cout << num_threads << endl;
	block_size = elements / num_threads;

}