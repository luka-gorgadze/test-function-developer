class ThreadPool {
public:
	ThreadPool(size_t);
	template<class F, class... Args>
	auto enqueue(F&& f, Args&&... args)
		->std::future<typename std::result_of<F(Args...)>::type>;
	~ThreadPool();
private:
	// need to keep track of threads so we can join them
	std::vector< std::thread > workers;
	// the task queue
	std::queue< std::function<void()> > tasks;

	// synchronization
	std::mutex queue_mutex;
	std::condition_variable condition;
	bool stop;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
	: stop(false)
{
	for (size_t i = 0; i < threads; ++i)
		workers.emplace_back(
			[this]
	{
		for (;;)
		{
			std::function<void()> task;

			{
				std::unique_lock<std::mutex> lock(this->queue_mutex);
				this->condition.wait(lock,
					[this] { return this->stop || !this->tasks.empty(); });
				if (this->stop && this->tasks.empty())
					return;
				task = std::move(this->tasks.front());
				this->tasks.pop();
			}

			task();
		}
	}
	);
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
-> std::future<typename std::result_of<F(Args...)>::type>
{
	using return_type = typename std::result_of<F(Args...)>::type;

	auto task = std::make_shared< std::packaged_task<return_type()> >(
		std::bind(std::forward<F>(f), std::forward<Args>(args)...)
		);

	std::future<return_type> res = task->get_future();
	{
		std::unique_lock<std::mutex> lock(queue_mutex);

		// don't allow enqueueing after stopping the pool
		if (stop)
			throw std::runtime_error("enqueue on stopped ThreadPool");

		tasks.emplace([task]() { (*task)(); });
	}
	condition.notify_one();
	return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
	{
		std::unique_lock<std::mutex> lock(queue_mutex);
		stop = true;
	}
	condition.notify_all();
	for (std::thread &worker : workers)
		worker.join();
}

ThreadPool *threadPool;

/*
In case when all threads must write gradient into one index
index - indicates index of congested gradient
n - number of items to sum, not number of elements in array
*/
double valgradPoolOneCongestion
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
		futures[i] = threadPool->enqueue(functionToCompute, x, g, &sums[i], index, block_start, block_end);
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
Accumulation with begin and end indices
begin - index where accumulation starts
end - index where accumulation ends
m - index shift, used by function
n - length of the array
*/
double valgradPoolRangedWithMWithN
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
		futures[i] = threadPool->enqueue(functionToCompute, x, g, block_start, block_end, m, n);
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
double valgradPoolRangedWithM
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
		futures[i] = threadPool->enqueue(functionToCompute, x, g, block_start, block_end, m, 0);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, block_start, end, m, 0);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}

/*
For simplest case of accumulation functions
n - number of items to sum, not number of elements in array
*/
double valgradPool
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
		futures[i] = threadPool->enqueue(functionToCompute, x, g, block_start, block_end);
		block_start = block_end;
	}
	double sum = functionToCompute(x, g, block_start, n);
	for (i = 0; i < futures.size(); i++) {
		sum += futures[i].get();
	}
	return sum;
}
