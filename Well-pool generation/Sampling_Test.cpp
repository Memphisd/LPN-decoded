#include "Sampling_Test.h"


uint64_t * right_secret;
extern std::vector<std::mt19937_64> random_engines;
extern std::vector<std::uniform_int_distribution<uint64_t>> dists;


uint64_t **results;

void init_list()
{

	uint64_t **first = (uint64_t **)malloc(sizeof *first * samp_amount);
	uint64_t *seconds = (uint64_t *)malloc(sizeof *seconds * samp_amount * int_per_row);


#pragma omp parallel for
	for (uint64_t i = 0; i < samp_amount; i++)
		first[i] = seconds + int_per_row * i;

	results = first;
}
void test_samples(int thread)
{
	uint64_t * g;
	uint64_t * save = new uint64_t[int_per_row];
	int act_arr = n % 64 == 0 ? int_per_row : int_per_row - 1;
	bool secret;
	bool zero = true;

	int counter_q = 0;
	uint64_t samp_amount_thread = samp_amount / threads;
	uint64_t index = samp_amount_thread*thread;
	while(counter_q<samp_amount_thread)
	{
		g = query(thread);
		secret = ((g[act_arr] >> (63 - (n % 64))) & 0x1);
		if (secret)
			g[act_arr] ^= (1ULL << (63 - (n % 64)));
		if ((n-x)%64!=0 && !(g[(n - x) / 64] & (1ULL << (64 - ((n - x) % 64)))-1) == 0)
			zero = false;
		for (int i = (n-x) % 64 != 0 ? (n - x) / 64 + 1 : (n - x) / 64;i < (n + 63) / 64 && zero;++i)
			if (g[i] != 0)
				zero = false;
		if (zero)
		{
			for (int i = 0;i < int_per_row;i++)
				results[index+counter_q][i] = g[i];
			counter_q++;
		}
		if (counter_q >= samp_amount)
			break;
		zero = true;
	}
}



int main(int argc, char *argv[])
{

	std::ofstream out(out_file, std::ofstream::app);
	std::streambuf *coutbuf = std::cout.rdbuf();
	std::cout.rdbuf(out.rdbuf());

	std::cout << "__________________________";

	uint64_t * det_secret = NULL;
	right_secret = initialize_Oracle(threads);

	omp_set_num_threads(threads);
	init_list();

	//start the test
	printf("\nstarted");
	fflush(stdout);

	const clock_t START = clock();
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
	for (int i = 0;i < threads;++i)
		test_samples(i);

	const double T_ELAPSED = (double)(clock() - START) / CLOCKS_PER_SEC;
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();


	std::cout << "\nDuration (seconds): " << ((double)duration);
	std::cout << "\nsampling took (cpu-time):" << T_ELAPSED << "\n";
	printf("\nfinished");

	std::cout.rdbuf(coutbuf);
	return 0;


}
