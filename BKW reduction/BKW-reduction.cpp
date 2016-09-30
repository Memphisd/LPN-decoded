
#include "BKW-reduction.h"

void show_bin(uint64_t num)
{
	uint64_t bit = 1ULL << 63;
	for (int i = 0; i <64;i++)
	{
		printf("%d", (bit&num) >> (63 - i));
		bit = bit >> 1;
	}
}

void l_shift(uint64_t * v, int shift_amount, int regs)
{
	//handle shifts greater then 64 via swaps of whole uints
	if (shift_amount > 63)
	{
		for (int i = 0;i + shift_amount / 64 < regs;++i)
		{
			v[i] = v[i + shift_amount / 64];
		}
		for (int i = regs - shift_amount / 64;i <regs;i++)
			v[i] = 0;
		shift_amount = shift_amount % 64;
	}

	//do shifts lower then 64 via carrys
	if (shift_amount != 0)
	{
		uint64_t zw;
		uint64_t bitmask = ((1ULL << shift_amount) - 1) << (64 - shift_amount);
		for (int i = 0;i<regs - 1;i++)
		{
			zw = v[i + 1] & bitmask;
			v[i] <<= shift_amount;
			v[i] ^= (zw >> (64 - shift_amount));
		}
		v[regs - 1] <<= shift_amount;
	}
}
uint64_t **bkw_list;// = new uint64_t *[start_amount];
int iteration = 0;

uint64_t get_label(uint64_t * g)
{

	uint64_t zw[int_per_row];
	for (int i = 0;i < int_per_row;++i)
		zw[i] = g[i];
	l_shift(zw, iteration*b, int_per_row);

	uint64_t bitmask = ((1ULL << b) - 1) << (64 - b);
	uint64_t label = (zw[0] & bitmask) >> (64 - b);

	return label;
}


void init_list()
{

	uint64_t **first = (uint64_t **)malloc(sizeof *first * start_amount);
	uint64_t *seconds = (uint64_t *)malloc(sizeof *seconds * start_amount * int_per_row);

	//for (int j = 0;j<start_amount/(1<<32))
#pragma omp parallel for
	for (uint64_t i = 0; i < start_amount; i++)
		first[i] = seconds + int_per_row * i;

	bkw_list = first;
	uint64_t * g;
#pragma omp parallel for
	for (uint64_t i = 0;i < start_amount;++i)
	{
		//bkw_list[i] = new uint64_t[int_per_row];
		g = query(omp_get_thread_num());
		for (int j = 0;j < int_per_row;++j)
			bkw_list[i][j] = g[j];

	}
}
bool sort_operator(uint64_t * g1, uint64_t * g2)
{
	return get_label(g1) > get_label(g2);
}

bool sort_operator_whole(uint64_t * g1, uint64_t * g2)
{
	for (int i = 0;i<int_per_row;++i)
		if (g1[i]>g2[i])
			return true;
		else if (g1[i] < g2[i])
			return false;
		else
			continue;
		return false;
}
void print_list(uint64_t size)
{
	for (int i = 0;i < size;++i)
	{
		for (int j = 0;j < int_per_row;++j)
			show_bin(bkw_list[i][j]);
		printf("\n");
	}
}
void add(uint64_t * g1, uint64_t * g2)
{
	for (int i = 0;i < int_per_row;++i)
	{
		g1[i] ^= g2[i];
	}
}

bool is_zero(uint64_t * g)
{
	for (int i = 0;i < (n + 63) / 64;++i)
		if (g[i] != 0)
			return false;
	return true;
}
int eliminate = 0;
bool zero = false;

uint64_t bkw_reduction()
{
	uint64_t size = start_amount;
	uint64_t new_size;
	for (int i = 0;i < a;++i)
	{
		iteration = i;
		__gnu_parallel::sort(bkw_list, bkw_list + size, sort_operator);

#pragma omp parallel for
		for (int t = 0;t < threads;++t)
		{
			uint64_t index = t*(size / threads);
			uint64_t stop = t == threads - 1 ? size : (t+1)*(size / threads);
			//find first label for this thread
			while (index < stop && index != 0 && get_label(bkw_list[index - 1]) == get_label(bkw_list[index]))
				index++;

			while (index < stop)
			{
				uint64_t act_label = get_label(bkw_list[index]);
				uint64_t iter = index + 1;
				while (iter < size && act_label == get_label(bkw_list[iter]))
				{
					add(bkw_list[iter], bkw_list[index]);
					++iter;
				}
				add(bkw_list[index], bkw_list[index]);
				index = iter;

			}
		}

		__gnu_parallel::sort(bkw_list, bkw_list + size, sort_operator_whole);
		for (new_size = size - 1;new_size > 0 && is_zero(bkw_list[new_size]);--new_size)
		{
		}
		size = new_size;
		if (size == 0)
		{
			zero = true;
			break;
		}

		printf("\n%d/%d Run finished\n", i + 1, a);


	}

	if (zero)
		printf("no vector survived");
	else
	{
		//print_list(size);
		printf("\n%d vectors suvived", size);
	}
	return size;
}

int main()
{
	std::ofstream out("out.txt", std::ofstream::app);
	std::streambuf *coutbuf = std::cout.rdbuf();
	std::cout.rdbuf(out.rdbuf());

	omp_set_num_threads(threads);
	initialize_Oracle(threads);
	init_list();
	printf("\n\n\n");

	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	uint64_t amount = bkw_reduction();
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	std::cout << "\n\n__________________________________\n";
	std::cout << "Duration (seconds): " << ((double)duration) / 1000 << "\n";
	std::cout << "Elements survived: " << amount << "\n";
	std::cout << "__________________________________\n";
	std::cout.rdbuf(coutbuf);
}