
#include "LPN_Oracle.h"

time_t t;
std::random_device rd;
std::mt19937_64 mt(static_cast<uint64_t>(rd()));
std::default_random_engine gen_ber(rd());
std::uniform_int_distribution<uint64_t> uint_dist;

std::vector<std::bernoulli_distribution> ber_dists;
std::vector<std::mt19937_64> random_engines;
std::vector<std::uniform_int_distribution<uint64_t>> dists;




std::bernoulli_distribution ber(tau);
uint64_t * secret = new uint64_t[int_per_row];
uint64_t ** res, ** vec;
int arr_for_secret = k/64;
extern int w;
int max_errors = 2000;
int errors = 0;
int querys = 0;
extern uint64_t * right_error;
//extern int k;

uint64_t* query(int thread,bool gen)
{
	if (gen && querys == 0)
		errors = 0;	

	for (int ipr = 0;ipr < int_per_row;ipr++)
		vec[thread][ipr] = dists[thread](random_engines[thread]);

	if((k%64) != 0)
		vec[thread][int_per_row-1]<<=(63 -((k-1)%64));

	bool l = scalar_secret(vec[thread],res[thread]);

	if (ber_dists[thread](random_engines[thread])&& errors<max_errors)
	{
		l ^= 0x1;
		if (gen)
		{
			right_error[querys / 64] ^= (1ULL << (63 - (querys % 64)));
			errors++;
		}
	}
	if (gen)
	{
		querys++;
		querys = querys %n;
	}
	
	

	if(arr_for_secret==int_per_row && k%64==0)
		vec[thread][int_per_row]=0;
	if(l)
		vec[thread][arr_for_secret]^=(1ULL<<(63-(k%64)));

	return vec[thread];
}



uint64_t * initialize_Oracle(int threads)
{

	vec = new uint64_t*[threads];
	if((k%64)!=0)
		for(int i =0;i<threads;i++)
			vec[i]=new uint64_t[int_per_row];
	else
		for(int i =0;i<threads;i++)
			vec[i]=new uint64_t[int_per_row+1];
	
	res = new uint64_t*[threads];
	for(int i = 0;i<threads;i++)
		res[i]=new uint64_t[int_per_row];

	for(int ipr=0;ipr<int_per_row;ipr++)
		secret[ipr]=uint_dist(mt);

	if((k%64) != 0)
		secret[int_per_row-1]<<=(63-((k-1)%64));

	for (int i = 0;i < threads;i++)
	{
		std::mt19937_64 mtt((rd()));
		std::uniform_int_distribution<uint64_t> uint_distri;
		std::bernoulli_distribution bern(tau);
		ber_dists.push_back(bern);
		random_engines.push_back(mtt);
		dists.push_back(uint_distri);
	}
	return secret;
}

bool scalar_secret(uint64_t * a,uint64_t * res)
{
	
	bool ret=0;
	uint64_t xor_sum=0;
	int c;
	for(int ipr=0;ipr<int_per_row;ipr++)
	{
		res[ipr]=(secret[ipr]&a[ipr]);
		xor_sum^=res[ipr];
	}
	for (c = 0;xor_sum;c++)
		xor_sum &= (xor_sum - 1);
	ret = c % 2;

	return ret;
}
