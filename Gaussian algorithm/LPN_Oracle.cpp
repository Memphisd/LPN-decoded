
#include "LPN_Oracle.h"

time_t t;
std::random_device rd;
std::mt19937_64 mt(static_cast<uint64_t>(rd()));
std::default_random_engine gen_ber(rd());
std::uniform_int_distribution<uint64_t> uint_dist;

std::vector<std::bernoulli_distribution> ber_dists;
std::vector<std::mt19937_64> random_engines;
std::vector<std::uniform_int_distribution<uint64_t>> dists;


std::bernoulli_distribution ber(p);
uint64_t * secret = new uint64_t[int_per_row];
uint64_t ** vec;
int arr_for_secret = n/64;

uint64_t* query(int thread)
{
	for (int ipr = 0;ipr < int_per_row;ipr++)
		vec[thread][ipr] = uint_dist((random_engines[thread]));

	if((n&bitmask_64) != 0)
		vec[thread][int_per_row-1]<<=(63-((n-1)&bitmask_64));

	bool l =scalar_secret(vec[thread]);

	if (ber_dists[thread](random_engines[thread]))
		l ^= 0x1;


	if(arr_for_secret==int_per_row && (n&bitmask_64)==0)
		vec[thread][int_per_row]=0;
	if(l)
		vec[thread][arr_for_secret]^=(1ULL<<(63-(n&bitmask_64)));

	return vec[thread];
}

uint64_t * initialize_Oracle(int threads)
{
	vec = new uint64_t*[threads];
	if((n&bitmask_64)!=0)
		for(int i =0;i<threads;i++)
			vec[i]=new uint64_t[int_per_row];
	else
		for(int i =0;i<threads;i++)
			vec[i]=new uint64_t[int_per_row+1];

	for(int ipr=0;ipr<int_per_row;ipr++)
		secret[ipr]=uint_dist(mt);

	if((n&bitmask_64) != 0)
		secret[int_per_row-1]<<=(63-((n-1)&bitmask_64));

	for (int i = 0;i < threads;i++)
	{
		ber_dists.push_back(std::bernoulli_distribution(p));
		random_engines.push_back(std::mt19937_64((rd())));
		dists.push_back(std::uniform_int_distribution<uint64_t>());
	}
	return secret;
}

bool scalar_secret(uint64_t * a)
{
	uint64_t xor_sum=0;
	for(int ipr=0;ipr<int_per_row;ipr++)
		xor_sum^= (secret[ipr] & a[ipr]);
	
	for (int i = 0;i < 6;++i)
		xor_sum = xor_sum ^ (xor_sum >> (1 << (5 - i)));

	return xor_sum & 0x1;
}