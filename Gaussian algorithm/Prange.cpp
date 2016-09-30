// bruteforce_linear_codes.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "Prange.h"

uint64_t *  querys = new uint64_t[threads];

bool debug = false;
bool * stop;
uint64_t * right_secret;
int vec_pool_size =  ceil(16*pow(n,2)/(1-p));

uint64_t ** vec_pool;
bool * vec_pool_sec;
extern std::vector<std::mt19937_64> random_engines;
extern std::vector<std::uniform_int_distribution<uint64_t>> dists;
double exec_time = 0;
double iter_per_sec = 0;
bool found_right;
uint64_t** premature_sec = new uint64_t*[threads];

uint64_t ** mat_test = new uint64_t*[test_amount];
bool * secret_test = new bool[test_amount];

inline void xor_entry(uint64_t** mat, int i, int j)
{
	for(int ipr=0;ipr<vec_sec_length;ipr++)
		mat[i][ipr]^=mat[j][ipr];
}

void gaussian_elimination(uint64_t **mat, int n)
{
	uint64_t bitmask = 1ULL<<63;
	int k;
	//convert matrix to echelon form
	for(int i = 0;i<n;i++)
	{
		for(k = i; k<n && ((bitmask & mat[k][i/64]) == 0);k++){}
		if(k!=n)
			if(k!=i)
				xor_entry(mat, i,k);
			for(int j_row=(k==i?i+1:k);j_row<n;j_row++)
				if((bitmask & mat[j_row][i/64]) > 0)
					xor_entry(mat, j_row,i);
		bitmask = bitmask>1?bitmask>>1:1ULL<<63;
	}

	//finish matrix diagnoalization
	bitmask = (n&bitmask_64)==0?1ULL:1ULL<<(64-(n&bitmask_64));
	for(int i = n-1;i>=0;i--)
	{
		if((mat[i][i/64]&bitmask)>0)
			for(int j=0;j<i;j++)
				if((mat[j][i/64]&bitmask)>0)
					xor_entry(mat,j,i);
		bitmask = (bitmask == (1ULL<<63))?1ULL:bitmask<<1;
	}

}

void gaussian_elimination_on_prep_mat(uint64_t **mat, int row, int row_start, int col, int col_start, int n)
{
	uint64_t bitmask = 1ULL << 63;
	int k;
	int act_arr;


	//eliminate ones below the pre_generated gaussian part
	for (int i = 0;i < pre_gen_mat_size;++i)
	{
		for (int j = pre_gen_mat_size;j < row;++j)
			if ((mat[j][i / 64] & bitmask)>0)
				xor_entry(mat, j, i);
		bitmask = bitmask>1 ? bitmask >> 1 : 1ULL << 63;
	}

	bitmask = 1ULL << (63 - (col_start &bitmask_64));
	//convert matrix to echelon form
	for (int i = row_start;i<n + row_start;i++)
	{
		act_arr = (i + col_start - row_start) / 64;
		for (k = i; k < n + row_start && ((bitmask & mat[k][act_arr]) == 0);k++) {}
		if (k != n + row_start)
			if (k != i)
				xor_entry(mat, i, k);

		for (int j_row = (k == i ? i + 1 : k);j_row < n + row_start;j_row++)
			if ((bitmask & mat[j_row][act_arr]) > 0)
				xor_entry(mat, j_row, i);
			
		
		bitmask = bitmask>1 ? bitmask >> 1 : 1ULL << 63;
	}

	//finish matrix diagonalization
	bitmask = ((col_start + n) &bitmask_64) == 0 ? 1ULL : 1ULL << (63 - ((n + col_start - 1) &bitmask_64));

	for (int i = row_start + n - 1;i >= row_start;i--)
	{
		int act_arr = (col_start + i - row_start) / 64;
		if ((mat[i][act_arr] & bitmask) > 0)
			for (int j = row_start;j < i;j++)
				if ((mat[j][act_arr] & bitmask)>0)
					xor_entry(mat, j, i);
		bitmask = (bitmask == (1ULL << 63)) ? 1ULL : bitmask << 1;
	}

	//eliminate ones above the just generated gaussian part
	bitmask = 1ULL << (63-(pre_gen_mat_size &bitmask_64));
	for (int i = pre_gen_mat_size;i < row;++i)
	{
		if(mat[i][i/64]>0)
			for (int j = 0;j < pre_gen_mat_size;++j)
				if ((mat[j][i / 64] & bitmask)>0)
					xor_entry(mat, j, i);
		bitmask = bitmask>1 ? bitmask >> 1 : 1ULL << 63;
	}
}
void set_secret(int setting,std::list<int>* zero_ind, uint64_t ** mat,uint64_t * guess)
{

	int i = 0;
	uint64_t bitmask;
	//set secret vector at zero-row positions to actual guess 
	for(iter_int it = zero_ind->begin(); it!=zero_ind->end() ;++it)
	{
		bitmask = 1ULL << (63 - ((*it) &bitmask_64));
		guess[(*it) / 64] ^= (guess[(*it) / 64] & bitmask) ^ (((setting >> i) & 0x1) << (63 - ((*it)&bitmask_64)));
		i++;
	}

	//calculate secretvector dependent on actual guess
	for(iter_int it = zero_ind->begin(); it!=zero_ind->end();++it)
	{
		bitmask = 1ULL << (63 - ((*it) &bitmask_64));
		if((guess[(*it) / 64] & bitmask)>0)
		{
			for (int j_row = 0;j_row < n;j_row++)
				if ((bitmask & mat[j_row][(*it) / 64]) > 0 && j_row != (*it))
					guess[(*it) / 64] ^= (1ULL << (63 - (j_row &bitmask_64)));
		}
	}
}

//determine secret from matrix and save zero-row positons, then guess these positions of the secret vector
void determine_secret(uint64_t ** mat,std::list<int>* zero_ind,std::list<uint64_t *>* secrets,int thread)
{
	//find zero-rows
	uint64_t bitmask;
	for(int i =n-1;i>=0;i--)
	{
		bitmask = 1ULL<<(63-(i&bitmask_64));
		if((mat[i][i/64]&bitmask)==0)
			zero_ind->push_back(i);
	}
	
	//determine premature secret from matrix 
	uint64_t bitmask_row = 1ULL << 63;
	bitmask = 1ULL << (63 - (n &bitmask_64));
	for (int i = 0;i < int_per_row;++i)
		premature_sec[thread][i] = 0;
	for (int i = 0;i < n;++i)
	{
		if ((mat[i][vec_sec_length - 1] & bitmask)>0)
			premature_sec[thread][i / 64] ^= bitmask_row;
		bitmask_row = bitmask_row == 1 ? 1ULL << 63 : bitmask_row>>1;
	}

	if(zero_ind->size()<=max_guesses)
	{
		iter_uint it = secrets->begin();
		for(int i = 0;i<(1<<zero_ind->size()) ;i++)
		{
			for (int i = 0;i < int_per_row;++i)
				(*it)[i] = premature_sec[thread][i];
			if(zero_ind->size()!=0)
				set_secret(i,zero_ind, mat, (*it));
			++it;
		}
	}
}

void set_test_matrix(uint64_t** mat, bool* secret, int rows, int thread)
{
	uint64_t * g;
	int act_arr = (n &bitmask_64) == 0 ? int_per_row : int_per_row - 1;;
	for(int i = 0; i < rows;i++)
	{
		g=query(thread);
		secret[i]=((g[act_arr]>>(63-(n&bitmask_64)))&0x1);
		if(secret[i])
			g[act_arr]^=(1ULL<<(63-(n&bitmask_64)));
		for(int ipr=0;ipr<int_per_row;ipr++)
			mat[i][ipr]=g[ipr];
	}
}

void set_matrix(uint64_t** mat, int rows, int thread, int start)
{
	uint64_t * g;
	for (int i = start; i < rows;i++)
	{
		g = query(thread);
		for (int ipr = 0;ipr<vec_sec_length;ipr++)
			mat[i][ipr] = g[ipr];
	}
}

bool scalar(uint64_t * a,uint64_t * b)
{
	uint64_t xor_sum=0;
	for(int ipr=0;ipr<int_per_row;ipr++)
		xor_sum^= (b[ipr] & a[ipr]);

	for (int i = 0;i < 6;++i)
		xor_sum = xor_sum ^ (xor_sum >> (1<<(5-i)));

	return xor_sum&0x1;
}

uint64_t* test_guess(std::list<uint64_t *>* sec,uint64_t** mat,bool* secret, int guess_count)
{
	int count;
	int gc=0;

	for(iter_uint it = sec->begin();it!=sec->end() && gc<guess_count;++it)
	{
		count=0;
		for(int i = 0;i<test_amount;i++)
		{
			if(scalar((*it),mat[i])!=secret[i])
				count++;
			if(count>=threshold)
				break;
		}
		if(count<threshold)
			return (*it);	
		gc++;
	}
	return NULL;
}

void set_mat_from_pool(uint64_t ** mat, int rows, int thread, bool half)
{
	for (int i = half ? 0 : pre_gen_mat_size;i < rows;++i)
	{
		int choose = dists[thread](random_engines[thread]) % vec_pool_size;
		for (int j = 0;j < vec_sec_length;++j)
			mat[i][j] = vec_pool[choose][j];
	}
}

bool check_gaussian(uint64_t ** mat, int row_start, int col_start, int n)
{
	int act_arr;
	uint64_t bitmask = 1ULL << (63 - ((col_start + n - 1) &bitmask_64));
	for (int i = row_start + n - 1;i >= row_start;i--)
	{
		act_arr = (i - row_start + col_start) / 64;
		if ((mat[i][act_arr] & bitmask) == 0)
			return false;
		bitmask = bitmask == 1ULL << 63 ? 1ULL : bitmask << 1;
	}
	return true;
}

void show_bin(uint64_t num)
{
	uint64_t bit = 1ULL << 63;
	for (int i = 0; i <64;i++)
	{
		printf("%d", (int)((bit&num) >> (63 - i)));
		bit = bit >> 1;
	}
}

//output matrix
void print_matrix(uint64_t **mat, int col, int row)
{
	int int_per_r = (col + 63) / 64;
	for (int i = 0; i<row;i++)
	{
		for (int ipr = 0;ipr < int_per_r;ipr++)
			show_bin(mat[i][ipr]);
		printf("\n");
	}
}

void start_bruteforce()
{
	//setup all necessary pointers for each thread
	uint64_t ** mat = new uint64_t*[n];
	for(int i =0;i<n;i++)
		mat[i]=new uint64_t[vec_sec_length];

	uint64_t ** pre_gen_mat = new uint64_t*[pre_gen_mat_size];

	for (int i = 0; i < pre_gen_mat_size;i++)
		pre_gen_mat[i] = new uint64_t[vec_sec_length];

	uint64_t * det_secret = NULL;

	std::list<uint64_t *>* sec= new std::list<uint64_t*>();
	for(int i= 0;i<(1<<max_guesses);i++)
	{
		uint64_t * g = new uint64_t[int_per_row];
		sec->push_back(g);
	}
	std::list<int>* zero_ind = new std::list<int>();
	int thread = omp_get_thread_num();
	
	while(det_secret == NULL)
	{	
		if (pre_gen_mat_size != 0)
		{
			do
			{
				if (prange)
					set_mat_from_pool(pre_gen_mat, pre_gen_mat_size, thread, true);
				else
					set_matrix(pre_gen_mat, pre_gen_mat_size, thread, 0);
				gaussian_elimination(pre_gen_mat, pre_gen_mat_size);
			} while (!(check_gaussian(pre_gen_mat, 0, 0, pre_gen_mat_size)));
		}
		for (uint64_t iter = 0; iter < pre_gen_iter && det_secret== NULL;++iter)
		{
			zero_ind->clear();
			//copy pregenerated matrix to working matrix
			for (int i = 0;i < pre_gen_mat_size;++i)
			{
				for (int ipr = 0;ipr < vec_sec_length;++ipr)
					mat[i][ipr] = pre_gen_mat[i][ipr];
			}
			if(prange)
				set_mat_from_pool(mat, n, thread, false);
			else
				set_matrix(mat, n, thread, pre_gen_mat_size);

			//do gaussian elimination
			gaussian_elimination_on_prep_mat(mat, n, pre_gen_mat_size, n, pre_gen_mat_size, n - pre_gen_mat_size);

			//determine secret from matrix (including all necessary guesses)
			determine_secret(mat, zero_ind, sec,thread);
		
			//test all the secrets
			det_secret = test_guess(sec, mat_test, secret_test, (1 << zero_ind->size()));

			querys[thread]++;

			if ((*stop))
				break;
		}
		if ((*stop))
			break;
	}

	//Output determined secret and runtime informations
	if(!(*stop))
	{
		*stop=true;
		found_right = true;
		for (int ipr = 0;ipr < int_per_row;ipr++)
			if (right_secret[ipr] != det_secret[ipr])
			{
				found_right = false;
				break;
			}

		if (found_right)
			std::cout << "Thread " << thread << " found the right secret\n";
		else
			std::cout << "Thread " << thread << " found a wrong secret\n";

	}


}
int main(int argc, char *argv[])
{

	std::ofstream out(out_file, std::ofstream::app);
	std::streambuf *coutbuf = std::cout.rdbuf();
	std::cout.rdbuf(out.rdbuf());
	stop = new bool;

	//initialize vectors and matrices
	vec_pool = new uint64_t *[vec_pool_size];
	vec_pool_sec = new bool[vec_pool_size];
	for (int i = 0;i < vec_pool_size;i++)
		vec_pool[i] = new uint64_t[vec_sec_length];
	for (int i = 0;i < threads;++i)
		premature_sec[i] = new uint64_t[int_per_row];
	for (int i = 0;i < test_amount;i++)
		mat_test[i] = new uint64_t[int_per_row];


	std::cout << "\n_____________Start with n=" << n << ", p=" << p << ", x="<<pre_gen_mat_size<<", delta="<<pre_gen_iter<<", threads=" << threads << "_____________\n";

	right_secret = initialize_Oracle(threads);

	
	omp_set_num_threads(threads);


	
	for (int i = 0;i < runs;i++)
	{
		set_test_matrix(mat_test, secret_test, test_amount, 0);

		for (int i = 0;i < threads;i++)
			querys[i] = 0;

		printf("Executing run %d of %d\n", i + 1, runs);

		set_matrix(vec_pool, vec_pool_size, 0,0);

		printf("started");
		*stop = false;
		const clock_t START = clock();
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
		for (int q = 0;q < omp_get_num_threads();q++)
			start_bruteforce();
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
		const double T_ELAPSED = (double)(clock() - START) / CLOCKS_PER_SEC;
		uint64_t iterations = 0;
		for (int i = 0;i < threads;i++)
			iterations += querys[i];
		exec_time += ((double)duration) / 1000;
		iter_per_sec += (double)iterations / T_ELAPSED;


			std::cout << "Duration (seconds): " << ((double)duration) / 1000 << "\n";
			std::cout << "CPU-Time (seconds): " << T_ELAPSED << "\n";
			std::cout << "Total iterations: " << (double)iterations << "\n";
			std::cout << "Iterations per second (real-time): " << (double)iterations / (((double)duration) / 1000) << "\n";
			std::cout << "Iterations per second (cpu -time): " << (double)iterations / T_ELAPSED << "\n";
		
	}
	if (runs>1)
	{
		std::cout << "Average iterations per second (cpu-time): " << iter_per_sec / runs << "\n\n\n";
		std::cout << "Average Execution Time: " << exec_time / runs << "\n\n\n";
	}

	printf("finished");
	std::cout.rdbuf(coutbuf);
	return 0;

}

