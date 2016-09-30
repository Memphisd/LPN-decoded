#include "Matrixoperations.h"

time_t ti;
std::mt19937_64 mtw(static_cast<uint64_t>(time(&ti)));
std::uniform_int_distribution<uint64_t> dist;

extern std::vector<std::vector<int>> precomp_bitpositions;
extern uint64_t bitmask_r1;
extern uint64_t * right_error;
extern std::vector<std::mt19937_64> random_engines;
extern std::vector<std::uniform_int_distribution<uint64_t>> dists;
extern std::vector<std::vector<int>> zero_rows;
//calculates res = A*B
void matrix_mult(uint64_t** A, uint64_t ** B, uint64_t ** res, int col_A,int row_A,int col_B)
{
	for(int i =0;i<row_A;i++)
		for(int j=0;j<(col_B+63)/64;j++)
			res[i][j]=0;
	
	
	uint64_t bitmask_A;
	uint64_t bitmask_B;
	for(int i = 0; i<row_A;i++)
	{
		bitmask_B=1ULL<<63;
		for(int j=0;j<col_B;j++)
		{
			bitmask_A=1ULL<<63;
			for(int k=0;k<col_A;k++)
			{
				if((A[i][k / 64]&bitmask_A)>0 && (B[k][j / 64]&bitmask_B)>0)
					res[i][j / 64] ^=bitmask_B;
				
				bitmask_A=bitmask_A>1?bitmask_A>>1:1ULL<<63;
			}
			bitmask_B=bitmask_B>1?bitmask_B>>1:1ULL<<63;
		}
	}
}
//A * v = res
void matrix_vector_mult(uint64_t ** A, uint64_t * v, uint64_t * res, int row_A, int col_A)
{
	uint64_t bitmask = 1ULL<<63, zw;
	uint8_t bit;
	int k;
	for (int i = 0;i < (row_A + 63) / 64;i++)
		res[i]=0;

	for(int i=0;i<row_A;i++)
	{
		bit = 0;
		for (int j = 0;j < (col_A + 63) / 64;j++)
		{
			zw=A[i][j]&v[j];
			for(k=0;zw;k++)
				zw &= (zw-1);
			bit += (k%2);
		}
		bit %= 2;
		if(bit)
			res[i / 64]^=bitmask;
		
		bitmask=bitmask>1?bitmask>>1:1ULL<<63;
	}
}
//fast matrix vector multiplication for sparse vector using precomputated list to determine set bits
//operating on column-major matrix (tranposed)
void matrix_vector_mult_trans(uint64_t ** A, uint64_t * v, uint64_t * res, int row_A, int col_A)
{
	for (int i = 0;i < (col_A + 63) / 64;++i)
		res[i] = 0;

	unsigned char byte;

	for (int i = 0;i < (row_A+63)/64;++i)
	{
		for (int j = 0;j < 8;++j)
		{
			byte = (v[i] >> (56-8 * j)) & 0xFF;
			for (int k = 0;k < (int) precomp_bitpositions[byte].size();++k)
				for (int o = 0;o<(col_A + 63) / 64;++o)
				res[o] ^= A[precomp_bitpositions[byte][k]+8*j+i*64][o];
		}
			
	}
}
void transpose_matrix(uint64_t ** A, uint64_t ** res, int row_start, int row_end, int col_start, int col_end)
{
	//set result matrix to zero
	for(int i =0;i<col_end-col_start;i++)
	{
		for(int j=0;j<(row_end-row_start+63)/64;j++)
			res[i][j]=0;
	}

	//Transpose given submatrix and save in res
	uint64_t bitmask_row, bitmask_col =  1ULL<<63;
	for(int i =row_start;i<row_end;i++)
	{
		bitmask_row = 1ULL << (63 - (col_start % 64));

		for(int j=col_start;j<col_end;j++)
		{
			if((A[i][j/64]&bitmask_row)>0)
				res[j-col_start][(i-row_start)/64]^=bitmask_col;
			bitmask_row = bitmask_row > 1 ? bitmask_row >> 1 : 1ULL << 63;
		}
		bitmask_col=bitmask_col>1?bitmask_col>>1:1ULL<<63;
	}
}
//set matrix to identity matrix
void set_to_identity(uint64_t ** mat, int rows)
{
	uint64_t bitmask = 1ULL << 63;
	for (int i = 0;i < rows;i++)
	{
		for (int ipr = 0;ipr < (rows+63)/64;ipr++)
			mat[i][ipr] = 0;
		mat[i][i / 64] ^= bitmask;
		bitmask = bitmask > 1 ? bitmask >> 1 : 1ULL << 63;
	}
}

void gaussian_elimination(uint64_t **mat, int row, int row_start, int col, int col_start, int n, uint64_t ** T)
{
	uint64_t bitmask = 1ULL << (63 - (col_start % 64));
	int k;
	int act_arr;
	int int_per_col = (col +63) / 64;

	//convert matrix to echelon form
	for (int i = row_start;i<n + row_start;i++)
	{
		act_arr = (i-row_start+col_start) / 64;
		
		for (k = i; k<(n + row_start) && ((bitmask & mat[k][act_arr]) == 0) ;k++) {}
		if (k != n + row_start && k!=i)
			for (int ipr = 0;ipr < int_per_col;ipr++)
			{
				mat[i][ipr] ^= mat[k][ipr];
				if (ipr<n_k_length)
					T[i][ipr] ^= T[k][ipr];
			}

		for (int j_row = (k <= i ? i + 1 : k);j_row<(n + row_start);j_row++)
			if ((bitmask & mat[j_row][act_arr]) > 0)
			{
				for (int ipr = 0;ipr < int_per_col;ipr++)
				{
					mat[j_row][ipr] ^= mat[i][ipr];
					if (ipr<n_k_length)
						T[j_row][ipr] ^= T[i][ipr];
				}
			}
		bitmask = bitmask>1 ? bitmask >> 1 : 1ULL << 63;
	}

	//finish matrix diagonalization
	bitmask = (col_start + n) % 64 == 0 ? 1ULL : 1ULL << (63 - ((n + col_start - 1) % 64));

	for (int i = row_start + n - 1;i >= row_start;i--)
	{
		int act_arr = (col_start + i - row_start) / 64;
		if ((mat[i][act_arr] & bitmask)>0)
			for (int j =0;j<i;j++)
				if ((mat[j][act_arr] & bitmask)>0)
					for (int ipr = 0;ipr < int_per_col;ipr++)
					{
						mat[j][ipr] ^= mat[i][ipr];
						if (ipr<n_k_length)
							T[j][ipr] ^= T[i][ipr];
					}
						
		
		bitmask = (bitmask == (1ULL << 63)) ? 1ULL : bitmask << 1;
	}
}
void gaussian_elimination(uint64_t **mat,int row, int row_start, int col, int col_start, int n)
{
	uint64_t bitmask = 1ULL<<(63-(col_start%64));
	int k;
	int act_arr;
	int int_per_col=(col+63)/64;

	//convert matrix to echelon form
	for(int i = row_start;i<n+row_start;i++)
	{
		act_arr=(i+col_start-row_start)/64;
		for(k = i; k<n+row_start && ((bitmask & mat[k][act_arr]) == 0);k++){}
		if(k!=n+row_start)
			if(k!=i)
				for(int ipr=0;ipr<int_per_col;ipr++)
					mat[i][ipr]^=mat[k][ipr];
			for(int j_row=(k==i?i+1:k);j_row<n+row_start;j_row++)
				if((bitmask & mat[j_row][act_arr]) > 0)
					for(int ipr=0;ipr<int_per_col;ipr++)
						mat[j_row][ipr]^=mat[i][ipr];
		bitmask = bitmask>1?bitmask>>1:1ULL<<63;
	}


	//finish matrix diagonalization
	bitmask = (col_start+n)%64==0?1ULL:1ULL<<(63-((n+col_start-1)%64));

	for(int i = row_start+n-1;i>=row_start;i--)
	{
		int act_arr = (col_start+i-row_start)/64;
		if((mat[i][act_arr]&bitmask)>0)
			for(int j= row_start;j<i;j++)
				if((mat[j][act_arr]&bitmask)>0)
					for(int ipr=0;ipr<int_per_col;ipr++)
						mat[j][ipr]^=mat[i][ipr];
		bitmask = (bitmask == (1ULL<<63))?1ULL:bitmask<<1;
	}
}

void swap_cols(uint64_t ** mat, int s1, int s2, int rows)
{
	uint64_t bitmask_s1 = 1ULL << (63 - (s1 % 64));
	uint64_t bitmask_s2 = 1ULL << (63 - (s2 % 64));
	bool save=false;

	for (int i = 0;i < rows;i++)
	{
		if ((mat[i][s1 / 64] & bitmask_s1)>0)
		{
			save = true;
			mat[i][s1 / 64] ^= bitmask_s1;
		}
		if ((mat[i][s2 / 64] & bitmask_s2)>0)
		{
			mat[i][s1 / 64] ^= bitmask_s1;
			mat[i][s2 / 64] ^= bitmask_s2;
		}

		if (save)
			mat[i][s2 / 64] ^= bitmask_s2;
		save = false;
	}
}

void copy_matrix(uint64_t ** Src, uint64_t ** Dest, int rows, int cols)
{
	for (int i = 0;i < rows;i++)
		for (int j = 0;j < (cols+63)/64;j++)
			Dest[i][j] =  Src[i][j];

}
void permute(uint64_t ** mat, uint64_t ** P, int rows, int cols, int thread)
{
	//permute cols of the matrix
	int swap1, swap2;
	for (int i = 0;i < cols;i++)
	{
		swap1 = i;
		swap2 = i+(dists[thread](random_engines[thread]) %(cols-i));
		if (swap1 != swap2)
		{
			swap_cols(mat, swap1, swap2, rows);
			swap_cols(P , swap1, swap2, cols);
		}
	}
}

void create_generator_matrix(uint64_t ** mat, uint64_t ** gen_mat, uint64_t * codew_errors, int dim, int thread)
{
	uint64_t * g;
	int act_arr;
	uint64_t bit;
	for (int ipr = 0;ipr < (dim+63)/64;ipr++)
		codew_errors[ipr] = 0;

	for (int i = 0;i < (n + 63) / 64;i++)
		right_error[i] = 0;

	act_arr = k % 64 == 0 ? int_per_row : int_per_row - 1;
	for (int i = 0; i < dim;i++)
	{
		g = query(thread,true);
		bit = ((g[act_arr] >> (63 - (k % 64))) & 0x1);
		if (bit)
		{
			g[act_arr] ^= (1ULL << (63 - (k % 64)));
			codew_errors[i / 64] ^= bit << (63 - (i % 64));
		}

		for (int ipr = 0;ipr<int_per_row;ipr++)
			mat[i][ipr] = g[ipr];
	}

	transpose_matrix(mat, gen_mat, 0, n, 0, k);

	//partial gaussian elimination on generator matrix
	gaussian_elimination(gen_mat, k, 0, n, 0, k);

}

void transform_to_parity_check(uint64_t ** gen_mat, uint64_t ** parity_mat)
{
	//transpose non identity part of generator matrix
	uint64_t ** transposed_gen = new uint64_t *[n-k];
	for (int i = 0;i <n-k;i++)
		transposed_gen[i] = new uint64_t[int_per_row];

	transpose_matrix(gen_mat, transposed_gen, 0, k, k, n);

	//save transposed matrix in parity check matrix and set remaining columns to zero
	for (int i = 0;i < n-k;i++)
	{
		for (int ipr = 0; ipr < int_per_row;ipr++)
			parity_mat[i][ipr] = transposed_gen[i][ipr];
		for (int ipr = int_per_row;ipr < (n+63)/64;ipr++)
			parity_mat[i][ipr] = 0;
	}
	//add identity part of parity check matrix
	uint64_t bitmask = 1ULL << (63 - (k % 64));
	for (int i = 0;i < (n-k);i++)
	{
		parity_mat[i][(i + k) / 64] ^= bitmask;
		bitmask = bitmask > 1 ? bitmask >> 1 : 1ULL << 63;
	}
}
//output uint64_t in binary format
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
	int int_per_r= (col+63) / 64;
	for (int i = 0; i<row;i++)
	{
		for (int ipr = 0;ipr < int_per_r;ipr++)
			show_bin(mat[i][ipr]);//,ipr!=(int_per_r-1) || col%64==0?64:col%64);
		printf("\n");
	}
}
void get_submatrix(uint64_t ** mat, uint64_t ** res, int row_end, int col_end)
{
	//copy necessary rows/columns
	int act_arr = (col_end+63)/64;
	for(int i = 0;i<row_end;i++)
		for(int j = 0;j<act_arr;j++)
			res[i][j]=mat[i][j];

	//erase not used columns
	if(col_end%64!=0)
		for(int i = 0;i<row_end;i++)
			res[i][act_arr-1]= (res[i][act_arr-1]>>(63-((col_end-1)%64)))<<(63-((col_end-1)%64));
}

bool check_gaussian(uint64_t ** mat,int row_start, int col_start, int n)
{
	int act_arr;
	uint64_t bitmask = 1ULL << (63 - ((col_start+n-1)% 64));
	for (int i = row_start+n-1;i >= row_start;i--)
	{
		act_arr = (i - row_start + col_start) / 64;
		if ((mat[i][act_arr] & bitmask) == 0)
			return false;
		bitmask = bitmask==1ULL<<63?1ULL:bitmask<<1;
	}
	return true;
}

int get_zero_rows(uint64_t ** mat, int row_start, int col_start, int n, int thread)
{
	int act_arr;
	uint64_t bitmask = 1ULL << (63 - ((col_start + n - 1) % 64));
	int count = 0;
	for (int i = row_start + n - 1;i >= row_start;i--)
	{
		act_arr = (i - row_start + col_start) / 64;
		if ((mat[i][act_arr] & bitmask) == 0)
		{
			zero_rows[thread][count] = i;
			count++;
		}
		bitmask = bitmask == 1ULL << 63 ? 1ULL : bitmask << 1;
	}
	return count;
}

bool swap_zero_cols(uint64_t ** mat, uint64_t ** perm_mat, uint64_t ** T, int amount, int row_start, int col_start, int thread, int rows, int cols, int n)
{
	uint64_t bitmask_zero, bitmask_row;
	int act_arr;
	int int_per_col = (cols + 63) / 64, j;

	for (int i = 0;i < amount;++i)
	{
		//check if zero_row still exist
		bitmask_zero = 1ULL << (63 - ((col_start + zero_rows[thread][i] - row_start)) % 64);
		act_arr = (zero_rows[thread][i] - row_start + col_start) / 64;
		if ((bitmask_zero & mat[zero_rows[thread][i]][act_arr]) == 0)
		{
			//search col with "1" in that row
			bitmask_row = 1ULL << 63;
			for (j = 0;j < col_start;++j)
			{
				if ((bitmask_row & mat[zero_rows[thread][i]][j / 64]) > 0)
				{
					//swap zero column with found column to eliminate zero row
					swap_cols(mat, j, zero_rows[thread][i]-row_start+col_start, rows);
					swap_cols(perm_mat, j, zero_rows[thread][i] - row_start + col_start, cols);
					//proceed with partial gaussian elimination
					for (int k = 0;k < rows;++k)
					{
						if(k!=zero_rows[thread][i] && (mat[k][act_arr]&bitmask_zero)>0)
							for (int ipr = 0;ipr < int_per_col;++ipr)
							{
								mat[k][ipr] ^= mat[zero_rows[thread][i]][ipr];
								if (ipr<n_k_length)
									T[k][ipr] ^= T[zero_rows[thread][i]][ipr];
							}
					}
					break;
				}
				//if no "1" is found in that row we need a whole new permutation
				if (j == col_start)
					return false;
			}
		}
	}
	return true;
}

void create_random_mat(uint64_t ** mat, int rows,int cols)
{
	for(int i=0;i<rows;i++)
		for(int j=0;j<(cols+63)/64;j++)
		{
			mat[i][j]=dist(mtw);
			if(j==cols/64 && cols %64 !=0)
				mat[i][j]= (mat[i][j]>>(64-(cols%64)))<<(64-(cols%64));
		}
}

void get_col(uint64_t ** mat, int col, int rows, uint64_t * vec)
{
	uint64_t bitmask_col = 1ULL << (63 - (col % 64));
	uint64_t bitmask_vec = 1ULL << 63;
	vec[(rows+63) / 64-1] = 0;

	for (int i = 0;i < rows;i++)
	{
		if ((mat[i][col / 64] & bitmask_col)>0)
			vec[i / 64] |= bitmask_vec;
		else if ((vec[i / 64] & bitmask_vec)>0)
			vec[i / 64] ^= bitmask_vec;
	
		bitmask_vec = bitmask_vec == 1 ? 1ULL << 63 : bitmask_vec >> 1;
	}
}

void set_matrix(uint64_t** mat, bool* secret, int rows, int thread)
{
	uint64_t * g;
	int act_arr = k % 64 == 0 ? int_per_row : int_per_row - 1;;
	for (int i = 0; i < rows;i++)
	{
		g = query(thread,false);
		secret[i] = ((g[act_arr] >> (63 - (k % 64))) & 0x1);
		if (secret[i])
			g[act_arr] ^= (1ULL << (63 - (k % 64)));
		for (int ipr = 0;ipr<int_per_row;ipr++)
			mat[i][ipr] = g[ipr];
	}
}

bool test_guess(uint64_t** mat_test, bool* secret_test, uint64_t * guess)
{
	int count = 0;
	for (int i = 0;i<m;i++)
	{
		if (scalar(guess, mat_test[i]) != secret_test[i])
			count++;
		if (count >= c)
			break;
	}
	if (count<c)
		return true;

	
	return false;
}