#include "BJMM-.h"

extern int errors;
extern std::mt19937_64 mt;
extern std::uniform_int_distribution<uint64_t> uint_dist;
std::vector<std::vector<int>> precomp_bitpositions(256, std::vector<int>());
std::vector<std::vector<int>> zero_rows(threads, std::vector<int>(20));

//basis vector information
static uint64_t permutations = choose((l + k) / 2, (p / 2 + eps) / 2);
std::vector<std::pair<int, int>> changes;

//flags
bool * stop = new bool;
bool * temp = new bool[threads];

//correct errror and secret (just for testing)
uint64_t * right_error = new uint64_t[(n + 63) / 64];
uint64_t * right_secret;


uint64_t * codeword_w_errors = new uint64_t[length_v_n];
uint64_t ** T_det_sec = new uint64_t *[k];
uint64_t ** codewords = new uint64_t *[threads];
uint64_t ** secrets = new uint64_t *[threads];

//hypothesis testing matrices and vectors 
uint64_t ** mat_test = new uint64_t*[m];
bool * secret_test = new bool[m];

void my_handler(int s)
{
	(*stop) = true;
}
//debug informations
uint64_t * iterations = new uint64_t[threads];
uint64_t * right_perms = new uint64_t[threads];

struct compare_l0_labels
{
	bool operator()(const list_entry& s1, const list_entry& s2) const
	{
		for (int i = 0;i<length_v_r1;i++)
		{
			if (s1.value[i]<s2.value[i])
				return true;
			else if (s1.value[i]>s2.value[i])
				return false;
		}
		return false;
	}
};

struct compare_l1_labels
{
	bool operator()(const list_entry& s1, const list_entry& s2) const
	{
		for (int i = 0;i<length_v_r2;i++)
		{
			if (s1.value[i]<s2.value[i])
				return true;
			else if (s1.value[i]>s2.value[i])
				return false;
		}
		return false;
	}
};

int choose(int n, int k)
{
	if (k == 0)
		return 1;
	return (n * choose(n - 1, k - 1)) / k;
}

void xor_regs(uint64_t * A, uint64_t * B, int regs)
{
	for (int i = 0;i<regs;i++)
		A[i] ^= B[i];
}

void xor_regs_save(uint64_t * A, uint64_t * B, uint64_t * C, int regs)
{
	for (int i = 0;i<regs;i++)
		A[i] = B[i]^C[i];
}

bool scalar(uint64_t * a, uint64_t * b)
{
	uint64_t xor_sum = 0;
	for (int ipr = 0;ipr<int_per_row;ipr++)
		xor_sum ^= b[ipr] & a[ipr];

	for (int i = 0;i < 6;++i)
		xor_sum = xor_sum ^ (xor_sum >> (1 << (5 - i)));

	return xor_sum & 0x1;
}


//right shift a uint64_t array bitwise
void r_shift(uint64_t * v, int shift_amount, int regs)
{
	//handle shifts greater then 64 via swaps of whole uints
	if (shift_amount > 63)
	{
		for (int i = regs - 1;i - shift_amount / 64 >= 0;i--)
		{
			v[i] = v[i - shift_amount / 64];
		}
		for (int i = 0;i < shift_amount / 64;i++)
			v[i] = 0;
		shift_amount = shift_amount % 64;
	}

	//do shifts lower then 64 via carrys
	if (shift_amount != 0)
	{
		uint64_t zw[2] = { 0,0 };
		uint64_t bitmask = (1ULL << shift_amount) - 1;
		for (int i = 0;i<regs;i++)
		{
			zw[i % 2] = v[i] & bitmask;
			v[i] >>= shift_amount;
			v[i] ^= (zw[(i + 1) % 2] << (64 - shift_amount));
		}
	}
}

//left shift a uint64_t array bitwise (shift_amount < 64)
void l_shift(uint64_t * v, int shift_amount, int regs)
{
	int i;
	for (i = 0; i < regs - 1; ++i)
	{
		v[i] = (v[i] << shift_amount) | ((v[i + 1] >> (64 - shift_amount)) & ((1 << shift_amount) - 1));
	}
	v[regs - 1] <<= shift_amount;
}

int weigth(uint64_t v)
{
	int w;
	for (w = 0;v;w++)
		v &= (v - 1);
	return w;
}

bool check_solution(int i1, int i2, int j1, int j2, list L_h, list L_l,uint64_t ** Q_trans, uint64_t * syn_tilde, uint64_t ** perm_mat,int thread)
{
	int we = 0;
	uint64_t Qe1[(n - k + 63) / 64];
	uint64_t e1[length_v_kl];
	uint64_t e2[length_v_n];
	uint64_t e_tilde[length_v_n];
	uint64_t error[length_v_n];

	//rebuild error based on given indices
	xor_regs_save(e1, (*L_l)[i1], (*L_l)[i2], length_v_kl);
	xor_regs(e1, (*L_h)[j1], length_v_kl);
	xor_regs(e1, (*L_h)[j2], length_v_kl);

	matrix_vector_mult_trans(Q_trans, e1, Qe1, k + l, n - k);

	//calculate we=weight(e2)
	for (int i = 0;i < (n - k + 63) / 64;i++)
		we += weigth(Qe1[i] ^ syn_tilde[i]);

	int we_e1 = 0;
	for (int i = 0;i < length_v_kl;i++)
		we_e1 += weigth(e1[i]);

	//if e2 has small weight
	if (we <= max_error_weight - we_e1)
	{
		//calculate corresponding secret
		for (int i = 0;i < length_v_n;i++)
			e2[i] = e_tilde[i] = 0;
		for (int i = 0;i < length_v_kl;++i)
			e_tilde[i] = e1[i];

		for (int i = 0;i < (n - k + 63) / 64;i++)
			e2[i] = Qe1[i] ^ syn_tilde[i];

		r_shift(e2, k, length_v_n);

		for (int i = 0;i < length_v_n;i++)
			e_tilde[i] ^= e2[i];
	
		matrix_vector_mult(perm_mat, e_tilde, error, n, n);

		for (int i = 0;i < (k + 63) / 64;++i)
			codewords[thread][i] = codeword_w_errors[i] ^ error[i];
		if (k % 64 != 0)
			codewords[thread][k / 64] = (codewords[thread][k / 64] >> (63 - (k - 1) % 64)) << (63 - (k - 1) % 64);
		matrix_vector_mult(T_det_sec, codewords[thread], secrets[thread], k, k);

		
		//do hypothesis test
		if (test_guess(mat_test, secret_test, secrets[thread]))
		{
			
			(*stop) = true;

			//check against correct secret (just for testing)
			bool right_sec = true;
			for (int i = 0;i < (k + 63) / 64;++i)
				if (right_secret[i] != secrets[thread][i])
					right_sec = false;
			if (right_sec)
				std::cout << "\ngot the RIGHT secret";
			else
				std::cout << "\ngot the WRONG secret";

			//output information
			uint64_t total_iter = 0;
			for (int i = 0;i < threads;i++)
				total_iter += iterations[i];
			std::cout << "\nneeded " << total_iter << " iterations";

			return true;
		}
		else //choose new hypothesis-test-matrix
			set_matrix(mat_test, secret_test, m, 0);
	}
	return false;
}

//Merge two lists via sorted lists
int merge_join(index_list L1, index_list L2, index_list merged_List, int indices, int L1_size, int L2_size, int thread, list L_h, list L_l, uint64_t ** Q_trans, uint64_t * syn_tilde, uint64_t ** perm_mat)
{
	list_entry new_entry;
	int i0, i1, j0, j1;
	uint64_t c = 0;
	int i = 0, j = 0;

	while (i < L1_size && j < L2_size)
	{
		if ((*L1)[i].value[0]<(*L2)[j].value[0])
			i++;
		else if ((*L1)[i].value[0]>(*L2)[j].value[0])
			j++;
		else if ((*L1)[i].value[0] == (*L2)[j].value[0])
		{
			i0 = i1 = i;
			j0 = j1 = j;
			while (i1 < L1_size && (*L1)[i0].value[0] == (*L1)[i1].value[0]) { i1++; }

			while (j1 < L2_size && (*L2)[j0].value[0] == (*L2)[j1].value[0]) { j1++; }

			for (i = i0;i < i1;i++)
				for (j = j0;j < j1;j++)
				{
					if (indices == 2)
					{
						if (c == merged_List->size())
						{
							new_entry.index = new int[indices];
							new_entry.value = new uint64_t[length_v_r2];
							merged_List->push_back(new_entry);
						}

						(*merged_List)[c].index[0] = (*L1)[i].index[0];
						(*merged_List)[c].index[1] = (*L2)[j].index[0];

					}
					else
						if (check_solution((*L1)[i].index[0], (*L2)[j].index[0], (*L1)[i].index[1], (*L2)[j].index[1], L_h, L_l, Q_trans, syn_tilde, perm_mat, thread))
						{
							temp[thread] = true;
							break;
						}
					c++;
				}
			if ((*stop) || temp[thread])
				break;
			i = i1;
			j = j1;
		}
	}

	return c;
}

void initialize_index_list(index_list L, int value_length, int amount_indices)
{
	for (index_list_iter it = L->begin();it != L->end();++it)
	{
		(*it).index = new int[amount_indices];
		(*it).value = new uint64_t[value_length];
	}
}

//get lower r1 bits of Q1*e1 out of List 
void get_lower_bits_as_list(index_list L, index_list L_low, uint64_t * value)
{
	int shift = (l - r1) % 64;
	for (uint64_t i = 0; i != L->size();i++)
	{
		for (int j = (l - r1) / 64;j<length_v_r1;j++)
			value[j - ((l - r1) / 64)] = (*L)[i].value[j];

		l_shift(value, shift, length_v_l);

		for (int j = 0;j<length_v_r1;j++)
			(*L_low)[i].value[j] = value[j];

		(*L_low)[i].index[0] = (*L)[i].index[0];
	}
}

//get lower r2 bits of Q1*e1 out of List 
void get_upper_bits_as_list(index_list L_l, index_list L_h, index_list level1_list, int size)
{
	for (int i = 0;i < size;i++)
	{
		for (int j = 0;j < length_v_r2;j++)
		{
			(*level1_list)[i].value[j] = (*L_l)[(*level1_list)[i].index[0]].value[j] ^ (*L_h)[(*level1_list)[i].index[1]].value[j];
		}
		if ((r2 % 64) != 0)
			(*level1_list)[i].value[length_v_r2 - 1] = (((*level1_list)[i].value[length_v_r2 - 1] >> (64 - (r2 % 64))) << (64 - (r2 % 64)));
	}
}

//add v to each entry in L
void add_to_list(index_list L, uint64_t * v, int regs, int size)
{
	for (int i = 0;i < size;++i)
		xor_regs ((*L)[i].value, v, regs);
}

//check if permutation achives right weight distribution (just for testing)
bool check_perm(uint64_t ** perm, uint64_t ** trans_perm, uint64_t * test_e_tilde, uint64_t * right_error, bool out)
{
	transpose_matrix(perm, trans_perm, 0, n, 0, n);

	matrix_vector_mult(trans_perm, right_error, test_e_tilde,n,n);
	int we = 0;

	for (int i = 0;i < (k + l + 63) / 64;++i)
	{
		if (i == ((k + l + 63) / 64 - 1))
			test_e_tilde[i] &= (((1ULL << ((k + l) % 64)) - 1)<<(64-((k + l) % 64)));
		we += weigth(test_e_tilde[i]);
	}

	int we2 = 0;
	for (int i = 0;i < ((k + l) / 2 + 63) / 64;++i)
	{
		if(i == ((k + l) / 2 + 63) / 64 - 1)
			test_e_tilde[i] &= (((1ULL << ((k + l)/2 % 64)) - 1) << (64 - ((k + l)/2 % 64)));
		we2 += weigth(test_e_tilde[i]);
	}

	if (we == p && we2==p/2)
	{
		return true;
	}
	else
		return false;

}

//iterative part of the BJMM-
void start_search(list L_h, list L_l, uint64_t ** parity_mat , uint64_t * syndrom, bool * stop)
{
	
	//------Listinitialization-------//
	#pragma region Lists


	//Initialize lists for Q*e for both permutation lists
	index_list Q_List_l = new std::vector<list_entry>(permutations);
	index_list Q_List_h = new std::vector<list_entry>(permutations);
	initialize_index_list(Q_List_l, length_v_l, 1);
	initialize_index_list(Q_List_h, length_v_l, 1);

	//create label lists of level 0 
	index_list l0_label_List1 = new std::vector<list_entry>(permutations);
	index_list l0_label_List2 = new std::vector<list_entry>(permutations);
	initialize_index_list(l0_label_List1, length_v_r1, 1);
	initialize_index_list(l0_label_List2, length_v_r1, 1);
	//create label lists of level 1
	index_list l1_label_List1 = new std::vector<list_entry>(permutations / (1 << r1)*permutations);
	index_list l1_label_List2 = new std::vector<list_entry>(permutations / (1 << r1)*permutations );
	initialize_index_list(l1_label_List1, length_v_r2, 2);
	initialize_index_list(l1_label_List2, length_v_r2, 2);
	index_list final_list;
	//index_list final_list = new std::vector<list_entry>(1);
	//initialize_index_list(l1_label_List2, 0, 1);
	#pragma endregion
	//------Listinitialization finished-------//

	 
	int thread = omp_get_thread_num();
	
	int l1_L1_size, l1_L2_size;
	l1_L1_size = l1_L2_size = permutations*permutations / (1 << r1);


	//initialize all necessary matrices and vectors
	#pragma region Declaration
	
	//matrices and vectors
	uint64_t * syn_tilde = new uint64_t[(n - k + 63) / 64];
	uint64_t * t1 = new uint64_t[length_v_r1];
	uint64_t * syn_r1 = new uint64_t[length_v_r1];
	uint64_t * syn_r2 = new uint64_t[length_v_r2];
	uint64_t * syn_l = new uint64_t[length_v_l];
	uint64_t ** T = new uint64_t*[n - k];
	uint64_t ** perm_mat = new uint64_t*[n];
	uint64_t ** Q1 = new uint64_t *[l];
	uint64_t ** Q1_trans = new uint64_t*[k + l];
	uint64_t ** Q = new uint64_t *[n - k];
	uint64_t ** Q_trans = new uint64_t *[k + l];
	uint64_t ** new_parity = new uint64_t *[n-k];
	uint64_t * syn_R = new uint64_t[(n - k + 63) / 64];
	uint64_t * value = new uint64_t[length_v_r1 + 1];
	uint64_t ** trans_perm = new uint64_t *[n];
	uint64_t * test_e_tilde = new uint64_t[(n + 63) / 64];
	#pragma endregion

	#pragma region Initializiation
	for (int i = 0;i < n - k;i++)
	{
		T[i] = new uint64_t[(n - k +63)/64];
		Q[i] = new uint64_t[(k  + l+63) / 64];
		new_parity[i] = new uint64_t[length_v_n];
	}
	for (int i = 0;i < k + l;++i)
	{
		Q1_trans[i] = new uint64_t[length_v_l];
		Q_trans[i] = new uint64_t[(n - k + 63) / 64];
	}

	for (int i = 0;i < n;++i)
		perm_mat[i] = new uint64_t[length_v_n];
	
	for (int i = 0;i<l;++i)
		Q1[i] = new uint64_t[(k + l+63) / 64];
		
	for (int i = 0;i < (n);++i)
		trans_perm[i] = new uint64_t[(n + 63) / 64];
#pragma endregion
	 //------------Initialization finished----------------//


	set_to_identity(perm_mat, n);
	set_to_identity(T, n - k);

	fflush(stdout);
	iterations[thread] = right_perms[thread] = 0;

	copy_matrix(parity_mat, new_parity, n - k, n);
	for (int i = 0;i < (n - k + 63) / 64;++i)
		syn_R[i] = syndrom[i];

	while (!(*stop))
	{
		iterations[thread]++;


		//---permute new_H and convert to quasi-symmetric form---//
		int perms = 0;
		int amount_zero_rows = 0;
		do
		{
			perms++;
			permute(new_parity, perm_mat, n - k, n, thread);
			gaussian_elimination(new_parity, n - k, l, n, k + l, n - k - l, T);
			amount_zero_rows = get_zero_rows(new_parity, l, k + l, n - k - l, thread);
		} while (amount_zero_rows != 0 && !(swap_zero_cols(new_parity, perm_mat, T, amount_zero_rows, l, k + l, thread, n - k, n, n - k - l)));
		//-----------quasi-symmetric form finished-------------//

		//check current permutation on solving-potential (just for testing)
		if (perm_check)
			if (check_perm(perm_mat, trans_perm, test_e_tilde, right_error, false))
				right_perms[thread]++;



		//calculate s~	
		matrix_vector_mult(T, syn_R, syn_tilde, n - k, n - k);

		//---get upper l bit of syndrom and splitt this l bits in lower r1 and upper r2 bit---//
		//get upper l bits of syndrom
		for (int i = 0;i < length_v_l;i++)
			syn_l[i] = syn_tilde[i];
		if (l % 64 != 0)
			syn_l[length_v_l - 1] = (syn_l[length_v_l - 1] >> (64 - (l % 64))) << (64 - (l % 64));
		//split syn_l in upper r2 and lower r1 bits
		for (int i = 0;i < length_v_r2;i++)
			syn_r2[i] = syn_l[i];
		if (r2 % 64 != 0)
		{
			syn_r2[length_v_r2 - 1] = (syn_r2[length_v_r2 - 1] >> (64 - (r2 % 64))) << (64 - (r2 % 64));
			l_shift(syn_l, (r2 % 64), length_v_l);
		}
		for (int i = r2 / 64;i < length_v_l;i++)
			syn_r1[i - (r2 / 64)] = syn_l[i];
		if (r1 % 64 != 0)
			syn_r1[length_v_r1 - 1] = (syn_r1[length_v_r1 - 1] >> (64 - (r1 % 64))) << (64 - (r1 % 64));

		//restore syn_l
		for (int i = 0;i < length_v_l;i++)
			syn_l[i] = 0;
		for (int i = 0;i < length_v_r1;i++)
			syn_l[i] = syn_r1[i];
		r_shift(syn_l, r2, length_v_l);
		for (int i = 0;i < length_v_r2;i++)
			syn_l[i] ^= syn_r2[i];
		//----------------------------splitting finished------------------------------//

		//---save Q*e1 for each e1 in L_l and L_h in the new lists Q_List_l and Q_List_h---//
		#pragma region Create Labellists

		//get both submatrices Q and Q1
		get_submatrix(new_parity, Q, n - k, k + l);
		get_submatrix(new_parity, Q1, l, k + l);
		transpose_matrix(Q1, Q1_trans, 0, l, 0, k + l);
		transpose_matrix(Q, Q_trans, 0, n - k, 0, k + l);

		index_list_iter it_h = Q_List_h->begin();
		index_list_iter it_l = Q_List_l->begin();

		//calculate first result via matrix vector multiplication
		int i = 0;
		matrix_vector_mult_trans(Q1_trans, (*L_l->begin()), (*it_l).value, l + k, l);
		matrix_vector_mult_trans(Q1_trans, (*L_h->begin()), (*it_h).value, l + k, l);
		(*it_h).index[0] = i;
		(*it_l).index[0] = i;
		++it_l, ++it_h, ++i;

		while (it_l != Q_List_l->end())
		{
			//compute Q_1*e1 for permutation list high via changelist (adding two columns of Q_1 to previous result)
			for (int o = 0;o < length_v_l;o++)
				(*it_h).value[o] = (*(it_h - 1)).value[o];
			xor_regs((*it_h).value, Q1_trans[changes[i - 1].first], length_v_l);
			xor_regs((*it_h).value, Q1_trans[changes[i - 1].second], length_v_l);
			(*it_h).index[0] = i;

			//compute Q_1*e1 for permutation list low via changelist (adding two columns of Q_1 to previous result)
			for (int o = 0;o < length_v_l;o++)
				(*it_l).value[o] = (*(it_l - 1)).value[o];
			xor_regs((*it_l).value, Q1_trans[changes[i - 1].first + (l + k) / 2], length_v_l);
			xor_regs((*it_l).value, Q1_trans[changes[i - 1].second + (l + k) / 2], length_v_l);
			(*it_l).index[0] = i;

			++it_h, ++it_l, ++i;
		}
		#pragma endregion
		//----------------------generating Q1e1 lists finished-------------------------//

		//------------Join to level 1-------------//
		//get lower r1 bits of Q1*e1 as list
		get_lower_bits_as_list(Q_List_l, l0_label_List1, value);
		get_lower_bits_as_list(Q_List_h, l0_label_List2, value);

		for (int e = 0;e < rep && !(*stop);++e)
		{
			if (rep != 1)
			{
				//choose t1 and add to label list1
				for (int i = 0;i < length_v_r1;i++)
					t1[i] = uint_dist(mt);
				if ((r1 % 64) != 0)
					t1[length_v_r1 - 1] = (t1[length_v_r1 - 1] >> (64 - (r1 % 64))) << (64 - (r1 % 64));
				add_to_list(l0_label_List1, t1, length_v_r1, l0_label_List1->size());
			}

			//merge join level 0 lists to obtain a list of indices for level 1
			std::sort(l0_label_List1->begin(), l0_label_List1->end(), compare_l0_labels());
			std::sort(l0_label_List2->begin(), l0_label_List2->end(), compare_l0_labels());
			l1_L1_size = merge_join(l0_label_List1, l0_label_List2, l1_label_List1, 2, l0_label_List1->size(), l0_label_List2->size(), thread, L_h, L_l, Q_trans, syn_tilde, perm_mat);

			//add upper bits of s~ to label list1
			add_to_list(l0_label_List1, syn_r1, length_v_r1, l0_label_List1->size());
			std::sort(l0_label_List1->begin(), l0_label_List1->end(), compare_l0_labels());

			//merge join again to obtain second list of indices 
			l1_L2_size = merge_join(l0_label_List1, l0_label_List2, l1_label_List2, 2, l0_label_List1->size(), l0_label_List2->size(), thread, L_h, L_l, Q_trans, syn_tilde, perm_mat);
			//--------------------Join to level 1 finished---------------------//

			//------------Join to final level--------------//
			//get level 1 labels
			get_upper_bits_as_list(Q_List_l, Q_List_h, l1_label_List1, l1_L1_size);
			get_upper_bits_as_list(Q_List_l, Q_List_h, l1_label_List2, l1_L2_size);

			//add syn_r2 to label list 2
			add_to_list(l1_label_List2, syn_r2, length_v_r2, l1_L2_size);

			
			//merge to final list (check elements on the fly)
			std::sort(l1_label_List1->begin(), l1_label_List1->begin() + l1_L1_size, compare_l1_labels());
			std::sort(l1_label_List2->begin(), l1_label_List2->begin() + l1_L2_size, compare_l1_labels());
			merge_join(l1_label_List1, l1_label_List2, final_list, 4, l1_L1_size, l1_L2_size, thread, L_h, L_l, Q_trans, syn_tilde, perm_mat);

			//------------Join to final level finished------------//

			if (temp[thread]||(*stop))
				break;
		}
	}
}

int main(int argc, char *argv[])
{
	uint64_t * syn_R = new uint64_t[(n - k + 63) / 64];
	uint64_t * syndrom = new uint64_t[(n - k + 63) / 64];
	uint64_t ** mat = new uint64_t *[n];
	uint64_t ** parity_mat = new uint64_t *[n - k];
	uint64_t ** T = new uint64_t*[n - k];
	uint64_t ** gen_mat = new uint64_t*[k];
	uint64_t ** R = new uint64_t*[n - k];
	uint64_t ** R_cpy = new uint64_t*[n - k];
	uint64_t ** new_parity = new uint64_t *[n - k];

	//Initialize list for all combinations of errorvectors
	list L_perm_l = new std::vector<uint64_t*>(permutations);
	list L_perm_r = new std::vector<uint64_t*>(permutations);

	list_iter it_r = L_perm_r->begin();
	for (list_iter it_l = L_perm_l->begin(); it_l != L_perm_l->end();++it_l)
	{
		(*it_l) = new uint64_t[length_v_kl];
		(*it_r) = new uint64_t[length_v_kl];
		++it_r;
	}
	signal(SIGINT, my_handler);

	std::ofstream out(out_file,std::ofstream::app);
	std::streambuf *coutbuf = std::cout.rdbuf(); 
	std::cout.rdbuf(out.rdbuf());

	for (int i = 0;i < n;i++)
		mat[i] = new uint64_t[int_per_row];
	for (int i = 0;i < m;++i)
		mat_test[i] = new uint64_t[(k + 63) / 64];
	for (int i = 0;i < n - k;i++)
	{
		T[i] = new uint64_t[(n - k + 63) / 64];
		parity_mat[i] = new uint64_t[length_v_n];
		R[i] = new uint64_t[(n - k + 63) / 64];
		R_cpy[i] = new uint64_t[(n - k + 63) / 64];
		new_parity[i] = new uint64_t[length_v_n];
	}
	for (int i = 0;i < k;++i)
	{
		T_det_sec[i] = new uint64_t[(k + 63) / 64];
		gen_mat[i] = new uint64_t[length_v_n];
	}
	for (int i = 0;i < threads;++i)
	{
		codewords[i] = new uint64_t[(k + 63) / 64];
		secrets[i] = new uint64_t[(k + 63) / 64];
	}


	right_secret = initialize_Oracle(threads);

	(*stop) = false;

	//create generator-, parity-check-matrix and syndrom
	do {
		create_generator_matrix(mat, gen_mat, codeword_w_errors, n, 0);
	} while (!(check_gaussian(gen_mat, 0, 0, k)) || errors > max_error_weight);
	int w = errors;
	transform_to_parity_check(gen_mat, parity_mat);
	matrix_vector_mult(parity_mat, codeword_w_errors, syndrom, n - k, n);

	//for hypothesis test
	set_matrix(mat_test, secret_test, m, 0);
	set_to_identity(T_det_sec, k);
	gaussian_elimination(mat, k, 0, k, 0, k,T_det_sec);

	//precompute set bit table:
	for (int i = 0;i < 256;++i)
		for (int j = 0;j < 8;++j)
			if ((i&(1 << j))>0)
				precomp_bitpositions[i].push_back(7-j);

	//randomize parity-check matrix
	do {
		create_random_mat(R, n - k, n - k);
		copy_matrix(R, R_cpy, n - k, n - k);
		gaussian_elimination(R_cpy, n - k, 0, n - k, 0, n - k);
	} while (!(check_gaussian(R_cpy, 0, 0, n - k)));
	//calculate syn_R = R*syndrom and new_H = R*H
	matrix_vector_mult(R, syndrom, syn_R, n - k, n - k);
	matrix_mult(R, parity_mat, new_parity, n - k, n - k, n);


	//enumerate error vectors
	create_permutation_lists(L_perm_l, L_perm_r);
	omp_set_num_threads(threads);


	std::cout << "\n-------n=" << n << " , k=" << k << " , l=" << l << ", eps=" << eps << " , w=" << w << ", r1=" << r1 << ", p=" << p << ", c=" << rep << "-------\n";
	printf("Error has weigth %d", errors);
	
	//timings
	auto start = std::chrono::high_resolution_clock::now();
	const clock_t START = clock();
	fflush(stdout);
	
#pragma omp parallel for
	for (int i = 0;i <  omp_get_num_threads();i++)
		start_search(L_perm_l, L_perm_r, new_parity, syn_R, stop);

	const double T_ELAPSED = (double)(clock() - START) / CLOCKS_PER_SEC;
	auto end = std::chrono::high_resolution_clock::now();
	
	//output informations
	uint64_t total_iter = 0;
	int total_r_perms = 0;
	for (int i = 0;i < threads;i++)
	{
		total_iter += iterations[i];
		total_r_perms += right_perms[i];
	}
	
	//to screen
	printf("\nwhole search took %d ms", (int) std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
	if(perm_check)
		printf("\nNeeded %d right permutations\n", total_r_perms);

	//to file
	std::cout << "\nwhole search took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";
	std::cout << "\nwhole search took (cpu-time)" << T_ELAPSED << "s";
	std::cout << "\ntotal iterations: " << total_iter;
	std::cout << "\nperms per second (cpu-time) " << total_iter / (T_ELAPSED);
	if (perm_check)
		std::cout << "\nright permutations: " << total_r_perms;

	std::cout.rdbuf(coutbuf);
	
}
