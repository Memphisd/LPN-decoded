#pragma once

#ifndef BJMM
#define BJMM
#include "Matrixoperations.h"
#include "LPN_Oracle.h"
#include "Permutation.h"

//LPN parameters
static const int k = 128;
static const double tau = 0.125;

//optimization parameters
static const int l = 24;
static const int n = 768;
static const int threads = 31;

//tree parameters
static const int r1 = 8;
static const int rep = 8;

//weight parameters
static const int eps = 0;
static const int p = 8;
static const int max_error_weight = 111;

//parameters for hypothesis-testing
static const int m = 92;
static const int c = 19;

static const char * out_file = "out.txt";

static const bool perm_check = false;

//arraylengths
static const int r2 = l - r1;
static const int length_v_l = (l + 63) / 64;
static const int length_v_r1 = (r1 + 63) / 64;
static const int length_v_r2 = (r2 + 63) / 64;
static const int length_v_n = (n + 63) / 64;
static const int length_v_kl = (k + l + 63) / 64;
static const int int_per_row = (k + 63) / 64;
static const int n_k_length = (n - k + 63) / 64;

struct list_entry
{
	int * index;
	uint64_t * value;
};


typedef std::vector<uint64_t*> * list;
typedef std::vector<uint64_t*>::iterator list_iter;
typedef std::vector<list_entry> * index_list;
typedef std::vector<list_entry>::iterator index_list_iter;

int choose(int n, int k);
void xor_regs(uint64_t * A, uint64_t * B, int regs);
void r_shift(uint64_t * v, int shift_amount, int regs);
void l_shift(uint64_t * v, int shift_amount, int regs);

int weigth(uint64_t v);
int merge_join(index_list L1, index_list L2, index_list merged_List, int indices, int L1_size, int L2_size, int thread, list L_h, list L_l, uint64_t ** Q_trans, uint64_t * syn_tilde, uint64_t ** perm_mat);
void initialize_index_list(index_list L, int value_length, int amount_indices);
void get_lower_bits_as_list(index_list L, index_list L_low, uint64_t * value);
void get_upper_bits_as_list(index_list L_l, index_list L_h, index_list level1_list, int size);
void start_search(list L_h, list L_l, uint64_t ** Q,uint64_t * syndrom,bool * stop);
bool scalar(uint64_t * a, uint64_t * b);
#endif