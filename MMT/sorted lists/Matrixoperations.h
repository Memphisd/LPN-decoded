#pragma once

#ifndef MatMult
#define MatMult

#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <random>
#include <chrono>
#include <omp.h>
#include <iostream>
#include <list>
#include <fstream>
#include <cstring>
#include <vector>
#include <algorithm>
#include "LPN_Oracle.h"
#include "MMT.h"
#include <signal.h>


void matrix_vector_mult(uint64_t ** A, uint64_t * v, uint64_t * res, int row_A, int col_A);
void matrix_mult(uint64_t** A, uint64_t ** B, uint64_t ** res, int col_A,int row_A,int col_B);
void transpose_matrix(uint64_t ** A, uint64_t ** res, int row_start, int row_end, int col_start, int col_end);
void set_to_identity(uint64_t ** mat, int rows);
void gaussian_elimination(uint64_t **mat, int row, int row_start, int col, int col_start, int n, uint64_t ** T);
void gaussian_elimination(uint64_t **mat, int row, int row_start, int col, int col_start, int n);
void create_generator_matrix(uint64_t ** mat, uint64_t ** gen_mat, uint64_t * codew_errors, int dim, int thread);
void transform_to_parity_check(uint64_t ** gen_mat, uint64_t ** parity_mat);
void get_submatrix(uint64_t ** mat, uint64_t ** res, int row_end, int col_end);
bool check_gaussian(uint64_t ** mat, int row_start, int col_start, int n);
void show_bin(uint64_t num);
void print_matrix(uint64_t **mat, int col, int row);
void permute(uint64_t ** mat, uint64_t ** P, int rows, int cols,int thread);
void copy_matrix(uint64_t ** Src, uint64_t ** Dest, int rows, int cols);
void create_random_mat(uint64_t ** mat, int rows,int cols);
void get_col(uint64_t ** mat, int col, int rows, uint64_t * vec);
void matrix_vector_mult_trans(uint64_t ** A, uint64_t * v, uint64_t * res, int row_A, int col_A);
bool swap_zero_cols(uint64_t ** mat, uint64_t ** perm_mat, uint64_t ** T, int amount, int row_start, int col_start, int thread, int rows, int cols, int n);
int get_zero_rows(uint64_t ** mat, int row_start, int col_start, int n, int thread);
void set_matrix(uint64_t** mat, bool* secret, int rows, int thread);
bool test_guess(uint64_t** mat_test, bool* secret_test, uint64_t * guess);

#endif 