#pragma once

#ifndef Brute
#define Brute

#include <stdio.h>
#include <stdint.h>
#include <stdio.h>
#include <random>
#include <chrono>
#include <omp.h>
#include <iostream>
#include <list>
#include <fstream>
#include <cstring>

#include "LPN_Oracle.h"
typedef std::list<int>::iterator iter_int;
typedef std::list<uint64_t *>::iterator iter_uint;

//use vector pool?
static const bool prange = false;

//LPN parameters
static const int n = 100;
static const double p = 0.125;

//parameter of fast-elimination [(0,1) = default]
static const int pre_gen_mat_size = 0;
static const uint64_t pre_gen_iter = 1;

//hypothesis-testing parameters
static const int test_amount = 92;
static const int threshold = 19;

//number of executions
static const int runs = 1;	

static const int threads = 6;				

static const char * out_file = "out.txt";
static const int bitmask_64 = (1 << 6) - 1;
static const uint64_t max_guesses = 10;				
static const int int_per_row = (n + 63) / 64;	
static const int vec_sec_length = n / 64 + 1;

inline void xor_entry(uint64_t** mat, int i, int j);

void gaussian_elimination(uint64_t **mat, int n);

void set_matrix(uint64_t** mat, int rows, int thread,int start);

void set_secret(int setting,std::list<int>* zero_ind, uint64_t ** mat,uint64_t * guess);

void determine_secret(uint64_t ** mat,std::list<int>* zero_ind,std::list<uint64_t *>* secrets,int thread);

uint64_t* test_guess(std::list<uint64_t *>* sec,uint64_t** mat,bool* secret,int guess_count);

bool scalar(uint64_t * a,uint64_t * b);

void print_matrix(uint64_t **mat, int col, int row);

void show_bin(uint64_t num);

#endif