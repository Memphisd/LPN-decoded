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
static const int x = 4;						//Number of bits that shall be zeroed by precomp
static const uint64_t samp_amount = 4000000;		//Number of vectors neeeded for generator matrix
static const int n = 128 + x;					//Bitlength of secret vector
static const double p = 0.125;					//error probability of the orcale
static const int threads = 6;					//Number of threads to use
static const char * out_file = "out.txt";
static const int int_per_row = (n + 63) / 64;	
#endif