#pragma once
#pragma once

#ifndef BKW
#define BKW
#define _GLIBCXX_PARALLEL

#include <chrono>
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
#include <omp.h>
#include <algorithm>
#include <signal.h>
#include "LPN_Oracle.h"
#include <parallel/algorithm>

using namespace std;

//LPN-Parameters
static const int n = 113;	
static const double tau = 0.25;

//d = bits to eliminate 
static const int d = 93;

//starting list size
static const uint64_t start_amount =6654011676;

//optimization parameters
static const int a = 3;
static const int b = 31;
static const int threads =64;



static const int int_per_row = n/ 64+1;


#endif
