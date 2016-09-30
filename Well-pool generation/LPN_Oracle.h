
#pragma once
#ifndef ORACLE
#define ORACLE
#include "Sampling_Test.h"

uint64_t* initialize_Oracle(int threads);
uint64_t * query(int thread);
bool scalar_secret(uint64_t * a,uint64_t * res);


#endif