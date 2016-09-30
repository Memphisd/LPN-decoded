
#pragma once
#ifndef ORACLE
#define ORACLE
#include "Matrixoperations.h"

uint64_t* initialize_Oracle(int threads);
uint64_t * query(int thread,bool gen);
bool scalar_secret(uint64_t * a,uint64_t * res);

#endif