#pragma once
#ifndef Perm	
#define Perm
#include "Matrixoperations.h"
#include "LPN_Oracle.h"
#include "MMT.h"
typedef std::vector<uint64_t*> * list;
void create_permutation_lists(list L_l, list L_r);

#endif