#include "Permutation.h"

extern std::vector<std::pair<int, int>> changes;
//extern int l, n, p, eps;

void get_diff(int * a, int * b, int weight, int n)
{
	//get bitposition difference between two consecutive permutations
	for (int i = 0;i < weight;i++)
	{
		std::pair<int, int> change;
		if (a[i] != b[i])
		{
			if (i == weight - 1 || a[i + 1] == b[i + 1])
			{
				change.first = a[i];
				change.second = b[i];
				changes.push_back(change);
				break;
			}
			else
			{
				if (a[i] == b[i + 1])
				{
					change.first = a[i + 1];
					change.second = b[i];
					changes.push_back(change);
					break;
				}
				else
				{
					change.first = a[i];
					change.second = b[i + 1];
					changes.push_back(change);
					break;
				}
			}
		}
	}
}

void copy(int*a, int*b, int weight)
{
	for (int i = 0;i < weight;i++)
		b[i] = a[i];
}

void save(int * bits, uint64_t * value_l, uint64_t *value_r, int weight)
{
	int amount_int = (int)ceil((l + k) / 64.0);
	uint64_t bitmask;

	for (int i = 0;i < amount_int;i++)
	{
		value_l[i] = 0;
		value_r[i] = 0;
	}

	for (int i = 0;i < weight;i++)
	{
		bitmask = 1ULL << (63 - (bits[i] % 64));
		value_l[bits[i] / 64] ^= bitmask;
		bitmask = 1ULL << (63 - ((bits[i] + (l + k) / 2) % 64));
		value_r[(bits[i] + (l + k) / 2) / 64] ^= bitmask;
	}

}

inline int save_entry(int * a, list L_l, list L_r,int weight,int * old, int c)
{
	get_diff(a, old, weight, (l + k) / 2);
	save(a, (*L_l)[c], (*L_r)[c], weight);
	copy(a, old, weight);
	return ++c;
}

void create_permutation_lists(list L_l, list L_r)
{
	int weight = (p / 2 + eps) / 2;
	int * a = new int[weight + 1];
	int * old = new int[weight];
	for (int i = 0;i < weight;i++)
		a[i] = i;
	a[weight] = (l + k) / 2;
	int j=0;
	int c = 0;

	save(a, (*L_l)[c], (*L_r)[c], weight);
	copy(a, old, weight);
	c++;
	bool c_flow = true;
	do
	{
		if (weight % 2 == 1)
			if (a[0] + 1 < a[1])
			{
				a[0] += 1;
				c = save_entry(a, L_l, L_r, weight, old, c);
				continue;
			}
			else if (weight != 1)
				j = 2;
			else
				break;
		else
			if (a[0] > 0)
			{
				a[0] -= 1;
				c = save_entry(a, L_l, L_r, weight, old, c);
				continue;
			}
			else
			{
				j = 2;
				c_flow = false;
			}
		do
		{
			if (c_flow)
			{
				if (a[j - 1] >= j)
				{
					a[j - 1] = a[j - 2];
					a[j - 2] = j - 2;
					c = save_entry(a, L_l, L_r, weight, old, c);
					break;
				}
				else
					j++;
			}
			c_flow = true;
			if ((a[j - 1] + 1) < a[j])
			{
				a[j - 2] = a[j - 1];
				a[j - 1] += 1;
				c = save_entry(a, L_l, L_r, weight, old, c);
				break;
			}
			else
				j++;
		} while (j <= weight);
	} while (j <= weight);
}
