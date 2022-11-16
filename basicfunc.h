#pragma once
#include<algorithm>

namespace myEC {
	const int EMPTYVALUE = -2139062144;
	const int MAX = 1e9;
	const int SINGLEOBJECTBESTNUM = 1;
	const double E_CONST = 2.718281828459045;
	const double PI_CONST = 3.1415926535624;

	double sigmoid(double value)
	{
		return 1 / (1 + pow(E_CONST, value));
	}

	//return a value between 0-1
	double rand01()
	{
		return double(rand()) / RAND_MAX;
	}

	//return a value between 0-1 but not be 0 or 1
	double rand01_()
	{
		return double(rand() + 1) / (RAND_MAX + 2);
	}

	double Eu_distance(double a[], double b[], size_t size = 2)
	{
		double back = 0;
		for (int i = 0; i < size; i++)
		{
			back += pow(a[i] - b[i], 2);
		}

		return sqrt(back);
	}

	template<class T>
	void sort(T* left, T* right)
	{
		T* p_buffer = (T*)malloc(sizeof(T));
		size_t ssize = right - left;

		//insert sort
		for (int i = 0; i < ssize; i++)
		{
			for (int j = 0; j < i; j++)
			{
				//insert individual
				if (left[i] < left[j])
				{
					*p_buffer = left[i];
					for (int j1 = i; j1 > j; j1--)
						left[j1] = left[j1 - 1];
					left[j] = *p_buffer;
					break;
				}
			}
		}

		free(p_buffer);
	}

	struct sortHelper {
		int id;
		double value;

		sortHelper() {}

		sortHelper(int id, double value)
		{
			this->id = id;
			this->value = value;
		}

		bool operator<(const sortHelper& a)const
		{
			return value < a.value;
		}

		bool operator>(const sortHelper& a)const
		{
			return value > a.value;
		}

	};
}
