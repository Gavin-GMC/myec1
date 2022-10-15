#define _CRT_SECURE_NO_WARNINGS
#pragma once
#include<iostream>
#include<fstream>
//#include"myec.h"

namespace basicTest1
{
	const int problemsize = 5;

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 0; i < problemsize; i++)
		{
			fitness[0] += pow(solution[i], 2);
		}
	}

	bool constrain(int demensionId, double value)
	{
		return value<10 && value>-10;
	}

	double repair(int demensionId, double value)
	{
		return rand() % 10 - 5;
	}
}

//A popular benchmark of multidimensional knapsack problem
namespace tspKroa100 {
	const int CITYNUMBER = 100;
	double connection[CITYNUMBER][CITYNUMBER];
	double citys[CITYNUMBER][2];
	bool visit[CITYNUMBER] = { 0 };

	int problemsize = CITYNUMBER;
	//int locate_city;

	double dis(double x[2], double y[2])
	{
		return pow(pow((x[0] - y[0]), 2) + pow((x[1] - y[1]), 2), 0.5);
	}

	void get_data()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = 0; j < CITYNUMBER; j++)
				connection[i][j] = -1;

		std::fstream data;
		data.open("data/kroA100.tsp", std::ios::in);
		int buffer;
		for (int i = 0; i < CITYNUMBER; i++)
			data >> buffer >> citys[i][0] >> citys[i][1];

		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = i; j < CITYNUMBER; j++)
				connection[i][j] = connection[j][i] = round(dis(citys[i], citys[j]));

		for (int i = 0; i < CITYNUMBER; i++)
			connection[i][i] = -1;
	}

	double heuristic(int id, double value)
	{
		return 1.0 / connection[id][int(value)];
	}

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 1; i < CITYNUMBER; i++)
		{
			fitness[0] += connection[int(solution[i - 1])][int(solution[i])];
		}
		fitness[0] += connection[int(solution[CITYNUMBER - 1])][int(solution[0])];
	}

	void ini()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			visit[i] = false;
		//locate_city = rand() % CITYNUMBER;
	}

	bool constrain(int demensionId, double value)
	{
		return !visit[int(value)] && connection[demensionId][int(value)] > 0;
	}

	void change(int demensionId, double value)
	{
		if (value < 0 || visit[int(value)])
			return;
		//locate_city = value;
		visit[int(value)] = true;
	 }

	double repair(int demensionId, double value)
	{
		int to = rand() % CITYNUMBER;
		while (visit[to])
		{
			to++;
			if (to == CITYNUMBER)
				to = 0;
		}
		return to;
	}

	void greedy(double solution[], size_t size)
	{
		myEC::sortHelper sortbuffer[CITYNUMBER];
		ini();
		int locate_city = rand() % CITYNUMBER;

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < CITYNUMBER; j++)
			{
				if (constrain(locate_city, j))
				{
					sortbuffer[j] = myEC::sortHelper(j, heuristic(locate_city, j));
				}
				else
					sortbuffer[j] = myEC::sortHelper(j, -10000);
			}
			sort(sortbuffer, sortbuffer + CITYNUMBER, std::greater<myEC::sortHelper>());
			solution[i] = sortbuffer[0].id;
			change(locate_city, solution[i]);
			locate_city = solution[i];
		}
	}
}

namespace tspEil51 {
	const int CITYNUMBER = 51;
	double connection[CITYNUMBER][CITYNUMBER];
	double citys[CITYNUMBER][2];
	bool visit[CITYNUMBER] = { 0 };

	int problemsize = CITYNUMBER;
	//int locate_city;

	double dis(double x[2], double y[2])
	{
		return pow(pow((x[0] - y[0]), 2) + pow((x[1] - y[1]), 2), 0.5);
	}

	void get_data()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = 0; j < CITYNUMBER; j++)
				connection[i][j] = -1;

		std::fstream data;
		data.open("data/eil51.tsp", std::ios::in);
		int buffer;
		for (int i = 0; i < CITYNUMBER; i++)
			data >> buffer >> citys[i][0] >> citys[i][1];

		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = i; j < CITYNUMBER; j++)
				connection[i][j] = connection[j][i] = round(dis(citys[i], citys[j]));
	}

	double heuristic(int id, double value)
	{
		return 1.0 / connection[id][int(value)];
	}

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 1; i < CITYNUMBER; i++)
		{
			fitness[0] += connection[int(solution[i - 1])][int(solution[i])];
		}
		fitness[0] += connection[int(solution[CITYNUMBER - 1])][int(solution[0])];
	}

	void ini()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			visit[i] = false;
	}

	bool constrain(int demensionId, double value)
	{
		return !visit[int(value)] && connection[demensionId][int(value)] > 0;
	}

	void change(int demensionId, double value)
	{
		if (value < 0 || visit[int(value)])
			return;
		visit[int(value)] = true;
	}

	double repair(int demensionId, double value)
	{
		int to = rand() % CITYNUMBER;
		while (visit[to])
		{
			to++;
			if (to == CITYNUMBER)
				to = 0;
		}
		return to;
	}

	void greedy(double solution[], size_t size)
	{
		myEC::sortHelper sortbuffer[CITYNUMBER];
		ini();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < CITYNUMBER; j++)
			{
				if (constrain(i, j))
				{
					sortbuffer[j] = myEC::sortHelper(j, heuristic(i, j));
				}
				else
					sortbuffer[j] = myEC::sortHelper(j, -10000);
			}
			sort(sortbuffer, sortbuffer + CITYNUMBER, std::greater<myEC::sortHelper>());
			solution[i] = sortbuffer[0].id;
			change(i, solution[i]);
		}
	}
}

namespace ACStspKroa100 {
	const int CITYNUMBER = 100;
	const int SWARMSIZE = 10;

	double connection[CITYNUMBER][CITYNUMBER];
	double citys[CITYNUMBER][2];
	bool visit[SWARMSIZE][CITYNUMBER] = { 0 };

	int problemsize = CITYNUMBER;
	int ant_now;
	//int locate_city;

	double dis(double x[2], double y[2])
	{
		return pow(pow((x[0] - y[0]), 2) + pow((x[1] - y[1]), 2), 0.5);
	}

	void get_data()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = 0; j < CITYNUMBER; j++)
				connection[i][j] = -1;

		std::fstream data;
		data.open("data/kroA100.tsp", std::ios::in);
		int buffer;
		for (int i = 0; i < CITYNUMBER; i++)
			data >> buffer >> citys[i][0] >> citys[i][1];

		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = i; j < CITYNUMBER; j++)
				connection[i][j] = connection[j][i] = round(dis(citys[i], citys[j]));
	}

	double heuristic(int id, double value)
	{
		return 1.0 / connection[id][int(value)];
	}

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 1; i < CITYNUMBER; i++)
		{
			fitness[0] += connection[int(solution[i - 1])][int(solution[i])];
		}
		fitness[0] += connection[int(solution[CITYNUMBER - 1])][int(solution[0])];
	}

	void ini()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			for (int a = 0; a < SWARMSIZE; a++)
				visit[a][i] = false;

		ant_now = 0;
	}

	bool constrain(int demensionId, double value)
	{
		return !visit[ant_now][int(value)] && connection[demensionId][int(value)] > 0;
	}

	void change(int demensionId, double value)
	{
		if (value < 0 || visit[ant_now][int(value)])
			return;

		//locate_city = value;
		visit[ant_now][int(value)] = true;
		ant_now++;

		if (ant_now == SWARMSIZE)
			ant_now = 0;
	}

	double repair(int demensionId, double value)
	{
		int to = rand() % CITYNUMBER;
		while (visit[ant_now][to])
		{
			to++;
			if (to == CITYNUMBER)
				to = 0;
		}
		return to;
	}

	void greedy(double solution[], size_t size)
	{
		myEC::sortHelper sortbuffer[CITYNUMBER];
		ini();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < CITYNUMBER; j++)
			{
				if (constrain(i, j))
				{
					sortbuffer[j] = myEC::sortHelper(j, heuristic(i, j));
				}
				else
					sortbuffer[j] = myEC::sortHelper(j, -10000);
			}
			sort(sortbuffer, sortbuffer + CITYNUMBER, std::greater<myEC::sortHelper>());
			solution[i] = sortbuffer[0].id;
			change(i, solution[i]);
		}
	}
}

//A popular benchmark of multidimensional knapsack problem
namespace mkpGK01 {
	int m, n, optimal;
	int* profit;
	int* resource;
	int* usedResource;
	int* capacity;
	double* h_matrix;

	int problemsize;
	
	void get_data()
	{
		std::fstream data;
		data.open("data/gk/gk01.dat", std::ios::in);

		data >> n >> m >> optimal;

		profit = new int[n];
		resource = new int[n * m];
		usedResource = new int[m];
		capacity = new int[m];
		h_matrix = new double[n];

		for (int i = 0; i < n; i++)
			data >> profit[i];

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				data >> resource[i * n + j];


		for (int i = 0; i < m; i++)
			data >> capacity[i];

		for (int i = 0; i < m; i++)
			usedResource[i] = 0;

		problemsize = n;

		double total;
		for (int i = 0; i < n; i++)
		{
			total = 0;
			for (int j = 0; j < m; j++)
				total += (double(resource[j * n + i]) / capacity[j]);
			h_matrix[i] = profit[i] / total;
		}
	}

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 0; i < n; i++)
			if (solution[i] != 0)
				fitness[0] += profit[i];
		fitness[0] *= -1;
	}

	void ini()
	{
		for (int i = 0; i < m; i++)
			usedResource[i] = 0;
	}

	bool constrain(int demensionId, double value)
	{
		if (value == 0)
			return true;
		for (int i = 0; i < m; i++)
			if (usedResource[i] + resource[i * n + demensionId] > capacity[i])
				return false;
		return true;
	}

	void change(int demensionId, double value)
	{
		if (value <= 0)
			return;
		for (int i = 0; i < m; i++)
			usedResource[i] += resource[i * n + demensionId];
	}

	double repair(int demensionId, double value)
	{
		return 0;
	}

	double heuristic(int id, double value)
	{
		if (value == 0)
			return 0;
		else return h_matrix[id];
	}

	void remove_data()
	{
		delete[] profit;
		delete[] resource;
		delete[] usedResource;
		delete[] capacity;
		delete[] h_matrix;
	}

	void greedy(double solution[], size_t size)
	{
		myEC::sortHelper sortbuffer[100];
		ini();

		for (int j = 0; j < 100; j++)
			sortbuffer[j] = myEC::sortHelper(j, h_matrix[j]);
		sort(sortbuffer, sortbuffer + 100, std::greater<myEC::sortHelper>());

		for (int i = 0; i < size; i++)
		{
			if (constrain(i, 1))
				solution[i] = 1;
			else solution[i] = 0;
			change(i, solution[i]);
		}
	}
}

//A benchmark of single object optiization
namespace CEC14{
	void cec14_test_func(double*, double*, int, int, int);

	double* OShift, * M, * y, * z, * x_bound;
	int ini_flag = 0, n_flag, func_flag, * SS;

	bool constrain(int demensionId, double value)
	{
		return value<100 && value>-100;
	}



	void evaluate(double solution[], double fitness[])
	{
		cec14_test_func(solution, fitness, n_flag, 1, func_flag);
	}

	double repair(int demensionId, double value)
	{
		return rand() % 200 - 100;
	}

	//d for 2,10,20,30,50,100. func_num for 1-30
	/*void test(int d, int func_num)
	{
		ini_flag = 0;
		n_flag = d;
		func_flag = func_num;

		myEC::OptimizerBuilder b(myEC::algorithm::PSO, d, evaluate, 100, d * 10000);
		b.setConstrainFunc(constrain);
		b.setRepaireFunc(repair);
		for (int i = 0; i < d; i++)
			b.setDemensionRange(i, -100, 100);
		//b.setGASelection(GA::selection::championship);
		b.setLogType(1);

		myEC::Optimizer* o = b.build();
		o->exe();


		myEC::Solution* s = o->getGbest();
		std::string ans = s->ansprint();

		std::cout << ans << std::endl;
	}*/
}
