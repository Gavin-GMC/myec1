 #pragma once
#include<ctime>
#include<cstring>
#include<string>
#include<algorithm>
#include<fstream>
#include<iostream>
#include"basicfunc.h"

namespace myEC {
	class Solution
	{
	protected:
		bool (*compare_func)(double f1[], double f2[], size_t size) = nullptr;
	public:

		int size = 0;
		double* result;
		double* fitness;
		int object_number = 0;

		Solution() {}

		Solution(int _size, int object_number = 1)
		{
			size = _size;
			result = new double[size];
			this->object_number = object_number;
			fitness = new double[object_number];
		}

		~Solution()
		{
			delete[] result;
			delete[] fitness;
		}

		virtual void setSize(int _size, int object_number = 1)
		{
			if (size != _size)
			{
				size = _size;
				if (!result)
					delete[] result;
				result = new double[size];
			}
			if (this->object_number != object_number)
			{
				this->object_number = object_number;
				if (!fitness)
					delete[] fitness;
				fitness = new double[object_number];
			}
		}

		void setCompareFunc(bool (*compare_func)(double f1[], double f2[], size_t size))
		{
			this->compare_func = compare_func;
		}

		void replace(Solution* object)
		{
			memcpy(fitness, object->fitness, object_number * sizeof(double));
			memcpy(result, object->result, size * sizeof(double));
		}

		void replace(double* result, double* fitness)
		{
			memcpy(this->fitness, fitness, object_number * sizeof(double));
			memcpy(this->result, result, size * sizeof(double));
		}

		std::string ansprint()
		{
			std::string back = "";
			for (int i = 0; i < object_number; i++)
				back += (std::to_string(i) + ":\t" + std::to_string(fitness[i]) + "\t");
			back += "\ndetail:\n";
			for (int i = 0; i < size; i++)
				back += (std::to_string(i) + ":\t" + std::to_string(result[i]) + "\n");
			return back;
		}

		std::string roughprint()
		{
			std::string back = "";
			for (int i = 0; i < object_number; i++)
				back += (std::to_string(i) + ":\t" + std::to_string(fitness[i]) + "\t");
			back += "\n";
			return back;
		}

		bool operator<(const Solution& a)const
		{
			if (compare_func != nullptr)
				return compare_func(fitness, a.fitness, object_number);
			else {
				for (int i = 0; i < object_number; i++)
					if (this->fitness[i] > a.fitness[i])
						return false;
				for (int i = 0; i < object_number; i++)
					if (this->fitness[i] < a.fitness[i])
						return true;
				return false;
			}
		}

		bool operator==(const Solution& a)const
		{
			return memcmp(this->result, a.result, size * sizeof(double)) == 0;
		}

	};

	class Evaluater
	{
	private:
		void (*evaluate_func)(double solution[], double fitness[]);
	public:
		Evaluater() {}

		Evaluater(void (*evaluate_func)(double solution[], double fitness[])) 
		{
			this->evaluate_func = evaluate_func;
		}

		~Evaluater() {}

		void setEvaluateFunc(void (*evaluate_func)(double solution[], double fitness[]))
		{
			this->evaluate_func = evaluate_func;
		}

		void operator() (double* s, double* f) const
		{
			evaluate_func(s, f);
		}

		void operator() (Solution* s) const
		{
			evaluate_func(s->result, s->fitness);
		}
	};

	class Optimizer
	{
	protected:
		int swarm_size;
		int problem_size;
		int objectnum;
		int generation;

		Solution* swarm;
		Solution* gbest;
		time_t exe_time;
		int gbestsize = 1;

		double* constrainRangeList;
		
		std::ofstream logFile;

		Evaluater evaluator;
		bool (*compare_func)(double f1[], double f2[], size_t size);
		void (*solution_ini_func)(double solution[], size_t size);
		void (*model_ini)(void);
		bool (*check_func)(int demensionId, double value);
		void (*model_change)(int demensionId, double value);
		double (*repair)(int demensionId, double value);
		double (*heuristic_func)(int demensionId, double value);

		bool roughLog = false;
		bool detailedLog = false;

		bool default_constrain(int demensionId, double value) const
		{
			return constrainRangeList[2 * demensionId] == EMPTYVALUE ||
				(constrainRangeList[2 * demensionId] <= value && value <= constrainRangeList[2 * demensionId + 1]);
		}

		bool constrain_check(int demensionId, double value)
		{
			if (check_func != nullptr)
				return check_func(demensionId, value);
			else return (demensionId < problem_size&& default_constrain(demensionId, value));
		}

		virtual void solutioncheck(Solution* s)
		{
			if (model_ini != nullptr)
				model_ini();
			for (int i = 0; i < problem_size; i++)
			{
				if (!constrain_check(i, s->result[i]))
					s->result[i] = repair(i, s->result[i]);
				if (model_change != nullptr)
					model_change(i, s->result[i]);
			}
		}
	public:

		Optimizer() {}

		Optimizer(int ps, int ss, int gen, void (*evaluate_func)(double solution[], double fitness[]), int objectNum,
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			int logType)
		{
			problem_size = ps;
			swarm_size = ss;
			generation = gen;
			this->objectnum = objectNum;

			evaluator = Evaluater(evaluate_func);

			this->solution_ini_func = solution_ini_func;
			this->model_ini = model_ini;
			this->check_func = constrain_check;
			this->model_change = model_change;
			this->repair = repair_func;

			constrainRangeList = new double[problem_size * 2];

			if(logType==1)
				this->roughLog = true;
			else if (logType==2)
				this->detailedLog = true;

			if (detailedLog)
				logFile.open("log/DLog_" + std::to_string(time(NULL)) + ".txt");
			else if (roughLog)
				logFile.open("log/RLog_" + std::to_string(time(NULL)) + ".txt");
		}

		~Optimizer()
		{
			delete[]constrainRangeList;
			if (logFile)
				logFile.close();
		}

		virtual void exe() = 0;

		Solution* getGbest()
		{
			return gbest;
		}

		Solution* getSwarm()
		{
			return swarm;
		}

		int getGbestSize()
		{
			return gbestsize;
		}

		time_t get_exe_time()
		{
			return exe_time;
		}

		void constructConstrainRangeList(double* list, size_t si)
		{
			memcpy(constrainRangeList, list, si * sizeof(double));
		}
	};
}

