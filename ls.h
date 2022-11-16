#pragma once

#include"optimizer.h"
#include"eaoperator.h"

namespace myEC {
	
	//Tabu Search
	class TabuSearch:public Optimizer
	{
	protected:
		int tabulength;
		Solution* tabulist;
		int* age;
		int searchtimes;
		Solution* xnew;
		Solution* lbest;
		Solution* neighbor;

		void (*search_func)(Solution* id, double* clist);

		static void default_search(Solution* id, double* clist)
		{
			int mid = rand() % id->size;
			EAOperator::bit_mutation(id, clist[2 * mid + 1], clist[2 * mid], mid, 0);
		}

		bool tabuCheck(Solution* s)
		{
			for (int i = 0; i < tabulength; i++)
			{
				if (age[i] != 0 && *s == tabulist[i])
					return false;
			}
			return true;
		}

		void updateTabuList(Solution* s)
		{
			for (int i = 0; i < tabulength; i++)
			{
				if (age[i] > 0)
					age[i]--;
			}
			for (int i = 0; i < tabulength; i++)
			{
				if (age[i] == 0)
				{
					tabulist[i].replace(s);
					age[i] = tabulength;
					return;
				}
			}
		}

		virtual void ini()
		{
			double valueBuffer;
			if (solution_ini_func != nullptr)
			{
				solution_ini_func(swarm->result, problem_size);
			}
			else {
				if (model_ini != nullptr)
					model_ini();
				for (int j = 0; j < problem_size; j++)
				{
					if (constrainRangeList[2 * j] != EMPTYVALUE)
					{
						valueBuffer = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
					}
					else {
						valueBuffer = rand01() * 20000 - 10000;
					}
					if (repair != nullptr && !constrain_check(j, valueBuffer))
						valueBuffer = repair(j, valueBuffer);
					if (model_change != nullptr)
						model_change(j, valueBuffer);
					swarm->result[j] = valueBuffer;
				}
			}
			evaluator(swarm);
			gbest->replace(swarm);
		}

	public:
		TabuSearch(int ps, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			int tabulength,
			void (*search_func)(Solution* id, double* clist), int searchtimes,
			int logType)
			:Optimizer(ps, 1, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->tabulength = tabulength;
			this->searchtimes = searchtimes;
			swarm = new Solution(problem_size, objectNumber);
			xnew = new Solution(problem_size, objectNumber);
			lbest = new Solution(problem_size, objectNumber);
			neighbor = new Solution(problem_size, objectNumber);
			gbest = new Solution(problem_size, objectNumber);
			
			
			tabulist = new Solution[tabulength];
			age = new int[tabulength];
			for (int i = 0; i < tabulength; i++)
			{
				tabulist[i].setSize(problem_size, objectNumber);
				age[i] = 0;
			}

			swarm->setCompareFunc(compare_func);
			xnew->setCompareFunc(compare_func);
			lbest->setCompareFunc(compare_func);
			neighbor->setCompareFunc(compare_func);

			if (search_func == nullptr)
				this->search_func = default_search;
			else this->search_func = search_func;
		}

		~TabuSearch()
		{
			delete swarm;
			delete xnew;
			delete lbest;
			delete neighbor;
			delete[] tabulist;
			delete age;
			delete gbest;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int tfes = 1;
			bool alltabu;
			while (tfes < generation)
			{
				//neighbor is the best solution not tabu and lbest is the best solution in this iteration
				alltabu = true;

				neighbor->replace(swarm);
				search_func(neighbor, constrainRangeList);
				evaluator(neighbor);
				if (tabuCheck(neighbor))
					alltabu = false;				
				else lbest->replace(neighbor);
								
				for (int i = 1; i < searchtimes; i++)
				{
					xnew->replace(swarm);
					search_func(xnew, constrainRangeList);
					evaluator(xnew);

					if (tabuCheck(xnew))
					{
						if (alltabu)
							neighbor->replace(xnew);
						alltabu = false;

						if(*xnew < *neighbor)
							neighbor->replace(xnew);
						
					}
					else if (alltabu && *xnew < *lbest)
						lbest->replace(xnew);
				}
				if (alltabu)
				{
					swarm->replace(lbest);
					updateTabuList(lbest);
				}
				else {
					swarm->replace(neighbor);
					updateTabuList(neighbor);
				}
	
				tfes += searchtimes;

				if (*swarm < *gbest)
				{
					gbest->replace(swarm);
					if (detailedLog)
						logFile << tfes << "\t" << gbest->ansprint() << "\n";
					else if (roughLog)
						logFile << tfes << "\t" << gbest->roughprint();
				}				
			}

			exe_time = time(NULL) - exe_time;
		}
	};

	//Simulated Annealing algorithm
	class SA :public Optimizer
	{
	public:
		enum class cooling_type { classic, quick };

	protected:
		double t0;
		double p0;
		double temperature;
		Solution* xnew;
		int searchtimes;
		int balanceIndex;
		double lastfitness;

		void (*search_func)(Solution* id, double* clist);

		static void default_search(Solution* id, double* clist)
		{
			int mid = rand() % id->size;
			EAOperator::bit_mutation(id, clist[2 * mid + 1], clist[2 * mid], mid, 0);
		}

		double accept()
		{
			if (*xnew < *swarm)
				return true;
			else return rand01() < pow(E_CONST, -1 * abs(xnew->fitness[0] - swarm->fitness[0]) / temperature);
		}

		double (*cooling)(int t0, int k);

		static double quickCooling(int t0, int k)
		{
			return t0 / (1 + k);
		}

		static double classicCooling(int t0, int k)
		{
			return t0 / log10(1 + k);
		}

		virtual void ini()
		{
			double valueBuffer;
			if (solution_ini_func != nullptr)
			{
				solution_ini_func(swarm->result, problem_size);
			}
			else {
				if (model_ini != nullptr)
					model_ini();
				for (int j = 0; j < problem_size; j++)
				{
					if (constrainRangeList[2 * j] != EMPTYVALUE)
					{
						valueBuffer = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
					}
					else {
						valueBuffer = rand01() * 20000 - 10000;
					}
					if (repair != nullptr && !constrain_check(j, valueBuffer))
						valueBuffer = repair(j, valueBuffer);
					if (model_change != nullptr)
						model_change(j, valueBuffer);
					swarm->result[j] = valueBuffer;
				}
			}
			evaluator(swarm);
			gbest->replace(swarm);

			balanceIndex = 0;
			lastfitness = swarm->fitness[0];

			if (t0 == EMPTYVALUE)
			{
				//sampling a group of solution
				double fitnessbuffer[10];
				for (int i = 0; i < 10; i++)
				{
					if (model_ini != nullptr)
						model_ini();
					for (int j = 0; j < problem_size; j++)
					{
						if (constrainRangeList[2 * j] != EMPTYVALUE)
						{
							valueBuffer = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							valueBuffer = rand01() * 20000 - 10000;
						}
						if (repair != nullptr && !constrain_check(j, valueBuffer))
							valueBuffer = repair(j, valueBuffer);
						if (model_change != nullptr)
							model_change(j, valueBuffer);
						xnew->result[j] = valueBuffer;
					}
					evaluator(xnew);
					fitnessbuffer[i] = xnew->fitness[0];
				}
				
				if (p0 == EMPTYVALUE)
				{
					double max = fitnessbuffer[0];
					double min = fitnessbuffer[0];
					for (int i = 1; i < 10; i++)
					{
						if (fitnessbuffer[i] > max)
							max = fitnessbuffer[i];
						else if (fitnessbuffer[i] < min)
							min = fitnessbuffer[i];
					}
					t0 = -1 * (max - min) / p0;
				}
				else {
					double mean = 0;
					double var = 0;
					for (int i = 0; i < 10; i++)
						mean += fitnessbuffer[i];
					mean /= 10;
					for (int i = 0; i < 10; i++)
						var += pow(fitnessbuffer[i] - mean, 2);
					t0 = var / 10;
				}
			}
		}

		virtual bool innerBalance()
		{
			if (searchtimes != EMPTYVALUE)
			{
				balanceIndex++;
				if (balanceIndex != searchtimes)
					return false;
				else {
					balanceIndex = 0;
					return true;
				}
			}
			else {
				if (abs((swarm->fitness[0] - lastfitness) / gbest->fitness[0]) < 1e4)
					balanceIndex++;
				else balanceIndex = 0;
				lastfitness = swarm->fitness[0];
				return balanceIndex > 9;
			}
		}

	public:
		SA(int ps, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double t0, cooling_type cool_func, double p0,
			void (*search_func)(Solution* id, double* clist), int searchtimes,
			int logType)
			:Optimizer(ps, 1, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->t0 = t0;
			this->searchtimes = searchtimes;
			this->p0 = p0;
			swarm = new Solution(problem_size, objectNumber);
			xnew = new Solution(problem_size, objectNumber);
			gbest = new Solution(problem_size, objectNumber);
			
			swarm->setCompareFunc(compare_func);
			xnew->setCompareFunc(compare_func);

			if (search_func == nullptr)
				this->search_func = default_search;
			else this->search_func = search_func;
			
			if (cool_func == cooling_type::classic)
				cooling = classicCooling;
			else this->cooling = quickCooling;
		}

		~SA()
		{
			delete swarm;
			delete xnew;
			delete gbest;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int tfes = 1;
			int k = 0;
			temperature = t0;
			while (tfes < generation)
			{
				xnew->replace(swarm);
				search_func(xnew, constrainRangeList);
				evaluator(xnew);

				if (accept())
					swarm->replace(xnew);

				if (innerBalance())
					temperature = cooling(t0, ++k);
				tfes++;

				if (*swarm < *gbest)
				{
					gbest->replace(swarm);
					if (detailedLog)
						logFile << tfes << "\t" << gbest->ansprint() << "\n";
					else if (roughLog)
						logFile << tfes << "\t" << gbest->roughprint();
				}
			}

			exe_time = time(NULL) - exe_time;
		}
	};
}

