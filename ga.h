#pragma once
#include"optimizer.h"
#include"multiobject.h"
#include"eaoperator.h"

namespace myEC {
	//Classic GA for single Object and discrete optimization
	class GA : public Optimizer
	{
	protected:
		Solution* offspring;
		bool elitist_strategy;

		Mutation mutation_func;
		Crossover crossover_func;
		Selection selection_func;

		virtual void ini()
		{
			double valueBuffer;
			for (int i = 0; i < swarm_size; i++)
			{
				if (solution_ini_func != nullptr)
				{
					solution_ini_func(swarm[i].result, problem_size);
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
						swarm[i].result[j] = valueBuffer;
					}
				}
				evaluator(&swarm[i]);
			}
			gbest->replace(&swarm[0]);
			for (int i = 1; i < swarm_size; i++)
			{
				if (swarm[i] < *gbest)
					gbest->replace(&swarm[i]);
			}
		}

	public:

		GA(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			crossoverParameter cp, mutationParameter mp,
			selectionParameter sp, bool elitist_strategy,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			swarm = new Solution[swarm_size * 2];
			gbest = new Solution(problem_size, objectNumber);

			offspring = swarm + swarm_size;

			for (int i = 0; i < 2 * swarm_size; i++)
			{
				swarm[i].setSize(ps, objectNumber);
				swarm[i].setCompareFunc(compare_func);;
			}

			crossover_func = Crossover(cp);
			mutation_func = Mutation(mp, constrainRangeList);
			selection_func = Selection(sp);

			this->elitist_strategy = elitist_strategy;
		}

		~GA()
		{
			if (offspring < swarm)
				swarm = offspring;
			delete[] swarm;
			delete[] gbest;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int fes = swarm_size;
			if (elitist_strategy)
				fes--;
			for (int i = swarm_size; i < generation; i += fes)
			{
				selection_func(swarm, swarm_size, offspring);
				std::swap(swarm, offspring);

				for (int j = 0; j < swarm_size - 1; j += 2)
				{
					crossover_func(&swarm[j], &swarm[j + 1]);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					mutation_func(&swarm[j]);
				}

				if (repair != nullptr)
				{
					for (int j = 0; j < swarm_size; j++)
						solutioncheck(&swarm[j]);
				}

				if (elitist_strategy)
				{
					swarm[fes].replace(gbest);
					for (int j = 0; j < fes; j++)
						evaluator(&swarm[j]);
				}
				else {
					for (int j = 0; j < swarm_size; j++)
						evaluator(&swarm[j]);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					if (swarm[j] < *gbest)
					{
						gbest->replace(&swarm[j]);
						if (detailedLog)
							logFile << i / swarm_size << "\t" << gbest->ansprint() << "\n";
						else if (roughLog)
							logFile << i / swarm_size << "\t" << gbest->roughprint();
					}
				}
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//A Fast and Elitist Multiobjective Genetic Algorithm: NSGA - II
	class NSGA2 : public GA
	{
	private:
		double* indicator;//1 for rank and 2 for the negative of crowding distance
		Solution* PF;

	public:
		NSGA2(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			crossoverParameter cp, mutationParameter mp,
			selectionParameter sp,
			int logType)
			:GA(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, cp, mp,
				sp, false, logType)
		{
			indicator = new double[2 * swarm_size * 2];
			PF = new Solution[PFSIZE];
			gbestsize = 0;

			for (int i = 0; i < PFSIZE; i++)
			{
				PF[i].setSize(problem_size, objectNumber);
			}

			delete gbest;
			gbest = PF;
		};

		~NSGA2()
		{
			delete[] indicator;
			PF = nullptr;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int fes = swarm_size;
			bool update = false;
			for (int i = swarm_size; i < generation; i += fes)
			{
				//selection_func(swarm, swarm_size, offspring, 2, 2 * swarm_size, indicator);
				selection_func(swarm, swarm_size, offspring);

				for (int j = 0; j < swarm_size - 1; j += 2)
				{
					crossover_func(&offspring[j], &offspring[j + 1]);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					mutation_func(&offspring[j]);
				}

				if (repair != nullptr)
				{
					for (int j = 0; j < swarm_size; j++)
						solutioncheck(&offspring[j]);
				}

				for (int j = 0; j < swarm_size; j++)
					evaluator(&offspring[j]);

				MO::fastNonDominatedSort<Solution>(swarm, swarm_size * 2, indicator);
				MO::buildPartialOrder<Solution>(swarm, swarm_size * 2, indicator, 2, 0);

				int beforeid = 0;
				int rank_now = 1;
				for (int j = 0; j < swarm_size * 2; j++)
				{
					if (indicator[j] != rank_now)
					{
						MO::normalizeCrowdDistance<Solution>(swarm + beforeid, j - beforeid, indicator + 2 * swarm_size + beforeid);
						beforeid = j;
						rank_now = indicator[j];
					}
				}
				MO::normalizeCrowdDistance<Solution>(swarm + beforeid, swarm_size * 2 - beforeid, indicator + 2 * swarm_size + beforeid);

				MO::buildPartialOrder(swarm, swarm_size * 2, indicator, 2);

				for (int j = 0; j < swarm_size * 2; j++)
				{
					if (indicator[j] != 1)
						break;
					if (MO::paretoFrontUpdate(&swarm[j], PF, gbestsize))
						update = true;
				}

				if ((i / fes) % 20 == 0 && update)
				{
					if (detailedLog)
					{
						for (int j = 0; j < gbestsize; j++)
						{
							logFile << i / swarm_size << "\t" << gbest[j].ansprint() << "\n";
						}
						logFile << "\n";
					}
					else if (roughLog)
					{
						for (int j = 0; j < gbestsize; j++)
						{
							logFile << i / swarm_size << "\t" << gbest[j].roughprint();
						}
						logFile << "\n";
					}
					update = false;
				}
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//A Multiobjective Evolutionary Algorithm Based on Decomposition: MOEA/D
	class MOEAD :public GA
	{
	public:
		enum class D_approach { weight_sum, tchebycheff, PBI };

	private:
		D_approach D_type;
		int neighbor_size;
		double cita = 1;

		double* lamda;
		int* neighborid;
		Solution* PF;
		double* z;
		
		void buildlamda()
		{
			if (objectnum < 2)
				return;

			int H, supportsize;
			//compute the exact H for problem decomposition
			H = 1;
			supportsize = objectnum;
			while (supportsize < swarm_size)
			{
				H++;
				supportsize = supportsize * (H + objectnum - 1) / H;
			}
			//the number of subproblems shouble be ignored
			int nig = supportsize - swarm_size;
			double pig = 1.0 / nig;

			int* lamda_H = new int[objectnum - 1];
			for (int i = 0; i < objectnum - 2; i++)
				lamda_H[i] = 0;
			lamda_H[objectnum - 2] = -1;
			int sid = 0;

			int index;

			for (int i = 0; i < supportsize; i++)
			{
				//find next decomposition in dictionary order
				index = objectnum - 2;
				lamda_H[index]++;
				while (lamda_H[index] == H + 1)
					lamda_H[--index]++;
				for (int i = index + 1; i < objectnum - 1; i++)
					lamda_H[i] = lamda_H[index];

				//ignore the needless subproblem based on probability
				if (nig > 0 && rand01() < pig)
				{
					nig--;
					continue;
				}
				//build the subproblem based on decomposition
				lamda[sid * objectnum] = double(lamda_H[0]) / H;
				for (int i = 1; i < objectnum - 1; i++)
				{
					lamda[sid * objectnum + i] = double(lamda_H[i] - lamda_H[i - 1]) / H;
				}
				lamda[sid * objectnum + objectnum - 1] = double(H - lamda_H[objectnum - 2]) / H;
				sid++;
				if (sid == swarm_size)
					break;
			}

			delete[] lamda_H;
		}

		void buildneighbor()
		{
			sortHelper* sortbuffer = new sortHelper[swarm_size];
			for (int i = 0; i < swarm_size; i++)
			{
				for (int j = 0; j < swarm_size; j++)
				{
					sortbuffer[j].id = j;
					if (i == j)
					{
						sortbuffer[j].value = MAX;
					}
					else {
						sortbuffer[j].value = Eu_distance(lamda + i * objectnum, lamda + j * objectnum, objectnum);
					}
				}
				sort(sortbuffer, sortbuffer + swarm_size);
				for (int j = 0; j < neighbor_size; j++)
					neighborid[i * neighbor_size + j] = sortbuffer[j].id;
			}
			delete[] sortbuffer;
		}

		void z_update(Solution* s, int ssize)
		{
			switch (D_type)
			{
			case D_approach::tchebycheff:
			case D_approach::PBI:
				for (int i = 0; i < ssize; i++)
				{
					for (int j = 0; j < objectnum; j++)
						if (s[i].fitness[j] < z[j])
							z[j] = s[i].fitness[j];
				}
				break;
			default:
				break;
			}
		}

		double sub_fitness(int id)
		{
			int lid =  id % swarm_size;
			double back = EMPTYVALUE;
			switch (D_type)
			{
			case D_approach::weight_sum:
				back = 0;
				for (int i = 0; i < objectnum; i++)
				{
					back += lamda[lid * objectnum + i] * swarm[id].fitness[i];
				}
				break;
			case D_approach::tchebycheff:
				back = -1 * MAX;
				double buf;
				for (int i = 0; i < objectnum; i++)
				{
					buf = lamda[lid * objectnum + i] * abs(swarm[id].fitness[i] - z[i]);
					if (buf > back)
						back = buf;
				}
				break;
			case D_approach::PBI:
			{
				double d1 = 0, d2 = 0;
				//calculate d1
				for (int i = 0; i < objectnum; i++)
				{
					d1 += lamda[lid * objectnum + i] * abs(swarm[id].fitness[i] - z[i]);
				}
				double dlamda = 0;
				for (int i = 0; i < objectnum; i++)
					dlamda += pow(lamda[lid * objectnum + i], 2);
				d1 /= sqrt(dlamda);
				//caculate d2
				for (int i = 0; i < objectnum; i++)
					d2 += pow(swarm[id].fitness[i] - (z[i] - d1 * lamda[lid * objectnum + i]), 2);
				d2 = sqrt(d2);
				back = d1 + cita * d2;
				break;
			}
			default:
				break;
			}
			return back;
		}

		virtual void ini()
		{
			GA::ini();
			buildlamda();
			buildneighbor();
			for (int i = 0; i < objectnum; i++)
				z[i] = MAX;
			z_update(swarm, swarm_size);
		}

	public:
		MOEAD(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			crossoverParameter cp, mutationParameter mp,
			D_approach decomposition_type, int neighbor_size,
			int logType)
			:GA(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, cp, mp,
				default_sp, false, logType)
		{
			D_type = decomposition_type;
			this->neighbor_size = neighbor_size;

			lamda = new double[swarm_size * objectnum];
			neighborid = new int[swarm_size * neighbor_size];
			PF = new Solution[PFSIZE];
			z = new double[objectnum];
			gbestsize = 0;

			for (int i = 0; i < PFSIZE; i++)
			{
				PF[i].setSize(problem_size, objectNumber);
			}

			delete gbest;
			gbest = PF;
		};

		~MOEAD()
		{
			delete[] lamda;
			delete[] neighborid;
			PF = nullptr;
			delete[] z;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int fes = swarm_size;
			bool update = false;
			Solution b1;
			b1.setSize(problem_size, objectnum);
			int bid0, bid1;
			//double sub_f;
			for (int i = swarm_size; i < generation; i += fes)
			{
				for (int j = 0; j < swarm_size; j++)
				{
					//selection
					bid0 = rand() % neighbor_size;
					do
					{
						bid1 = rand() % neighbor_size;
					} while (bid1 == bid0);
					offspring[j].replace(&swarm[neighborid[j * neighbor_size + bid0]]);
					b1.replace(&swarm[neighborid[j * neighbor_size + bid1]]);

					//crossover
					crossover_func(&offspring[j], &b1);

					//mutation
					mutation_func(&offspring[j]);
				}

				if (repair != nullptr)
				{
					for (int j = 0; j < swarm_size; j++)
						solutioncheck(&offspring[j]);
				}

				for (int j = 0; j < swarm_size; j++)
					evaluator(&offspring[j]);

				z_update(offspring, swarm_size);

				for (int j = 0; j < swarm_size; j++)
				{
					//sub_f = sub_fitness(j + swarm_size);
					if (sub_fitness(j + swarm_size) < sub_fitness(j))
					{
						double sub_ff = sub_fitness(j);
						double sub_fo = sub_fitness(j + swarm_size);

						swarm[j].replace(&swarm[j + swarm_size]);
						if (MO::paretoFrontUpdate(&offspring[j], PF, gbestsize))
							update = true;
					}
				}

				if ((i / fes) % 20 == 0 && update)
				{
					if (detailedLog)
					{
						logFile << i / swarm_size << "\n";
						for (int j = 0; j < gbestsize; j++)
						{
							logFile << i / swarm_size << "\t" << gbest->ansprint() << "\n";
						}
						logFile << "\n";
					}
					else if (roughLog)
					{
						logFile << i / swarm_size << "\n";
						for (int j = 0; j < gbestsize; j++)
						{
							logFile << i / swarm_size << "\t" << gbest->roughprint();
						}
						logFile << "\n";
					}
					update = false;
				}
			}
			exe_time = time(NULL) - exe_time;
		}
	};

}