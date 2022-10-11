#pragma once
#include"optimizer.h"
#include"multiobject.h"
#include"eaoperator.h"

namespace myEC {
	//Classic GA for single Object and discrete optimization
	class GA : public Optimizer
	{
	protected:
		double p_m, p_c;
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
			double pc, double pm, EAOperator::crossover crossover_type, EAOperator::mutation mutation_type, int mutation_times,
			EAOperator::selection selection_type, bool elitist_strategy,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->p_m = pm;
			this->p_c = pc;
			swarm = new Solution[swarm_size * 2];
			gbest = new Solution(problem_size, objectNumber);

			offspring = swarm + swarm_size;

			for (int i = 0; i < 2 * swarm_size; i++)
			{
				swarm[i].setSize(ps, objectNumber);
				swarm[i].setCompareFunc(compare_func);;
			}

			crossover_func = Crossover(crossover_type);
			mutation_func = Mutation(mutation_type, mutation_times, constrainRangeList);
			selection_func = Selection(selection_type);

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
					if (rand01() < p_c)
						crossover_func(&swarm[j], &swarm[j + 1]);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					if (rand01() < p_m)
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
			double pc, double pm, EAOperator::crossover crossover_type, EAOperator::mutation mutation_type, int mutation_times,
			EAOperator::selection selection_type,
			int logType)
			:GA(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, pc, pm, crossover_type, mutation_type, mutation_times,
				selection_type, false, logType)
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
				selection_func(swarm, swarm_size, offspring);

				for (int j = 0; j < swarm_size - 1; j += 2)
				{
					if (myEC::rand01() < p_c)
						crossover_func(&offspring[j], &offspring[j + 1]);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					if (myEC::rand01() < p_m)
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
				MO::buildPartialOrder<Solution>(swarm, swarm_size * 2, indicator, 1);

				int beforeid = 0;
				int rank_now = 1;
				for (int j = 0; j < swarm_size * 2; j++)
				{
					if (indicator[j] != rank_now)
					{
						MO::CrowdDistance<Solution>(swarm + beforeid, j - beforeid, indicator + swarm_size + beforeid);
						beforeid = j;
						rank_now = indicator[j];
					}
				}

				for (int j = 0; j < swarm_size * 2; j++)
				{
					if (indicator[j] != 1)
						break;
					if (MO::paretoFrontUpdate(&offspring[j], PF, gbestsize))
						update = true;
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