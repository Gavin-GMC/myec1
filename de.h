#pragma once
#include"optimizer.h"

namespace myEC {

	//Differential Evolution Algorithm for continues optimization
	class DE :public Optimizer
	{
	protected:
		Solution* offspring;
		double f;
		double cr;
		double fitness_max;
		double fitness_min;
		int bestid;
		int worstid;

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

			bestid = worstid = 0;
			for (int i = 1; i < swarm_size; i++)
			{
				if (swarm[i] < swarm[bestid])
					bestid = i;
				else if (swarm[worstid] < swarm[i])
					worstid = i;
			}
			fitness_min = swarm[bestid].fitness[0];
			fitness_max = swarm[worstid].fitness[0];

			gbest->replace(&swarm[bestid]);
		}

		virtual void crossover(int id)
		{
			if (cr == EMPTYVALUE)
			{
				double cri;
				cri = 0.1 + 0.5 * (swarm[id].fitness[0] - fitness_min) / (fitness_max - fitness_min);

				for (int i = 0; i < problem_size; i++)
					if (rand01() > cri)
						offspring[id].result[i] = swarm[id].result[i];
			}
			else {
				for (int i = 0; i < problem_size; i++)
					if (rand01() > cr)
						offspring[id].result[i] = swarm[id].result[i];
			}
		}

		virtual void mutation(int id)
		{
			int p1, p2, p3;

			//make sure donnot choose the same individual
			p1 = rand() % swarm_size;
			p2 = rand() % swarm_size;
			p3 = rand() % swarm_size;
			while (p3 == p1 || p3 == p2)
				p3 = rand() % swarm_size;
			while (p2 == p1 || p2 == p3)
				p2 = rand() % swarm_size;

			if (f == EMPTYVALUE)
			{
				if (swarm[p3] < swarm[p2])
					std::swap(p2, p3);
				if (swarm[p2] < swarm[p1])
				{
					std::swap(p1, p2);
					if (swarm[p3] < swarm[p2])
						std::swap(p2, p3);
				}

				//Fi = Fl+(Fu-Fl)*(fm-fb)/(fw-fb)
				for (int i = 0; i < problem_size; i++)
					offspring[id].result[i] = swarm[p1].result[i] +
					(0.1 + 0.8 * (swarm[p2].fitness[0] - swarm[p1].fitness[0]) / (swarm[p3].fitness[0] - swarm[p1].fitness[0]))
					* (swarm[p2].result[i] - swarm[p3].result[i]);
			}
			else {
				for (int i = 0; i < problem_size; i++)
					offspring[id].result[i] = swarm[p1].result[i] +
					f * (swarm[p2].result[i] - swarm[p3].result[i]);
			}
		}

	public:
		DE(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double f, double cr,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->f = f;
			this->cr = cr;
			swarm = new Solution[2 * swarm_size];
			gbest = new Solution(problem_size, objectNumber);

			offspring = swarm + swarm_size;

			for (int i = 0; i < 2 * swarm_size; i++)
			{
				swarm[i].setSize(ps, objectNumber);
				swarm[i].setCompareFunc(compare_func);;
			}
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int fes = swarm_size;
			for (int i = swarm_size; i < generation; i += fes)
			{
				for (int j = 0; j < swarm_size; j++)
				{
					mutation(j);
					crossover(j);
				}

				for (int j = 0; j < swarm_size; j++)
				{
					evaluator(&offspring[j]);
					if (offspring[j] < swarm[j])
						swarm[j].replace(&offspring[j]);
				}


				bestid = worstid = 0;
				for (int i = 1; i < swarm_size; i++)
				{
					if (swarm[i] < swarm[bestid])
						bestid = i;
					else if (swarm[worstid] < swarm[i])
						worstid = i;
				}
				fitness_min = swarm[bestid].fitness[0];
				fitness_max = swarm[worstid].fitness[0];

				if (swarm[bestid] < *gbest)
				{
					gbest->replace(&swarm[bestid]);
					if (detailedLog)
						logFile << i / swarm_size << "\t" << gbest->ansprint() << "\n";
					else if (roughLog)
						logFile << i / swarm_size << "\t" << gbest->roughprint();
				}
			}
			exe_time = time(NULL) - exe_time;
		}

	};
}