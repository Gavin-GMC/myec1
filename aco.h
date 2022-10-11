#pragma once
#include"optimizer.h"

namespace myEC {
	//classic ACO for sigle object optimization
	class AS : public Optimizer
	{
	protected:
		double* pheromone;
		double alpha = 1;
		double belta = 2;
		double rho = 0.5;
		double t0 = -1;
		int choicenumber;
		bool is_direct;
		bool is_related;
		sortHelper* sortbuffer;
		int pre_choice;

		virtual void global_update(Solution* best)
		{
			for (int j = 0; j < problem_size * choicenumber; j++)
				pheromone[j] *= (1 - rho);

			for (int j = 0; j < swarm_size; j++)
			{
				if (is_related)
				{
					for (int i = 1; i < best[j].size; i++)
					{
						pheromone[(int)best[j].result[i - 1] * choicenumber + (int)best[j].result[i]] += (1 / best[j].fitness[0]);
					}
					pheromone[(int)best[j].result[best[j].size - 1] * choicenumber + (int)best[j].result[0]] += (1 / best[j].fitness[0]);

					if (!is_direct)
					{
						for (int i = 1; i < best[j].size; i++)
						{
							pheromone[(int)best[j].result[i] * choicenumber + (int)best[j].result[i - 1]] =
								pheromone[(int)best[j].result[i - 1] * choicenumber + (int)best[j].result[i]];
						}
						pheromone[(int)best[j].result[0] * choicenumber + (int)best[j].result[best[j].size - 1]] =
							pheromone[(int)best[j].result[best[j].size - 1] * choicenumber + (int)best[j].result[0]];
					}
				}
				else {
					for (int i = 0; i < best[j].size; i++)
						pheromone[i * choicenumber + (int)best[j].result[i]] += (1 / best[j].fitness[0]);
				}
			}			
		}

		virtual double priority(int dimension, int choice)
		{
			double h;
			if (heuristic_func != nullptr)
				h = heuristic_func(dimension, choice);
			else h = rand01();
			return pow(pheromone[dimension * choicenumber + choice], alpha) * pow(h, belta);
		}

		virtual int find_next(int dimension)
		{
			int next = 0;
			double total = 0;
			for (int j = 0; j < choicenumber; j++)
			{
				if (constrain_check(dimension, j))
					sortbuffer[j] = sortHelper(j, priority(dimension, j));
				else
					sortbuffer[j] = sortHelper(j, 0);
					
				total += sortbuffer[j].value;
			}
			if (total == 0)
			{
				std::cerr << "No feasible aolution!\n";
				return EMPTYVALUE;
			}

			for (int j = 0; j < choicenumber; j++)
				sortbuffer[j].value /= total;
			total = rand01();

			while (true)
			{
				total -= sortbuffer[next].value;
				if (total < 0)
					break;
				next++;
				if (next == choicenumber-1)
					break;
			}

			while (sortbuffer[next].value == 0)
				next--;

			return sortbuffer[next].id;
		}

		virtual void ini()
		{
			Solution greedy(problem_size, objectnum);
			int state;

			if (model_ini != nullptr)
				model_ini();
			pre_choice = rand() % problem_size;
			for (int i = 0; i < problem_size; i++)
			{
				for (int j = 0; j < choicenumber; j++)
				{
					if (is_related)
						state = pre_choice;
					else state = i;

					if (constrain_check(state, j))
					{
						if (heuristic_func != nullptr)
							sortbuffer[j] = sortHelper(j, heuristic_func(state, j));
						else
							sortbuffer[j] = sortHelper(j, rand01());
					}
					else
						sortbuffer[j] = sortHelper(j, -1 * MAX);
				}
				sort(sortbuffer, sortbuffer + choicenumber, std::greater<sortHelper>());
				pre_choice = greedy.result[i] = sortbuffer[0].id;
				if (model_change != nullptr)
					model_change(i, greedy.result[i]);
			}
			evaluator(&greedy);

			if (t0 < 0)
				t0 = 1 / (swarm_size * greedy.fitness[0]);

			if (is_related)
			{
				for (int i = 0; i < choicenumber * choicenumber; i++)
					pheromone[i] = t0;
			}
			else {
				for (int i = 0; i < problem_size * choicenumber; i++)
					pheromone[i] = t0;
			}

			gbest->replace(&greedy);
		}

	public:

		AS(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), double objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double (*heuristic_func)(int demensionId, double value),
			double alpha, double belta, double rho,
			double t0, int choicenumber, bool is_direct, bool is_related,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			swarm = new Solution[swarm_size];
			gbest = new Solution(problem_size, objectNumber);
			sortbuffer = new sortHelper[choicenumber];

			for (int i = 0; i < swarm_size; i++)
			{
				swarm[i].setSize(ps, objectNumber);
				swarm[i].setCompareFunc(compare_func);
			}

			this->heuristic_func = heuristic_func;
			this->choicenumber = choicenumber;
			this->is_direct = is_direct;
			this->is_related = is_related;
			if (is_related)
				pheromone = new double[choicenumber * choicenumber];
			else
				pheromone = new double[ps * choicenumber];

			this->t0 = t0;
			this->alpha = alpha;
			this->belta = belta;
			this->rho = rho;
		}

		~AS()
		{
			delete[] swarm;
			delete[] gbest;
			delete[] sortbuffer;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = 0; i < generation; i += swarm_size)
			{
				for (int j = 0; j < swarm_size; j++)
				{
					//build an ant
					if (model_ini != nullptr)
						model_ini();
					pre_choice = rand() % choicenumber;
					for (int id = 0; id < problem_size; id++)
					{
						if (is_related)
							pre_choice = swarm[j].result[id] = find_next(pre_choice);
						else
							swarm[j].result[id] = find_next(id);

						if (model_change != nullptr)
							model_change(id, swarm[j].result[id]);
					}
					evaluator(&swarm[j]);
				}

				global_update(swarm);

				for (int j = 0; j < swarm_size; j++)
				{
					if (swarm[j] < *gbest)
					{
						gbest->replace(&swarm[j]);
						if (detailedLog)
							logFile << i / swarm_size << "\t" << gbest->ansprint();
						else if (roughLog)
							logFile << i / swarm_size << "\t" << gbest->roughprint();
					}
				}
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//Ant Conoly System
	//User should build a external counter for the prallel solution building
	class ACS : public AS
	{
	protected:
		Solution* lbest;
		int* pre_choice;
		double p0;

		virtual void global_update(Solution* best)
		{
			if (is_related)
			{
				for (int i = 1; i < best->size; i++)
				{
					pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]] =
						(1 - rho) * pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]] + rho * (1 / best->fitness[0]);
				}
				pheromone[(int)best->result[best->size - 1] * choicenumber + (int)best->result[0]] += (1 / best->fitness[0]);

				if (!is_direct)
				{
					for (int i = 1; i < best->size; i++)
					{
						pheromone[(int)best->result[i] * choicenumber + (int)best->result[i - 1]] =
							pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]];
					}
					pheromone[(int)best->result[0] * choicenumber + (int)best->result[best->size - 1]] =
						pheromone[(int)best->result[best->size - 1] * choicenumber + (int)best->result[0]];
				}
			}
			else {
				for (int i = 0; i < best->size; i++)
					pheromone[i * choicenumber + (int)best->result[i]] = 
					(1 - rho) * pheromone[i * choicenumber + (int)best->result[i]] + rho * (1 / best->fitness[0]);
			}
		}

		virtual void local_update(int i, int j)
		{
			if (!is_direct)
			{
				pheromone[i * choicenumber + j] = pheromone[j * choicenumber + i] = (1 - rho) * pheromone[i * choicenumber + j] + rho * t0;
			}
			else pheromone[i * choicenumber + j] = (1 - rho) * pheromone[i * choicenumber + j] + rho * t0;
		}

		virtual int find_next(int dimension)
		{
			bool explore = (rand01() > p0);

			if (explore)
			{
				int next = 0;
				double total = 0;

				for (int j = 0; j < choicenumber; j++)
				{
					if (constrain_check(dimension, j))
						sortbuffer[j] = sortHelper(j, priority(dimension, j));
					else
						sortbuffer[j] = sortHelper(j, 0);

					total += sortbuffer[j].value;
				}
				if (total == 0)
				{
					std::cerr << "No feasible aolution!\n";
					return EMPTYVALUE;
				}

				for (int j = 0; j < choicenumber; j++)
					sortbuffer[j].value /= total;
				total = rand01();
				while (true)
				{
					total -= sortbuffer[next].value;
					if (total < 0)
						break;
					next++;
					if (next == choicenumber - 1)
						break;
				}
				while (sortbuffer[next].value == 0)
					next--;

				return sortbuffer[next].id;				
			}
			else {
				double max = -1 * MAX;
				int back = 0;
				double pribuffer;
				for (int j = 0; j < choicenumber; j++)
				{
					if (!constrain_check(dimension, j))
						continue;
					pribuffer = priority(dimension, j);
					if (pribuffer > max)
					{
						max = pribuffer;
						back = j;
					}
				}
				if (max == -1 * MAX)
				{
					std::cerr << "No feasible aolution!\n";
					return EMPTYVALUE;
				}
				return back;
			}		
		}

		virtual void ini()
		{
			Solution greedy(problem_size, objectnum);
			int state;

			if (model_ini != nullptr)
				model_ini();
			pre_choice[0] = rand() % problem_size;
			
			for (int i = 0; i < problem_size; i++)
			{
				for (int counter = 0; counter < swarm_size; counter++)
				{
					if (counter == 0)
					{
						for (int j = 0; j < choicenumber; j++)
						{
							if (is_related)
								state = pre_choice[0];
							else state = i;

							if (constrain_check(state, j))
							{
								if (heuristic_func != nullptr)
									sortbuffer[j] = sortHelper(j, heuristic_func(state, j));
								else
									sortbuffer[j] = sortHelper(j, rand01());
							}
							else
								sortbuffer[j] = sortHelper(j, -1 * MAX);
						}
						sort(sortbuffer, sortbuffer + choicenumber, std::greater<sortHelper>());
						pre_choice[0] = greedy.result[i] = sortbuffer[0].id;
						if (model_change != nullptr)
							model_change(state, greedy.result[i]);
					}
					else
						model_change(counter, i);
				}
			}
			evaluator(&greedy);

			if (t0 < 0)
				t0 = 1 / (problem_size * greedy.fitness[0]);				

			if (is_related)
			{
				for (int i = 0; i < choicenumber * choicenumber; i++)
					pheromone[i] = t0;
			}
			else {
				for (int i = 0; i < problem_size * choicenumber; i++)
					pheromone[i] = t0;
			}

			gbest->replace(&greedy);
		}
	public:
		ACS(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), double objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double (*heuristic_func)(int demensionId, double value),
			double alpha, double belta, double rho, double p0,
			double t0, int choicenumber, bool is_direct, bool is_related,
			int logType)
			:AS(ps, ss, gen, evaluate_func, objectNumber, compare_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
				alpha, belta, rho, t0, choicenumber, is_direct, is_related, logType)
		{
			pre_choice = new int[ss];
			this->p0 = p0;
		}

		~ACS()
		{
			delete[] pre_choice;
			lbest = nullptr;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = 0; i < generation; i += swarm_size)
			{
				//build swarm 
				if (model_ini != nullptr)
					model_ini();
				for (int j = 0; j < swarm_size; j++)
					pre_choice[j] = rand() % choicenumber;
				for (int id = 0; id < problem_size; id++)
				{
					for (int j = 0; j < swarm_size; j++)
					{
						if (is_related)
						{
							swarm[j].result[id] = find_next(pre_choice[j]);
							if (id != 0)
								local_update(pre_choice[j], swarm[j].result[id]);
							pre_choice[j] = swarm[j].result[id];
						}
						else
						{
							swarm[j].result[id] = find_next(id);
							local_update(id, swarm[j].result[id]);
						}
						if (model_change != nullptr)
							model_change(id, swarm[j].result[id]);
					}
				}
				for (int j = 0; j < swarm_size; j++)
					local_update(swarm[j].result[problem_size - 1], swarm[j].result[0]);

				for (int j = 0; j < swarm_size; j++)
					evaluator(&swarm[j]);

				//find local best in this iteration
				lbest = swarm;
				for (int j = 1; j < swarm_size; j++)
				{
					if (swarm[j] < *lbest)
						lbest = swarm + j;
				}
				if (*lbest < *gbest)
				{
					gbest->replace(lbest);
					if (detailedLog)
						logFile << i / swarm_size << "\t" << gbest->ansprint();
					else if (roughLog)
						logFile << i / swarm_size << "\t" << gbest->roughprint();
				}

				global_update(lbest);
				global_update(gbest);
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//Max-Min Ant System
	/*
	class MMAS : public AS
	{
	protected:
		int age;
		double t_max;
		double t_min;

		virtual void global_update(Solution* best)
		{
			if (is_related)
			{
				for (int i = 1; i < best->size; i++)
				{
					pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]] =
						rho * pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]] + (1 / best->fitness[0]);
				}
				pheromone[(int)best->result[best->size - 1] * choicenumber + (int)best->result[0]] += (1 / best->fitness[0]);

				if (!is_direct)
				{
					for (int i = 1; i < best->size; i++)
					{
						pheromone[(int)best->result[i] * choicenumber + (int)best->result[i - 1]] =
							pheromone[(int)best->result[i - 1] * choicenumber + (int)best->result[i]];
					}
					pheromone[(int)best->result[0] * choicenumber + (int)best->result[best->size - 1]] =
						pheromone[(int)best->result[best->size - 1] * choicenumber + (int)best->result[0]];
				}
			}
			else {
				for (int i = 0; i < best->size; i++)
					pheromone[i * choicenumber + (int)best->result[i]] =
					rho * pheromone[i * choicenumber + (int)best->result[i]] + (1 / best->fitness[0]);
			}
		}

	public:
		MMAS(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), double objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double (*heuristic_func)(int demensionId, double value),
			double alpha, double belta, double rho,
			double t0, int choicenumber, bool is_direct, bool is_related,
			int logType)
			:AS(ps, ss, gen, evaluate_func, objectNumber, compare_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
				alpha, belta, rho, t0, choicenumber, is_direct, is_related, logType)
		{}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = 0; i < generation; i += swarm_size)
			{
				for (int j = 0; j < swarm_size; j++)
				{
					//build an ant
					if (model_ini != nullptr)
						model_ini();
					pre_choice = rand() % problem_size;
					for (int id = 0; id < problem_size; id++)
					{
						if (is_related)
							pre_choice = swarm[j].result[id] = find_next(pre_choice);
						else
							swarm[j].result[id] = find_next(id);

						if (model_change != nullptr)
							model_change(id, swarm[j].result[id]);
					}
					evaluator(&swarm[j]);
				}

				global_update(swarm);

				for (int j = 0; j < swarm_size; j++)
				{
					if (swarm[j] < *gbest)
					{
						gbest->replace(&swarm[j]);
						if (detailedLog)
							logFile << i / swarm_size << "\t" << gbest->ansprint();
						else if (roughLog)
							logFile << i / swarm_size << "\t" << gbest->roughprint();
					}
				}
			}
			exe_time = time(NULL) - exe_time;
		}

	};*/
}