#pragma once
#include"optimizer.h"

namespace myEC {
	class Set_Particle: public Solution
	{
	public:
		Solution* pbest;
		double* velocity;
		int velocitySize;

		Set_Particle() {}

		Set_Particle(int size)
		{
			pbest = new Solution(size);
			velocity = new double[size];
		}

		~Set_Particle()
		{
			delete pbest;
			delete[]velocity;
		}

		void setSize(int _size, int velocitySize, double objectNumber = 1)
		{
			if (size != _size)
			{
				size = _size;
				if (!result)
					delete[] result;
				result = new double[size];
			}
			if (this->velocitySize != velocitySize)
			{
				this->velocitySize = velocitySize;
				if (!velocity)
					delete[] velocity;
				velocity = new double[velocitySize];
			}
			if (this->object_number != objectNumber)
			{
				this->object_number = objectNumber;
				if (!fitness)
					delete[] fitness;
				fitness = new double[object_number];
			}
			if (!pbest)
				delete pbest;
			pbest = new Solution(size, object_number);
		}
	};

	//Set-based Particle Swarm Optimization for discrete optimization
	class SPSO :public Optimizer
	{
	protected:
		bool is_direct;
		bool is_related;
		int pre_choice;
		int choicenumber;
		int v_size;

		bool v_heuristic;
		bool f_heuristic;

		double c1, c2;
		double w, r1, r2;
		Set_Particle* pswarm;
		sortHelper* sortbuffer;

		virtual void velocityUpdate(int id)
		{
			int to;
			double va;
			for (int i = 0; i < v_size; i++)
				pswarm[id].velocity[i] = w * pswarm[id].velocity[i];
			if (is_related)
			{
				int* example = new int[problem_size * 2];

				//learn from pbest
				for (int i = 1; i < problem_size; i++)
				{
					example[2 * i] = pswarm[id].pbest->result[i - 1];
					example[2 * i + 1] = pswarm[id].pbest->result[i];
				}
				example[0] = pswarm[id].pbest->result[problem_size - 1];
				example[1] = pswarm[id].pbest->result[0];

				//difference set build
				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] == pswarm[id].result[problem_size - 1] && example[2 * i + 1] == pswarm[id].result[0])
					{
						example[2 * i] = EMPTYVALUE;
						continue;
					}
					for (int j = 1; j < problem_size; j++)
					{
						if (example[2 * i] == pswarm[id].result[j - 1] && example[2 * i + 1] == pswarm[id].result[j])
						{
							example[2 * i] = EMPTYVALUE;
							break;
						}
					}					
				}			

				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] != EMPTYVALUE)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]])
							pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]] = va;
					}
				}
				if (!is_direct)
				{
					for (int i = 0; i < problem_size; i++)
					{
						if (example[2 * i] != EMPTYVALUE)
							pswarm[id].velocity[example[2 * i + 1] * choicenumber + example[2 * i]] = pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]];
					}
				}

				//learn from gbest
				for (int i = 1; i < problem_size; i++)
				{
					example[2 * i] = gbest->result[i - 1];
					example[2 * i + 1] = gbest->result[i];
				}
				example[0] = gbest->result[problem_size - 1];
				example[1] = gbest->result[0];

				//difference set build
				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] == pswarm[id].result[problem_size - 1] && example[2 * i + 1] == pswarm[id].result[0])
					{
						example[2 * i] = EMPTYVALUE;
						continue;
					}
					for (int j = 1; j < problem_size; j++)
					{
						if (example[2 * i] == pswarm[id].result[j - 1] && example[2 * i + 1] == pswarm[id].result[j])
						{
							example[2 * i] = EMPTYVALUE;
							break;
						}
					}
				}

				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] != EMPTYVALUE)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]])
							pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]] = va;
					}
				}
				if (!is_direct)
				{
					for (int i = 0; i < problem_size; i++)
					{
						if (example[2 * i] != EMPTYVALUE)
							pswarm[id].velocity[example[2 * i + 1] * choicenumber + example[2 * i]] = pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]];
					}
				}

				delete[] example;
			}
			else {
				for (int i = 0; i < problem_size; i++)//task
				{
					to = pswarm[id].pbest->result[i];
					if (to != EMPTYVALUE && pswarm[id].result[i] != to)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[i * choicenumber + to])
							pswarm[id].velocity[i * choicenumber + to] = va;
					}

					to = gbest->result[i];
					if (to != EMPTYVALUE && swarm[id].result[i] != to)
					{
						va = c2 * r2;
						if (va > pswarm[id].velocity[i * choicenumber + to])
							pswarm[id].velocity[i * choicenumber + to] = va;
					}
				}
			}
		}

		virtual void positionUpdate(int id)
		{
			Solution newSolution(problem_size, objectnum);

			for (int i = 0; i < problem_size; i++)
				newSolution.result[i] = EMPTYVALUE;

			if (model_ini != nullptr)
				model_ini();
			int counter;

			if (is_related)
			{
				pre_choice = rand() % choicenumber;

				int* buffer1 = new int[choicenumber];
				memset(buffer1, 128, choicenumber * sizeof(int));
				for (int i = 1; i < problem_size; i++)
					buffer1[(int)pswarm[id].result[i - 1]] = (int)pswarm[id].result[i];
				buffer1[(int)pswarm[id].result[problem_size - 1]] = (int)pswarm[id].result[0];

				for (int p = 0; p < problem_size; p++)
				{
					counter = 0;
					//velocity crisp
					for (int i = 0; i < choicenumber; i++)
					{
						if (rand01() < pswarm[id].velocity[pre_choice * choicenumber + i])
						{
							if (v_heuristic && heuristic_func != nullptr)
								sortbuffer[counter++] = sortHelper(i, heuristic_func(pre_choice, i));
							else
								sortbuffer[counter++] = sortHelper(i, rand01());
						}
					}
					std::sort(sortbuffer, sortbuffer + counter, std::greater<sortHelper>());

					for (int c = 0; c < counter; c++)
					{
						if (constrain_check(pre_choice, sortbuffer[c].id))
						{
							newSolution.result[p] = sortbuffer[c].id;
							model_change(pre_choice, newSolution.result[p]);
							pre_choice = newSolution.result[p];
							break;
						}
					}

					if (newSolution.result[p] != EMPTYVALUE)
						continue;

					//position crisp
					if (buffer1[pre_choice] != EMPTYVALUE)
					{
						if (constrain_check(pre_choice, buffer1[pre_choice]))
						{
							newSolution.result[p] = buffer1[pre_choice];
							model_change(pre_choice, newSolution.result[p]);
							pre_choice = newSolution.result[p];
							continue;
						}
					}

					//full crisp
					counter = 0;
					for (int i = 0; i < choicenumber; i++)
					{
						if (f_heuristic && heuristic_func != nullptr)
							sortbuffer[counter++] = sortHelper(i, heuristic_func(pre_choice, i));
						else
							sortbuffer[counter++] = sortHelper(i, rand01());
					}
					std::sort(sortbuffer, sortbuffer + counter, std::greater<sortHelper>());

					for (int j = 0; j < choicenumber; j++)
					{
						if (constrain_check(pre_choice, sortbuffer[j].id))
						{
							newSolution.result[p] = sortbuffer[j].id;
							model_change(pre_choice, newSolution.result[p]);
							pre_choice = newSolution.result[p];
							break;
						}
					}
				}

				delete[] buffer1;
			}
			else {
				for (int p = 0; p < problem_size; p++)
				{
					counter = 0;
					//velocity crisp
					for (int i = 0; i < choicenumber; i++)
					{
						if (rand01() < pswarm[id].velocity[p * choicenumber + i])
						{
							if (v_heuristic && heuristic_func != nullptr)
								sortbuffer[counter++] = sortHelper(i, heuristic_func(p, i));
							else
								sortbuffer[counter++] = sortHelper(i, rand01());
						}
					}
					std::sort(sortbuffer, sortbuffer + counter, std::greater<sortHelper>());

					for (int c = 0; c < counter; c++)
					{
						if (constrain_check(p, sortbuffer[c].id))
						{
							newSolution.result[p] = sortbuffer[c].id;
							model_change(p, newSolution.result[p]);
							break;
						}
					}

					if (newSolution.result[p] >= 0)
						continue;

					//position crisp
					if (0 <= pswarm[id].result[p] && pswarm[id].result[p] < choicenumber)
					{
						if (constrain_check(p, pswarm[id].result[p]))
						{
							newSolution.result[p] = pswarm[id].result[p];
							model_change(p, newSolution.result[p]);
							continue;
						}
					}

					//full crisp
					counter = 0;
					for (int i = 0; i < choicenumber; i++)
					{
						if (f_heuristic && heuristic_func != nullptr)
							sortbuffer[counter++] = sortHelper(i, heuristic_func(p, i));
						else
							sortbuffer[counter++] = sortHelper(i, rand01());
					}
					std::sort(sortbuffer, sortbuffer + counter, std::greater<sortHelper>());

					for (int j = 0; j < choicenumber; j++)
					{
						if (constrain_check(p, sortbuffer[j].id))
						{
							newSolution.result[p] = sortbuffer[j].id;
							model_change(p, newSolution.result[p]);
							break;
						}
					}
				}
			}

			evaluator(&newSolution);
			pswarm[id].replace(&newSolution);
		}

		virtual void ini()
		{
			double valueBuffer;

			if (is_related)
			{
				int to;
				for (int i = 0; i < swarm_size; i++)
				{
					//velocity initial
					memset(pswarm[i].velocity, 0, v_size * sizeof(double));
					for (int j = 0; j < choicenumber; j++)
					{
						to = rand() % choicenumber;
						pswarm[i].velocity[j * choicenumber + to] = 1;
						if (!is_direct)
						{
							pswarm[i].velocity[to * choicenumber + j] = pswarm[i].velocity[j * choicenumber + to];
						}
					}

					if (solution_ini_func != nullptr)
					{
						solution_ini_func(pswarm[i].result, problem_size);
					}
					else {
						if (model_ini != nullptr)
							model_ini();
						pre_choice = rand() % choicenumber;

						for (int j = 0; j < choicenumber; j++)
						{
							valueBuffer = rand() % choicenumber;
							for (int count = 0; count < choicenumber; count++)
							{
								if (!constrain_check(pre_choice, valueBuffer))
								{
									valueBuffer++;
									if (valueBuffer == choicenumber)
										valueBuffer = 0;
								}
								else break;
							}

							if (model_change != nullptr)
								model_change(pre_choice, valueBuffer);
							pre_choice = pswarm[i].result[j] = valueBuffer;
						}
					}

					evaluator(&pswarm[i]);
					pswarm[i].pbest->replace(&pswarm[i]);
				}
			}
			else {
				for (int i = 0; i < swarm_size; i++)
				{
					memset(pswarm[i].velocity, 0, v_size * sizeof(double));
					for (int j = 0; j < problem_size; j++)
					{
						pswarm[i].velocity[j * choicenumber + (rand() % choicenumber)] = 1;
					}

					if (solution_ini_func != nullptr)
					{
						solution_ini_func(pswarm[i].result, problem_size);	
					}
					else {
						if (model_ini != nullptr)
							model_ini();
						for (int j = 0; j < problem_size; j++)
						{														
							pswarm[i].result[j] = EMPTYVALUE;
							valueBuffer = rand() % choicenumber;

							for (int count = 0; count < choicenumber; count++)
							{
								if (!constrain_check(j, valueBuffer))
								{
									valueBuffer++;
									if (valueBuffer == choicenumber)
										valueBuffer = 0;
								}
								else
								{
									pswarm[i].result[j] = valueBuffer;
									break;
								}
							}
							if (model_change != nullptr)
								model_change(j, valueBuffer);
						}
					}
					evaluator(&pswarm[i]);
					pswarm[i].pbest->replace(&pswarm[i]);
				}
			}

			gbest->replace(pswarm[0].pbest);
			for (int i = 1; i < swarm_size; i++)
			{
				if (*pswarm[i].pbest < *gbest)
					gbest->replace(pswarm[i].pbest);
			}
		}

		virtual void gbestUpdate(int generation)
		{
			for (int j = 0; j < swarm_size; j++)
			{
				if (pswarm[j] < *pswarm[j].pbest)
				{
					pswarm[j].pbest->replace(&pswarm[j]);
					if (*pswarm[j].pbest < *gbest)
					{
						gbest->replace(pswarm[j].pbest);
						if (detailedLog)
							logFile << generation / swarm_size << "\t" << gbest->ansprint() << "\n";
						else if (roughLog)
							logFile << generation / swarm_size << "\t" << gbest->roughprint();
					}
				}
			}
		}

	public:
		SPSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double (*heuristic_func)(int demensionId, double value),
			double c1, double c2, bool v_heuristic, bool f_heuristic,
			int choicenumber, bool is_direct, bool is_related,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->c1 = c1;
			this->c2 = c2;
			swarm = new Set_Particle[swarm_size];
			pswarm = (Set_Particle*)swarm;
			gbest = new Solution(problem_size, objectNumber);

			this->v_heuristic = v_heuristic;
			this->f_heuristic = f_heuristic;
			this->heuristic_func = heuristic_func;

			this->choicenumber = choicenumber;
			this->is_direct = is_direct;
			this->is_related = is_related;
			sortbuffer = new sortHelper[choicenumber];

			if (is_related)
				v_size = choicenumber * choicenumber;
			else v_size = problem_size * choicenumber;

			for (int i = 0; i < swarm_size; i++)
			{
				pswarm[i].setSize(ps, v_size, objectNumber);
				pswarm[i].setCompareFunc(compare_func);
				pswarm[i].pbest->setCompareFunc(compare_func);
			}
		}

		~SPSO()
		{
			pswarm = nullptr;
			delete[] swarm;
			delete[] sortbuffer;
			delete gbest;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			for (int i = swarm_size; i < generation; i += swarm_size)
			{
				r1 = rand01();
				r2 = rand01();
				w = 0.9 - (double)i / generation * 0.5;

				for (int j = 0; j < swarm_size; j++)
				{
					velocityUpdate(j);
					positionUpdate(j);
				}
				gbestUpdate(i);
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//Set-based Comprehensive Learning Particle Swarm Optimizer
	class S_CLPSO :public SPSO
	{
	protected:
		double* pc;
		int* fid;
		bool* fi;

		int m;
		int* age;

		void fiUpdate(int id)
		{
			int id1, id2;
			bool allid = true;

			do {
				id1 = rand() % swarm_size;
			} while (id1 == id);

			do {
				id2 = rand() % swarm_size;
			} while (id2 == id || id2 == id1);

			if (*pswarm[id2].pbest < *pswarm[id1].pbest)
				id1 = id2;
			fid[id] = id1;

			for (int i = 0; i < problem_size; i++)
			{
				fi[id * problem_size + i] = rand01() < pc[id];
			}
			for (int i = 0; i < problem_size; i++)
			{
				if (fi[id * problem_size + i])
				{
					allid = false;
					break;
				}
			}

			if (allid)
				fi[id * problem_size + rand() % problem_size] = true;

			age[id] = 0;
		}

		virtual void ini()
		{
			SPSO::ini();

			for (int i = 0; i < swarm_size; i++)
				fiUpdate(i);
		}

		virtual void velocityUpdate(int id)
		{
			int to;
			double va;
			for (int i = 0; i < v_size; i++)
				pswarm[id].velocity[i] = w * pswarm[id].velocity[i];

			if (is_related)
			{
				int* example = new int[problem_size * 2];

				//build learn example
				for (int i = 1; i < problem_size; i++)
				{
					if (fi[id * problem_size + i])
					{
						example[2 * i] = pswarm[fid[id]].pbest->result[i - 1];
						example[2 * i + 1] = pswarm[fid[id]].pbest->result[i];
					}
					else {
						example[2 * i] = pswarm[id].pbest->result[i - 1];
						example[2 * i + 1] = pswarm[id].pbest->result[i];
					}
					
				}
				if (fi[id * problem_size])
				{
					example[0] = pswarm[fid[id]].pbest->result[problem_size - 1];
					example[1] = pswarm[fid[id]].pbest->result[0];
				}
				else {
					example[0] = pswarm[id].pbest->result[problem_size - 1];
					example[1] = pswarm[id].pbest->result[0];
				}
				
				//difference set build
				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] == pswarm[id].result[problem_size - 1] && example[2 * i + 1] == pswarm[id].result[0])
					{
						example[2 * i] = EMPTYVALUE;
						continue;
					}
					for (int j = 1; j < problem_size; j++)
					{
						if (example[2 * i] == pswarm[id].result[j - 1] && example[2 * i + 1] == pswarm[id].result[j])
						{
							example[2 * i] = EMPTYVALUE;
							break;
						}
					}
				}

				for (int i = 0; i < problem_size; i++)
				{
					if (example[2 * i] != EMPTYVALUE)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]])
							pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]] = va;
					}
				}
				if (!is_direct)
				{
					for (int i = 0; i < problem_size; i++)
					{
						if (example[2 * i] != EMPTYVALUE)
							pswarm[id].velocity[example[2 * i + 1] * choicenumber + example[2 * i]] = pswarm[id].velocity[example[2 * i] * choicenumber + example[2 * i + 1]];
					}
				}

				delete[] example;
			}
			else {
				for (int i = 0; i < problem_size; i++)//task
				{
					if (fi[id * problem_size + i])
						to = pswarm[fid[id]].pbest->result[i];
					else 
						to = pswarm[id].pbest->result[i];
					
					if (to != EMPTYVALUE && pswarm[id].result[i] != to)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[i * choicenumber + to])
							pswarm[id].velocity[i * choicenumber + to] = va;
					}
				}
			}
		}

		virtual void gbestUpdate(int generation)
		{
			for (int j = 0; j < swarm_size; j++)
			{
				if (pswarm[j] < *pswarm[j].pbest)
				{
					pswarm[j].pbest->replace(&pswarm[j]);
					if (*pswarm[j].pbest < *gbest)
					{
						gbest->replace(pswarm[j].pbest);
						if (detailedLog)
							logFile << generation / swarm_size << "\t" << gbest->ansprint() << "\n";
						else if (roughLog)
							logFile << generation / swarm_size << "\t" << gbest->roughprint();
					}
				}
				else {
					age[j]++;
					if (age[j] >= m)
						fiUpdate(j);
				}
			}
		}

	public:
		S_CLPSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double (*heuristic_func)(int demensionId, double value),
			double c, int m,
			bool v_heuristic, bool f_heuristic,
			int choicenumber, bool is_direct, bool is_related,
			int logType)
			:SPSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func,
				heuristic_func, c, 0, v_heuristic, f_heuristic, choicenumber, is_direct, is_related, logType)
		{
			pc = new double[swarm_size];
			age = new int[swarm_size];
			fid = new int[swarm_size];
			fi = new bool[swarm_size * problem_size];

			this->m = m;

			for (int i = 0; i < swarm_size; i++)
				pc[i] = 0.05 + 0.45 * (pow(E_CONST, (10 * i / (swarm_size - 1))) - 1) / (pow(E_CONST, 10) - 1);
		}

		~S_CLPSO()
		{
			delete[] pc;
			delete[] age;
			delete[] fid;
			delete[] fi;
		}
	};
}

