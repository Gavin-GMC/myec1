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

	void checking(Solution* s)
	{
		bool v[100];
		for (int i = 0; i < 100; i++)
			v[i] = false;

		for (int i = 0; i < 100; i++)
			if (!v[(int)s->result[i]])
				v[(int)s->result[i]] = true;
			else
				std::cerr << "mistake!!";
	}

	class SPSO :public Optimizer
	{
	protected:
		bool is_direct;
		bool is_related;
		int pre_choice;
		int choicenumber;
		int v_size;

		bool v_heuristic;
		bool f_heuristic = true;

		double c1, c2;
		double w, r1, r2;
		Set_Particle* pswarm;
		sortHelper* sortbuffer;

		virtual void velocityUpdate(int id)
		{
			int to;
			int va;
			for (int i = 0; i < v_size; i++)
				pswarm[id].velocity[i] = w * pswarm[id].velocity[i];

			if (is_related)
			{
				int* buffer1 = new int[choicenumber];
				int* buffer2 = new int[choicenumber];

				memset(buffer1, 128, choicenumber * sizeof(int));
				for (int i = 1; i < problem_size; i++)
					buffer1[(int)pswarm[id].result[i - 1]] = (int)pswarm[id].result[i];
				buffer1[(int)pswarm[id].result[problem_size - 1]] = (int)pswarm[id].result[0];

				memset(buffer2, 128, choicenumber * sizeof(int));
				for (int i = 1; i < problem_size; i++)
					buffer2[(int)pswarm[id].pbest->result[i - 1]] = (int)pswarm[id].pbest->result[i];
				buffer2[(int)pswarm[id].pbest->result[problem_size - 1]] = (int)pswarm[id].pbest->result[0];

				for (int i = 0; i < choicenumber; i++)
				{
					if (buffer2[i] != buffer1[i] && buffer2[i] != EMPTYVALUE)
					{
						va = c1 * r1;
						if (va > pswarm[id].velocity[i * choicenumber + buffer2[i]])
							pswarm[id].velocity[i * choicenumber + buffer2[i]] = va;
					}
				}
				if (!is_direct)
				{
					for (int i = 0; i < choicenumber; i++)
					{
						if (buffer2[i] != EMPTYVALUE)
							pswarm[id].velocity[buffer2[i] * choicenumber + i] = pswarm[id].velocity[i * choicenumber + buffer2[i]];
					}
				}

				memset(buffer2, 128, choicenumber * sizeof(int));
				for (int i = 1; i < problem_size; i++)
					buffer2[(int)gbest->result[i - 1]] = (int)gbest->result[i];
				buffer2[(int)gbest->result[problem_size - 1]] = (int)gbest->result[0];

				for (int i = 0; i < choicenumber; i++)
				{
					if (buffer2[i] != buffer1[i] && buffer2[i] != EMPTYVALUE)
					{
						va = c2 * r2;
						if (va > pswarm[id].velocity[i * choicenumber + buffer2[i]])
							pswarm[id].velocity[i * choicenumber + buffer2[i]] = va;
					}
				}
				if (!is_direct)
				{
					for (int i = 0; i < choicenumber; i++)
					{
						if (buffer2[i] != EMPTYVALUE)
							pswarm[id].velocity[buffer2[i] * choicenumber + i] = pswarm[id].velocity[i * choicenumber + buffer2[i]];
					}
				}

				delete[] buffer1;
				delete[] buffer2;
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
					std::sort(sortbuffer, sortbuffer + counter, std::less<sortHelper>());

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
					if (0 <= pswarm[id].result[p] && buffer1[pre_choice] < choicenumber)
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
					std::sort(sortbuffer, sortbuffer + counter, std::less<sortHelper>());

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
					std::sort(sortbuffer, sortbuffer + counter, std::less<sortHelper>());

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
					std::sort(sortbuffer, sortbuffer + counter, std::less<sortHelper>());

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
			checking(&newSolution);
		}

		virtual void ini()
		{
			double valueBuffer;

			if (is_related)
			{
				int to;
				for (int i = 0; i < swarm_size; i++)
				{
					memset(pswarm[i].velocity, 0, v_size * sizeof(double));
					if (solution_ini_func != nullptr)
					{
						solution_ini_func(pswarm[i].result, problem_size);

						for (int j = 0; j < choicenumber; j++)
						{
							to = rand() % choicenumber;
							pswarm[i].velocity[j * choicenumber + to] = rand01();
							if (!is_direct)
							{
								pswarm[i].velocity[j * choicenumber + to] = pswarm[i].velocity[to * choicenumber + j];
							}
						}
					}
					else {
						if (model_ini != nullptr)
							model_ini();
						pre_choice = rand() % choicenumber;

						for (int j = 0; j < choicenumber; j++)
						{
							to = rand() % choicenumber;
							pswarm[i].velocity[j * choicenumber + to] = rand01();
							if (!is_direct)
								pswarm[i].velocity[to * choicenumber + j] = pswarm[i].velocity[j * choicenumber + to];

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
					if (solution_ini_func != nullptr)
					{
						solution_ini_func(pswarm[i].result, problem_size);

						for (int j = 0; j < problem_size; j++)
						{
							pswarm[i].velocity[j * choicenumber + (rand() % choicenumber)] = rand01();
						}
					}
					else {
						if (model_ini != nullptr)
							model_ini();
						for (int j = 0; j < problem_size; j++)
						{
							valueBuffer = rand() % choicenumber;
							pswarm[i].velocity[j * choicenumber + (rand() % choicenumber)] = rand01();

							if (repair != nullptr && !constrain_check(j, valueBuffer))
								valueBuffer = repair(j, valueBuffer);
							if (model_change != nullptr)
								model_change(j, valueBuffer);
							pswarm[i].result[j] = valueBuffer;
						}
					}

					evaluator(&pswarm[i]);
					pswarm[i].pbest->replace(&pswarm[i]);
				}
			}

			gbest = pswarm[0].pbest;
			for (int i = 1; i < swarm_size; i++)
			{
				if (*pswarm[i].pbest < *gbest)
					gbest = pswarm[i].pbest;
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

			this->v_heuristic = v_heuristic;
			this->f_heuristic = f_heuristic;

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
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = 0; i < swarm_size; i++)
				checking(&pswarm[i]);


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
				for (int j = 0; j < swarm_size; j++)
				{
					if (pswarm[j] < *pswarm[j].pbest)
					{
						pswarm[j].pbest->replace(&pswarm[j]);
						if (*pswarm[j].pbest < *gbest)
						{
							gbest = pswarm[j].pbest;
							if (detailedLog)
								logFile << i / swarm_size << "\t" << gbest->ansprint() << "\n";
							else if (roughLog)
								logFile << i / swarm_size << "\t" << gbest->roughprint();
						}
					}
				}
			}
			exe_time = time(NULL) - exe_time;
		}
	};
}

