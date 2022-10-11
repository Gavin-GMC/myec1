#pragma once
#include"optimizer.h"
#include"multiobject.h"

namespace myEC {
	class Particle :public Solution
	{
	public:
		Solution* pbest;
		double* velocity;

		Particle() {}

		Particle(int size)
		{
			pbest = new Solution(size);
			velocity = new double[size];
		}

		~Particle()
		{
			delete pbest;
			delete[]velocity;
		}

		void setSize(int _size, int object_number = 1)
		{
			if (size != _size)
			{
				size = _size;
				if (!result)
					delete[] result;
				result = new double[size];
				if (!velocity)
					delete[] velocity;
				velocity = new double[size];

			}
			if (this->object_number != object_number)
			{
				this->object_number = object_number;
				if (!fitness)
					delete[] fitness;
				fitness = new double[object_number];
			}
			if (!pbest)
				delete pbest;
			pbest = new Solution(size, object_number);
		}
	};

	//Classic PSO for single Object
	class PSO : public Optimizer
	{
	protected:
		double c1, c2;
		double w, r1, r2;
		Particle* pswarm;

		virtual void velocityUpdate(int id)
		{
			for (int i = 0; i < problem_size; i++)
				pswarm[id].velocity[i] = w * pswarm[id].velocity[i]
				+ c1 * r1 * (gbest->result[i] - pswarm[id].result[i])
				+ c2 * r2 * (pswarm[id].pbest->result[i] - pswarm[id].result[i]);
		}

		virtual void positionUpdate(int id)
		{
			double valueBuffer;
			if (model_ini != nullptr)
				model_ini();
			for (int i = 0; i < problem_size; i++)
			{
				valueBuffer = pswarm[id].result[i] + pswarm[id].velocity[i];
				if (repair != nullptr && !constrain_check(i, valueBuffer))
					valueBuffer = repair(i, valueBuffer);
				if (model_change != nullptr)
					model_change(i, valueBuffer);
				pswarm[id].result[i] = valueBuffer;
			}
			evaluator(&pswarm[id]);
		}

		virtual void ini()
		{
			double valueBuffer;
			for (int i = 0; i < swarm_size; i++)
			{
				if (solution_ini_func != nullptr)
				{
					solution_ini_func(pswarm[i].result, problem_size);

					for (int j = 0; j < problem_size; j++)
					{
						if (constrainRangeList[2 * j] != myEC::EMPTYVALUE)
						{
							pswarm[i].velocity[j] = myEC::rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							pswarm[i].velocity[j] = myEC::rand01() * 20000 - 10000;
						}
					}
				}
				else {
					if (model_ini != nullptr)
						model_ini();
					for (int j = 0; j < problem_size; j++)
					{
						if (constrainRangeList[2 * j] != myEC::EMPTYVALUE)
						{
							valueBuffer = myEC::rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
							pswarm[i].velocity[j] = myEC::rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							valueBuffer = myEC::rand01() * 20000 - 10000;
							pswarm[i].velocity[j] = myEC::rand01() * 20000 - 10000;
						}
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
			gbest->replace(pswarm[0].pbest);
			for (int i = 1; i < swarm_size; i++)
			{
				if (*pswarm[i].pbest < *gbest)
					gbest->replace(pswarm[i].pbest);
			}
		}

	public:
		PSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), double objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double c1, double c2,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			this->c1 = c1;
			this->c2 = c2;
			swarm = new Particle[swarm_size];
			gbest = new Solution(problem_size, objectNumber);
			pswarm = (Particle*)swarm;

			for (int i = 0; i < swarm_size; i++)
			{
				pswarm[i].setSize(ps, objectNumber);
				pswarm[i].pbest->setCompareFunc(compare_func);
			}
		}

		~PSO()
		{
			delete gbest;
			swarm = nullptr;
			delete[] pswarm;
		}

		virtual void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = swarm_size; i < generation; i += swarm_size)
			{
				r1 = myEC::rand01();
				r2 = myEC::rand01();
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
							gbest->replace(pswarm[j].pbest);
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


	//Competitive Swarm Optimizer for single object optimization proposed by R. Chen
	class CSO : public PSO
	{
	protected:
		double phi;
		
		Solution* meanSolution;
		bool* update_finished;

		void meanSolution_build()
		{
			if (phi == 0)
				return;

			for (int i = 0; i < problem_size; i++)
			{
				meanSolution->result[i] = 0;
				for (int j = 0; j < swarm_size; j++)
					meanSolution->result[i] += pswarm[j].result[i];
				meanSolution->result[i] /= swarm_size;
			}
		}

		void velocityUpdate(int lid, int wid)
		{
			for (int i = 0; i < problem_size; i++)
				pswarm[lid].velocity[i] = w * pswarm[lid].velocity[i]
				+ c1 * r1 * (pswarm[wid].result[i] - pswarm[lid].result[i])
				+ phi * c2 * r2 * (meanSolution->result[i] - pswarm[lid].result[i]);
		}
	public:
		CSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double phi,
			double c1, double c2,
			int logType)
			:PSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, c1, c2, logType)
		{
			this->phi = phi;
			meanSolution = new Solution(ps, objectNumber);
			update_finished = new bool[swarm_size];
		}

		~CSO()
		{
			delete meanSolution;
			delete[] update_finished;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int fes = swarm_size / 2;
			int id1, id2;
			for (int i = swarm_size; i < generation; i += fes)
			{			
				meanSolution_build();

				for (int j = 0; j < swarm_size; j++)
					update_finished[j] = false;

				for (int j = 1; j < swarm_size; j+=2)
				{
					w = myEC::rand01();
					r1 = myEC::rand01();
					r2 = myEC::rand01();

					//choice two un-update particle
					id1 = rand() % swarm_size;
					while (update_finished[id1])
					{
						id1++;
						if (id1 == swarm_size)
							id1 = 0;
					}
					update_finished[id1] = true;
					id2 = rand() % swarm_size;
					while (update_finished[id2])
					{
						id2++;
						if (id2 == swarm_size)
							id2 = 0;
					}
					update_finished[id2] = true;

					if (pswarm[id1] < pswarm[id2])
						std::swap(id1, id2);

					velocityUpdate(id1, id2);
					positionUpdate(id1);
				}
				for (int j = 0; j < swarm_size; j++)
				{
					if (pswarm[j] < *pswarm[j].pbest)
					{
						pswarm[j].pbest->replace(&pswarm[j]);

						if (*pswarm[j].pbest < *gbest)
						{
							gbest->replace(pswarm[j].pbest);

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

	//Level-based Learning Swarm Optimizer for single object optimization proposed by Q. Yang
	class LLSO : public PSO
	{
	protected:
		int level_number;
		int level_size;
		bool betterReplace;

		Solution* solutionBuffer;

		void level_build()
		{
			sort(pswarm, pswarm + swarm_size);
		}

		void velocityUpdate(int id, int id1, int id2)
		{
			for (int i = 0; i < problem_size; i++)
				pswarm[id].velocity[i] = w * pswarm[id].velocity[i]
				+ c1 * r1 * (pswarm[id1].result[i] - pswarm[id].result[i])
				+ c2 * r2 * (pswarm[id2].result[i] - pswarm[id].result[i]);
		}

		void positionUpdate(int id)
		{
			double valueBuffer;
			if (model_ini != nullptr)
				model_ini();
			for (int i = 0; i < problem_size; i++)
			{
				valueBuffer = pswarm[id].result[i] + pswarm[id].velocity[i];
				if (repair != nullptr && !constrain_check(i, valueBuffer))
					valueBuffer = repair(i, valueBuffer);
				if (model_change != nullptr)
					model_change(i, valueBuffer);
				solutionBuffer->result[i] = valueBuffer;
			}
			evaluator(solutionBuffer);

			if (betterReplace && !(*solutionBuffer < pswarm[id]))
				return;
			pswarm[id].replace(solutionBuffer);
		}

	public:
		LLSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			int level_number, bool betterRelpace,
			double c1, double c2,
			int logType)
			:PSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, c1, c2, logType)
		{
			this->level_number = level_number;
			level_size = ss / level_number;
			if (ss % level_number != 0)
				level_size++;
			solutionBuffer = new Solution(ps, objectNumber);

			this->betterReplace = betterRelpace;
		}

		~LLSO()
		{
			delete solutionBuffer;
		}

		void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();

			int l1, l2;
			int fes = swarm_size - level_size;
			for (int i = swarm_size; i < generation; i += fes)
			{
				level_build();

				for (int l = 1; l < level_number; l++)
				{
					for (int j = 0; j < level_size && l * level_size + j < swarm_size; j++)
					{
						w = myEC::rand01();
						r1 = myEC::rand01();
						r2 = myEC::rand01();

						l1 = rand() % l;
						l2 = rand() % l;

						if (l1 > l2)
							std::swap(l1, l2);
						l1 = l1 * level_size + rand() % level_size;
						l2 = l2 * level_size + rand() % level_size;

						velocityUpdate(l * level_size + j, l1, l2);
						positionUpdate(l * level_size + j);
					}
				}
				for (int j = 0; j < swarm_size; j++)
				{
					if (pswarm[j] < *pswarm[j].pbest)
					{
						pswarm[j].pbest->replace(&pswarm[j]);

						if (*pswarm[j].pbest < *gbest)
						{
							gbest->replace(pswarm[j].pbest);
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
