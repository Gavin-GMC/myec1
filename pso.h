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
						if (constrainRangeList[2 * j] != EMPTYVALUE)
						{
							pswarm[i].velocity[j] = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							pswarm[i].velocity[j] = rand01() * 20000 - 10000;
						}
					}
				}
				else {
					if (model_ini != nullptr)
						model_ini();
					for (int j = 0; j < problem_size; j++)
					{
						if (constrainRangeList[2 * j] != EMPTYVALUE)
						{
							valueBuffer = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
							pswarm[i].velocity[j] = rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							valueBuffer = rand01() * 20000 - 10000;
							pswarm[i].velocity[j] = rand01() * 20000 - 10000;
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
					w = rand01();
					r1 = rand01();
					r2 = rand01();

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
				
				gbestUpdate(i);
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
						w = rand01();
						r1 = rand01();
						r2 = rand01();

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
				
				gbestUpdate(i);
			}
			exe_time = time(NULL) - exe_time;
		}
	};

	//Binary Particle Swarm Optimization for binary optimization
	class BPSO : public PSO
	{
	protected:
		virtual void positionUpdate(int id)
		{
			double valueBuffer;
			if (model_ini != nullptr)
				model_ini();
			for (int i = 0; i < problem_size; i++)
			{
				if (rand01() < sigmoid(pswarm[id].velocity[i]))
					valueBuffer = 1;
				else valueBuffer = 0;

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
						pswarm[i].velocity[j] = rand01() * (-6) + 3;
					}
				}
				else {
					if (model_ini != nullptr)
						model_ini();
					for (int j = 0; j < problem_size; j++)
					{
						if (rand01() < 0.5)
							valueBuffer = 0;
						else valueBuffer = 1;
						pswarm[i].velocity[j] = rand01() * (-6) + 3;

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
		BPSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double c1, double c2,
			int logType)
			:PSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, c1, c2, logType)
		{
		}

		~BPSO() { }
	};

	//Sticky Binary Particle Swarm Optimization for binary optimization
	class SBPSO : public BPSO
	{
	protected:
		bool dynamic_modal;
		double is_L;
		double is_U;
		double ustkS_L;
		double ustkS_U;
		double alpha;

		double* stk;
		double is;
		double ip;
		double ig;
		double ustkS;

		virtual void velocityUpdate(int id)
		{
			for (int i = 0; i < problem_size; i++)
				pswarm[id].velocity[i] = is * (1 - stk[id * problem_size + i])
				+ ip * abs(pswarm[id].pbest->result[i] - pswarm[id].result[i])
				+ ig * abs(gbest->result[i] - pswarm[id].result[i]);
		}

		virtual void positionUpdate(int id)
		{
			double valueBuffer;
			if (model_ini != nullptr)
				model_ini();
			for (int i = 0; i < problem_size; i++)
			{
				if (rand01() < pswarm[id].velocity[i])
					valueBuffer = 1 - pswarm[id].result[i];
				else
					valueBuffer = pswarm[id].result[i];

				if (repair != nullptr && !constrain_check(i, valueBuffer))
					valueBuffer = repair(i, valueBuffer);
				if (model_change != nullptr)
					model_change(i, valueBuffer);

				if (pswarm[id].result[i] != valueBuffer)
					stk[id * problem_size + i] = 1;
				else
				{
					stk[id * problem_size + i] -= 1 / ustkS;
					if (stk[id * problem_size + i] < 0)
						stk[id * problem_size + i] = 0;
				}

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
						pswarm[i].velocity[j] = rand01();
					}
				}
				else {
					if (model_ini != nullptr)
						model_ini();
					for (int j = 0; j < problem_size; j++)
					{
						if (rand01() < 0.5)
							valueBuffer = 0;
						else valueBuffer = 1;
						pswarm[i].velocity[j] = rand01();

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

			for (int i = 0; i < swarm_size; i++)
				for (int j = 0; j < problem_size; j++)
					stk[i * problem_size + j] = rand01();
		}

	public:
		SBPSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			int ustkS, double is, double alpha,
			int logType)
			:BPSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, 0, 0, logType)
		{
			if (alpha == EMPTYVALUE)
				alpha = 2;
			this->alpha = alpha;

			if (ustkS == EMPTYVALUE && is == EMPTYVALUE)
			{
				dynamic_modal = true;

				is_L = 1.0 / problem_size;
				is_U = 10.0 / problem_size;
				ustkS_L = double(generation) / swarm_size / 100;
				ustkS_U = double(generation) / swarm_size / 10;
			}
			else {
				dynamic_modal = false;

				if (ustkS == EMPTYVALUE)
					ustkS = 8 * generation / swarm_size / 100;
				this->ustkS = ustkS;

				if (is == EMPTYVALUE)
					is = 4.0 / problem_size;
				this->is = is;

				this->ig = (1 - this->is) / (1 + alpha);
				this->ip = alpha * this->ig;
			}

			stk = new double[swarm_size * problem_size];
		}

		~SBPSO()
		{
			delete[] stk;
		}

		virtual void exe()
		{
			srand(time(NULL));
			exe_time = time(NULL);
			ini();
			for (int i = swarm_size; i < generation; i += swarm_size)
			{
				if (dynamic_modal)
				{
					ustkS = ustkS_L + (ustkS_U - ustkS_L) * (double(i) / generation);
					is = is_U - (is_U - is_L) * (double(i) / generation);
					ig = (1 - is) / (1 + alpha);
					ip = ig * alpha;
				}

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

	//Comprehensive Learning Particle Swarm Optimizer
	class CLPSO : public PSO
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
			PSO::ini();

			for (int i = 0; i < swarm_size; i++)
				fiUpdate(i);
		}

		virtual void velocityUpdate(int id)
		{
			for (int i = 0; i < problem_size; i++)
			{
				if (fi[id * problem_size + i]) 
				{
					pswarm[id].velocity[i] = w * pswarm[id].velocity[i]
						+ c1 * r1 * (pswarm[fid[id]].pbest->result[i] - pswarm[id].result[i]);
				}
				else
				{
					pswarm[id].velocity[i] = w * pswarm[id].velocity[i]
						+ c1 * r1 * (pswarm[id].pbest->result[i] - pswarm[id].result[i]);
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
		CLPSO(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), int objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			double c, int m,
			int logType)
			:PSO(ps, ss, gen, evaluate_func, objectNumber, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, c, 0, logType)
		{
			pc = new double[swarm_size];
			age = new int[swarm_size];
			fid = new int[swarm_size];
			fi = new bool[swarm_size * problem_size];

			this->m = m;

			for (int i = 0; i < swarm_size; i++)
			{
				pc[i] = 0.05 + 0.45 * (pow(E_CONST, (10.0 * i / (swarm_size - 1))) - 1) / (pow(E_CONST, 10) - 1);
			}
			
		}

		~CLPSO() 
		{
			delete[] pc;
			delete[] age;
			delete[] fid;
			delete[] fi;
		}
	};

	//Fuzzy Discrete Particle Swarm Optimization
	class FDPSO : public PSO
	{

	};
}
