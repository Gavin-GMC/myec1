#pragma once
#include"optimizer.h"
#include"eaoperator.h"

namespace myEC {
	class PDA
	{
	public:
		PDA() {}

		~PDA() {}

		virtual void ini(Solution* sset, size_t size) = 0;
		virtual void distribution_estimation(Solution* sset, size_t size, double alpha) = 0;
		virtual void solution_build(Solution* sset) = 0;
	};

	class Gaussian_PDA :public PDA
	{
	private:
		double* mean;
		double* stdv;
		int problem_size;
		double guassion(double mean = 0, double stdv = 1)
		{
			return mean + stdv * (
				sqrt((-2) * log(myEC::rand01()))
				* sin(2 * 3.1415926 * myEC::rand01()));
		}

	public:
		Gaussian_PDA(int problem_size) :PDA()
		{
			mean = new double[problem_size];
			stdv = new double[problem_size];
			this->problem_size = problem_size;
		}

		~Gaussian_PDA()
		{
			delete[] mean;
			delete[] stdv;
		}

		void ini(Solution* sset, size_t size)
		{
			double e;
			double v;

			for (int i = 0; i < problem_size; i++)
			{
				e = 0;
				for (int j = 0; j < size; j++)
					e += sset[j].result[i];
				e /= size;
				mean[i] = e;

				v = 0;
				for (int j = 0; j < size; j++)
					v += pow(sset[j].result[i] - e, 2);
				stdv[i] = sqrt(v / (size - 1));

			}
		}

		void distribution_estimation(Solution* sset, size_t size, double alpha = 0)
		{
			double e;
			double v;

			for (int i = 0; i < problem_size; i++)
			{
				e = 0;
				for (int j = 0; j < size; j++)
					e += sset[j].result[i];
				e /= size;
				mean[i] = e;

				v = 0;
				for (int j = 0; j < size; j++)
					v += pow(sset[j].result[i] - e, 2);
				stdv[i] = sqrt(v / (size - 1));

			}
		}

		void solution_build(Solution* solu)
		{
			for (int i = 0; i < problem_size; i++)
				solu->result[i] = guassion(mean[i], stdv[i]);
		}
	};

	// Classic EDA
	class EDA : public Optimizer
	{
	public:
		static const enum class pda { gaussian, gaussian_PBILc };

	protected:
		Solution* offspring;
		pda pda_type;
		PDA* estimator;
		int good_number;
		double alpha;
		bool elitist_strategy;

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
						if (constrainRangeList[2 * j] != myEC::EMPTYVALUE)
						{
							valueBuffer = myEC::rand01() * (constrainRangeList[2 * j + 1] - constrainRangeList[2 * j]) + constrainRangeList[2 * j];
						}
						else {
							valueBuffer = myEC::rand01() * 20000 - 10000;
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

			estimator->ini(swarm, swarm_size);
		}

	public:
		EDA(int ps, int ss, int gen,
			void (*evaluate_func)(double solution[], double fitness[]), double objectNumber,
			bool (*compare_func)(double f1[], double f2[], size_t size),
			void (*solution_ini_func)(double solution[], size_t size),
			void (*model_ini)(void),
			bool (*constrain_check)(int demensionId, double value),
			void (*model_change)(int demensionId, double value),
			double (*repair_func)(int demensionId, double value),
			pda pda_type, int good_number, double alpha,
			EAOperator::selection selection_type, bool elitist_strategy,
			int logType)
			:Optimizer(ps, ss, gen, evaluate_func, objectNumber, solution_ini_func, model_ini, constrain_check, model_change, repair_func, logType)
		{
			swarm = new Solution[swarm_size * 2];
			gbest = new Solution(problem_size);
			offspring = swarm + swarm_size;

			for (int i = 0; i < (2 * swarm_size); i++)
			{
				swarm[i].setSize(ps, objectNumber);
				swarm[i].setCompareFunc(compare_func);
			}

			this->pda_type = pda_type;
			switch (pda_type)
			{
			case pda::gaussian:
				estimator = new Gaussian_PDA(problem_size);
				break;
			default:
				break;
			}

			this->good_number = good_number;
			this->alpha = alpha;
			selection_func = Selection(selection_type);
			this->elitist_strategy = elitist_strategy;
		}

		~EDA()
		{
			if (offspring < swarm)
				swarm = offspring;
			offspring = nullptr;
			delete[] swarm;
			delete[] gbest;
			delete estimator;
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
				
				sort(swarm, swarm + swarm_size);
				estimator->distribution_estimation(swarm, good_number, alpha);

				for (int j = 0; j < fes; j++)
				{
					estimator->solution_build(&swarm[j]);
					if (repair != nullptr)
					{
						solutioncheck(&swarm[j]);
					}
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
}
