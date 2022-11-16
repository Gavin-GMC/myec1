#pragma once

#include<vector>
#include<iostream>
#include"pso.h"
#include"ga.h"
#include"aco.h"
#include"eda.h"
#include"spso.h"
#include"de.h"
#include"ls.h"
#include"gp.h"

namespace myEC {
	static const enum class algorithm { PSO, CSO, LLSO, GA, NSGA2, AS, ACS, EDA, SPSO, DE, TabuSearch, SA, BPSO, SBPSO, CLPSO, S_CLPSO, MOEAD};

	class OptimizerBuilder
	{
	private:

		int problemsize;
		int swarmsize;
		int generation;
		algorithm type;

		double algorithm_parameter[10];
		double* constrainRangeList;
		int objectNum;

		int logType;

		void (*evaluate_func)(double solution[], double fitness[]);
		bool (*compare_func)(double f1[], double f2[], size_t size);
		void (*model_ini)(void);
		bool (*constrain_check)(int demensionId, double value);
		void (*model_change)(int demensionId, double value);
		double (*repair_func)(int demensionId, double value);
		double (*heuristic_func)(int demensionId, double value);
		void (*solution_ini_func)(double solution[], size_t size);
		void (*search_func)(Solution* id, double* clist);

		selectionParameter sp;
		crossoverParameter cp;
		mutationParameter mp;


	public:
		OptimizerBuilder(algorithm type, int size, void (*evaluate_func)(double solution[], double fitness[]),
			int swarmsize = 100, int generation = 1e5, int objectNumer = 1, int logType = 0,
			bool (*compare_func)(double f1[], double f2[], size_t size) = nullptr, double (*heuristic_func)(int demensionId, double value) = nullptr,
			void (*solution_ini_func)(double solution[], size_t size) = nullptr,
			void(*model_ini)(void) = nullptr, bool (*constrain_check)(int demensionId, double value) = nullptr,
			void (*model_change)(int demensionId, double value) = nullptr, double (*repair_func)(int demensionId, double value) = nullptr)
		{
			this->type = type;
			problemsize = size;
			this->swarmsize = swarmsize;
			this->generation = generation;
			this->evaluate_func = evaluate_func;

			objectNum = objectNumer;

			this->logType = logType;


			this->compare_func = compare_func;
			this->heuristic_func = heuristic_func;
			this->solution_ini_func = solution_ini_func;

			this->model_ini = model_ini;
			this->constrain_check = constrain_check;
			this->model_change = model_change;
			this->repair_func = repair_func;

			constrainRangeList = new double[problemsize * 2];
			for (int i = 0; i < problemsize; i++)
			{
				constrainRangeList[i * 2] = EMPTYVALUE;
			}

			if (type == algorithm::PSO || type == algorithm::BPSO)
			{
				algorithm_parameter[0] = algorithm_parameter[1] = 2;
			}
			else if (type == algorithm::CLPSO)
			{
				algorithm_parameter[0] = 2;
				algorithm_parameter[1] = 5;
			}
			else if (type == algorithm::GA || type == algorithm::NSGA2)
			{
				algorithm_parameter[0] = 0;
			}
			else if (type == algorithm::AS)
			{
				algorithm_parameter[0] = 1;
				algorithm_parameter[1] = 2;
				algorithm_parameter[2] = 0.5;
				algorithm_parameter[3] = -1;
				algorithm_parameter[4] = EMPTYVALUE;
				algorithm_parameter[5] = 0;
				algorithm_parameter[6] = 1;
			}
			else if (type == algorithm::ACS)
			{
				algorithm_parameter[0] = 1;
				algorithm_parameter[1] = 2;
				algorithm_parameter[2] = 0.1;
				algorithm_parameter[3] = -1;
				algorithm_parameter[4] = EMPTYVALUE;//choice number
				algorithm_parameter[5] = 0;//is_direct
				algorithm_parameter[6] = 1;//is_related
				algorithm_parameter[7] = 0.9;
			}
			else if (type == algorithm::LLSO)
			{
				algorithm_parameter[0] = algorithm_parameter[1] = 2;
				algorithm_parameter[2] = 4;
				algorithm_parameter[3] = true;
			}
			else if (type == algorithm::EDA)
			{
				algorithm_parameter[0] = 0;
				algorithm_parameter[1] = 0.5 * swarmsize;
				algorithm_parameter[2] = 0.01;
				algorithm_parameter[3] = 0;
				algorithm_parameter[4] = 0;
			}
			else if (type == algorithm::SPSO)
			{
				algorithm_parameter[0] = 2;
				algorithm_parameter[1] = 2;
				algorithm_parameter[2] = 0;
				algorithm_parameter[3] = 1;
				algorithm_parameter[4] = EMPTYVALUE;
				algorithm_parameter[5] = 1;
				algorithm_parameter[6] = 6;
			}
			else if (type == algorithm::S_CLPSO)
			{
				algorithm_parameter[0] = 2;
				algorithm_parameter[1] = 5;
				algorithm_parameter[2] = 0;
				algorithm_parameter[3] = 1;
				algorithm_parameter[4] = EMPTYVALUE;
				algorithm_parameter[5] = 1;
				algorithm_parameter[6] = 6;
			}
			else if (type == algorithm::CSO)
			{
				algorithm_parameter[0] = algorithm_parameter[1] = 2;
				algorithm_parameter[2] = 0;
			}
			else if (type == algorithm::DE)
			{
				algorithm_parameter[0] = algorithm_parameter[1] = EMPTYVALUE;
			}
			else if (type == algorithm::TabuSearch)
			{
				algorithm_parameter[0] = problemsize;
				algorithm_parameter[1] = problemsize;
				this->search_func = nullptr;
			}
			else if (type == algorithm::SA)
			{
				algorithm_parameter[0] = EMPTYVALUE;
				algorithm_parameter[1] = 0;
				algorithm_parameter[2] = EMPTYVALUE;
				algorithm_parameter[3] = problemsize;
				this->search_func = nullptr;
			}
			else if (type == algorithm::SBPSO)
			{
				algorithm_parameter[0] = algorithm_parameter[1] = EMPTYVALUE;
				algorithm_parameter[2] = 2;
			}
		}

		~OptimizerBuilder() {
			delete[]constrainRangeList;
		}

		void setSwarmSize(int ss)
		{
			this->swarmsize = ss;
		}

		void setGeneration(int gen)
		{
			this->generation = gen;
		}

		void setObjectNumber(int number)
		{
			this->objectNum = number;
		}

		void setCompareFunc(bool (*compare_func)(double f1[], double f2[], size_t size))
		{
			this->compare_func = compare_func;
		}

		void setSolutionIniFunc(void (*solution_ini_func)(double solution[], size_t size))
		{
			this->solution_ini_func = solution_ini_func;
		}

		void setModelIniFunc(void (*model_ini)(void))
		{
			this->model_ini = model_ini;
		}

		void setConstrainFunc(bool (*constrain_check)(int demensionId, double value))
		{
			this->constrain_check = constrain_check;
		}

		void setModelChangFunc(void (*model_change)(int demensionId, double value))
		{
			this->model_change = model_change;
		}

		void setRepaireFunc(double (*repair_func)(int demensionId, double value))
		{
			this->repair_func = repair_func;
		}

		void setHeuristicFunc(double (*heuristic_func)(int demensionId, double value))
		{
			this->heuristic_func = heuristic_func;
		}

		void setDemensionRange(int demensionID, double minV, double maxV)
		{
			if (demensionID >= problemsize)
			{
				std::cerr << "demension id overlap the boundary\n";
			}
			constrainRangeList[demensionID * 2] = minV;
			constrainRangeList[demensionID * 2 + 1] = maxV;
		}

		void setChoiceNumber(int number)
		{
			for (int j = 0; j < problemsize; j++)
			{
				constrainRangeList[j * 2] = 0;
				constrainRangeList[j * 2 + 1] = number;
			}
		}

		//0 for no log, 1 for rough log, 2 for detailed log
		void setLogType(int type)
		{
			logType = type;
		}

		void setPSOparameter(double c1 = 2, double c2 = 2)
		{
			if (type !=algorithm::PSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c1;
			algorithm_parameter[1] = c2;
		}

		void setBPSOparameter(double c1 = 2, double c2 = 2)
		{
			if (type != algorithm::BPSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c1;
			algorithm_parameter[1] = c2;
		}

		//ustkS & is both EMPTYVALUE for adaptive setting, ustkS: 8/iteration  is:4/problem_size, alpha: 2
		void setSBPSOparameter(int ustkS = EMPTYVALUE, double is=EMPTYVALUE, double alpha = 2)
		{
			if (type != algorithm::SBPSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = ustkS;
			algorithm_parameter[1] = is;
			algorithm_parameter[2] = alpha;
		}

		//phi|swarmsize : 0|ss<100, 0-0.1|ss=200, 0.1-0.2|ss=400-600, 0.1-0.3|ss=1000
		void setCSOparameter(double c1 = 2, double c2 = 2, double phi = 0)
		{
			if (type != algorithm::CSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c1;
			algorithm_parameter[1] = c2;
			algorithm_parameter[2] = phi;
		}

		//pc usually between 0.4-0.99
		void setGACorssover(double pc = 0.9, EAOperator::crossover crossover_type = EAOperator::crossover::singlepoint, double extra_parameter = 1)
		{
			if (type != algorithm::GA && type != algorithm::NSGA2)
			{
				std::cerr << "wrong setting of crossover operation\n";
				return;
			}
			cp.pc = pc;
			cp.type = crossover_type;
			cp.additional_parameter = extra_parameter;
		}

		//pm usually between 0.001-0.1
		void setGAMutation(double pm = 0.01, EAOperator::mutation mutation_type = EAOperator::mutation::bit, double extra_parameter = 1)
		{
			if (type != algorithm::GA && type != algorithm::NSGA2)
			{
				std::cerr << "wrong setting of mutation operation\n";
				return;
			}
			mp.pm = pm;
			mp.type = mutation_type;
			mp.additional_parameter = extra_parameter;
		}

		void setGASelection(EAOperator::selection selection_type = EAOperator::selection::roulette, bool elitist_strategy = false, double extra_parameter = 1)
		{
			if (type != algorithm::GA && type != algorithm::NSGA2)
			{
				std::cerr << "wrong setting of selection operation\n";
				return;
			}
			algorithm_parameter[0] = elitist_strategy;
			sp.type = selection_type;
			sp.additional_parameter = extra_parameter;
		}

		void setASparameter(double alpha = 1, double belta = 2, double rho = 0.5, double t0 = -1)
		{
			if (type != algorithm::AS && type!=algorithm::ACS)
			{
				std::cerr << "wrong setting of algorithm type\n";
				return;
			}
			algorithm_parameter[0] = alpha;
			algorithm_parameter[1] = belta;
			algorithm_parameter[2] = rho;
			algorithm_parameter[3] = t0;
		}

		void setACSparameter(double alpha = 1, double belta = 2, double rho = 0.1, double p0 = 0.9, double t0 = -1)
		{
			if (type != algorithm::AS && type != algorithm::ACS)
			{
				std::cerr << "wrong setting of algorithm type\n";
				return;
			}
			algorithm_parameter[0] = alpha;
			algorithm_parameter[1] = belta;
			algorithm_parameter[2] = rho;
			algorithm_parameter[3] = t0;
			algorithm_parameter[7] = p0;
		}

		void setDiscreteProblemType(bool is_direct = false, bool is_related = true)
		{
			if (type != algorithm::AS && type != algorithm::ACS 
				&& type != algorithm::SPSO && type != algorithm::S_CLPSO)
			{
				std::cerr << "wrong setting of algorithm type\n";
				return;
			}

			algorithm_parameter[5] = is_direct;
			algorithm_parameter[6] = is_related;
		}

		void setLLSOparameter(double c1 = 2, double c2 = 2, int levelnumber = 4, bool betterReplace = true)
		{
			if (type != algorithm::LLSO)
			{
				std::cerr << "wrong setting of algorithm type\n";
				return;
			}
			algorithm_parameter[0] = c1;
			algorithm_parameter[1] = c2;
			algorithm_parameter[2] = levelnumber;
			algorithm_parameter[3] = betterReplace;
		}

		void setCLPSOparameter(double c = 2, double m = 5)
		{
			if (type != algorithm::CLPSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c;
			algorithm_parameter[1] = m;
		}

		void setEDASelection(EAOperator::selection selection_type = EAOperator::selection::roulette, bool elitist_strategy = false)
		{
			if (type != algorithm::EDA)
			{
				std::cerr << "wrong setting of selection operation\n";
				return;
			}
			algorithm_parameter[3] = double(selection_type);
			algorithm_parameter[4] = elitist_strategy;
		}

		void setEDAEstimatorType(EDA::pda pda_type = EDA::pda::gaussian, int good_number = 50, double attenuation = 0.01)
		{
			if (type != algorithm::EDA)
			{
				std::cerr << "wrong setting of estimator\n";
				return;
			}
			algorithm_parameter[0] = double(pda_type);
			algorithm_parameter[1] = good_number;
			algorithm_parameter[2] = attenuation;
		}

		void setSPSOparameter(double c1 = 2, double c2 = 2, bool v_heuristic = false, bool f_heuristic = true)
		{
			if (type != algorithm::SPSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c1;
			algorithm_parameter[1] = c2;

			algorithm_parameter[2] = v_heuristic;
			algorithm_parameter[3] = f_heuristic;
		}

		void setS_CLPSOparameter(double c = 2, int m = 5, bool v_heuristic = false, bool f_heuristic = true)
		{
			if (type != algorithm::S_CLPSO)
			{
				std::cerr << "wrong setting of parameter c\n";
				return;
			}
			algorithm_parameter[0] = c;
			algorithm_parameter[1] = m;
			algorithm_parameter[2] = v_heuristic;
			algorithm_parameter[3] = f_heuristic;
		}

		//f usually between 0-2, cr between 0.1-0.6, and EMPTYVALUE for adaptive setting
		void setDEparameter(double f = 0.5, double cr = EMPTYVALUE)
		{
			algorithm_parameter[0] = f;
			algorithm_parameter[1] = cr;
		}

		void setTabuSearchparameter(int tabulength = 20, void (*search_func)(Solution* id, double* clist) = nullptr, int searchtimes = 20)
		{
			algorithm_parameter[0] = tabulength;
			algorithm_parameter[1] = searchtimes;
			this->search_func = search_func;
		}

		void setSAparameter(double t0 = EMPTYVALUE, SA::cooling_type cool_func=SA::cooling_type::classic, double p0=EMPTYVALUE,
			void (*search_func)(Solution* id, double* clist)=nullptr, int searchtimes=20)
		{
			algorithm_parameter[0] = t0;
			algorithm_parameter[1] = double(cool_func);
			algorithm_parameter[2] = p0;
			algorithm_parameter[3] = searchtimes;
			this->search_func = search_func;
		}

		//pc usually between 0.4-0.99
		void setMOEADCorssover(double pc = 0.9, EAOperator::crossover crossover_type = EAOperator::crossover::singlepoint, double extra_parameter = 1)
		{
			if (type != algorithm::MOEAD)
			{
				std::cerr << "wrong setting of crossover operation\n";
				return;
			}
			cp.pc = pc;
			cp.type = crossover_type;
			cp.additional_parameter = extra_parameter;
		}

		//pm usually between 0.001-0.1
		void setMOEADMutation(double pm = 0.01, EAOperator::mutation mutation_type = EAOperator::mutation::bit, double extra_parameter = 1)
		{
			if (type != algorithm::MOEAD)
			{
				std::cerr << "wrong setting of mutation operation\n";
				return;
			}
			mp.pm = pm;
			mp.type = mutation_type;
			mp.additional_parameter = extra_parameter;
		}

		void setMOEADDecomponent(int neighborsize=10, MOEAD::D_approach decomponent_type = MOEAD::D_approach::tchebycheff)
		{
			algorithm_parameter[0] = int(decomponent_type);
			algorithm_parameter[1] = neighborsize;
		}

		Optimizer* build()
		{
			Optimizer* back;
			if (type == algorithm::PSO)
			{
				back = new PSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::CSO){
				back = new CSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[2], algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::LLSO){
				back = new LLSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func,
					algorithm_parameter[2], algorithm_parameter[3], algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::GA) {
				back = new GA(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, cp, mp, sp, algorithm_parameter[0], logType);
			}
			else if (type == algorithm::NSGA2){
				back = new NSGA2(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, cp, mp, sp, logType);
			}
			else if (type == algorithm::AS) {
				for (int i = 0; i < problemsize; i++)
				{
					if (constrainRangeList[2 * i] == EMPTYVALUE)
					{
						std::cerr << "failing optimizer build for not set choice range\n";
						return nullptr;
					}
					else if (constrainRangeList[2 * i + 1] > algorithm_parameter[4])
						algorithm_parameter[4] = constrainRangeList[2 * i + 1];
				}

				back = new AS(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
					algorithm_parameter[0], algorithm_parameter[1], algorithm_parameter[2], algorithm_parameter[3], int(algorithm_parameter[4]),
					algorithm_parameter[5], algorithm_parameter[6], logType);
			}
			else if (type == algorithm::ACS) {
				for (int i = 0; i < problemsize; i++)
				{
					if (constrainRangeList[2 * i] == EMPTYVALUE)
					{
						std::cerr << "failing optimizer build for not set choice range\n";
						return nullptr;
					}
					else if (constrainRangeList[2 * i + 1] > algorithm_parameter[4])
						algorithm_parameter[4] = constrainRangeList[2 * i + 1];
				}

				back = new ACS(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
					algorithm_parameter[0], algorithm_parameter[1], algorithm_parameter[2], algorithm_parameter[7], algorithm_parameter[3], int(algorithm_parameter[4]),
					algorithm_parameter[5], algorithm_parameter[6], logType);
			}
			else if (type == algorithm::EDA){
				back = new EDA(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, EDA::pda(algorithm_parameter[0]), algorithm_parameter[1],
					algorithm_parameter[2], sp, algorithm_parameter[4], logType);
			}
			else if (type == algorithm::SPSO){
				for (int i = 0; i < problemsize; i++)
				{
					if (constrainRangeList[2 * i] == EMPTYVALUE)
					{
						std::cerr << "failing optimizer build for not set choice range\n";
						return nullptr;
					}
					else if (constrainRangeList[2 * i + 1] > algorithm_parameter[2])
						algorithm_parameter[4] = constrainRangeList[2 * i + 1];
				}

				back = new SPSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
					algorithm_parameter[0], algorithm_parameter[1], bool(algorithm_parameter[2]), bool(algorithm_parameter[3]), int(algorithm_parameter[4]), bool(algorithm_parameter[5]), bool(algorithm_parameter[6]), logType);
			}
			else if (type == algorithm::S_CLPSO) {
				for (int i = 0; i < problemsize; i++)
				{
					if (constrainRangeList[2 * i] == EMPTYVALUE)
					{
						std::cerr << "failing optimizer build for not set choice range\n";
						return nullptr;
					}
					else if (constrainRangeList[2 * i + 1] > algorithm_parameter[2])
						algorithm_parameter[4] = constrainRangeList[2 * i + 1];
				}

				back = new S_CLPSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, heuristic_func,
					algorithm_parameter[0], algorithm_parameter[1], bool(algorithm_parameter[2]), bool(algorithm_parameter[3]), int(algorithm_parameter[4]), bool(algorithm_parameter[5]), bool(algorithm_parameter[6]), logType);
			}
			else if (type == algorithm::DE) {
				back = new DE(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::TabuSearch) {
				back = new TabuSearch(problemsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0],search_func, algorithm_parameter[1], logType);
			}
			else if (type == algorithm::SA) {
				back = new SA(problemsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func,
					algorithm_parameter[0],SA::cooling_type(algorithm_parameter[1]),algorithm_parameter[2], search_func, algorithm_parameter[3], logType);
			}
			else if (type == algorithm::BPSO)
			{
				back = new BPSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::SBPSO)
			{
				back = new SBPSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0], algorithm_parameter[1], algorithm_parameter[2], logType);
			}
			else if (type == algorithm::CLPSO) {
				back = new CLPSO(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, algorithm_parameter[0], algorithm_parameter[1], logType);
			}
			else if (type == algorithm::MOEAD) {
				back = new MOEAD(problemsize, swarmsize, generation, evaluate_func, objectNum, compare_func, solution_ini_func, model_ini, constrain_check, model_change, repair_func, cp, mp, MOEAD::D_approach(algorithm_parameter[0]), algorithm_parameter[1], logType);
			}

			else {
				std::cerr << "failing optimizer build\n";
				return nullptr;
			}

			back->constructConstrainRangeList(constrainRangeList, 2 * problemsize);
			return back;
		}
	};
}
