#pragma once
#include"optimizer.h"

namespace myEC {
	class EAOperator
	{
	public:
		static const enum class crossover { singlepoint, doublepoint };
		static const enum class selection { roulette, championship };
		static const enum class mutation { bit, exchange, overturn };

		static void singlepoint_crossover(Solution* id1, Solution* id2)
		{
			int point = rand() % id1->size;
			double* buffer = new double[point + 1];
			memcpy(buffer, id1->result, point * sizeof(double));
			memcpy(id1->result, id2->result, point * sizeof(double));
			memcpy(id2->result, buffer, point * sizeof(double));
			delete[] buffer;
		}

		static void doublepoint_crossover(Solution* id1, Solution* id2)
		{
			int point = rand() % id1->size;
			double* buffer = new double[point + 1];
			memcpy(buffer, id1->result, point * sizeof(double));
			memcpy(id1->result, id2->result, point * sizeof(double));
			memcpy(id2->result, buffer, point * sizeof(double));

			point = rand() % id1->size;
			buffer = new double[point + 1];
			memcpy(buffer, id1->result, point * sizeof(double));
			memcpy(id1->result, id2->result, point * sizeof(double));
			memcpy(id2->result, buffer, point * sizeof(double));
			delete[] buffer;

		}

		static void championship_selection(Solution* fswarm, size_t ss, Solution* oswarm)
		{
			int id1, id2;
			for (int i = 0; i < ss; i++)
			{
				id1 = rand() % ss;
				id2 = rand() % ss;
				if (fswarm[id1] < fswarm[id2])
					oswarm[i].replace(&fswarm[id1]);
				else oswarm[i].replace(&fswarm[id2]);
			}
		}

		static void roulette_selection(Solution* fswarm, size_t ss, Solution* oswarm)
		{
			double total = 0;
			double* p = new double[ss];
			double max = -1;
			double min = MAX;

			for (int i = 0; i < ss; i++)
			{
				if (max < fswarm[i].fitness[0])
					max = fswarm[i].fitness[0];
				if (min > fswarm[i].fitness[0])
					min = fswarm[i].fitness[0];
			}

			max *= 1.001;
			for (int i = 0; i < ss; i++)
			{
				p[i] = (max - fswarm[i].fitness[0]) / (max - min);
				total += p[i];
			}

			for (int i = 0; i < ss; i++)
				p[i] /= total;
			double p0 = rand01();
			int id = 0;
			int j = 0;
			while (j < ss)
			{
				while (p0 > 0)
				{
					id++;
					if (id == ss)
						id = 0;
					p0 -= p[id];
				}
				oswarm[j].replace(&fswarm[id]);
				j++;
				p0 += rand01();
			}
			delete[] p;
		}

		static void bit_mutation(Solution* id, int times, double* clist)
		{
			int rid;
			for (int i = 0; i < times; i++)
			{
				rid = rand() % id->size;
				if (clist[2 * rid] != EMPTYVALUE)
				{
					id->result[rid] = rand01() * (clist[2 * rid + 1] - clist[2 * rid]) + clist[2 * rid];
				}
				else {
					id->result[rid] = rand01() * 20000 - 10000;
				}
			}

		}

		static void exchange_mutation(Solution* id, int times, double* clist)
		{
			int eid1;
			int eid2;
			for (int i = 0; i < times; i++)
			{
				eid1 = rand() % id->size;
				eid2 = rand() % id->size;
				std::swap(id->result[eid1], id->result[eid2]);
			}
		}

		static void overturn_mutation(Solution* id, int times, double* clist)
		{
			int oid1;
			int oid2;

			for (int i = 0; i < times; i++)
			{
				oid1 = rand() % id->size;
				oid2 = rand() % id->size;
				if (oid1 > oid2)
					std::swap(oid1, oid2);
				while (oid1 < oid2)
					std::swap(id->result[oid1++], id->result[oid2--]);
			}
		}
	};

	class Mutation
	{
	private:
		int mtimes;
		double* clist;
		void (*mutation_func)(Solution* id, int times, double* clist);
	public:
		Mutation() {}

		Mutation(EAOperator::mutation m)
		{
			switch (m)
			{
			case EAOperator::mutation::bit:
				mutation_func = EAOperator::bit_mutation;
				break;
			case EAOperator::mutation::exchange:
				mutation_func = EAOperator::exchange_mutation;
				break;
			case EAOperator::mutation::overturn:
				mutation_func = EAOperator::overturn_mutation;
				break;
			default:
				break;
			}
		}

		Mutation(EAOperator::mutation m, int times)
		{
			switch (m)
			{
			case EAOperator::mutation::bit:
				mutation_func = EAOperator::bit_mutation;
				break;
			case EAOperator::mutation::exchange:
				mutation_func = EAOperator::exchange_mutation;
				break;
			case EAOperator::mutation::overturn:
				mutation_func = EAOperator::overturn_mutation;
				break;
			default:
				break;
			}
			this->mtimes = times;
		}

		Mutation(EAOperator::mutation m, int times, double* clist)
		{
			switch (m)
			{
			case EAOperator::mutation::bit:
				mutation_func = EAOperator::bit_mutation;
				break;
			case EAOperator::mutation::exchange:
				mutation_func = EAOperator::exchange_mutation;
				break;
			case EAOperator::mutation::overturn:
				mutation_func = EAOperator::overturn_mutation;
				break;
			default:
				break;
			}
			this->mtimes = times;
			this->clist = clist;
		}

		~Mutation()
		{
			clist = nullptr;
		}

		void operator()(Solution* id, int times, double* clist) const
		{
			mutation_func(id, times, clist);
		}

		void operator()(Solution* id, double* clist) const
		{
			mutation_func(id, mtimes, clist);
		}

		void operator()(Solution* id) const
		{
			mutation_func(id, mtimes, clist);
		}
	};

	class Crossover
	{
	private:
		void (*crossover_func)(Solution* id1, Solution* id2);
	public:
		Crossover() {}
		Crossover(EAOperator::crossover c)
		{
			switch (c)
			{
			case EAOperator::crossover::singlepoint:
				crossover_func = EAOperator::singlepoint_crossover;
				break;
			case EAOperator::crossover::doublepoint:
				crossover_func = EAOperator::doublepoint_crossover;
				break;
			default:
				break;
			}
		}

		void operator()(Solution* id1, Solution* id2)
		{
			crossover_func(id1, id2);
		}
	};

	class Selection
	{
	private:
		void (*selection_func)(Solution* fswarm, size_t ss, Solution* oswarm);
	public:
		Selection() {}
		Selection(EAOperator::selection s)
		{
			switch (s)
			{
			case EAOperator::selection::championship:
				selection_func = EAOperator::championship_selection;
				break;
			case EAOperator::selection::roulette:
				selection_func = EAOperator::roulette_selection;
				break;
			default:
				break;
			}
		}

		void operator()(Solution* fswarm, size_t ss, Solution* oswarm) const
		{
			selection_func(fswarm, ss, oswarm);
		}
	};
}
