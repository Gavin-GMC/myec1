#pragma once
#include"optimizer.h"

namespace myEC {
	class EAOperator
	{
	public:
		enum class crossover { singlepoint, doublepoint, SBX, uniform };
		enum class selection { roulette, championship };
		enum class mutation { bit, exchange, overturn, PM, binary };

		static void singlepoint_crossover(Solution* id1, Solution* id2, double exparameter)
		{
			int point = rand() % id1->size;
			double* buffer = new double[point + 1];
			memcpy(buffer, id1->result, point * sizeof(double));
			memcpy(id1->result, id2->result, point * sizeof(double));
			memcpy(id2->result, buffer, point * sizeof(double));
			delete[] buffer;
		}

		static void doublepoint_crossover(Solution* id1, Solution* id2, double exparameter)
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

		static void SBX_crossover(Solution* id1, Solution* id2, double exparameter)
		{
			double eta = exparameter;
			double r;
			double belta;
			double c1, c2;

			for (int i = 0; i < id1->size; i++)
			{
				if (rand01() > 0.5)
					continue;
				else {
					r = rand01_();
					if (r > 0.5)
						belta = pow(2 - r * 2, 1 / (1 + eta));
					else belta = pow(r * 2, 1 / (1 + eta));

					if (rand01() > 0.5)
						belta *= -1;

					c1 = 0.5 * ((1 + belta) * id1->result[i] + (1 - belta) * id2->result[i]);
					c2 = 0.5 * ((1 - belta) * id1->result[i] + (1 + belta) * id2->result[i]);
					id1->result[i] = c1;
					id2->result[i] = c2;
				}				
			}
		}

		static void uniform_crossover(Solution* id1, Solution* id2, double exparameter)
		{
			for (int i = 0; i < id1->size; i++)
			{
				if (rand01() < 0.5)
					std::swap(id1->result[i], id2->result[i]);
			}
		}

		static void championship_selection(Solution* fswarm, size_t ss, Solution* oswarm, int index_num, int offset_length, double* index)
		{
			int id1, id2;
			if (index_num == 0)
			{
				for (int i = 0; i < ss; i++)
				{
					id1 = rand() % ss;
					id2 = rand() % ss;
					if (fswarm[id1] < fswarm[id2])
						oswarm[i].replace(&fswarm[id1]);
					else oswarm[i].replace(&fswarm[id2]);
				}
			}
			else {
				bool better;
				for (int i = 0; i < ss; i++)
				{
					better = false;
					id1 = rand() % ss;
					id2 = rand() % ss;
					
					for (int j = 0; j < index_num; j++)
					{
						if (index[j * offset_length + id1] != index[j * offset_length + id2])
						{
							better = (index[j * offset_length + id1] < index[j * offset_length + id2]);
							break;
						}
					}

					if (better)
						oswarm[i].replace(&fswarm[id1]);
					else oswarm[i].replace(&fswarm[id2]);
				}
			}		
		}

		static void roulette_selection(Solution* fswarm, size_t ss, Solution* oswarm, int index_num, int offset_length, double* index)
		{
			double total = 0;
			double* p = new double[ss];
			double max = -1;
			double min = MAX;

			if (index_num == 0)
			{
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
			}
			else {
				for (int i = 0; i < ss; i++)
				{
					if (max < index[i])
						max = index[i];
					if (min > index[i])
						min = index[i];
				}

				max *= 1.001;
				for (int i = 0; i < ss; i++)
				{
					p[i] = (max - index[i]) / (max - min);
					total += p[i];
				}				
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

		static void bit_mutation(Solution* id, double upb, double lowb, int geneid, double exparameter)
		{
			if (lowb == EMPTYVALUE)
			{
				lowb = -10000;
				upb = 10000;
			}
			id->result[geneid] = rand01() * (upb - lowb) + lowb;
		}

		static void exchange_mutation(Solution* id, double upb, double lowb, int geneid, double exparameter)
		{
			int eid1;
			int eid2;
			//exparameter is the times of exchange
			for (int i = 0; i < exparameter; i++)
			{
				eid1 = rand() % id->size;
				eid2 = rand() % id->size;
				std::swap(id->result[eid1], id->result[eid2]);
			}
		}

		static void overturn_mutation(Solution* id, double upb, double lowb, int geneid, double exparameter)
		{
			int oid1;
			int oid2;

			//exparameter is the times of overturn
			for (int i = 0; i < exparameter; i++)
			{
				oid1 = rand() % id->size;
				oid2 = rand() % id->size;
				if (oid1 > oid2)
					std::swap(oid1, oid2);
				while (oid1 < oid2)
					std::swap(id->result[oid1++], id->result[oid2--]);
			}
		}

		static void PM_mutation(Solution* id, double upb, double lowb, int geneid, double exparameter)
		{
			double eta = exparameter;
			double r;
			double sigma;

			if (lowb == EMPTYVALUE)
			{
				upb = 10000;
				lowb = -10000;
			}

			if (id->result[geneid] > upb)
				id->result[geneid] = upb;
			if (id->result[geneid] < lowb)
				id->result[geneid] = lowb;

			r = rand01_();
			if (r > 0.5)
				sigma = 1 - pow(2 * (1 - r) + 2 * (r - 0.5) * pow(1 - (upb - id->result[geneid]) / (upb - lowb), eta + 1), 1 / (eta + 1));
			else
				sigma = pow(2 * r + (1 - 2 * r) * pow(1 - (id->result[geneid] - lowb) / (upb - lowb), eta + 1), 1 / (eta + 1)) - 1;

			id->result[geneid] += sigma * (upb - lowb);
		}

		static void binary_muattion(Solution* id, double upb, double lowb, int geneid, double exparameter)
		{
			if (id->result[geneid] > 0.5)
				id->result[geneid] = 0;
			else id->result[geneid] = 1;
		}
	};

	struct mutationParameter
	{
		EAOperator::mutation type = EAOperator::mutation::bit;
		double pm = 0.01;
		double additional_parameter = 1;
	};
	class Mutation
	{
	private:
		double pm;
		double* clist;
		double exparameter;
		void (*mutation_func)(Solution* id, double upb, double lowb, int geneid, double exparameter);
		bool mutationtype;

		void setMutationFunc(EAOperator::mutation m)
		{
			switch (m)
			{
			case EAOperator::mutation::bit:
				mutation_func = EAOperator::bit_mutation;
				mutationtype = true;
				break;
			case EAOperator::mutation::exchange:
				mutation_func = EAOperator::exchange_mutation;
				mutationtype = false;
				break;
			case EAOperator::mutation::overturn:
				mutation_func = EAOperator::overturn_mutation;
				mutationtype = false;
				break;
			case EAOperator::mutation::PM:
				mutation_func = EAOperator::PM_mutation;
				mutationtype = true;
				break;
			case EAOperator::mutation::binary:
				mutation_func = EAOperator::binary_muattion;
				mutationtype = true;
				break;
			default:
				break;
			}
		}

		void bychromosome(Solution* id, double pm, double exparameter) const
		{
			if (rand01() < pm)
				mutation_func(id, 0, 0, 0, exparameter);
		}

		void bygene(Solution* id, double pm, double* clist, double exparameter) const
		{
			for (int i = 0; i < id->size; i++)
			{
				if (rand01() < pm)
				{
					mutation_func(id, clist[2 * i + 1], clist[2 * i], i, exparameter);
				}
			}
		}

	public:
		Mutation() {}

		Mutation(EAOperator::mutation m, double pm = 0.01, double* clist = nullptr, double additional_parameter = EMPTYVALUE)
		{
			setMutationFunc(m);
			this->pm = pm;
			this->clist = clist;
			this->exparameter = additional_parameter;
		}

		Mutation(mutationParameter m, double* clist = nullptr)
		{
			setMutationFunc(m.type);
			this->pm = m.pm;
			this->clist = clist;
			this->exparameter = m.additional_parameter;
		}

		~Mutation()
		{
			clist = nullptr;
		}

		void operator()(Solution* id, double pm, double exparameter) const
		{
			if (mutationtype)
				bygene(id, pm, clist, exparameter);
			else bychromosome(id, pm, exparameter);
		}

		void operator()(Solution* id, double pm) const
		{
			if (mutationtype)
				bygene(id, pm, clist, exparameter);
			else bychromosome(id, pm, exparameter);
		}

		void operator()(Solution* id) const
		{
			if (mutationtype)
				bygene(id, pm, clist, exparameter);
			else bychromosome(id, pm, exparameter);
		}
	};

	struct crossoverParameter
	{
		EAOperator::crossover type = EAOperator::crossover::singlepoint;
		double pc = 0.01;
		double additional_parameter = 1;
	};
	class Crossover
	{
	private:
		double exparameter;
		double pc;
		void (*crossover_func)(Solution* id1, Solution* id2, double exparameter);

		void setCrossoverFunc(EAOperator::crossover c)
		{
			switch (c)
			{
			case EAOperator::crossover::singlepoint:
				crossover_func = EAOperator::singlepoint_crossover;
				break;
			case EAOperator::crossover::doublepoint:
				crossover_func = EAOperator::doublepoint_crossover;
				break;
			case EAOperator::crossover::SBX:
				crossover_func = EAOperator::SBX_crossover;
				break;
			case EAOperator::crossover::uniform:
				crossover_func = EAOperator::uniform_crossover;
				break;
			default:
				break;
			}
		}
	public:
		Crossover() {}
		Crossover(EAOperator::crossover c, double pc = 0.9, double exparameter = EMPTYVALUE)
		{
			setCrossoverFunc(c);
			this->pc = 0.9;
			this->exparameter = exparameter;
		}
		Crossover(crossoverParameter c)
		{
			setCrossoverFunc(c.type);
			this->pc = c.pc;
			this->exparameter = c.additional_parameter;
		}

		void operator()(Solution* id1, Solution* id2)
		{
			if (rand01() < pc)
				crossover_func(id1, id2, exparameter);
		}

		void operator()(Solution* id1, Solution* id2, double pc)
		{
			if (rand01() < pc)
				crossover_func(id1, id2, exparameter);
		}

		void operator()(Solution* id1, Solution* id2, double pc,double exparameter)
		{
			if (rand01() < pc)
				crossover_func(id1, id2, exparameter);
		}
	};

	struct selectionParameter
	{
		EAOperator::selection type = EAOperator::selection::championship;
		double additional_parameter = 1;
	};
	const selectionParameter default_sp;
	class Selection
	{
	private:
		double exparameter;
		void (*selection_func)(Solution* fswarm, size_t ss, Solution* oswarm, int index_num, int offset_length, double* index);

		void setSelectionFunc(EAOperator::selection s)
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
	public:
		Selection() {}
		Selection(EAOperator::selection s, double exparameter = EMPTYVALUE)
		{
			setSelectionFunc(s);
			this->exparameter = exparameter;
		}
		Selection(selectionParameter s)
		{
			setSelectionFunc(s.type);
			this->exparameter = s.additional_parameter;
		}

		void operator()(Solution* fswarm, size_t ss, Solution* oswarm) const
		{
			selection_func(fswarm, ss, oswarm, 0, 0, nullptr);
		}

		void operator()(Solution* fswarm, size_t ss, Solution* oswarm, int index_num, int offset_length, double* index) const
		{
			selection_func(fswarm, ss, oswarm, index_num, offset_length, index);
		}
	};
}
