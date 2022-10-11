#pragma once
#include<vector>
#include"optimizer.h"

namespace myEC {
	const int PFSIZE = 200;

	double getSlope(double e1, double f1, double e2, double f2)
	{
		if (f1 == f2)
			std::cerr << "error: The slope doesn't exist\n";

		return atan((e1 - e2) / (f1 - f2));
	}

	//can it set as a private static parameter of MO? 
	sortHelper* MOsortbuffer;
	int MObuffersize;

	class MO
	{
	public:
		template <class T>
		static void fastNonDominatedSort(T swarm[], size_t ss, double* rank_back)
		{
			std::vector<int> n;
			int np, ra;
			std::vector<std::vector<int>> S;
			S.resize(ss);
			std::vector<int>F;
			std::vector<int>Q;

			for (int i = 0; i < ss; i++)
			{
				np = 0;
				for (int j = 0; j < ss; j++)
				{
					if (i != j)
					{
						if (swarm[i] < swarm[j])
							S[i].push_back(j);
						else if (swarm[j] < swarm[i])
							np += 1;
					}
				}
				n.push_back(np);
				if (np == 0)
				{
					rank_back[i] = 1;
					F.push_back(i);
				}
			}

			ra = 1;
			while (!F.empty())
			{
				Q.clear();
				for (int i = 0; i < F.size(); i++)
				{
					for (int j = 0; j < S[F[i]].size(); j++)
					{
						n[S[F[i]][j]] -= 1;
						if (n[S[F[i]][j]] == 0)
						{
							rank_back[S[F[i]][j]] = ra + 1;
							Q.push_back(S[F[i]][j]);
						}
					}
				}
				ra += 1;
				F = Q;
			}
		}

		template <class T = Solution>
		static void CrowdDistance(T swarm[], size_t ss, double* distance_back)
		{
			if (MObuffersize != ss)
			{
				delete[] MOsortbuffer;
				MOsortbuffer = new sortHelper[ss];
				MObuffersize = ss;
			}

			for (int i = 0; i < ss; i++)
				distance_back[i] = 0;

			for (int o = 0; o < swarm->object_number; o++)
			{
				for (int i = 0; i < ss; i++)
				{
					MOsortbuffer[i].id = i;
					MOsortbuffer[i].value = swarm[i].fitness[o];
				}
				std::sort(MOsortbuffer, MOsortbuffer + ss);

				for (int i = 1; i < ss - 1; i++)
				{
					distance_back[MOsortbuffer[i].id] -=
						MOsortbuffer[i + 1].value - MOsortbuffer[i - 1].value;
				}
				distance_back[MOsortbuffer[0].id] = distance_back[MOsortbuffer[ss - 1].id] = -1 * MAX;
			}
		}

		template <class T = Solution>
		static void normalizeCrowdDistance(T swarm[], size_t ss, double* distance_back)
		{
			if (MObuffersize != ss)
			{
				delete[] MOsortbuffer;
				MOsortbuffer = new sortHelper[ss];
				MObuffersize = ss;
			}

			for (int i = 0; i < ss; i++)
				distance_back[i] = 0;

			for (int o = 0; o < swarm->object_number; o++)
			{
				for (int i = 0; i < ss; i++)
				{
					MOsortbuffer[i].id = i;
					MOsortbuffer[i].value = swarm[i].fitness[o];
				}
				std::sort(MOsortbuffer, MOsortbuffer + ss);

				for (int i = 1; i < ss - 1; i++)
				{
					distance_back[MOsortbuffer[i].id] -=
						(MOsortbuffer[i + 1].value - MOsortbuffer[i - 1].value) / (MOsortbuffer[ss - 1].value - MOsortbuffer[0].value);
				}
				distance_back[MOsortbuffer[0].id] = distance_back[MOsortbuffer[ss - 1].id] = -1 * MAX;
			}
		}

		template <class T = Solution>
		static void referenceDRank(T swarm[], size_t ss, double*referencePoint, int objective_number, int rank_number, double* rank_back)
		{
			if (MObuffersize != ss)
			{
				delete[] MOsortbuffer;
				MOsortbuffer = new sortHelper[ss];
				MObuffersize = ss;
			}

			int rank_size = ss / rank_number;
			if (ss % rank_number != 0)
				rank_size++;

			for (int i = 0; i < ss; i++)
			{
				MOsortbuffer[i].id = i;
				MOsortbuffer[i].value = Eu_distance(swarm[i].fitness, referencePoint, objective_number);
			}

			std::sort(MOsortbuffer, MOsortbuffer + ss);
			for (int i = 0; i < ss; i++)
			{
				rank_back[MOsortbuffer[i].id] = i / rank_size;
			}
		}

		//only for bi-objective optimization
		template <class T = Solution>
		static void CrowdAngle(T swarm[], size_t ss, double* angle_back)
		{
			if (MObuffersize != ss)
			{
				delete[] MOsortbuffer;
				MOsortbuffer = new sortHelper[ss];
				MObuffersize = ss;
			}

			double f1range, f2range;
			int length = 0;

			for (int i = 0; i < ss; i++)
			{
				MOsortbuffer[i].id = i;
				MOsortbuffer[i].value = swarm[i].fitness[0];
			}
			std::sort(MOsortbuffer, MOsortbuffer + ss);

			f1range = abs(swarm[MOsortbuffer[0].id].fitness[0] - swarm[MOsortbuffer[ss - 1].id].fitness[0]);
			f2range = abs(swarm[MOsortbuffer[0].id].fitness[1] - swarm[MOsortbuffer[ss - 1].id].fitness[1]);

			//remove the coincident points
			for (int i = 0; i < ss - 1; i++)
			{
				if (MOsortbuffer[i].value == MOsortbuffer[i + 1].value)
				{
					angle_back[MOsortbuffer[i].id] = 1;
					MOsortbuffer[i].id = EMPTYVALUE;
				}
			}
			for (int i = 0; i < ss; i++)
			{
				if (MOsortbuffer[i].id != EMPTYVALUE)
					MOsortbuffer[length++].id = MOsortbuffer[i].id;
			}

			for (int i = 1; i < length - 1; i++)
			{
				angle_back[MOsortbuffer[i].id] = -1 * abs(
					getSlope(swarm[MOsortbuffer[i - 1]].fitness[0] / f1range, swarm[MOsortbuffer[i - 1]].fitness[1] / f2range,
						swarm[MOsortbuffer[i].id].fitness[0] / f1range, swarm[MOsortbuffer[i].id].fitness[1] / f2range)
					- getSlope(swarm[MOsortbuffer[i].id].fitness[0] / f1range, swarm[MOsortbuffer[i].id].fitness[1] / f2range,
						swarm[MOsortbuffer[i + 1].id].fitness[0] / f1range, swarm[MOsortbuffer[i + 1].id].fitness[1] / f2range));
			}
			angle_back[MOsortbuffer[0].id] = angle_back[MOsortbuffer[length - 1].id] = 0;
		}

		template <class T>
		static void buildPartialOrder(T swarm[], size_t ss, double* indicator, int indicator_number)
		{
			bool better;
			double indi_buffer;
			T* t_buffer = (T*)malloc(sizeof(T));

			//insert sort
			for (int i = 0; i < ss; i++)
			{
				for (int j = 0; j < i; j++)
				{
					//compare individual
					better = false;
					for (int k = 0; k < indicator_number; k++)
					{
						if (indicator[k * ss + i] != indicator[k * ss + j])
						{
							better = (indicator[k * ss + i] < indicator[k * ss + j]);
							break;
						}
					}

					//insert individual
					if (better)
					{
						*t_buffer = swarm[i];
						for (int j1 = i; j1 > j; j1--)
							swarm[j1] = swarm[j1 - 1];
						swarm[j] = *t_buffer;

						for (int k = 0; k < indicator_number; k++)
						{
							indi_buffer = indicator[k * ss + i];
							for (int j1 = i; j1 > j; j1--)
								indicator[k * ss + j1] = indicator[k * ss + j1 - 1];
							indicator[k * ss + j] = indi_buffer;
						}
						break;
					}	
				}
			}

			free(t_buffer);
		}

		static bool paretoFrontUpdate(Solution* candidate, Solution pf[], int& pf_size)
		{
			for (int i = 0; i < pf_size; i++)
			{
				if (pf[i] < *candidate)
					return false;
				if (*candidate < pf[i])
				{
					std::swap(pf[i], pf[pf_size - 1]);
					i--;
					pf_size--;
				}
			}
			if (pf_size == (PFSIZE - 1))
				return false;

			pf[pf_size++].replace(candidate);
			return true;
		}
	};
}