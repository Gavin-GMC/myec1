# Name: myec

## Intro: 这是一个尚不成熟的元启发式算法的综合型算法库，旨在提供一个解决复杂非凸优化问题的友好的算法库，或在进行对比试验的时候免去复现他人算法的麻烦。其中封装了遗传算法、分布估计算法、粒子群优化算法、蚁群优化算法和它们的一些变体以及一些迭代式局部搜索算法。用户可以通过指定格式的模型定义将问题传入优化器中并获取最终优化结果，并支持粗细两种不同粒度的优化过程打印。

##Attention: 本算法库中的所有算法目前只进行了有限的测试，不保证在所有情况下的准确性，请谨慎使用运行结果。同时欢迎测试与建议

##Requirements: 目前只支持c++，其他平台暂未列入计划。

##Usage: 用户通过构造器类 "OptimizerBuilder" 进行优化器的生成。在设置好参数后运行优化器类的 build( ) 函数可以返回一个生成好的优化器类“Optimizer”的指针。执行Optimizer对象的运行函数 exe( )进行优化，并通过优化器类的getGbest( )，get_exe_time( )等方法获取最终的优化信息。本算法库所有方法与对象均处于命名空间myEC中，所有的有关调用请加上命名空间。

构造器类 "OptimizerBuilder" ：
必选参数：
type：算法的种类，所有支持的算法种类可以通过myEC::algorithm查看。
size：问题的维度。
evaluate_func：解的评估函数，接收两个参数，第一个参数类型为double*为待评估的解，第二个参数类型为double*为存放评估适应度值的容器。
可选参数：
swarmsize：种群大小，默认100
generation：生成、评估的次数，默认100000
objectNumber：目标数，默认1。如果设置大于1请使用多目标算法或提供适当的比较函数
compare_func：解的比较函数，默认升序且越靠前的目标优先级越高，接收三个参数，前两个参数类型为double*，为待比较的解，第三个参数类型为size_t为解的规模。
solution_ini_func：初始解的生成函数，默认随机生成，接收两个参数，第一个参数类型为double*为储存初始解的容器，第二个参数类型为size_t为解的规模。
model_ini：问题模型初始化函数，默认无，不接收参数。如TSP问题将所有城市设为未造访，背包问题将所有背包置空
constrain_check：问题模型约束函数，默认无，接收两个参数，第一个参数类型为int为问题的维度，第二个参数类型为double为待检测的值。需要注意的是，对于某些算法当is_related设为真时，第一个参数应为上一个的选择如TSP问题中上一步到访的城市
model_change：问题模型的改变函数，表示当前选择对模型的影响，默认无，接收两个参数，第一个参数类型为int为问题的维度，第二个参数类型为double为当前执行的值。如TSP问题将本次到访的城市设为已到访
repair_func：解的修复函数，指定当某一维违反约束后该如何进行修复，默认无，接收两个参数，第一个参数类型为int为问题的维度，第二个参数类型为double为待修复的值。
heuristic_func：启发式函数，默认无，接收两个参数，第一个参数类型为int为问题的维度，第二个参数类型为double为待检测的值。对于某些可利用启发式信息的优化器提供启发式信息，如不设置则使用随机值
logType：日志种类，默认0不打印日志，1表示粗粒度日志，2表示细粒度日志（请谨慎使用细粒度日志，将极大增大优化时间）
其他算法相关参数：可使用set+算法名称检索对于的设置函数。

对于必选参数必须在生成构造器时传入，可选参数可以在生成时传入或生成后使用相应的set函数进行设置

示范样例1：
#include"myec.h"

namespace basicTest1
{
	const int problemsize = 5;

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 0; i < problemsize; i++)
		{
			fitness[0] += pow(solution[i], 2);
		}
	}

	bool constrain(int demensionId, double value)
	{
		return value<10 && value>-10;
	}

	double repair(int demensionId, double value)
	{
		return rand() % 10 - 5;
	}
}

int main()
{
	myEC::OptimizerBuilder b(myEC::algorithm::PSO, basicTest1::problemsize, basicTest1::evaluate, 100, 1e5);
	b.setConstrainFunc(basicTest1::constrain);
	b.setRepaireFunc(basicTest1::repair);
	for (int i = 0; i < basicTest1::problemsize; i++)
		b.setDemensionRange(i, -10, 10);
	//b.setGASelection(myEC::EAOperator::selection::championship);
	b.setLogType(1);

	myEC::Optimizer* o = b.build();
	o->exe();


	myEC::Solution* s = o->getGbest();
	std::string ans = s->ansprint();

	std::cout << ans << std::endl;

	return 0;
}

示范样例2：
#include"myec.h"

namespace tspKroa100 {
	const int CITYNUMBER = 100;
	double connection[CITYNUMBER][CITYNUMBER];
	double citys[CITYNUMBER][2];
	bool visit[CITYNUMBER] = { 0 };
	//int locate_city;

	double dis(double x[2], double y[2])
	{
		return pow(pow((x[0] - y[0]), 2) + pow((x[1] - y[1]), 2), 0.5);
	}

	void get_data()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = 0; j < CITYNUMBER; j++)
				connection[i][j] = -1;

		std::fstream data;
		data.open("kroA100.tsp", std::ios::in);
		int buffer;
		for (int i = 0; i < CITYNUMBER; i++)
			data >> buffer >> citys[i][0] >> citys[i][1];

		for (int i = 0; i < CITYNUMBER; i++)
			for (int j = i; j < CITYNUMBER; j++)
				connection[i][j] = connection[j][i] = round(dis(citys[i], citys[j]));
	}

	double heuristic(int id, double value)
	{
		return 1.0 / connection[id][int(value)];
	}

	void evaluate(double solution[], double fitness[])
	{
		fitness[0] = 0;
		for (int i = 1; i < CITYNUMBER; i++)
		{
			fitness[0] += connection[int(solution[i - 1])][int(solution[i])];
		}
		fitness[0] += connection[int(solution[CITYNUMBER - 1])][int(solution[0])];
	}

	void ini()
	{
		for (int i = 0; i < CITYNUMBER; i++)
			visit[i] = false;
		//locate_city = rand() % CITYNUMBER;
	}

	bool constrain(int demensionId, double value)
	{
		return !visit[int(value)] && connection[demensionId][int(value)] > 0;
	}

	void change(int demensionId, double value)
	{
		//locate_city = value;
		visit[int(value)] = true;
	 }

	double repair(int demensionId, double value)
	{
		int to = rand() % CITYNUMBER;
		while (visit[to])
		{
			to++;
			if (to == CITYNUMBER)
				to = 0;
		}
		return to;
	}

	void greedy(double solution[], size_t size)
	{
		myEC::sortHelper sortbuffer[CITYNUMBER];
		ini();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < CITYNUMBER; j++)
			{
				if (constrain(i, j))
				{
					sortbuffer[j] = myEC::sortHelper(j, heuristic(i, j));
				}
				else
					sortbuffer[j] = myEC::sortHelper(j, -10000);
			}
			sort(sortbuffer, sortbuffer + CITYNUMBER, std::greater<myEC::sortHelper>());
			solution[i] = sortbuffer[0].id;
			change(i, solution[i]);
		}
	}
}

int main()
{
	tspKroa100::get_data();

	myEC::OptimizerBuilder b(myEC::algorithm::AS, CITYNUMBER, evaluate, 100, 1e5);
	b.setConstrainFunc(tspKroa100::constrain);
	b.setRepaireFunc(tspKroa100::repair);
	b.setModelChangFunc(tspKroa100::change);
	b.setChoiceNumber(CITYNUMBER);
	b.setModelIniFunc(tspKroa100::ini);
	b.setHeuristicFunc(tspKroa100::heuristic);
	b.setLogType(1);

	myEC::Optimizer* o = b.build();
	o->exe();


	myEC::Solution* s = o->getGbest();
	std::string ans = s->ansprint();

	std::cout << ans << std::endl;

	return 0;
}

##Contact    csmc_geng@mail.scut.edu.cn

##Authors    Mingcan Geng


