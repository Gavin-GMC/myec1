#include"myec.h"
#include"testset.h"


int main()
{
	/*myEC::OptimizerBuilder b(myEC::algorithm::CSO, basicTest1::problemsize, basicTest1::evaluate, 50, 1e4);
	b.setConstrainFunc(basicTest1::constrain);
	b.setRepaireFunc(basicTest1::repair);
	for (int i = 0; i < basicTest1::problemsize; i++)
		b.setDemensionRange(i, -10, 10);
	//b.setGASelection(myEC::EAOperator::selection::championship);
	b.setLogType(2);*/

	/*tspKroa100::get_data();

	myEC::OptimizerBuilder b(myEC::algorithm::AS, tspKroa100::problemsize, tspKroa100::evaluate, 20, 20 * 1e3);
	b.setConstrainFunc(tspKroa100::constrain);
	b.setRepaireFunc(tspKroa100::repair);
	b.setModelChangFunc(tspKroa100::change);
	b.setChoiceNumber(tspKroa100::CITYNUMBER);
	b.setModelIniFunc(tspKroa100::ini);
	b.setHeuristicFunc(tspKroa100::heuristic);
	b.setLogType(1);*/

	
	ACStspKroa100::get_data();

	myEC::OptimizerBuilder b(myEC::algorithm::ACS, ACStspKroa100::problemsize, ACStspKroa100::evaluate, ACStspKroa100::SWARMSIZE, ACStspKroa100::SWARMSIZE * 1e3);
	b.setConstrainFunc(ACStspKroa100::constrain);
	b.setRepaireFunc(ACStspKroa100::repair);
	b.setModelChangFunc(ACStspKroa100::change);
	b.setChoiceNumber(ACStspKroa100::CITYNUMBER);
	b.setModelIniFunc(ACStspKroa100::ini);
	b.setHeuristicFunc(ACStspKroa100::heuristic);
	b.setACSparameter();
	//b.setSPSOProblemType(false, true);
	//b.setSPSOparameter(2, 2, 1, 1);
	b.setLogType(1);

	myEC::Optimizer* o = b.build();
	o->exe();


	myEC::Solution* s = o->getGbest();
	std::string ans = s->ansprint();

	std::cout << ans << std::endl;

	return 0;
}