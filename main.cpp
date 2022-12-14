#include"myec.h"
#include"testset.h"


int main()
{
	/*myEC::OptimizerBuilder b(myEC::algorithm::GA, basicTest1::problemsize, basicTest1::evaluate, 50, 1e4);
	b.setConstrainFunc(basicTest1::constrain);
	b.setRepaireFunc(basicTest1::repair);
	for (int i = 0; i < basicTest1::problemsize; i++)
		b.setDemensionRange(i, -10, 10);
	//b.setGASelection(myEC::EAOperator::selection::championship);
	//b.setGAMutation(0.1);
	b.setLogType(2);*/

	/*myEC::OptimizerBuilder b(myEC::algorithm::GA, basicTest2::problemsize, basicTest2::evaluate, 50, 1e4);
	b.setConstrainFunc(basicTest2::constrain);
	b.setRepaireFunc(basicTest2::repair);
	for (int i = 0; i < basicTest2::problemsize; i++)
		b.setDemensionRange(i, -10, 10);
	b.setCompareFunc(basicTest2::compare);
	//b.setGASelection(myEC::EAOperator::selection::championship);
	b.setGAMutation(0.1);
	b.setLogType(2);*/

	/**/tspKroa100::get_data();

	myEC::OptimizerBuilder b(myEC::algorithm::SPSO, tspKroa100::problemsize, tspKroa100::evaluate, 30, 300 * 100);
	b.setConstrainFunc(tspKroa100::constrain);
	b.setRepaireFunc(tspKroa100::repair);
	b.setModelChangFunc(tspKroa100::change);
	b.setChoiceNumber(tspKroa100::CITYNUMBER);
	b.setModelIniFunc(tspKroa100::ini);
	b.setHeuristicFunc(tspKroa100::heuristic);
	b.setDiscreteProblemType(0, 1);
	b.setSolutionIniFunc(tspKroa100::greedy);
	b.setSPSOparameter(2, 2, 1, 1);
	//b.setS_CLPSOparameter(2, 1, 1);
	b.setLogType(1);
	
	/*
	ACStspKroa100::get_data();

	myEC::OptimizerBuilder b(myEC::algorithm::ACS, ACStspKroa100::problemsize, ACStspKroa100::evaluate, ACStspKroa100::SWARMSIZE, ACStspKroa100::SWARMSIZE * 1e3);
	b.setConstrainFunc(ACStspKroa100::constrain);
	b.setRepaireFunc(ACStspKroa100::repair);
	b.setModelChangFunc(ACStspKroa100::change);
	b.setChoiceNumber(ACStspKroa100::CITYNUMBER);
	b.setModelIniFunc(ACStspKroa100::ini);
	b.setHeuristicFunc(ACStspKroa100::heuristic);
	b.setACSparameter();
	b.setLogType(1);*/

	/*mkpGK01::get_data();
	myEC::OptimizerBuilder b(myEC::algorithm::SPSO, mkpGK01::problemsize, mkpGK01::evaluate, 100, 100 * 1000);
	//b.setGACorssover(0.9, myEC::EAOperator::crossover::uniform);
	//b.setGAMutation(1.0 / mkpGK01::problemsize, myEC::EAOperator::mutation::binary);
	b.setDiscreteProblemType(false, false);
	b.setSPSOparameter(2, 2, 1, 1);
	b.setChoiceNumber(2);
	b.setConstrainFunc(mkpGK01::constrain);
	b.setRepaireFunc(mkpGK01::repair);
	b.setModelChangFunc(mkpGK01::change);
	b.setModelIniFunc(mkpGK01::ini);
	b.setHeuristicFunc(mkpGK01::heuristic);
	b.setLogType(1); */

	/*
	myEC::OptimizerBuilder b(myEC::algorithm::MOEAD, ZDT::ZDT1::problemsize, ZDT::ZDT1::evaluate, 100, 250*100, 2);
	b.setConstrainFunc(ZDT::ZDT1::constrain);
	b.setRepaireFunc(ZDT::ZDT1::repair);
	for (int i = 0; i < ZDT::ZDT1::problemsize; i++)
		b.setDemensionRange(i, 0, 1);
	//b.setGACorssover(0.9, myEC::EAOperator::crossover::SBX, 20);
	//b.setGAMutation(1.0 / ZDT::ZDT1::problemsize, myEC::EAOperator::mutation::PM, 20);
	//b.setGASelection(myEC::EAOperator::selection::championship);
	b.setMOEADCorssover(0.9, myEC::EAOperator::crossover::SBX, 20);
	b.setMOEADMutation(1.0 / ZDT::ZDT1::problemsize, myEC::EAOperator::mutation::PM, 20);
	b.setMOEADDecomponent(10, myEC::MOEAD::D_approach::tchebycheff);*/

	b.setLogType(1);

	myEC::Optimizer* o = b.build();
	o->exe();


	myEC::Solution* s = o->getGbest();

	int gsize = o->getGbestSize();
	for (int i = 0; i < gsize; i++)
	{
		std::string ans = s[i].roughprint();
		//std::string ans = s[i].ansprint();
		std::cout << ans << std::endl;
	}

	std::cout << o->get_exe_time() << std::endl;

	return 0;
}