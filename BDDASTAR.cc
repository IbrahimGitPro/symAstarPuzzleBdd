////////////////////////////////////////////////////////
// File  : main.cc   contains the main function
////////////////////////////////////////////////////////

#include <stdio.h> 
#include <stdlib.h>
#include "BDDASTAR.h"
#include <time.h>

void  printUsage();
int main(int argc,char **argv)  {

	if(argc !=2){
	printUsage();
	exit(1);	
	}
	
	int initBddNodes        = 14000000;
	int initCache           = 140000;
	int bddVarNum = -1; // will be instantiated later;
	//Time measure
	clock_t start, finish;
	double runtime;
	int size;
	vector<string> pvars;
	vector<int> ivars;
	vector<int> gvars;
	int bits;
	int fmax;

	start = clock();
	parser(argv[1],size, pvars, ivars,gvars, bits,fmax);	

	// 1) Init Problem
	ProbDescription Problem(pvars, bits, ivars, gvars);	
	bddVarNum = Problem.initialize(fmax); 
	//Problem.print();
		
	// 2) Init BuDDy package
	bdd_init(initBddNodes,initCache);
	bdd_setvarnum(bddVarNum);
	bdd_setmaxincrease(1000000);
	//bdd_setcacheratio(64);
	bdd_setminfreenodes(2);
	
		
	// 3) Build Transition Relation + h and formula structures
	bdd  Trans;
	bdd heuristicBdd;
	bdd init;
	bdd goal;
	FHinfo fhInfo;
	// build transition relation
	Trans = make_transition_relation(Problem, sqrt(Problem.NumVars));
	cout<< "Transition function size: "<< bdd_nodecount(Trans)<<endl;
		
	//make init & goal
	init = make_init(Problem);
	goal = make_goal(Problem);
	
	//init fhInfo
	fhInfo.initialize(Problem);
		
	//build heuristic
	heuristicBdd = SymbolicManhattan(Problem,fhInfo);
	
	//init BDDAinfo
	BDDAinfo bddaInfo;
	bddaInfo.initialize(Problem);
	
	//Call A*
	
	BDDAstar(Problem,fhInfo,bddaInfo,heuristicBdd,Trans,init,goal);
	cout << "BddVarNum=           " << bddVarNum << endl;
	cout << "InitBddNodes=        " << initBddNodes << endl;
	cout << "InitCache=           " << initCache << endl;
	
	finish = clock();
	runtime = (double)(finish-start)/CLOCKS_PER_SEC;
	cout<< "Time: "<< runtime<< " s"<<endl;
	bdd_done();
	return 0; 
}




void  printUsage() {
	cout<< "USAGE:  ./BDDASTAR Problem.txt" <<endl;
	cout<< "Example: ./BDDASTAR Problems/15puzzle16.txt"<<endl;
}
