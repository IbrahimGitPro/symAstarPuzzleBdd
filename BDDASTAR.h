#ifndef BDDASTAR_H
#define BDDASTAR_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <bdd.h>
#include <bvec.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <stdio.h>




using namespace std;


struct  varDef {
	int        length;        // num of variables
	int        *current;      // mapping to current state bddvars 
	int        *next;         // mapping to next state bddvars
	varDef() {}
	varDef(int l,int *c, int *n) 
	{ length = l; current = c; next = n; }
	void print();
};


// structure  for representing binary number expresions
struct binary {
	int  digitNum;
	bdd *digit;
	binary(int dn, bdd *d)
    { digitNum = dn; digit = d; }
};




struct ProbDescription {
	vector<string> variablesNames;
	int NumVars;       //number of variables in the problem
	int Tile_Length;   //number of bit to encode a tile
	int fhMax;          // maxvalue of f and h
	varDef fVars;       // f-value bit vector representation
	varDef hVars;       // h-value bit vector representation
	vector<varDef> States;       // States definition
	vector<int> init; //start state
	vector<int> goal; //goal state
	ProbDescription(vector<string> varnames, int tilelen, vector<int> Init, vector<int> Goal)
	{variablesNames = varnames; NumVars = varnames.size();
		Tile_Length = tilelen; init = Init; goal = Goal;
	}
	int initialize(int fhMax);
	void print();
};


// BDDA* 
// structure holding quantification info
// of F and H variables
struct FHinfo {
	bdd Fcurrent;
	bdd Fnext;    
	bddPair *Fcurrent2next;
	bddPair *Fnext2current;  
	bdd Hcurrent;
	bdd Hnext;
	bddPair *Hcurrent2next;
	bddPair *Hnext2current;
	bddPair *Hcurrent2Fcurrent;
	
	void initialize(ProbDescription& Problem);
};


// BDDA*
// structure holding quantification info
// needed by BDDA*
struct BDDAinfo {
	set< pair<int,int> > stateVarPairs;
	// currentVar nextVar pairs
	bdd currentVars;        // current BDD vars
	bdd nextVars;           // next BDD vars
	bddPair *current2next;  // current to next substitution set
	bddPair *next2current;  // next to current substitution set
	vector<bdd> formulaBackward;    // f = (fmin + 1) + h - h'
	vector<bdd> formulaForward;     // f = (fmin - 1) + h - h'    
	void initialize(ProbDescription& Problem);
};


/*returns varnum*/
int ProbDescription::initialize(int fhMAX) {
	
	cout << "Intialize Problem..."<<endl;
	varDef vtmp;
	int bddVarNum = 0;
	fhMax = fhMAX;
	// 1) assign the first BDD variables to the f-encoding and h-encoding
	fVars.length = int(ceil(log(fhMax)/log(2)));
	fVars.current = new int[fVars.length];
	fVars.next = new int[fVars.length];
	for (int i = 0; i < fVars.length; i++)
		fVars.current[i] = bddVarNum++;
	for (int i = 0; i < fVars.length; i++)
		fVars.next[i] = bddVarNum++;
    
	hVars.length = int(ceil(log(fhMax)/log(2)));
	hVars.current = new int[hVars.length];
	hVars.next = new int[hVars.length];
	for (int i = 0; i < hVars.length; i++)
		hVars.current[i] = bddVarNum++;
	for (int i = 0; i < hVars.length; i++)
		hVars.next[i] = bddVarNum++;

	for (int i= 0; i<NumVars; i++) {
		 
		int *current = new int[Tile_Length];
		int *next = new int[Tile_Length];
		vtmp = varDef(Tile_Length,current, next);
		States.push_back(vtmp);
	}
	
	for (int j= 0; j<States.size(); j++) {
		for (int i = 0; i < States[j].length; i++)
			States[j].current[i] = bddVarNum++;
		for (int i = 0; i < States[j].length; i++)
			States[j].next[i] = bddVarNum++;
	}
	
	
	
	
	
	return bddVarNum;
	
}
//conversion of int to bdd representation
bdd int2bdd(int *var,int length,int n)  {
	bdd res;
	int j;  
	
	res = bddtrue;
	for (j=length-1; j >= 0; j--)
    {
		if (n % 2) 
			res &= bdd_ithvar(var[j]); 
		else
			res &= bdd_nithvar(var[j]);
		n /= 2;
    }
	return res;
}

void FHinfo::initialize(ProbDescription& Problem) {
	cout << "Intialize FHinfo..."<<endl;
	Fcurrent = bddtrue;
	Fnext = bddtrue;
	Fcurrent2next = bdd_newpair();
	Fnext2current = bdd_newpair();
	
	for (int i = 0; i < Problem.fVars.length; i++)
	{
		Fcurrent &= bdd_ithvar(Problem.fVars.current[i]);
		Fnext    &= bdd_ithvar(Problem.fVars.next[i]);
		bdd_setpair(Fcurrent2next,Problem.fVars.current[i],Problem.fVars.next[i]);
		bdd_setpair(Fnext2current,Problem.fVars.next[i],Problem.fVars.current[i]);
	}
	
	Hcurrent = bddtrue;
	Hnext    = bddtrue;
	Hcurrent2next = bdd_newpair();
	Hnext2current = bdd_newpair();
	
	for (int i = 0; i < Problem.hVars.length; i++)
	{
		Hcurrent &= bdd_ithvar(Problem.hVars.current[i]);
		Hnext    &= bdd_ithvar(Problem.hVars.next[i]);
		bdd_setpair(Hcurrent2next,Problem.hVars.current[i],Problem.hVars.next[i]);
		bdd_setpair(Hnext2current,Problem.hVars.next[i],Problem.hVars.current[i]);
	}
	
	Hcurrent2Fcurrent = bdd_newpair();
	for (int i = 0; i < Problem.hVars.length; i++)
		bdd_setpair(Hcurrent2Fcurrent,Problem.hVars.current[i],Problem.fVars.current[i]);	    	    	 	    
}

//
void BDDAinfo::initialize(ProbDescription& Problem) {
	
	cout << "Intialize BDDAinfo..."<<endl;
	int totalFormulaSize =0;
	
	// make current and next var sets and
	// next state substitution sets
	currentVars = bddtrue;
	nextVars    = bddtrue;
	current2next = bdd_newpair();
	next2current = bdd_newpair();
	
	for (int i = 0; i < Problem.NumVars; i++) 
		for (int j = 0; j < Problem.Tile_Length; j++)
		{
			pair<int,int> elem(Problem.States[i].current[j],Problem.States[i].next[j]);
			stateVarPairs.insert(elem);	 
			currentVars &= bdd_ithvar(Problem.States[i].current[j]);
			nextVars    &= bdd_ithvar(Problem.States[i].next[j]);
			bdd_setpair(current2next,Problem.States[i].current[j],Problem.States[i].next[j]);
			bdd_setpair(next2current,Problem.States[i].next[j],Problem.States[i].current[j]);
		}

	
	// generate BDDA* formulas for finding new f-values 
	bvec h  = bvec_varvec(Problem.hVars.length,Problem.hVars.current);
	bvec hm = bvec_varvec(Problem.hVars.length,Problem.hVars.next);
	bvec f  = bvec_varvec(Problem.fVars.length,Problem.fVars.current);
	formulaBackward.resize(Problem.fhMax);
	// make a formula for each possible fmin values
	for (int i = 1; i < Problem.fhMax; i++)
    {
		bvec sum = bvec(Problem.hVars.length,i-1);
		// make f = (fmin - 1) + h - h'
		sum = bvec_add(sum,h);
		sum = bvec_sub(sum,hm);
		
		formulaBackward[i] = (f == sum);
		totalFormulaSize += bdd_nodecount(formulaBackward[i]);
		
    }
	
	formulaForward.resize(Problem.fhMax);
	// make a formula for each possible fmin values
	for (int i = 0; i < Problem.fhMax; i++)
    {
		bvec sum = bvec(Problem.hVars.length,i+1);
		// make f = (fmin + 1) + h' - h 
		sum = bvec_add(sum,hm);
		sum = bvec_sub(sum,h);
		
		formulaForward[i] = (f == sum);
		totalFormulaSize += bdd_nodecount(formulaForward[i]);
    }
	
	cout << "TotalSize of formulas=" << totalFormulaSize << endl;
}

//
int anAbs(int x) {
	if (x > 0) 
		return x;
	else
		return -x;
}


// Compute Manhattan distance for the heuristic
int manhattanDist(int x,int y,int j,int posNum) {
	
	int dim = int(sqrt(double(posNum)));
	
	// find goal x and y position
	int xg = (j-1) % dim;
	int yg = (dim - 1) - ((j-1) / dim);
	
	return ( anAbs(x - xg) + anAbs(y - yg) );
} 

// Building Manhattan heuristic
bdd SymbolicManhattan(ProbDescription& Problem, FHinfo& fhInfo) {
	
	//     initially exp(h,state) = (h = 0)
	//     for each position for each tile
	//             exp(h,h',state) := exp(h,state) /\ if (state(position)) = tile then (h' = h + dist(position,tile)) else (h' = h) 
	//             exp(h',state) := exist h
	//             exp(h,state) := exp[h'/h] 

		cout << "Building Manhattan heuristic..."<<endl;
	
	
	bvec h     = bvec_varvec(Problem.hVars.length,Problem.hVars.current);
	bvec hm    = bvec_varvec(Problem.hVars.length,Problem.hVars.next);
	
	
	// init exp
	bdd exp = (bvec_varvec(Problem.hVars.length,Problem.hVars.current) == bvec(Problem.hVars.length,0));
	
	// find number of positions in problem
	// (= tile num + 1 since we also track the empty pos)
	int posNum =  Problem.NumVars;
	
	
	// iterate over positions
	for (int i = 0; i < posNum; i++)
    {
		// find x and y of position
		char buf[32];
		strcpy(buf,(Problem.variablesNames[i]).c_str());
		
		int x = int(buf[1] - '0');
		int y = int(buf[2] - '0');
		
		// iterate over tile numbers
		for (int j = 1; j < posNum; j++)
		{
			bvec hPlusd = bvec_add(bvec_varvec(Problem.hVars.length,Problem.hVars.current),
								   bvec(Problem.hVars.length,manhattanDist(x,y,j,posNum)));
			
			bdd r1 = hm == hPlusd;
			bdd r2 = hm == h;
			bdd l =  int2bdd(Problem.States[i].current,Problem.States[i].length,j);
			exp = exp & bdd_ite(l,r1,r2);
			exp = bdd_exist(exp,fhInfo.Hcurrent);
			exp = bdd_replace(exp,fhInfo.Hnext2current);
		}
    }
	
	
	return exp;
}


// Return minimum f-value
int fMin(ProbDescription& Problem,bdd& min,bdd& open) {
	
	int res = 0;
	min = open;
	
	for (int i = Problem.fVars.length - 1; i >= 0; i--)
    {
		bdd b = min & bdd_nithvar(Problem.fVars.current[i]);
		if (b != bddfalse)
		{
			min = b;
			res = 2*res;
		}
		else
		{
			b = min & bdd_ithvar(Problem.fVars.current[i]);
			if (b != bddfalse)
			{
				min = b;
				res = 2*res + 1;
			}
			else
				return -1;
		}
    } 	       
	return res;
}


// compute next successors
bdd imageOpen(BDDAinfo& bddaInfo, FHinfo& fhInfo, bdd& T, bdd h, bdd hm, bdd min, int fMin) {
	
	bdd newOpen = bddfalse;  // set of (f,s) pairs of the open queue of BDDA*
	
	
	// increment newOpen for each partition in turn
	// make a forward step 
	bdd currentNextVarPairs = min & T;
	// add heuristic of current vars and quantify current variables
	bdd nextVarsH = bdd_exist(currentNextVarPairs & h,bddaInfo.currentVars);
	// add heuristic of next vars
	bdd nextVarsHmH = nextVarsH & hm;
	
	// apply formula to get new f and quantify the two h values
	bdd openFrac = bddaInfo.formulaForward[fMin] & nextVarsHmH;
	openFrac = bdd_exist(openFrac,fhInfo.Hcurrent);
	openFrac = bdd_exist(openFrac,fhInfo.Hnext);
	newOpen |= openFrac;      
	
	// change newOpen to current variables (we are searching forward)
	newOpen = bdd_replace(newOpen,bddaInfo.next2current);
	
	return newOpen;
}

//IN
// bddaInfo   : varaible info
// fhInfo     : BDD variable position of f- and h-values
// T          : transition relation
// h          : BDD encoding of heuristic (h current state vars)
// hm         : BDD encoding of heuristic (h next state vars)
// min        : fmin * states  (f, states pair in open with min f value) (both in current vars)
//OUT
// states (current vars) reached backward from min, paired with their f-value (in current vars)
bdd preImageOpen(BDDAinfo& bddaInfo, FHinfo& fhInfo, bdd& T, bdd h, bdd hm, bdd min, int fMin) {
	
	
	bdd newOpen = bddfalse;  // set of (f,s) pairs of the open queue of BDDA*
	
	// change min to next variables (we are searching backwards)
	min = bdd_replace(min,bddaInfo.current2next);
	
	// increment newOpen for each partition in turn
	
		// make a backward step 
			
	bdd currentNextVarPairs = T & min;
	// add heuristic of next vars and quantify next variables
	bdd currentVarsHm = bdd_exist(currentNextVarPairs & hm,bddaInfo.nextVars);
	// add heuristic of current vars
	bdd currentVarsHmH = currentVarsHm & h;
	
	// apply formula to get new f and quantify the two h values
	bdd openFrac = bddaInfo.formulaBackward[fMin] & currentVarsHmH;
	openFrac = bdd_exist(openFrac,fhInfo.Hcurrent);
	openFrac = bdd_exist(openFrac,fhInfo.Hnext);
	newOpen |= openFrac;      
    
	
	return newOpen;
}



//Make initial state bdd
bdd make_init(ProbDescription &Problem){
	cout << "make init..."<<endl;
	bdd b = bddtrue;
	for (int j= 0; j<Problem.init.size(); j++) {
	 b &= int2bdd(Problem.States[j].current, Problem.States[j].length, Problem.init[j]);
	}
	return b;

}

//make goal state bdd
bdd make_goal(ProbDescription &Problem){
	cout << "make goal..."<< endl;
	bdd b = bddtrue;
	for (int j= 0; j<Problem.goal.size(); j++) {
	 b &= int2bdd(Problem.States[j].current, Problem.States[j].length, Problem.goal[j]);
		
	}
	return b;
	
}

// Used to indicate not moving tiles
bdd all_other_idle(ProbDescription &Problem, int src, int dst)
{
	bdd idle = bddtrue;
	int i;
	for(i=0; i<Problem.NumVars; i++)
	{
		if(i != src && i != dst){
			for (int j =0; j<Problem.States[i].length; j++) {
				idle &= bdd_biimp(bdd_ithvar(Problem.States[i].current[j]),bdd_ithvar(Problem.States[i].next[j]));
			}
		}	
	}
	return idle;
}

// make move bdd
bdd make_move(ProbDescription &Problem, int src, int dst, int num)
{
	bdd move = bddtrue;
	move &= int2bdd(Problem.States[src].current, Problem.States[src].length, 0);
	move &= int2bdd(Problem.States[dst].current, Problem.States[dst].length, num);
	move &= int2bdd(Problem.States[src].next, Problem.States[src].length, num);
	move &= int2bdd(Problem.States[dst].next, Problem.States[dst].length, 0);
	move &= all_other_idle(Problem, src, dst);
	return move;
}

// Build transition relation
bdd make_transition_relation(ProbDescription &Problem, int X)
{
	bdd T = bddfalse;
	int blank,dst,num;
	int SIZE = Problem.NumVars;
	for(blank=0;blank<SIZE;blank++)
	{
		//move up
		if(blank>(X-1))
		{
			dst = blank - X;
			for(num = 1; num < SIZE; num++)
			{
				T |= make_move(Problem, blank, dst, num);
			}
		}
		//move down
		if(blank<(X*(X-1)))
		{
			dst = blank + X;
			for(num = 1; num < SIZE; num++)
			{
				T |= make_move(Problem, blank, dst, num);
			}
		}
		//move left
		if(blank%X > 0)
		{
			dst = blank - 1;
			for(num = 1; num < SIZE; num++)
			{
				T |= make_move(Problem, blank, dst, num);
			}
		}
		//move right
		if(blank%X < (X-1))
		{
			dst = blank + 1;
			for(num = 1; num < SIZE; num++)
			{
				T |= make_move(Problem, blank, dst, num);
			}
		}
	}
	return T;
}

// the value of var, if it exist. Otherwise -1. Notice that this function
// non-deterministic will choose one of the values if several are possible  
int bvecEcoding2int(int varLength, int* var,bdd state) {
	
	int res = 0;
	bdd c = state;
	
	for (int i = varLength - 1; i >= 0; i--)
    {
		bdd b = c & bdd_nithvar(var[i]);
		if (b != bddfalse)
		{
			c = b;
			res = 2*res;
		}
		else
		{
			b = c & bdd_ithvar(var[i]);
			if (b != bddfalse)
			{
				c = b;
				res = 2*res + 1;
			}
			else
				return -1;
		}
    } 	       
	return res;
}

// the value of var, if it exist. Otherwise -1. Notice that this function
// non-deterministic will choose one of the values if several are possible  
int bvec2int(int varLength, int* var,bdd state) {
	
	int res = 0;
	bdd c = state;
	
	for (int i = 0;i < varLength; i++)
    {
		bdd b = c & bdd_nithvar(var[i]);
		if (b != bddfalse)
		{
			c = b;
			res = 2*res;
		}
		else
		{
			b = c & bdd_ithvar(var[i]);
			if (b != bddfalse)
			{
				c = b;
				res = 2*res + 1;
			}
			else
				return -1;
		}
    } 	       
	return res;
}

/***************Parse the input file: Read the following info*******************/
/*
 int size;
 string *vars;
 int *Initarr;
 int *Goalarr;
 int Bits;
 int Fmax;
 */

void parser(char *inputfile, int &size, vector<string> &vars, vector<int> &Initarr, vector<int> &Goalarr, int &Bits, int &Fmax)
{

FILE *fp;
if(!(fp = fopen(inputfile,"r")))
{
cout<<"File read error!"<<endl;
return;
}
char s[20];


fscanf(fp, "%s", s);
if(strcmp(s,"VarNum")==0)
{
fscanf(fp, "%*s%d",&size);
}
else if(strcmp(s, "VarNum:")==0)
fscanf(fp, "%d", &size);


fscanf(fp,"%s",s);
if(strcmp(s,"Vars")==0)
fscanf(fp,"%*s");
for(int i=0; i<size; i++)
{
fscanf(fp,"%s",s);
vars.push_back(s);
}



fscanf(fp,"%s",s);
if(strcmp(s,"Init")==0)
fscanf(fp,"%*s");
for(int i=0; i<size;i++)
{
fscanf(fp, "%s", s);
Initarr.push_back(atoi(s));
}


fscanf(fp,"%s",s);
if(strcmp(s,"Goal")==0)
fscanf(fp,"%*s");
for(int i=0;i<size;i++)
{
fscanf(fp,"%s",s);
Goalarr.push_back(atoi(s));
}


fscanf(fp,"%s", s);
if(strcmp(s,"BITS")==0)
fscanf(fp,"%*s");
fscanf(fp,"%d",&Bits);
fscanf(fp,"%s",s);

if(strcmp(s,"FMAX")==0)
fscanf(fp,"%*s");
fscanf(fp,"%d",&Fmax);
fclose(fp);
}


// Print the solution path
void print_solution(ProbDescription &Problem, vector<bdd> sol){
	int X = sqrt(Problem.NumVars);
	cout << "Solution path: "<<endl;
	for (int i = sol.size() -1; i>= 0; i--) {
		
		for (int j =0; j<Problem.NumVars; j++) {
			cout << bvec2int(Problem.States[j].length,Problem.States[j].current,sol[i])<< " ";
			if (j%X == (X-1)) {
				cout << endl;
			}
		}
		cout<< "; "<<endl;
	}
	


}






bdd previousState(BDDAinfo& bddaInfo, FHinfo& fhInfo, bdd& T, bdd h, bdd hm, bdd min, int fMin) {
	
	
	bdd newOpen = bddfalse;  // set of (f,s) pairs of the open queue of BDDA*
	
	// change min to next variables (we are searching backwards)
	min = bdd_replace(min,bddaInfo.current2next);
	
	// make a backward step 
	bdd currentNextVarPairs = T & min;
	// add heuristic of next vars and quantify next variables
	bdd currentVarsHm = bdd_exist(currentNextVarPairs & hm,bddaInfo.nextVars);
	// add heuristic of current vars
	bdd currentVarsHmH = currentVarsHm & h;
	
	// apply formula to get new f and quantify the two h values
	bdd openFrac = bddaInfo.formulaBackward[fMin] & currentVarsHmH;
	openFrac = bdd_exist(openFrac,fhInfo.Hcurrent);
	openFrac = bdd_exist(openFrac,fhInfo.Hnext);
	newOpen |= openFrac;      
	
	
	return newOpen;
}

// Symbolic A* implementation
void BDDAstar(ProbDescription& Problem,FHinfo& fhInfo,BDDAinfo& bddaInfo,
		      bdd heuristic,bdd &Trans,bdd init, bdd goal) {
	
	bdd open;       // search queue of BDDA*
	bdd Succ;    // Successors states f pairs
	bdd min;        // part of search queue with minimum f-value
	bdd heuristicm; // heuristic in next H variables	
	vector<bdd> fringe; // vector of reached states
	vector<int> fringeSize;
	int it = 0;
	double aveFringeSize;
	vector<bdd> solution;
	// initialize open, reach, fringe and heursiticm
	open = bdd_replace(heuristic,fhInfo.Hcurrent2Fcurrent) & init;
	fringe.push_back(open);
	// construct h' 
	heuristicm = bdd_replace(heuristic,fhInfo.Hcurrent2next);
	heuristicm = bdd_replace(heuristicm,bddaInfo.current2next);
	
	
	// main loop
	while (open != bddfalse)
    {      
		// find hmin set
		int fmin = fMin(Problem,min,open);
		
		// check if goal reached
		if ( (min & goal) != bddfalse )
			break; // break out of while loop
		
		it++;
		
		// update fringe
		fringe.push_back(min); 
		
		
		// subtract min from open
		open &= !min;
		
		// abstract f values from min (we know it equals fmin for all states)
		min = bdd_exist(min,fhInfo.Fcurrent);

		// find next states from min 
		Succ = imageOpen(bddaInfo,fhInfo,Trans,heuristic,heuristicm,min, fmin);
		
		
		// update open
		open |= Succ;
		
		//cout << "Step=" << it << " Fmin=" << fmin << " minSize= " << bdd_nodecount(min) << endl;
		fringeSize.push_back(bdd_nodecount(min));
		
    }


	//Solution extraction (if there is one):
	bdd ints = min & goal;
	if (ints != bddfalse) // look for a solution only if there is one....
    {
		bdd currentState = bdd_satone(ints);          
		bdd nextState;
		bdd nextStates; 
		int fringeNo = fringe.size() - 1;
		//push back the goal first
		solution.push_back(goal);
		
		//select one state in intersection
		// Extract solution from backward traversal of images
		while ( (currentState & init) == bddfalse )
		{
			// 1) find the first fringe overlapping with the next states of the current state
			int fVal = bvecEcoding2int(Problem.fVars.length,Problem.fVars.current,currentState);
			currentState = bdd_exist(currentState,fhInfo.Fcurrent);
			nextStates = preImageOpen(bddaInfo,fhInfo,Trans,heuristic,heuristicm,currentState,fVal);
			while (  ( fringe[fringeNo] & nextStates ) ==  bddfalse )
				fringeNo--;
			
			// 2) find one of the states in the overlap
			nextState = bddfalse;
			while (nextState == bddfalse)
			{
				nextState = previousState(bddaInfo,fhInfo,Trans,heuristic,heuristicm,currentState,fVal);
				nextState &= fringe[fringeNo];
			}
			currentState = nextState;
			solution.push_back(nextState);
			fringeNo--; // look at next fringe
		}
    }
	cout<< "Solution length: " <<(solution.size()- 1)<<endl;
	print_solution(Problem,solution);
		
	cout << "Number of Iterations = "<< it<<endl;
	double aveFringe = 0.0;
	for (int i = 0; i < fringeSize.size(); i++)
		aveFringe += fringeSize[i];
	aveFringeSize = aveFringe / double(fringeSize.size()); 
	cout << "Average number of nodes in fringe: " 
	<< aveFringeSize << endl;
    
	
}




void varDef::print() {
	
	cout << "    length : " << length << endl;
	cout << "    current : ";
	for (int i=0; i < length; i++)
		cout << current[i] << " ";
	cout << endl;
	cout << "    next : ";
	for (int i=0; i < length; i++)
		cout << next[i] << " ";
	cout << endl;
}



void ProbDescription::print() {
	cout << "\nProbDescription\n";
	cout << "States:\n\n";
	for (int i = 0; i < States.size(); i++)
    {
		States[i].print();
		cout << "\n\n";
    }
	
	cout << "\nvarNames : ";
	for (int i = 0; i < variablesNames.size(); i++)
		cout << variablesNames[i] << " ";
	cout << "\n\n";
	
	cout << " f encoding\n";
	fVars.print();
	cout << endl;
	cout << " h encoding\n";
	hVars.print();
	
}  

#endif
