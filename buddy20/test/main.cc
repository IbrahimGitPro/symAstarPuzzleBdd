/****************************************************** 
 * File  : main.cc   
 * Desc. : test 1 of action split implementation 
 * Author: Rune M. Jensen
 * Date  : 8/20/01
 ******************************************************/

#include <stream.h>
#include <string.h>
#include <stdio.h> 
#include <stdlib.h>
#include <bdd.h>

int main(void)  {


  // initialize bdd package
  bdd_init(1000,100);
  bdd_setvarnum(6);

  bdd x,xp,y,yp,z,zp,A,G,Ap1,Ap2,Am1,A0,B,B0;

  x  = bdd_ithvar(0);
  xp = bdd_ithvar(1); 
  y  = bdd_ithvar(2);
  yp = bdd_ithvar(3);
  z  = bdd_ithvar(4);
  zp = bdd_ithvar(5);
  
  A  = !x & (yp >> ! y) & (yp << !y) & (zp >>  z) & (zp << z);
  A |= x & !y & !z & yp & zp;

  B = (xp >> !x) & (xp << !x);
  
  G = y & z;

  cout << "A: " << A << endl;

  Ap1 = bdd_split(A,G,1);

  cout << "A+1: " << Ap1 << endl;

  Ap2 = bdd_split(A,G,2);

  cout << "A+2: " << Ap2 << endl;
 
  Am1 = bdd_split(A,G,-1);

  cout << "A-1: " << Am1 << endl;

  A0 = bdd_split(A,G,0);

  cout << "A0: " << A0 << endl;

  cout << "B: " << B << endl;

  B0 = bdd_split(B,G,0);

  cout << "B0: " << B0 << endl;

  

 bdd_done();
}

 


 





