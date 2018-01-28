/****************************************************** 
 * File  : main.cc   
 * Desc. : test 1 of SubBDD implementation 
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
  bdd_setvarnum(3);

  bdd x,y,z,q,g;

  x = bdd_ithvar(0);
  y = bdd_ithvar(1);
  z = bdd_ithvar(2);

  q = y;
  g = !y & x;

  int lower,upper;

  bdd_findrange(&lower,&upper,q,g);


  cout << lower << " : " << lower << endl;
  cout << upper << " : " << upper << endl;

 bdd_done();
}










