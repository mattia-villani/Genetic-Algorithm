/**
 * Author: Mattia Villani
 * Email : mattiavillani94@gmail.com
 * Github: https://github.com/mattia-villani
 */
#include "genetic.h"
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<sys/types.h>

#define abs(A) ((A)>=0?(A):-(A))

#ifndef GENERATIONS 
#define GENERATIONS 200
#endif 

#ifndef POP_SIZE
#define POP_SIZE 10
#endif

double function( double x ){
	return pow(x,3)+27;
}

double fitness(void* descriptor, void* args, void* chromosome){
	double v = 1. / (0.01 + function( *(double*) chromosome ));
	return v*v;
}
char terminationCondition(void* descriptor, void* args, void* chromosome){
	return fabs( function( *(double*) chromosome ) ) < DOUBLE_TOLLERANCE*10;
}
int main(){
	int l = 16;
	int gens;
	double result;
	srand(time(0));

	result = simpleDoubleGeneticProblem(-50, 50, GENERATIONS, POP_SIZE, NULL, fitness, terminationCondition, &gens);

	if ( result != NAN ){
		int i;
		printf("RESULT (after %d generations) : " FLOAT_FORMATTER "\n", gens, result);
	}else {
		printf("RESULT UNFOUND\n");
	}

	return 0;
}
