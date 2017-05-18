/**
 * Author: Mattia Villani
 * Email : mattiavillani94@gmail.com
 * Github: https://github.com/mattia-villani
 */
#include "genetic.h"
#include "genetic.c"
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<sys/types.h>

#define abs(A) ((A)>=0?(A):-(A))

#ifndef GENERATIONS 
#define GENERATIONS -1
#endif 

#ifndef POP_SIZE
#define POP_SIZE 100
#endif

double fitness_close_to_1994(void* descriptor, void* args, void* chromosome){
	int n = 1994;
	int big =*(int*)descriptor;
	int i;
	int count = 0;
	for ( i=big-1; i>=0; i--, n=n/2 )
		if ( n % 2 == ( (char*)chromosome )[i] )
			count++;
	return (double)count;
}
char termination_score_is_max(void* descriptor, void* args, void* chromosome){
	return fitness_close_to_1994( descriptor, args, chromosome ) == *(int*)descriptor;
}
int main(){
	int l = 16;
	int gens;
	char* result;
	srand(time(0));
	result = simpleBinaryGeneticProblem( l, GENERATIONS, POP_SIZE, NULL, fitness_close_to_1994, termination_score_is_max, &gens );
	if ( result ){
		int i;
		printf("RESULT (after %d generations) : ", gens);
		singleChromosomeDumper_binary(&l, result);
		printf("\n");
	}else {
		printf("RESULT UNFOUND\n");
	}

	return 0;
}