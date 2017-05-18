/**
 * Author: Mattia Villani
 * Email : mattiavillani94@gmail.com
 * Github: https://github.com/mattia-villani
 */
#include "genetic.h"
#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<math.h>

inline void* _realloc( void* old, ssize_t size ){
	void *ret ;
	assert( size );
	ret = realloc(old, size);
	assert(ret);
	return ret;
}

inline void* _malloc ( ssize_t size ){
	void* ret;
	assert( size );
	ret = malloc(size);
	assert( ret );
	return ret;
}


// utility functions
inline double _rand(){
    return ((double)RAND_FUNCTION_APPLIED) / (double)(RAND_MAX + 1.);
}
inline int _intRand(int max){
	return (RAND_FUNCTION_APPLIED%max);
}


int compare( const void*_a, const void*_b ){
	chromosome_t *a = (chromosome_t*)_a;
	chromosome_t *b = (chromosome_t*)_b;
	if ( a->isValid && !b->isValid ) return -1;
	if ( !a->isValid && b->isValid ) return 1;
	if ( b->score - a->score ) return - (a->score - b->score);
	return ((a->chromosome) < (b->chromosome))? 1 : -1;
}

void printChromosomePopulation(void*descriptor, chromosome_t*chromosomes, int size, void (*singleChromosomeDumper)(void*, void*)){
	int i;
	printf("Dump of the population (%d entities)\n", size);
	if ( size > OUTPUT_MAX_SIZE ){ 
		printf("\t\t  [output hidden because lareeger than %d]\n", OUTPUT_MAX_SIZE);
		return;
	}
	assert(singleChromosomeDumper);
	for(i=0;i<size;i++){
		printf("\t%5d) %8s | score = " FLOAT_FORMATTER " | value = ", i, chromosomes[i].isValid?"valid":"INVALID" ,chromosomes[i].score);
		singleChromosomeDumper(descriptor, chromosomes[i].chromosome);
		printf("\n");
	}
}
void printJustGenes(void*descriptor, void**chromosomes, int size, void (*singleChromosomeDumper)(void*, void*)){
	int i;
	printf("Dump of the genes (%d genes)\n", size);
	if ( size > OUTPUT_MAX_SIZE ){ 
		printf("\t\t  [output hidden because lareeger than %d]\n", OUTPUT_MAX_SIZE);
		return;
	}
	assert(singleChromosomeDumper);
	for(i=0;i<size;i++){
		printf("\t\t  %5d) ",i);
		singleChromosomeDumper(descriptor, chromosomes[i]);
		printf("\n");
	}
}


/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						GENERIC SOLVER - START
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/

void* genericGeneticProblem(
	void* descriptor, // whatever may be usefull in the other function, this descriptor is only passed, never used.
	int max_generation, 
	void* populationSizeArgs,
	unsigned int (*populationSize)(void*, void*, int), // descriptor + args + given the generation returns the population size
	void (*populationInitializer)(void*, void*[], int), // descriptor + fills the array (2°arg) of size 3° arg with an initialization
	void (*chromosomeFree)(void*, void*), // descriptor + chromosome, free the memory of the chromosome
	char (*isChromosomeValid)(void*, void*), // descriptor + chromosome returns true iff the chromosome is valid (e.g. not outside domain)
	void* terminationConditionArgs,
	char (*terminationCondition)(void*, void*, void*), // descriptor + args + chromosome returns true iff chromosome satisfy the problem
	void* fitnessArgs,
	double (*fitness)(void*, void*, void*), // descriptor + args + chromosome returns the fitness value
	void* howManyCrossOverArgs,
	unsigned int (*howManyCrossOver)(int, void*), // returns how many crossOver have to be done for the generation given as first arg and then its args
	char canFathersSurvive,
	char canFathersBeReElected,
	unsigned int crossOverArity,
	unsigned int howManyChildrenCrossOverGenerates,
	void* crossOverArgs, 
	void (*crossOver)(void*, void*, void*[], int, void*[], int), // descriptor + args + fathers + his size + array to be filled out of children + its length 
	void* parentSelectionArgs, // argument to be passed to parentSelection function
	void (*parentSelection)(void*, void*, int, chromosome_t[], int, chromosome_t*[], int), // descriptor + arguments + generation + poulation + its size + selected + howManay to select (NOTE, if you what the fathers not to survive, invalidate them here)
	void* mutateArgs,
	void (*mutate)(void*, void*, void*), // descriptor + argument + chromosome : mutates the chromosome
	void* reorderArgs,
	void (*reorder)(void*, void*, void*), // descriptor + argument + chromosome : reorder the chromosome
	void* dumperArgs,
	void (*dumper)(void*, void*,int, void*[], int), // descriptor + args + chromosomes + length. Called each cycle, do what you want: e.g print
	void (*singleChromosomeDumper)(void*, void*), // descriptor + chromosome, it prints out it
	char (*comparer)(void*, void*, void*), // descriptro + chromosome1 + chromosome2 returns 0 iff chromosome1 != chromosome2
	int* generationCreated 		// output
){
	int i;
	int generation = 0;
	unsigned int population_size = 0;
	chromosome_t* chromosomes = NULL;
	int bufferSize = 0;
	void* buffer = NULL;

	int terminationReached = 0;

	int childrenBufferSize = 0;
	void** childrenBuffer = NULL;
	int fathersBufferSize = 0;
	chromosome_t** fathersBuffer = NULL;

	___VERBOSE_OUTPUT___(printf("Genetic: initializing\n"))
	// init the population 
	population_size = populationSize(descriptor, populationSizeArgs, generation);
	___VERBOSE_OUTPUT___(printf("Genetic: got population size of %d\n", population_size))
	buffer = _malloc( population_size * sizeof(void*) );
	populationInitializer(descriptor, (void**)buffer, population_size);
	___VERBOSE_OUTPUT___(printf("Genetic: got population initialized by population initializer\n"))	
	chromosomes = (chromosome_t*)_malloc ( population_size * sizeof(chromosome_t) );
	___VERBOSE_OUTPUT___(printf("Genetic: validating...\n"))	
	for ( i = 0; i<population_size; i++ ){
		___VERBOSE_OUTPUT___(printf("\tElement %5d, ", i))
		chromosomes[i].chromosome = ((void**)buffer)[i];
		___VERBOSE_OUTPUT___(singleChromosomeDumper(descriptor, chromosomes[i].chromosome))
		chromosomes[i].isValid = !isChromosomeValid || isChromosomeValid(descriptor, chromosomes[i].chromosome);
		___VERBOSE_OUTPUT___(printf(", valid=%d, ",(int)chromosomes[i].isValid));
		assert(chromosomes[i].isValid); // i expect that the init is valid
		chromosomes[i].score = fitness(descriptor, fitnessArgs, chromosomes[i].chromosome);
		___VERBOSE_OUTPUT___(printf("score=" FLOAT_FORMATTER "\n",chromosomes[i].score));
	}	
	___VERBOSE_OUTPUT___(printf("Genetic: chromosomes evaluated with fitness\n"))	
	qsort(chromosomes,population_size,sizeof(chromosome_t),compare);
	___VERBOSE_OUTPUT___(printf("Genetic: chromosomes sorted by fitness\n"))	
	free(buffer);
	buffer = NULL;
	// end init : array chromosomes is filled with initial population and ordered by score

	// loop
	// invariant : chromosomes is sorted by score
	___VERBOSE_OUTPUT___(printf("Genetic: starting event loop\n"))	
	for ( ; (generation < max_generation || max_generation < 0) && !terminationReached && population_size>0 ; generation++ ){
		int newPopSize;
		int howManyCrossOverThisGen = howManyCrossOver(generation,howManyCrossOverArgs);

		terminationReached = ( terminationCondition && terminationCondition(descriptor, terminationConditionArgs, chromosomes->chromosome) );
		if ( terminationReached ) continue;

		// resizing the children buffer size in case it changed
		if ( howManyCrossOverThisGen*howManyChildrenCrossOverGenerates != childrenBufferSize )
			childrenBuffer = (void**)_malloc(sizeof(void*)*( childrenBufferSize = howManyChildrenCrossOverGenerates * howManyCrossOverThisGen ));
		// resizing the all fathers buffer size in case it changed
		if ( howManyCrossOverThisGen*crossOverArity != fathersBufferSize )
			fathersBuffer = (chromosome_t**)_malloc(sizeof(chromosome_t*)*( fathersBufferSize = crossOverArity * howManyCrossOverThisGen ));

		___VERBOSE_OUTPUT___(printf("Genetic: loop infos\n"))
		___VERBOSE_OUTPUT___(printf("\tgeneration: %d\n",generation))
		___VERBOSE_OUTPUT___(printf("\tTermination reached: %s\n",terminationReached?"true":"false"))
		___VERBOSE_OUTPUT___(printf("\tcrossOver to perform: %d\n",howManyCrossOverThisGen))
		___VERBOSE_OUTPUT___(printf("\tchildren for each crossover:   %d => children buffer: %d\n", howManyChildrenCrossOverGenerates, childrenBufferSize))
		___VERBOSE_OUTPUT___(printf("\tfathers f.e. crossover(arity): %d => fathers buffer : %d\n", crossOverArity, fathersBufferSize))
		___VERBOSE_OUTPUT___(printf("\t");printChromosomePopulation(descriptor, chromosomes, population_size, singleChromosomeDumper))		
		
		// performing crossovers
		for ( i=0; i<howManyCrossOverThisGen; i++ ){
			int j;
			___VERBOSE_OUTPUT___(printf("\tCrossOvering loop %d\n", i))
			// fathers selections
			parentSelection(descriptor, parentSelectionArgs, generation, chromosomes, population_size, fathersBuffer+i*crossOverArity, crossOverArity);
			___VERBOSE_OUTPUT___(
				printf("\t\tFathersSelected: \n");
				for ( j=0; j<crossOverArity; j++ ){
					int dif = -(int)( chromosomes - (fathersBuffer+i*crossOverArity)[j] );
					printf("\t\t\t%2d) referitng to %4d-th (%7s) chromosome : ",j, dif, (fathersBuffer+i*crossOverArity)[j]?"valid":"INVALID" );
					singleChromosomeDumper(descriptor, (fathersBuffer+i*crossOverArity)[j]->chromosome );
					printf("\n");
				})
			// performing the crossover 
			// but first prepare fathers 
			if ( bufferSize < crossOverArity ){	
				if ( buffer ) buffer = (void**)_realloc(buffer, sizeof(void*) * (bufferSize = crossOverArity) );
				else buffer = (void**)_malloc(sizeof(void*) * (bufferSize = crossOverArity) );
			}
			// copy fathers in the buffer
			for(j=0; j<crossOverArity; j++ )
				((void**)buffer)[j] = fathersBuffer[j]->chromosome;
			// actuall crossover
			___VERBOSE_OUTPUT___(printf("\t\tChildrenProduced: \n"))
			crossOver(descriptor, crossOverArgs, (void**)buffer, crossOverArity, (void**)childrenBuffer+i*howManyChildrenCrossOverGenerates, howManyChildrenCrossOverGenerates );
			___VERBOSE_OUTPUT___(
				for ( j=0; j<howManyChildrenCrossOverGenerates; j++ ){
					printf("\t\t\t%2d) ",j);
					singleChromosomeDumper(descriptor, (childrenBuffer+i*howManyChildrenCrossOverGenerates)[j]);
					printf("\n");
				})
			// mutating and reordening
			for ( j = i*howManyChildrenCrossOverGenerates ; j < (i+1)*howManyChildrenCrossOverGenerates ; j++ ){
				if ( mutate ){ 
					___VERBOSE_OUTPUT___(printf("\t\t\tMutating -> ");singleChromosomeDumper(descriptor,childrenBuffer[j]))
					mutate( descriptor, mutateArgs, childrenBuffer[j] );
					___VERBOSE_OUTPUT___(printf(" into -> ");singleChromosomeDumper(descriptor,childrenBuffer[j]); printf("\n")) 
				}
				if ( reorder ){ 
					___VERBOSE_OUTPUT___(printf("\t\t\tReordering -> ");singleChromosomeDumper(descriptor,childrenBuffer[j])) 
					reorder(descriptor, reorderArgs, (void**)childrenBuffer+j );
					___VERBOSE_OUTPUT___(printf(" into -> ");singleChromosomeDumper(descriptor,childrenBuffer[j]); printf("\n")) 
				}
			}
			// invalidating fathers if they can't be re-elected
			___VERBOSE_OUTPUT___(printf("\t\tInvalidating or restoring fathers because of re-election\n"))
			for ( j=i*crossOverArity; j<(i+1)*crossOverArity; j++ )
				if ( fathersBuffer[j]->isValid != canFathersBeReElected ){
					___VERBOSE_OUTPUT___({
						int dif = -(int)( chromosomes - fathersBuffer[j] );
						printf("\t\t\t%15s of the %3d-th father : ", canFathersBeReElected ? "Enabling" : "Invalidating",dif);
						singleChromosomeDumper(descriptor, fathersBuffer[j]->chromosome);
						printf("\n");
					})
					fathersBuffer[j]->isValid = canFathersBeReElected;
				}
		}
		//___VERBOSE_OUTPUT___(printf("\tPost procration : ");printChromosomePopulation(descriptor, chromosomes, population_size, singleChromosomeDumper))		

		if ( !canFathersBeReElected ) // they were invalidated before, re-enabeling
			for(i=0;i<fathersBufferSize;i++)
				fathersBuffer[i]->isValid = 1;
		if ( !canFathersSurvive ){ // invalidate fathers
			___VERBOSE_OUTPUT___(printf("\tInvalidating fathers because they can't survive\n"))
			for(i=0;i<fathersBufferSize;i++)
				if ( fathersBuffer[i] && fathersBuffer[i]->isValid && fathersBuffer[i]->chromosome ) {
					___VERBOSE_OUTPUT___({
							int dif = -(int)( chromosomes - fathersBuffer[i] );
							printf("\t\tinvalidation of the %3d-th father : ", dif);
							singleChromosomeDumper(descriptor, fathersBuffer[i]->chromosome);
							printf(".\n");
						})
					fathersBuffer[i]->isValid = 0;
					if ( fathersBuffer[i]->chromosome )
						chromosomeFree( descriptor, fathersBuffer[i]->chromosome );
					fathersBuffer[i]->chromosome = NULL;
				}
		}
		
		// some fathers may have been invalidated, i will reorder the population so that they goes at the end
		qsort(chromosomes, population_size, sizeof(chromosome_t), compare);
		___VERBOSE_OUTPUT___(printf("\tResorting all the population : ");printChromosomePopulation(descriptor, chromosomes, population_size, singleChromosomeDumper))		
		// selecting the survivals 
		// seeing how big will be the new generation
		newPopSize = populationSize(descriptor, populationSizeArgs, generation+1);
		___VERBOSE_OUTPUT___(printf("\tNew population size : %d. (the old one was %d)\n", newPopSize, population_size))				
		if ( newPopSize > population_size ){ // new spots needed
			chromosomes = (chromosome_t*)_realloc(chromosomes, newPopSize);
			for ( i = population_size; i<newPopSize; i++ )
				chromosomes[i].isValid = 0;
			population_size = newPopSize;
		}else if ( newPopSize < population_size){ // the lasts are out of the game
			for ( i = newPopSize; i<population_size; i++ )
				if ( chromosomes[i].chromosome ){
					chromosomes[i].isValid = 0;
					chromosomeFree( descriptor, chromosomes[i].chromosome );
				}
			chromosomes = (chromosome_t*)_realloc(chromosomes, newPopSize );
			population_size = newPopSize;
		}
		___VERBOSE_OUTPUT___(printf("\tIntegrating children in the population...\n"))						
		// the new population size is ready!
		// now childrenBuffer is full of childrens
		// let's selection insert them
		for ( i=0; i<childrenBufferSize; i++){
			void *child = childrenBuffer[i];
			___VERBOSE_OUTPUT___(printf("\t\tChild %3d : ", i); singleChromosomeDumper(descriptor, child));				
			if ( !isChromosomeValid || isChromosomeValid(descriptor, child) ){
				int firstMove = 0;
				int j;
				double score = fitness(descriptor, fitnessArgs, child);
				___VERBOSE_OUTPUT___(printf(", valid, score( " FLOAT_FORMATTER " )", score));
				for( j=0; j<population_size && child; j++)
					if ( chromosomes[j].isValid == 0 ){
						chromosomes[j].isValid = 1;
						chromosomes[j].score = score;
						chromosomes[j].chromosome = child;
						child = NULL;
						___VERBOSE_OUTPUT___(printf(", sitting at %3d", j));
					}else if ( score > chromosomes[j].score ){
						// swapping
						double t = score;
						void *el = child;
						score = chromosomes[j].score;
						child = chromosomes[j].chromosome;
						chromosomes[j].score = t;
						chromosomes[j].chromosome = el;
						___VERBOSE_OUTPUT___(if ( !firstMove ) printf(", insert at  %3d", j));
						firstMove = 1;
					}else if ( comparer && score == chromosomes[j].score ){
						int equals = comparer( descriptor, child, chromosomes[j].chromosome );
						if ( equals ){
							___VERBOSE_OUTPUT___( printf(" ALREADY PRESENT "));
							break;
						}
					}
				// the last is out of population
				if ( child ){ 
					___VERBOSE_OUTPUT___({
						printf(", killing ");
						singleChromosomeDumper(descriptor, child);
						printf(" of score " FLOAT_FORMATTER, score);
					})
					if ( child )
						chromosomeFree(descriptor, child);
				}
			}else { ___VERBOSE_OUTPUT___(printf("INVALID => DROPPED.")); }
			___VERBOSE_OUTPUT___(printf("\n"));
		}
		// reduce pop size if there are invalid at the end ( in case the array was enlarged )
		___VERBOSE_OUTPUT___( printf("\tCleaning population\n") )
		for ( ; population_size>0 && !chromosomes[population_size-1].isValid ; population_size-- ){
			___VERBOSE_OUTPUT___( printf("\t\tKilling element at position %d\n", population_size-1 ) );
		}
		// chromosomes is ordered based on score and it has the sons joint if they have a better score
		// if requested, print 
		if ( dumper ){ 
			int j;
			// but first prepare fathers 
			if ( bufferSize < population_size ){	
				if ( buffer ) buffer = (void**)_realloc(buffer, sizeof(void*) * (bufferSize = population_size) );
				else buffer = (void**)_malloc(sizeof(void*) * (bufferSize = population_size) );
			}
			// copy fathers in the buffer
			for(j=0; j<population_size; j++ )
				((void**)buffer)[j] = chromosomes[j].chromosome;
			dumper(descriptor, dumperArgs, generation, (void**)buffer, population_size );
		}
	}

	if ( generationCreated ) *generationCreated = generation;
	// returns the best chromosome: the first of chromosomes 
	if (fathersBuffer ) free(fathersBuffer);
	if ( childrenBuffer) free(childrenBuffer);
	if ( buffer ) free(buffer); 
	for ( i=1; i<population_size; i++ )
		if ( chromosomes[i].chromosome )
			chromosomeFree( descriptor, chromosomes[i].chromosome );
	buffer = chromosomes ? chromosomes->chromosome : NULL;
	free(chromosomes);
	___VERBOSE_OUTPUT___(if ( buffer ) { printf("SOLUTION got at generation (%d) : ",generation); singleChromosomeDumper(descriptor, buffer); printf("\n"); })
	return buffer;
}
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						GENERIC SOLVER - END
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/

/**
 * POPULATION SIZE UTILITY: CONSTANT VALUE
 */
// argument args is the pointer to the returned value
unsigned int populationSize_constant( void* descriptor , void* args, int generation ){
	assert ( args );
	return *(int*) args;
}

/**
 * HOW MANY CROSS OVERS: CONSTANT
 */
unsigned int howManyCrossOver_constant(int generation, void* args){
	return args ? *(int*) args: 1;
}


/**
 ** PARENT SELECTION ALGORITHMS
 **
 **/

/** PROPORCIONA SELECTION **/
void parentSelection_roulette( void*descriptor, void* arg, int generation,
	 chromosome_t pop[], int popSize,chromosome_t* fathers[], int fatherSize){
	int f, i;
	double scores = 0;
	int reelection = arg? *(char*)arg : 0;
	for ( i=0; i<popSize; i++ ) 
		if ( pop[i].isValid ) 
			scores += pop[i].score;
	for ( i=0; i<fatherSize; i++ ){
		double searched = scores*_rand();
		int f=0;
		___VERBOSE_OUTPUT___(printf("\t\t\tRouletting between 0 and " FLOAT_FORMATTER ". Got " FLOAT_FORMATTER "\n", scores, searched));
		while( searched > 0 && f<popSize ){
			int j;
			int alreadyUsed = ! pop[i].isValid ;
			for ( j=0; j<i && ! alreadyUsed ; j++ )
				alreadyUsed = fathers[j] == pop+f; 
			if ( ! alreadyUsed )
				searched-=pop[f].score;
			if ( searched > 0 )
				f++;
		}
		if ( f>=popSize )
			for ( f=0; f<popSize && !pop[f].isValid; f++ );
		f = (f<popSize ? f : popSize-1);
		fathers[i] = pop + f;
		if ( !reelection ){
			pop[f].isValid = 0;
			scores = scores - pop[f].score;
		}
	}		
}
/** TOURNAMENT SELECTION **/
/*
 * Wants arg as an int array of three positions. The first is the presure, the second is the determinism
 *	and the thirth indicates if the father repetition should be allowed
 */
void parentSelection_tournament( void*descriptor, void* arg, int generation,
	 chromosome_t pop[], int popSize, chromosome_t* fathers[], int fatherSize){
	int i;
	int presure = arg?*(int*)arg:2;
	int deterministic = arg?((int*)arg)[1]:1;
	int repetitionAllowed = arg?(((int*)arg)[2]):0;
	int countValid = 0;
	assert( presure > 0);
	for ( i=0; i<popSize; i++ ) countValid += (int)(pop[i].isValid);
	for ( i=0; i<fatherSize; i++ ){
		int minIndex = -1; 
		int j;
		for ( j=0; j<presure; j++ ){
			int r = _intRand(countValid);
			int k;
			for (k=0;k<popSize && k<r; k++)
				if ( ! pop[k].isValid ) r++;
			if ( deterministic ){
				if ( minIndex==-1 || r < minIndex )
					minIndex = r;
			}else {
				if ( minIndex==-1 || _rand()<1./(double)presure )
					minIndex = r;
			}
			if ( minIndex<0 || minIndex>popSize )
				for ( minIndex=0; minIndex<popSize && !pop[minIndex].isValid; minIndex++ );
			assert ( minIndex>=0 && minIndex<popSize );
		}
		fathers[i] = pop + minIndex; 
		if ( ! repetitionAllowed ){
			pop[minIndex].isValid = 0;
			countValid --;
		}
	}		
 }


/******************************************************************************
 ******************************************************************************
 * 							BINARY RAPRESENTATION
 ******************************************************************************
 ******************************************************************************/
void populationInitializer_binary_random(void* descriptor, void* toFill[] , int size){
	// descriptor is int* and contains the length of the binary array
	int big = *(int*)descriptor;
	int i;
	for (i=0;i<size;i++){
		int j;
		char* chr = (char*)_malloc(sizeof(char)*big);
		for ( j=0;j<big;j++ )
			chr[j] = _rand()>0.5;
		toFill[i] = chr;
	}
}

void chromosomeFree_binary(void* descriptor, void* chromosome){
	if (chromosome)
		free(chromosome);
}

// set arity y child as 2 -> 2 
void crossOver_binary_onePoint(void* descriptor, void* args, void* fathers[], int fsize, void* childs[], int csize){
	// descriptor is int* and contains the length of the binary array
	int big = *(int*)descriptor;
	char* ch1 = (char*)_malloc(sizeof(char)*big);
	char* ch2 = (char*)_malloc(sizeof(char)*big);
	int split;
	int i;
	do split=_intRand(big); while ( big>2 && ( !split || split>=big ) );
	assert ( csize == 2 && fsize == 2 );
	___VERBOSE_OUTPUT___(printf("\t\t\t\t(one point cross over at index %d (over %d bits)) \n", split, big))
	for (i=0; i<big; i++){
		if ( split == i ){
			char*t = ch1;
			ch1 = ch2;
			ch2 = t;
		}
		ch1[i] = ((char**)fathers)[0][i];
		ch2[i] = ((char**)fathers)[1][i];
	}
	childs[0] = ch1;
	childs[1] = ch2;
}

//mutation
// args has to be the probability 
void mutation_binary_canonical( void* descriptor, void*args, void* chromosome ){
	// descriptor is int* and contains the length of the binary array
	int big = *(int*)descriptor;
	double prob = args?*(double*)args:1./1000000;
	int i;
	
	for(i=0;i<big;i++)
		if(_rand()<prob)
			((char*)chromosome)[i] = ! ((char*)chromosome)[i];
}	

/*
 * It waits as args an array of two pos: the fitness args and the fitness function
 */
void singleChromosomeDumper_binary(void* descriptor, void* chromosome){
	int big = *(int*)descriptor;
	int j=0;
	long val=0;
	int do_val = big < sizeof(long)*8-1;
	for(j=0;j<big;j++){
		printf("%d", (int)(((char*)(chromosome))[j]));
		if ( do_val )
			val = val*2 + (int)(((char*)(chromosome))[j]);
	}
	if ( do_val ) printf("(%10ld)", val);
}
void dumper_print_binary(void* descriptor, void* args, int generation, void* pop [] , int size ){
	int i,j;
	void* fitArgs = (((void**)args)[0]);
	double (*fitness)(void*, void*, void*) = (double (*)(void*, void*, void*))(((void**)args)[1]);
	printf("Generation %d\n",generation);
	printf("\t index     score        value \n");
	for (i=0;i<size;i++){
		printf("\t%7d|" FLOAT_FORMATTER "| ", i, fitness(descriptor, fitArgs, pop[i]));
		singleChromosomeDumper_binary(descriptor, pop[i]);
		printf("\n");
	}
}
char comparer_binary(void* descriptor, void* ch1, void* ch2){
	int big = *(int*) descriptor;
	int i;
	for (i=0;i<big;i++)
		if ( ((char*)ch1)[i] != ((char*)ch2)[i] ) return 0;
	return 1; 
}


char* simpleBinaryGeneticProblem(int byteArraySize, int max_gen, int pop_size, void* fitnessArgs, double (*fitness)(void*, void*, void*), char (*terminationCondition)(void*, void*, void*), int *gen){
	void* dumper_args[] = {fitnessArgs, (void*)fitness};
	int howManyCross = pop_size>4 ?pop_size/4: 1 ;
	char var = 0;
	assert(pop_size >= 2);
	return (char*)genericGeneticProblem(
		&byteArraySize, 								// void* descriptor, // whatever may be usefull in the other function, this descriptor is only passed, never used.
		max_gen, 										// int max_generation,
		&pop_size, 				 						// void* populationSizeArgs,
		populationSize_constant, 						// unsigned int (populationSize*)(void*, int), // descriptor + given the generation returns the population size
		populationInitializer_binary_random, 			// void (populationInitializer*)(void*, void*[], int), // descriptor + fills the array (2°arg) of size 3° arg with an initialization
		chromosomeFree_binary, 							// void (chromosomeFree*)(void*, void*), // descriptor + chromosome, free the memory of the chromosome
		NULL,					 						// char (isChromosomeValid*)(void*, void*), // descriptor + chromosome returns true iff the chromosome is valid (e.g. not outside domain)
		NULL,											// void* terminationConditionArgs
		terminationCondition,     						// char (terminationCondition*)(void*, void*, void*), // descriptor + args + chromosome returns true iff chromosome satisfy the problem
		fitnessArgs, 									// void* fitnessArgs,
		fitness, 										// double (fitness*)(void*, void*, void*), // descriptor + args + chromosome returns the fitness value
		&howManyCross, 									// void* howManyCrossOverArgs,
		howManyCrossOver_constant,						// unsigned int (howManyCrossOver*)(int), // returns how many crossOver have to be done for the generation given as first arg
		1,												// char canFathersSurvive,
		var, 												// char canFathersBeReElected,
		2, 												// unsigned int crossOverArity,
		2, 												// unsigned int howManyChildrenCrossOverGenerates,
		NULL, 											// void* crossOverArgs,						 
		crossOver_binary_onePoint,						// void (crossOver*)(void*, void*, void*[], int, void*[], int), // descriptor + args + fathers + his size + array to be filled out of children + its length 
		&var, 											// void* parentSelectionArgs, // argument to be passed to parentSelection function
		parentSelection_roulette, 						// void (parentSelection*)(void*, void*, int, chromosome_t[], int, chromosome_t[], int), // descriptor + arguments + generation + poulation + its size + selected + howManay to select 
		NULL, 											// void* mutateArgs,
		mutation_binary_canonical,						// void (mutate*)(void*, void*, void*), // descriptor + argument + chromosome : mutates the chromosome
		NULL, 											// void* reorderArgs,
		NULL, 											// void (reorder*)(void*, void*, void*), // descriptor + argument + chromosome : reorder the chromosome
		NULL,//dumper_args,								// void* dumperArgs,
		NULL,//dumper_print_binary, 					// void (dumper*)(void*, void*,int, void*[], int) // descriptor + args + chromosomes + length. Called each cycle, do what you want: e.g print
		singleChromosomeDumper_binary,					// void singleChromosomeDumper(void*, void*) // descriptor + chromosome, it prints out it 
		comparer_binary,								// char (*comparer)(void*, void*, void*)
		gen
	);
}
























/******************************************************************************
 ******************************************************************************
 * 							double RAPRESENTATION
 ******************************************************************************
 ******************************************************************************/
void populationInitializer_double_random(void* descriptor, void* toFill[] , int size){
	// descriptor is int* and contains the length of the binary array
	int i;
	double min = *(double*)descriptor;
	double max = *((double*)descriptor+1);
	for (i=0;i<size;i++){
		toFill[i] = _malloc(sizeof(double));
		*(double*)(toFill[i]) = min + (max-min)*(_rand());
	}
}

void chromosomeFree_double(void* descriptor, void* chromosome){
	if (chromosome)
		free(chromosome);
}

// set arity y child as 2 -> 1 
void crossOver_double_intermediate(void* descriptor, void* args, void* fathers[], int fsize, void* childs[], int csize){
	double alpha = _rand();
	assert(fsize == 2 && csize == 1);
	childs[0] = _malloc(sizeof(double));
	*(double*)(childs[0]) = alpha**(double*)(fathers[0]) + (1-alpha)**(double*)(fathers[1]);
}

//mutation
// args has to be the probability 
void mutation_double_non_uniform( void* descriptor, void*args, void* chromosome ){
	double min = *(double*)descriptor;
	double max = *(((double*)descriptor)+1);
	double prob = *(double*)args;
	double b = *(((double*)args)+1);
	int T = (int) *(((double*)args)+2);
	int t = (int) *(((double*)args)+3);
	double vk = *(double*)chromosome;

	double vk1 = vk + (max-vk)*( 1 - pow( _rand(), pow(1-t/T,b) ) );
	double vk2 = vk + (vk-min)*( 1 - pow( _rand(), pow(1-t/T,b) ) );

	if ( _rand() <= prob )	
		*(double*)chromosome = _rand()<0.5 ? vk1 : vk2;
}	

/*
 * It waits as args an array of two pos: the fitness args and the fitness function
 */
void singleChromosomeDumper_double(void* descriptor, void* chromosome){
	printf(FLOAT_FORMATTER, *(double*)chromosome);
}

char comparer_double(void* descriptor, void* ch1, void* ch2){
	if ( *(double*) ch1 == 0 || *(double*) ch2 == 0 )
		return fabs(*(double*)ch1) < DOUBLE_TOLLERANCE && fabs(*(double*)ch2) < DOUBLE_TOLLERANCE;
	return fabs(*(double*) ch1 / *(double*) ch2)<DOUBLE_TOLLERANCE;
}

void dumperIncreaseCounter(void* descriptor, void* args,int generation, void* chromosomes[], int length) {
	*(((double*)args)+3)+=1;
}

double simpleDoubleGeneticProblem(double min, double max, int max_gen, int pop_size, void* fitnessArgs, double (*fitness)(void*, void*, void*), char (*terminationCondition)(void*, void*, void*), int *gen){
	double descriptor[]={min, max};
	int howManyCross = pop_size>4 ?pop_size/4: 1 ;
	double mutationArgs[]={1./1000., 3., (double)max_gen, 1};
	double* val;
	assert(max_gen > 0);
	val = (double*)genericGeneticProblem(
		descriptor,	 								// void* descriptor, // whatever may be usefull in the other function, this descriptor is only passed, never used.
		max_gen, 										// int max_generation,
		&pop_size, 				 						// void* populationSizeArgs,
		populationSize_constant, 						// unsigned int (populationSize*)(void*, int), // descriptor + given the generation returns the population size
		populationInitializer_double_random, 			// void (populationInitializer*)(void*, void*[], int), // descriptor + fills the array (2°arg) of size 3° arg with an initialization
		chromosomeFree_double, 							// void (chromosomeFree*)(void*, void*), // descriptor + chromosome, free the memory of the chromosome
		NULL,					 						// char (isChromosomeValid*)(void*, void*), // descriptor + chromosome returns true iff the chromosome is valid (e.g. not outside domain)
		NULL,											// void* terminationConditionArgs
		terminationCondition,     						// char (terminationCondition*)(void*, void*, void*), // descriptor + args + chromosome returns true iff chromosome satisfy the problem
		fitnessArgs, 									// void* fitnessArgs,
		fitness, 										// double (fitness*)(void*, void*, void*), // descriptor + args + chromosome returns the fitness value
		&howManyCross, 									// void* howManyCrossOverArgs,
		howManyCrossOver_constant,						// unsigned int (howManyCrossOver*)(int), // returns how many crossOver have to be done for the generation given as first arg
		1,												// char canFathersSurvive,
		0, 												// char canFathersBeReElected,
		2, 												// unsigned int crossOverArity,
		1, 												// unsigned int howManyChildrenCrossOverGenerates,
		NULL, 											// void* crossOverArgs,						 
		crossOver_double_intermediate,					// void (crossOver*)(void*, void*, void*[], int, void*[], int), // descriptor + args + fathers + his size + array to be filled out of children + its length 
		NULL, 											// void* parentSelectionArgs, // argument to be passed to parentSelection function
		parentSelection_tournament, 					// void (parentSelection*)(void*, void*, int, chromosome_t[], int, chromosome_t[], int), // descriptor + arguments + generation + poulation + its size + selected + howManay to select 
		mutationArgs,									// void* mutateArgs,
		mutation_double_non_uniform,					// void (mutate*)(void*, void*, void*), // descriptor + argument + chromosome : mutates the chromosome
		NULL, 											// void* reorderArgs,
		NULL, 											// void (reorder*)(void*, void*, void*), // descriptor + argument + chromosome : reorder the chromosome
		mutationArgs,//dumper_args,								// void* dumperArgs,
		dumperIncreaseCounter,//dumper_print_binary, 					// void (dumper*)(void*, void*,int, void*[], int) // descriptor + args + chromosomes + length. Called each cycle, do what you want: e.g print
		singleChromosomeDumper_double,					// void singleChromosomeDumper(void*, void*) // descriptor + chromosome, it prints out it 
		comparer_double,								// char (*comparer)(void*, void*, void*)
		gen
	);
	if ( val ){
		double v = *val;
		free(val);
		return v;
	}else return NAN;
}

