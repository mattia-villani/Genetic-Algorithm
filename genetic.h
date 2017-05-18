/**
 * Author: Mattia Villani
 * Email : mattiavillani94@gmail.com
 * Github: https://github.com/mattia-villani
 */
#ifndef GENETIC
#define GENETIC

/* defines how the score is formatted when printted */
#define FLOAT_FORMATTER "%10.5f"
/* defines the max population size to be printed ( if bigger, it won't be printted ) */
#define OUTPUT_MAX_SIZE 30
/* application of the rand used. if thread safeness is wanted, modify this */
#define RAND_FUNCTION_APPLIED (rand())
/* precision. two doubles with a ratio difference smaller than this are considered iguals */
#define DOUBLE_TOLLERANCE 0.0000001
/*
 * If compiled with VERBOSE defined, lot of debug output will be provided 
 */
#ifdef VERBOSE
#define ___VERBOSE_OUTPUT___(A) ({A;});
#else
#define ___VERBOSE_OUTPUT___(A) ({});
#endif

/*
 * 	CHROMOSOMEs are the basic entities manipulated and its rapresentation is transparent. 
 * 	The idea of this genetic algorithm implementation is to be as generic as possible
 * for this reason, no explicit rapresentation of a chromosome is provide, but some
 * frequent rapresentation will be provieded. For this reason the chromosomes will be 
 * rapresented as void*. Anyway, since often the rapresentation will need a descriptor, 
 * such as the dimension or a domain, each time a chromosome will be "read", it will be
 * passed his descriptor to the function using it as well. This allows an easier rapresentation
 * of some chromosomes. For instance, a chromosome may be an array and his size passed as 
 * descriptor. For more complicated datastructior anyway, it will be necessary to create a custom
 * struct containing the data, at this point the data structur may contain the chromosome and 
 * be considered itself a chromosome, or have a separetad descriptor.
 * For example, for the Array example, if the size is constant, it is sufficient to pass the 
 * size as descriptor. But if the size is dynamic and changes between a chromosome and another,
 * then the chromosome will have to be a structure containing its size and the actual array.
 */

/*
 * 	This is the basic data structure used (actually the only one).
 * It will contain 3 fields:
 *		chromosome : this is the pointer to the actual chromosome, whatever it is.
 *		score	   : this is the value that the fitness function has given to the chromosome
 *		isValid	   : this is a boolean rapresenting the validity of the chromosome. 
 *					Exemples of invalid chromosomes are the out domain ones 
 *					or the dying-due-to-selection ones.
 */
typedef struct{
	double score;
	char isValid;
	void* chromosome;	
} chromosome_t;



/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						GENERIC GENETIC PROBLEM SOLVER
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
/*
 * 	As highlighted before, the idea is to make this implementation as generic as possible
 * For this reason, the behaviour is wholly determinated by the parameters passed which will allow a huge 
 * customization. 
 *
 *	The algorithm is thread safe! Anyway,some of the predefined function, provided later,uses the rand
 * function. This makes the whole application thread unsafe, but a custom implementation may be provided 
 * as argument to the genericGeneticProblem, using rand_r. 
 * In order to keep the application as thread safe as possible,no global or static variables were used. 
 * Each time a function needs arguments unknown to the genericGeneticProblem, 
 * they will be passed as args to genericGeneticProblem and passed to the function interested.
 * For instance: the probability for a specific crossOver function may be passed as crossOverArgs
 *
 * 	Moreover, for a transparent chromosome rapresentation, each time a chromosome is passed to a function,
 * its descriptor will be as well. In no part of the algorithm, the rapresentation of the chromosome was used.	
 *
 *	Some of the arguments are optional and they will be commented as NULLABLE,
 * some others are not, so they will be indicated as NOT-nullable
 *	
 *	RETURNS: 
 * 		it returns a chromosome (which will have to be freed) 
 *		that betterly suits the problem with the args given.
 */
void* genericGeneticProblem(
	
	/* 1st arg
	 * DESCRIPTOR : this is the chromosome descriptor that the user defined.
	 * 				It is supposed to be unique and general for all the chromosomes.
	 *				Transparent for the algorithm
	 *				NULLABLE
	 */
	void* descriptor,

	/* 2nd arg
	 * MAX_GENERATION : It is the top generation reachable. 
	 *				If the solution wasn't found at MAX_GENERATION-th generation, the algorithm stops
	 *				and the best result is returned (but it will not be the solution of the problem)			
	 *				A value less than 0 rapresent an infinite number of generations.
	 */
	int max_generation, 

	/* 3th arg
	 *	POPULATION SIZE ARGS : It is the argument passed to population size
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* populationSizeArgs,

	/* 4th arg
	 *	POPULATION SIZE function : returns the wished population size
	 *							This function may be used to change the population size depending on the generation
	 *				1st arg: chromosome' descriptor
	 *				2nd arg: populationSizeArgs ( the one passed previously to genericGeneticProblem )
	 *				3th arg: generation index
	 *				NOTE: if the population size is constant, then, no malloc are done each loop and
	 *					a better performance is achieved.
	 *				NOT-nullble
	 */
	unsigned int (*populationSize)(void*, void*, int), 

	/* 5th arg
	 *	POPULATION INITIALIZER : provides the initial population
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: array of chromosomes to be filled (with pointers to the chromosomes)
	 *				3th arg: size of the array passed as 2nd arg. (the population size)
	 *				NOTE: the provided population has to be valid!
	 *				NOT-nullable
	 */
	void (*populationInitializer)(void*, void*[], int), 

	/* 6th arg
	 *	CHROMOSOME FREE : offers the service of freeing from the memory the chromosome
	 *					It may be a good idea that the users caches the momory locations and then reuses them
	 *					when it comes to create new chromosomes ( in the crossOver application ).
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: chromosome to be freed
	 *				NOT-nullable
	 */
	void (*chromosomeFree)(void*, void*), 

	/* 7th arg
	 *	IS CHROMOSOME VALID : it test the chromosome and return 0 iff the chromosome 
	 *						does not belog to the application domain 
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: chromosome to be tested
	 *				NULLABLE (rapresent a total application)
	 */
	char (*isChromosomeValid)(void*, void*), 

	/* 8th arg
	 *	TERMINATION CONDITION ARGS : arguments to be passed to terminationCondition (9th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* terminationConditionArgs,

	/* 9th arg
	 * 	TERMINATION CONDITION : it test the population champion and 
	 *							return 0 iff it does not satisfy the application
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: the termination condition args (the genericGeneticProblem's 8-th arg)
	 *				3th arg: the population champion (chromosome with the best score)
	 *				NULLABLE (the algorithm does not stop until the max_generation is reached
	 *							this option may be used when a terminal condition is unknown)
	 *				NOTE: if this argument is NULL and the MAX_GENERATION is -1, the algorithm will loop.
	 */
	char (*terminationCondition)(void*, void*, void*), 

	/* 10th arg
	 *	FITNESS ARGS : arguments to be passed to fitness (11th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* fitnessArgs,

	/* 11th arg
	 *	FITNESS : function that evaluates the fitness (score) of a chromosome.
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: fitnessArgs (10th arg)
	 *				3th arg: chromosome to be mesured
	 *				NOT-nullable
	 */
	double (*fitness)(void*, void*, void*),

	/* 12th arg
	 *	HOW MANY CROSS OVER ARGS : arguments to be passed to howManyCrossOver (13th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* howManyCrossOverArgs,

	/* 13th arg
	 *	HOW MANY CROSS OVER : Function that tells how many cross over to perform this generation
	 *				1st arg: the generation
	 *				2nd arg: howManyCrossOverArgs (12th arg)
	 *				NOTE: the user will have to be sure that return value is feasable
	 *					For example, if the father repetition is not allowed, then a too high
	 *					return value may lead to unexpected behaviours.
	 *					If father repetitio is not allowed, then the return value should be 
	 *					at most population_size / crossOverArity
	 *				NOT-nullable
	 */	
	unsigned int (*howManyCrossOver)(int, void*), 

	/* 14th arg
	 *	CAN FATHERS SURVIVE : 0 iff user does not want the father to survive in the population 
	 *						after the procreation. If set to true, then, after the procreation,
	 *						fathers will be set as invalid and evenctually replaced in the population.
	 *						BOOLEAN
	 */
	char canFathersSurvive,

	/* 15th arg
	 *	CAN FATHERS BE RE-ELECTED : 0 iff user does not want the fathers to be elected in more than
	 *						one crossover. If set to true, the fathers will be marked as invalid after
	 *						have beeing used by a crossOver operation. Anyway, at the end of the copulation
	 *						period, they will be reintroduced in the population.
	 *						BOOLEAN
	 */
	char canFathersBeReElected,

	/* 16th arg
	 *	CROSS OVER ARITY : indicates the arity of the cross over function
	 */
	unsigned int crossOverArity,

	/* 17th arg
	 *	HOW MANY CHILDRE CROSS OVER GENERATES : indicates how many children are 'returned' by cross over
	 */	
	unsigned int howManyChildrenCrossOverGenerates,

	/* 18th arg
	 *	CROSS OVER ARGS : arguments to be passed to crossOverArgs (19th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* crossOverArgs, 

	/* 19th arg
	 *	CROSS OVER : function that, given the parents, creates the childrens.
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: crossOVerArgs (18th arg)
	 *				3th arg: array of parent chromosomes (provided)
	 *				4th arg: parent array's size (=crossOverArity (16th arg))
	 *				5th arg: array of just born children (to be filled)
	 *				6th arg: children array's size (=howManyChildrenCrossOverGenerates (17th arg))
	 *				NOT-nullable
	 *				NOTE: if canFathersSurvive(14th arg) or canFathersBeReElected(15th arg)
	 *						then cross over MUST NOT modify the parents. 
	 *						else if the cross over funtion modify the parents, be sure to invalidate
	 *							the .chromosome reference or it will be freed. 
	 *				NOTE: the same problem may occour if the parent selection function (21th arg)
	 *						allows the fathers to be selected several times for the same crossover, 
	 *						in this case different fathers may refer to the same parent chromosome
	 */
	void (*crossOver)(void*, void*, void*[], int, void*[], int), 
	
	/* 20th arg
	 *	PARENT SELECTION ARGS : arguments to be passed to parentSelection (21th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* parentSelectionArgs, // argument to be passed to parentSelection function

	/* 21th arg
	 *	PARENT SELECTION : function that choses the parents which to perform a cross over with.
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: parentSelectionArgs (20th arg)
	 *				3th arg: generation
	 *				4th arg: the whole population (contains element marked as invalid) (provided)
	 *				5th arg: population size (size of 4th arg)
	 *				6th arg: parent refereces (array of reference to the population) (to be filled)
	 *				7th arg: parent references array size ( = crossOverArity (16th arg))
	 *				NOT-nullable
	 */
	void (*parentSelection)(void*, void*, int, chromosome_t[], int, chromosome_t*[], int), 

	/* 22th arg
	 *	MUTATE ARGS : arguments to be passed to mutate (23th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* mutateArgs,

	/* 23th arg
	 *	MUTATE : function that given a chromosome, it modifies it 
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: mutateArgs (22th arg)
	 *				3th arg: chromosome to be modified
	 *				NULLABLE (no modification at all)
	 */
	void (*mutate)(void*, void*, void*), 

	/* 24th arg
	 *	REORDER ARGS : arguments to be passed to reorder (25th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* reorderArgs,

	/* 25th arg
	 *	REORDER : function that given a chromosome, it reorder it 
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: reorderArgs (25th arg)
	 *				3th arg: chromosome to be reorder
	 *				NULLABLE (no reordering at all)
	 */
	void (*reorder)(void*, void*, void*), 

	/* 26th arg
	 *	DUMPER ARGS : arguments to be passed to dumper (27th arg)
	 *				Transparent for the algorith
	 *				NULLABLE
	 */
	void* dumperArgs,

	/* 27th arg
	 *	DUMPER : function for user porpouse called at the end of each loop cycle. (e.g. print)
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: dumperArgs (26th arg)
	 *				3th arg: generation
	 *				4th arg: population (chromosome array)
	 *				5th arg: population size
	 *				NULLABLE
	 */
	void (*dumper)(void*, void*,int, void*[], int), 

	/* 28th arg
	 *	SINGLE CHROMOSOME DUMPER : function that give a chromosome, it prints it. 
	 *							used for debugging when VERBOSE is defined
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: the chromosome
	 *				NULLABLE
	 */
	void (*singleChromosomeDumper)(void*, void*), 

	/* 29th arg
	 *	COMPARER : compares two chromosomes and returns 0 iff they are different
	 *							used to create a population of unique chromosomes.
	 *							If comparer is not null, then children will be inserted in 
	 *							the population only if no equals chromosomes are present already
	 *				1st arg: chromosome's descriptor
	 *				2nd arg: chromosome 1
	 *				3th arg: chromosome 2
	 *				NULLABLE 
	 *				NOTE: if nullable, the population may be full of clones.
	 */
	char (*comparer)(void*, void*, void*), 
	
	/* 30th arg
	 *	GENERATION CREATED : output of the number of the generation created before stopping.
	 *				reference to a int to be written
	 *				NULLABLE
	 */
	int* generationCreated 	
);

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						GENERIC GENETIC PROBLEM SOLVER ENDE
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						GENERAL UTILITY FUNCTIONS
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/


/**
 * POPULATION SIZE : CONSTANT VALUE
 *					function to be used as population size paramater (4th arg).
 *					When genericGeneticProblem is called, pass as 3th argument the reference to the 
 *					wished population size. which will be passed to this function as second argument,
 *					and then returned
 *			1st arg: chromosome's descriptor
 *			2nd arg: pointer to the value to return
 *			3th arg: generation 
 *			RETURNS: *(int*)(2nd arg)
 */
unsigned int populationSize_constant( void* descriptor , void* args, int generation );

/**
 * HOW MANY CROSS OVER : CONSTANT VALUE
 *					function to be used as howManyCrossOVer paramater (13th arg).
 *					When genericGeneticProblem is called, pass as 12th argument the reference to the 
 *					wished amount of cross over. which will be passed to this function as second argument,
 *					and then returned
 *			1st arg: generation
 *			2nd arg: pointer to the value to return. or null
 *			RETURNS: *(int*)(2nd arg) or 1 id 2nd arg is null
 */
unsigned int howManyCrossOver_constant(int generation, void* args);

/**
 * PARENT SELECTION
 * function of the type 
 * 	void (parentSelection*)(void* ag1, void*, int, chromosome_t[] ag2, int ag3, chromosome_t*[] ag4, int ag5)
 * WHERE
 *  arg1 is the descriptor of the chromosome, mostly useless here
 *	arg2 is the parentSelectionArgs (20th arg)
 *  arg3 is the generation
 *  arg4 is the input array of chromosomes
 *  arg5 is arg3's size.
 *  arg6 is the array containing the selected fathers (pointers to the original population)
 *  arg7 is arg5's size. 
 *
 *	arg6 is the output.
 *	arg7 is the crossOverArity provided to the genericGeneticProblem
 *
 * Behaviour expected:
 * a subset of fathers is selected in base of custom algorithms.
 *
 * NOTE: 
 *  1) if you don't want the same father to be more than one parent for a cross over, 
 *		be sure not to allow this in the parentSelection function implementation
 *	2) for each call, the rapresentation invariant is that the chromosomes are sorted
 *  3) do srand or not before calling these functions. they will not call srand
 */
 /*************************************************************************************/
 /*************************************************************************************/
 /*************************************************************************************/
/** 
 * 	PROPORCIONA SELECTION 
 *	ROULETE: the higher the score, the probabler is for a chromosome to become father.
 *	Argument arg (the 2nd arg): can father be re-elected for the same crossover
 *			If not null, it has to be a pointer to a char used as boolean
 *				If &0 then the output array will be filled with unique references,
 *				else if !=&0, the output array may have the same father repetided
 *			else if null, no re-election fot the father is allowed  

 *	NOTE: this implementation ensures that 
 *		1) only the valid chromosomes are selected (the population may contain invalid ones)
 *		2) usable either if the fathers can or cannot be re elected ( depending on args )
 *		3) does not fill fatherSize with repeted values
 **/
void parentSelection_roulette(void*,void*,int,chromosome_t[],int,chromosome_t*[],int);

/** 
 *	TOURNAMENT SELECTION: the higher the score, the probabler is for a chromosome to become father.
 *	Argument arg (the 2nd arg): an array of 3 ints
 *			(2nd arg)[0] = (int) presure (between how many the tournament is)
 *			(2nd arg)[1] = (int) is deterministic. if 1 the champion is chosen by score, otherwise randomly
 *			(2nd arg)[2] = (int) is father allowed to be re-elected

 *	NOTE: this implementation ensures that 
 *		1) only the valid chromosomes are selected (the population may contain invalid ones)
 *		2) usable either if the fathers can or cannot be re elected ( depending on args )
 *		3) does not fill fatherSize with repeted values
 **/
void parentSelection_tournament(void*,void*,int,chromosome_t[],int,chromosome_t*[],int);

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						BINARY rapresentation:
/*	This rapresentation rapresent a chromosome as a sequences of 0s and 1s.
 *  They are contained in an array of char and the size of the array is constant.
 * 	The size is pointed by the descriptor.
 */
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
/**
 * POPULATION INITIALIZER RANDOM:
 *						Initialize the population with a random siquence of booleans
 *						rapresented as char in a array of length *(int*)2nd arg.
 *		1st arg: chromosome's descriptor
 *		2nd arg: array to be fill with the population
 *		3th arg: size of the population
 */
void populationInitializer_binary_random(void*, void*[], int);

/**
 * CHROMOSOME FREE BINARY : 
 *						Free the array of chars.
 *		1st arg: chromosome's descriptor
 *		2nd arg: chromosome to be freed
 */
void chromosomeFree_binary(void*, void*);

/**
 *	CROSS OVER BINARY ONE POINT:
 *			Function that goes from two parent chromosome to two children chromosome.
 *			A randon number is sorted and the children will have a first part 
 *			(whose boundry is the number) of genes from a parent, and the other part
 *			from the other.
 *
 *			NOTE: crossOverArity and howManyChildrenCrossOverGenerates must be of value 2.
 *		
 *		1st arg: chromosome's descriptro
 *		2nd arg: cross over arguments 
 *		3th arg: array of chromosomes that have the role of parent.
 *		4th arg: the size. it has to be two
 *		5th arg: array to be filled with childs
 *		6th arg: the size. it has to be two
 */
void crossOver_binary_onePoint(void*, void*, void*[], int, void*[], int);

/*
 *	MUTATION BINARY CANONICAL:
 *				Given a chromosome, it flip each of its bit with probability provided
 *		1st arg: chromosome's descriptor
 *		2nd arg: mutation arg. It has to be a pointer to a double indicating the 
 *			mutation probability. If it is null, the probability is 0.000001.
 *		3th arg: the chromosome to modify
 */
void mutation_binary_canonical(void*, void*, void*);

/*
 * SINGLE CHROMOSOME DUMPER BINARY:
 *				it just print the sequence of bit rapresenting the chromosome
 *				and if the sequence is small enoguh (actually quite big), 
 *				his numeric value as well
 *		1st arg: chromosome's descriptor
 *		2nd arg: the chromosome to print
 */
void singleChromosomeDumper_binary(void* descriptor, void* chromosome);

/*
 * COMPARER BINARY:
 *				returns 0 if at leat one of the bit of the two chromosomes are different
 *		1st arg: chromosome's descriptor
 *		2nd arg: one chromosome
 *		3th arg: another one.
 */
char comparer_binary(void*, void*, void*);

/*
 * SIMPLE BINARY GENETIC PROBLEM:
 *				it wraps the generic genetic problem with default decisions
 *		1st arg: the size of the bit sequence 
 *		2nd arg: max generation
 *		3th arg: population size
 *		4th arg: fitness arguments
 *		5th arg: fitness function
 *		6th arg: termination test function
 *		7th arg: pointer to the int to be filled with the generation reached
 *		
 *		NOTE: the default settings are:
 *			pointer to sequence size as chromosome's descriptro
 *			constant population size
 *			population initialized randomly
 *			a fourth of the population is the number of crossover done
 *				since arity is two and repetition is not allowe, half of the population will be involved
 *			constant number of cross overs
 *			parents can survive 
 *			parents can't be re - elected for more than one cross over each generation
 *			arity = 2
 *			number of children generated by cross over = 2
 *			cross over used: one point
 *			parent selection used : roulete
 *				with no repetition of father
 *			canonical mutation
 *			comparer bit a bit.
 */
char* simpleBinaryGeneticProblem(int,int,int,void*,double(*)(void*,void*,void*),char(*)(void*,void*,void*),int*);



/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
//						DOUBLE rapresentation:
/*	This rapresentation rapresent a chromosome as a pointer to a double.
 * 	The descriptor is an array of two doubles, where the first value is the min and the second is the max.
 */
/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
/**
 * POPULATION INITIALIZER RANDOM:
 *						Initialize the population with a random siquence doubles
 *		1st arg: chromosome's descriptor
 *		2nd arg: array to be fill with the population
 *		3th arg: size of the population
 */
void populationInitializer_double_random(void*, void*[], int);

/**
 * CHROMOSOME FREE DOUBLE : 
 *						Free the double pointer.
 *		1st arg: chromosome's descriptor
 *		2nd arg: chromosome to be freed
 */
void chromosomeFree_double(void*, void*);

/**
 *	CROSS OVER DOUBLE INTERMEDIATE:
 *			Function that goes from two parent chromosome to one children chromosome.
 *			A randon number is sorted and the children will be the weigthed mean
 * 			of the two fathers where the sorted nomber is the weight for the first one 
 *			and the to one complement is the weight for the second
 *
 *			NOTE: crossOverArity and howManyChildrenCrossOverGenerates must be of value 2 and 1.
 *		
 *		1st arg: chromosome's descriptro
 *		2nd arg: cross over arguments 
 *		3th arg: array of chromosomes that have the role of parent.
 *		4th arg: the size. it has to be two
 *		5th arg: array to be filled with childs
 *		6th arg: the size. it has to be one
 */
void crossOver_double_intermediate(void*, void*, void*[], int, void*[], int);

/*
 *	MUTATION DOUBLE NON UNIFORM FOR REALs:
 *				With aprobability prob, the double is mutated in the following way.
 *				With a probability of 0.5, the chromosomes takes the value vk1 or vk2
 *				where :
 *					vk1 = chromosome + delta( t, max - chromosome )
 *					vk2 = chromosome + delta( t, chromosome - min )
 *					delta( t, v ) = v * ( 1 - r^( 1 - (t/T)^b ) )
 *				where : 
 *					t = current generation
 *					T = max generation
 *					r = random between 0 y 1
 *					min, max passed as argument
 *					prob, b are params.
 *		1st arg: chromosome's descriptor 
 *				which has to be an array of double[2] = { min, max }
 *		2nd arg: mutation arg. 
 *				it has to be the following array of double
 *				double[4] = { prob, b, T, t }
 *		3th arg: the chromosome to modify
 */
void mutation_double_non_uniform(void*, void*, void*);


/*
 * SINGLE CHROMOSOME DUMPER DOUBLE:
 *				it just print the double pointed
 *		1st arg: chromosome's descriptor
 *		2nd arg: the chromosome to print
 */
void singleChromosomeDumper_double(void* descriptor, void* chromosome);

/*
 * COMPARER DOUBLE:
 *				returns 0 if the ratio between the two chromosome is bigger than DOUBLE_TOLLERANCE
 *		1st arg: chromosome's descriptor
 *		2nd arg: one chromosome
 *		3th arg: another one.
 */
char comparer_double(void*, void*, void*);

/*
 * SIMPLE DOUBLE GENETIC PROBLEM:
 *				it wraps the generic genetic problem with default decisions
 *		1st arg: min value
 *		2nd arg: max value
 *		3th arg: max generation
 *		4th arg: population size
 *		5th arg: fitness arguments
 *		6th arg: fitness function
 *		7th arg: termination test function
 *		8th arg: pointer to the int to be filled with the generation reached
 *		RETURNS the best double or NAN
 *		
 *		NOTE: the default settings are:
 *			double[2]={min, max} as chromosome's descriptro
 *			constant population size
 *			population initialized randomly
 *			a fourth of the population is the number of crossover done
 *				since arity is two and repetition is not allowe, half of the population will be involved
 *			constant number of cross overs
 *			parents can survive 
 *			parents can't be re - elected for more than one cross over each generation
 *			arity = 2
 *			number of children generated by cross over = 1
 *			cross over used: intermediate for doubles
 *			parent selection used : tournament
 *				presure = 2,
 *				deterministic
 *				father can't be re-elected
 *			non uniform mutation for double
 *			comparer if ratio smaller than DOUBLE_TOLLERANCE.
 */
double simpleDoubleGeneticProblem(double,double,int,int,void*,double(*)(void*,void*,void*),char(*)(void*,void*,void*),int*);

#endif