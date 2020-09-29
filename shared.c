
//CM30225 Coursework 1 - Shared memory programming
//Candidate Number: 11066

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>


//Function prototypes
void *calcMatrix(void *argsStruct);
int* allocateSections (int cores, int arrayRows);
void createMatrix(int borderValue, int arrayRows);
void printMatrix(int arrayRows);
int arrayEmpty(int sections);
void calcResult(double **matrix, double precision, int sections, int arrayRows);
void freeArrays();

//Global arrays
int *borders; // Internal borders for each thread
double** myMatrix; //Matrix to be processed

// Array to store whether precision has been met for each thread
int* precisionMet; 

//Parallel variables (thread / mutex / barrier)
pthread_t *myThreads;
pthread_mutex_t *myMutex;
pthread_barrier_t myBarrier;

//Struct to be passed into each thread
struct argumentsForFunct 
{ 
    double **myMatrix;
    int section;
    int sections;
    int startPointCol;
    int endPointCol;
    int startPointRow;
    int endPointRow;
    double threadPrecision;
    int arrayRows;
}; 


//main
// Proc: Sets up conditions for calcResult function and calls calcResult
//       Creates and frees arrays used
int main() {
    //Main function variables
    
    // Timing arrays for benchmarking
    time_t begin[10];
    time_t end[10];   
    
    // Adjustable variables
    int arrayRows = 15; // Size of array
    double precision = 0.0001; // Precision (delta) to work towards
    int sections = 2; // Threads to utilise
    
    
    // Prevent more threads being requested than rows in matrix
    if(sections>arrayRows)
    {
        sections = arrayRows;
        printf("Using maximum threads: %d\n\n", arrayRows);
    }
    
    //*****Create process and print array*****
    printf("Using a %d square array with precision %lf\n", arrayRows, precision);
    createMatrix(10, arrayRows); // Create 2D array for processing
    printMatrix(arrayRows); 
    //Calculate solution
    begin[0] = time(NULL); //Begin timer
    calcResult(myMatrix, precision, sections, arrayRows); //Process matrix
    end[0] = time(NULL); //End timer
    printMatrix(arrayRows); 
    freeArrays();
    //****************************************
    
    //Print time taken
    printf("Time taken on %d threads = %d sec\n", sections, end[0] - begin[0]);
    
    return (EXIT_SUCCESS);
}


//calcResult (matrix solver)
//INPUT: matrix to be processed (matrix), accuracy to work towards (precision)
//       number of threads to be allocated to the task (sections)
//       number of rows (and columns) in the array (arrayRows)
// PROC: Initialises mutex and barriers
//       Creates and joins threads for calcMatrix and sets up their parameters
//       Frees malloced arrays created
//OUTPUT: N/A (array is now processed)
void calcResult(double **matrix, double precision, int sections, int arrayRows)
{
    int i = 0;
    
    //Malloc arrays
    myThreads = malloc(sizeof(pthread_t)*sections);
    myMutex = malloc((sections+1) * sizeof(pthread_mutex_t));
    precisionMet = calloc(sections, sizeof(int));
    
    //Intialise parallel variables
    for (i=0;i<sections+1;i++) 
    {
        if((pthread_mutex_init(&myMutex[i], NULL)) == 0)
            printf("Initialised mutex %d successfully\n", i);
        pthread_mutex_unlock(&myMutex[i]);
    }
    
    pthread_barrier_init(&myBarrier, NULL, sections); 
    
    //Allocate boundaries of each threads processing area
    borders = allocateSections(sections, arrayRows);
    
    
    //Dynamically generate argument structs for each thread
    struct argumentsForFunct *allArguments = 
    (struct argumentsForFunct *)malloc(sections*sizeof(struct argumentsForFunct));
    for(i=0; i<sections; i++){
        (allArguments+i)->section = i;
        (allArguments+i)->sections = sections;
        (allArguments+i)->myMatrix = myMatrix;
        (allArguments+i)->startPointCol = borders[(allArguments+i)->section];
        (allArguments+i)->endPointCol = borders[(allArguments+ i)->section+1];
        (allArguments+i)->threadPrecision = precision;
        (allArguments+i)->arrayRows = arrayRows;
    }
    
    
    //Create threads for program
    for(i = 0; i < sections; i++)
            pthread_create(&myThreads[i], NULL, &calcMatrix, allArguments+i);
    
    //Join threads for program
    for(i = 0; i < sections; i++)
        pthread_join(myThreads[i], NULL);
    
    free(allArguments);
}


//calcMatrix (thread function)
//INPUT: struct of arguments unique to the thread being used (argsStruct)
//PROC:  Calculates the value of a given segment of the matrix
//OUT:   N/A (the required section of the matrix is processed)
void *calcMatrix (void *argsStruct)
{
    // Counter variables for function
    int a = 0;
    int b = 0;  
    
    //Separate out struct elements
    int startPoint = ((struct argumentsForFunct*)argsStruct)->startPointCol;
    int endPoint = ((struct argumentsForFunct*)argsStruct)->endPointCol;
    double **myMatrix = ((struct argumentsForFunct*)argsStruct)->myMatrix;
    int section = ((struct argumentsForFunct*)argsStruct)->section;
    int sections = ((struct argumentsForFunct*)argsStruct)->sections; 
    double threadPrecision = ((struct argumentsForFunct*)argsStruct)->threadPrecision;
    int arrayRows = ((struct argumentsForFunct*)argsStruct)->arrayRows;
    printf("\nThread %d has started and has the following properties \n "
                "Section: %d\nstartPoint: %d\nendPoint: %d\nprecision: %lf\n"
                "array total rows: %d\ntotal threads: %d\n\n", section, section,
                startPoint, endPoint, threadPrecision, arrayRows, sections);
    do
    {
        pthread_barrier_wait(&myBarrier);
        
        precisionMet[section] = 0;
    // Iterate through each element in allocated area
    for(a = startPoint; a<endPoint; a++)
    { 
        //check if calc uses the bottom boundary and not the actual border
        if((a == borders[section]) & (section!= 0)) 
        {
            
            pthread_mutex_lock(&myMutex[section]);
            //printf("Thread %d has locked row %d\n", section, a); 
            
            for(b = 1; b<arrayRows-1; b++){
                double oldValue = myMatrix[a][b];
                myMatrix[a][b] = (myMatrix[a-1][b] + myMatrix[a][b-1]
                        + myMatrix[a+1][b] + myMatrix[a][b+1]) / 4;
                //printf("Thread %d, has accessed row %d col %d\n", section, a, b);
                if((fabs(oldValue - myMatrix[a][b]) > threadPrecision))
                    precisionMet[section] = 1;
            }
            pthread_mutex_unlock(&myMutex[section]);
            //printf("Thread %d has unlocked row %d\n", section, a); 
        }
        //check if calc uses the bottom boundary and not the actual border
        else if((a+1 == borders[section+1]) & (section!= arrayRows)) 
        {
            //printf("Thread %d, wants to access row %d\n", section, a+1);
            pthread_mutex_lock(&myMutex[section+1]);
            //printf("Thread %d has locked row %d\n", section, a+1); 
            
            for(b = 1; b<arrayRows-1; b++){
                double oldValue = myMatrix[a][b];
                myMatrix[a][b] = (myMatrix[a-1][b] + myMatrix[a][b-1]
                        + myMatrix[a+1][b] + myMatrix[a][b+1]) / 4;
                //printf("Thread %d, has accessed row %d col %d\n", section, a, b);
                if((fabs(oldValue - myMatrix[a][b]) > threadPrecision))
                    precisionMet[section] = 1;
            }
            
            pthread_mutex_unlock(&myMutex[section+1]);
            //printf("Thread %d has unlocked row %d\n", section, a+1); 
        }
        else
        {
            for(b = 1; b<arrayRows-1; b++){
                double oldValue = myMatrix[a][b];
                myMatrix[a][b] = (myMatrix[a-1][b] + myMatrix[a][b-1]
                        + myMatrix[a+1][b] + myMatrix[a][b+1]) / 4;
                //printf("Thread %d, has accessed r %d c %d\n", section, a, b);
                if((fabs(oldValue - myMatrix[a][b]) > threadPrecision))
                    precisionMet[section] = 1;
            }
        }      
    }
    pthread_barrier_wait(&myBarrier);
    }while(arrayEmpty(sections) == 0); //While thread precision is too low
    return NULL;
}


//createMatrix
//INPUT: Value to be placed on the boundaries of the matrix (borderValue)
//       Size of the matrix (arrayRows)
//PROC:  Creates a dynamically sized array with boundary values (borderValue)
//OUT:   The created array exists in memory
void createMatrix(int borderValue, int arrayRows)
{
    int i = 0;
    int j = 0;
    //Dynamically generate array for values
    myMatrix = malloc(arrayRows* sizeof(double *));
    for(i=0; i<arrayRows; i++)
    {
        myMatrix[i] = malloc(arrayRows * sizeof(double *));
    }
    //Fill outside of array with values
    for(i=0; i<arrayRows;i++)
    {
        for(j=0; j<arrayRows; j++)
        {
            if( j==0 || i==0 || j == arrayRows-1 || i == arrayRows-1)
                myMatrix[i][j] = borderValue;
        }
    }
}


//printMatrix
// Prints the content of the main matrix
// INPUT: arrayRows (size of matrix)
// PROC: Prints the content of the matrix
// OUT: (N/A)
void printMatrix(int arrayRows)
{
    //Function variables
    int i = 0;
    int j = 0;
    
    //Print preprocessed array
    for (i = 0; i<arrayRows; i++)
    {
        for (j = 0; j<arrayRows; j++)
        {
            printf("%f  ", myMatrix[i][j]);
        }
        printf("\n");
    }
}


//allocateSections
//INPUT:   sections (number of sections for array to be divided in)
//         arrayRows (size of array to be divided)
//PROC: Take matrix dimensions and return array of row borders for each thread
//OUTPUT: borders (array of row borders for each thread)
int* allocateSections (int sections, int arrayRows)
{
    int *borders = malloc((sections+1) * sizeof(int *));
    
    //Calculate quotient and extra values
    int quot = arrayRows / sections;
    int extra = arrayRows % sections;
    int j = 0;
    
    borders[0] = 1;
    borders[1] = quot;
    
    for (j = 1; j < sections; j++)
    {
        // Place next border at least quot ahead of current border
        borders[j+1] = borders[j] + quot;
        
        //If remainder values that need to be shared across border distance
        //still exist
        if (j < extra) 
            borders[j+1] ++; //Increment next border value
    }
    // Avoid placing border value on boundary
    if(borders[sections] == arrayRows)
        borders[sections] --;
                  
    //Test code print border ranges
    for (j = 0; j < sections+1; j++)
        printf("\n\nBORDERS: %d\n", borders[j]);
    return borders;
}


//arrayEmpty
//INPUT: sections (number of elements in precisionMet array)
//PROC: If all array elements are 0 return 1, else return 0
//OUT int determining array state
int arrayEmpty(int sections)
{
    //printf("I'm running arrayEmpty\n\n");
    int i = 0;
    for(i = 0; i<sections; i++)
    {
        //printf("PRECISION MET %d\n\n", precisionMet[i]);
        if(precisionMet[i]!= 0)
            return 0;
    }
    return 1;
}


//freeArrays
//PROC: Frees any remaining malloced arrays
void freeArrays()
{
    //Free malloced arrays
    free(myThreads);
    free(myMutex);
    free(precisionMet);
    free(myMatrix);
    free(borders);
}