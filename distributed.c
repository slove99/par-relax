// Matrix relaxation using MPI
// Candidate Number: 11066

// Included libraries
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Function prototypes
void allocateSections (int cores, int arrayRows, int* borders);
double** createMatrix(int borderValue, int arrayRows);
void printMatrix(int arrayRows, double** myMatrix);
void calcMatrix(int startPoint, int endPoint, double** myMatrix, int world_size,
        int world_rank, double precision, int arrayRows);

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    
    // Find thread rank
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // Find world size
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get processor name
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Timing variables
    time_t begin;
    time_t end;
    
    // ENTER MATRIX PARAMETERS
    int arrayRows = 500; // Width and height of matrix
    double precision = 0.001; // Precision of matrix
  
    double **myMatrix = createMatrix(12, arrayRows);
    int *borderArray = malloc((world_size+1) * sizeof(int *));
    allocateSections(world_size, arrayRows, borderArray);
    //printf("BORDER ARRAY: %d", borderArray[world_rank]);
    
    begin = time(NULL);
    calcMatrix(borderArray[world_rank], borderArray[world_rank+1], 
            myMatrix,world_size, world_rank, precision, arrayRows);
    int j = 0;
    if(world_rank == 0)
    {
        //Recv sections from threads
        for(j=1;j<world_size;j++){
            int blockSize = borderArray[j+1] - borderArray[j];
            //Store recv data in correct section of matrix
            MPI_Recv(&myMatrix[borderArray[j]][0], blockSize*arrayRows
                    + (2*(blockSize-1)), MPI_DOUBLE, j, 0, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);
        }
        //printMatrix(arrayRows, myMatrix);
        end = time(NULL);
        printf("\nCompleted processing array of %dx%d elements with precision "
                "%f\n",arrayRows, arrayRows, precision);
        printf("\nComputation took %d seconds", end-begin);
        printf("\nComputation used %d threads", world_size);

    }
    MPI_Finalize();
}

//allocateSections
//INPUT:   sections (number of sections for array to be divided in)
//         arrayRows (size of array to be divided)
//PROC: Take matrix dimensions and return array of row borders for each thread
//OUTPUT: borders (array of row borders for each thread)
void allocateSections (int sections, int arrayRows, int* borders)
{   
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
    //for (j = 0; j < sections+1; j++)
        //printf("\n\nBORDERS: %d\n", borders[j]);
}

//createMatrix
//INPUT: Value to be placed on the boundaries of the matrix (borderValue)
//       Size of the matrix (arrayRows)
//PROC:  Creates a dynamically sized array with boundary values (borderValue)
//OUT:   The created array exists in memory
double** createMatrix(int borderValue, int arrayRows)
{
    int i = 0;
    int j = 0;
    double **myMatrix;
    //Dynamically generate array for values
    myMatrix = calloc(arrayRows, sizeof(double *));
    for(i=0; i<arrayRows; i++)
    {
        myMatrix[i] = calloc(arrayRows, sizeof(double *));
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
    return myMatrix;
}


//printMatrix
// Prints the content of the main matrix
// INPUT: arrayRows (size of matrix)
// PROC: Prints the content of the matrix
// OUT: (N/A)
void printMatrix(int arrayRows, double** myMatrix)
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

// calcMatrix
// Main function to perform matrix calculation
// INPUT: start (start row of allocated section) 
//        end (end row of allocated section)
//        myMatrix (threads local matrix to calculate values with)
//        world_size (number of threads in the world)
//        world_rank (thread rank in the world)
//        precision (precision / accuracy to work towards)
//        arrayRows (number of row and columns in the entire matrix)
// PROC:  Calculates the allocated section of the matrix
// OUT:   N/A (The processed section of the matrix to be sent to root thread)
void calcMatrix(int start, int end, double** myMatrix, int world_size, 
        int world_rank, double precision, int arrayRows)
{
    printf("\nThread %d has started and has the following properties \n "
                "Section: %d\nstartPoint: %d\nendPoint: %d\nprecision: %lf\n"
                "array total rows: %d\ntotal threads: %d\n\n", world_rank,
            world_rank, start, end, precision, arrayRows, world_size);
    
    int a = 0;
    int b = 0;
    int c = 0;
    int f = 0;
    int globalPrecisionNotMet; //Sum of all thread precisionNotMet status
    double *recRow = malloc(arrayRows * sizeof(double *));
    
    do{
        int precisionNotMet = 0; //Individual thread precisionNotMet status
        globalPrecisionNotMet = 0;
        for(a = start; a<end-1; a++) //Iterate through all allocated rows
        {
            for(b = 1; b<arrayRows-1; b++){ //Iterate through all columns
                double oldValue = myMatrix[a][b];
                myMatrix[a][b] = (myMatrix[a-1][b] + myMatrix[a][b-1]
                            + myMatrix[a+1][b] + myMatrix[a][b+1]) / 4;
                if((fabs(oldValue - myMatrix[a][b]) > precision))
                    precisionNotMet = 1;
            }
       }
        // Code for exchange

        // Thread ahead sends computed top row so previous thread can use it 
        // to calc bottom row 
        if(world_rank != 0) // Thread 0 computes top section, can't send data up
            MPI_Ssend(&myMatrix[start][0], arrayRows, MPI_DOUBLE, 
                    world_rank-1, 0, MPI_COMM_WORLD);

        if(world_rank!= (world_size-1)){ //End thread has no data to recv
                                         //from below
            //Recv computed top row from prev thread
            MPI_Recv(recRow, arrayRows, MPI_DOUBLE, world_rank+1, 0, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(c=1; c<arrayRows-1; c++) //Write recv data to local matrix
                myMatrix[end][c] = recRow[c];
            }
        // Calculate final allocated row using received values
        for(c=1; c<arrayRows-1; c++){
            double oldValue = myMatrix[a][b];
            myMatrix[end-1][c] = (myMatrix[end-2][c] + myMatrix[end-1][c-1]
            + myMatrix[end][c] + myMatrix[end-1][c+1]) / 4;
            if((fabs(oldValue - myMatrix[a][b]) > precision))
                precisionNotMet = 1;
        }
        // Thread behind sends back computed bottom row so next thread 
        // can use it to calc top row
        if(world_rank!= (world_size-1))
            MPI_Ssend(&myMatrix[end-1][0], arrayRows, MPI_DOUBLE, 
                    world_rank+1, 0, MPI_COMM_WORLD);
        if(world_rank!=0){
            //Recv computed top row from prev thread
            MPI_Recv(recRow, arrayRows, MPI_DOUBLE, world_rank-1, 0, 
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(c=1; c<arrayRows-1; c++) //Write recv data to local matrix
                myMatrix[start-1][c] = recRow[c];
        }
        //if(precisionNotMet == 0)
        //    printf("\nPrecision MET on thread %d\n", world_rank);
        MPI_Barrier(MPI_COMM_WORLD); //Ensure all threads are on same iteration
        
        //Sum all precisionNotMet together - store in globalPrecisionNotMet
        //When globalPrecisionNotMet is 0, all threads have reached precision
        MPI_Allreduce(&precisionNotMet, &globalPrecisionNotMet, 1, MPI_INT, 
                MPI_SUM, MPI_COMM_WORLD);
        
  }while(globalPrecisionNotMet != 0);
  printf("\nPRECISION MET ON ALL THREADS\n");
  
  //Send computed section back to main thread (thread 0)
  int blockSize = end - start;
  MPI_Barrier(MPI_COMM_WORLD); //Ensure all threads have exited while loop
  if(world_rank != 0) //Thread 0 doesn't need to recv its own data
      //Offset of (2*(blockSize-1)) added to compensate for pairs of 0's
      //at the end of each sent row
    MPI_Ssend(&myMatrix[start][0], blockSize*arrayRows + (2*(blockSize-1)), 
            MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}
