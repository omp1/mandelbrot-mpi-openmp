#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

extern void matToImage(char* filename, int* mat, int* dims);
extern void matToImageColor(char* filename, int* mat, int* dims);


int main(int argc, char** argv){

    int rank;
    int numranks;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);

    omp_set_num_threads(6);

    //new
    

    int nx=30000;
    int ny=20000;
    // int* matrix=NULL;
    // if (rank==0){
        
    // }
    int* matrix = NULL;
    if (rank==0){
        matrix=(int*)malloc(nx*ny*sizeof(int));
        for (int i = 0; i < nx * ny; i++) {
            matrix[i] = 0;
        }
    }
    
    
    int maxIter=255;
    double xStart=-2;
    double xEnd=1;
    double yStart=-1;
    double yEnd=1;

    //int n =20;
    //int n=10000000;
    
    int work = 300; //numrows per task


    //master assigns chunks of work and workers keep coming back for more work once they finish their chunk
    if (rank==0){
        
       
        int current = 0;
        int doneWorkers = 0;

        

        //master here
        //assign current to -1 when no more work left
        
        for (int i=1; i<numranks; i++){
            int range[] = {0, 0}; 
            if (current <=ny){
                range[0] = current;
                range[1] = current+work-1;

                current+=work;
                if (range[1] > ny){
                    range[1] = ny;
                }

                MPI_Send(range, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            else if (current>ny){
                //no work left
                //invalid range
                range[0] = -1;
                MPI_Send(range, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
                doneWorkers++;
    

            }
        }

        //manager here
        //when workers report back after finishing work
        while (doneWorkers < numranks-1){
            
            MPI_Status stat;
            int rowIndex;
            int range[] = {0, 0}; 
            
            //recv from workers
            MPI_Recv(range, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
            
            int startingRow = range[0];
            int numRows = range[1]-range[0]+1;
            int* gatheringBuffer = (int*)malloc(numRows*nx*sizeof(int));

            MPI_Recv(gatheringBuffer, numRows*nx, MPI_INT, stat.MPI_SOURCE, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int i=0; i<numRows; i++){
                for (int j=0; j<nx; j++){
                    matrix[(startingRow+i)*nx+j] = gatheringBuffer[i*nx+j];
                }
            }
            free(gatheringBuffer);
            
            if (current>ny){
                range[0] = -1;
                doneWorkers++;
                MPI_Send(range, 2, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
            }
            else{
                //same as the else in upper code for when workers come back
                range[0] = current;
                range[1] = current+work-1;

                current+=work;

                if (range[1] >= ny){
                    range[1] = ny-1;
                }

                MPI_Send(range, 2, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
            }




        }
        int dims[2] = {ny, nx};
        matToImageColor("mandelbrot.jpg", matrix, dims);
        
        
        

    } //workers here
    else { //rank!=0
        
        double mpiStart = MPI_Wtime();
        int range[] = {0,0};
        
        //new
        // int numThreads;
        // double* threadTimes = NULL;
        int numThreads = omp_get_max_threads();
        double* threadTimes = (double*)malloc(numThreads*sizeof(double));
        for (int i=0; i<numThreads; i++){
            threadTimes[i]=0.0;
        }


        while (1){ // no true keyword in c
            MPI_Recv(range, 2, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (range[0] == -1){
                break;
            }

            

            int rows = range[1]-range[0]+1;
            int* localBuffer = (int*)malloc(nx*rows*sizeof(int));


           
            #pragma omp parallel
            {
                
                
                

                int threadID = omp_get_thread_num();
                double startone = omp_get_wtime(); 

                #pragma omp for schedule(dynamic, 10), nowait
                for (int i=range[0]; i<=range[1]; i++){
                    for(int j=0;j<nx;j++){
                        int index=i*nx+j;
                        int iter=0;
                        //Z=x*iy->0 to start
                        double x=0;
                        double y=0;
                        //C=x0+iy0 -> based on the pixel
                        double x0=xStart+(1.0*j/nx)*(xEnd-xStart);
                        double y0=yStart+(1.0*i/ny)*(yEnd-yStart);
                        while(iter<maxIter){
                            iter++;
                            double temp=x*x-y*y+x0;
                            y=2*x*y+y0;
                            x=temp;
                            if(x*x+y*y>4){
                                break;
                            }
                        }
                      localBuffer[(i-range[0])*nx+j]=iter;
                    }   
                    
                    

                    
                
                }

                double endone = omp_get_wtime();
                threadTimes[threadID] += (endone-startone);
            
                
                
            }
            
            
            MPI_Send(range,2,MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(localBuffer, rows*nx, MPI_INT, 0, 2, MPI_COMM_WORLD);
            if (localBuffer){
                free(localBuffer);
            }
            
        
            

           
            
        }
        double mpiend = MPI_Wtime();
        double elapsedTime = mpiend-mpiStart;
        for (int i = 0; i < numThreads; i++) {
            printf("Rank %d Thread %d total time: %f seconds\n", rank, i, threadTimes[i]);
        }
        
        printf("Rank %d Time %f \n", rank, elapsedTime);
        
        

        free(threadTimes);

    }
    

    if (matrix){
        free(matrix);
    }
    
    



    


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;

}






