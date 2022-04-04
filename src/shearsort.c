#define _XOPEN_SOURCE
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void quicksort(double *data, int left, int right, int n, int sign);
int partition(double *data, int left, int right, int pivotIndex, int sign);
void printMatrix(double *data, int n);

int main(int argc, char **argv) {
    if (4 != argc) {
        printf("Usage: rows/columns input-file output-file\n");
        return 0;
    } 

    /* Declaration of variables */ 
    int i,j,k,p,d;
    int rank, nproc, send_batch;
    double *data, *row_data, *column_data;
    double startTime, stopTime,sortTime,distTime,gatherTime,totalTime;

    int n = atoi(argv[1]);
    int N = n*n;

    /* Setting up MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    send_batch = n/nproc;
    row_data=(double *)calloc(sizeof(double),N/nproc);
    column_data=(double *)calloc(sizeof(double),N/nproc);

    /* Read initial data from file */
    if (rank ==0) {
        char *file_read = argv[2];
        char *output = argv[3];
        char *scan = (char*)malloc(sizeof(char)*8);

        FILE *input_file = fopen(file_read, "r");
        (void)!fscanf(input_file, "%s", scan);
        n = atoi(scan);
        N=n*n;

        data=(double *)malloc(sizeof(double)*N);

        for (i=0; i<N; i++){
            (void)!fscanf(input_file, "%s", scan);
            data[i] = strtof(scan,NULL);
        }
        free(scan);
        fclose(input_file);
    }

    /* Create initial data */
    /*if (rank == 0) {
        srand48(2);
        data=(double *)malloc(sizeof(double)*N);
        for (i=0;i<N;i++) {
            // Uniform random numbers
            data[i]=drand48();
        }
    }*/

    /*if (rank==0) {
        printf("Initial data:\n"); 
        printMatrix(data,n);
    }*/

    /* Creating new data types, sending and reciving type */
    MPI_Datatype MPI_temp,MPI_send_type,MPI_recv_type;
    MPI_Type_vector(send_batch,1,n,MPI_DOUBLE,&MPI_temp);
    MPI_Type_commit(&MPI_temp);

    MPI_Type_create_resized(MPI_temp, 0, 1*sizeof(double), &MPI_send_type);
    MPI_Type_commit(&MPI_send_type);

    MPI_Type_vector(send_batch,send_batch,n,MPI_DOUBLE,&MPI_temp);
    MPI_Type_commit(&MPI_temp);

    MPI_Type_create_resized(MPI_temp, 0, send_batch*sizeof(double), &MPI_recv_type);
    MPI_Type_commit(&MPI_recv_type);

    startTime = MPI_Wtime();
    /* Distribution of data */
    MPI_Scatter(data,N/nproc,MPI_DOUBLE,row_data,N/nproc,MPI_DOUBLE,0,MPI_COMM_WORLD);
    stopTime = MPI_Wtime();
    distTime = stopTime - startTime;

    d = ceil(log2(n));
    startTime = MPI_Wtime();
    /* Sorting loop */
    for (i=0; i<=d+1; i++) {
        for (k=0; k<send_batch; k++) {
            quicksort(row_data,k*n,(k+1)*n-1,n,1-2*((k+(rank%2)*(send_batch%2))%2)); //row sort ascending/descending order (even/odd)
        }
        /* Distribution of data */
        MPI_Alltoall(row_data,send_batch,MPI_send_type,column_data,1,MPI_recv_type,MPI_COMM_WORLD);

        if (i<=d) {
            for (k=0; k<send_batch; k++) {
                quicksort(column_data,k*n,(k+1)*n-1,n,1); //column sort ascending
            }
            /* Distribution of data */
            MPI_Alltoall(column_data,send_batch,MPI_send_type,row_data,1,MPI_recv_type,MPI_COMM_WORLD);
        }
    } 
    stopTime = MPI_Wtime();
    sortTime = stopTime-startTime;

    startTime = MPI_Wtime();
    MPI_Gather(row_data,send_batch*n,MPI_DOUBLE,data,send_batch*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    stopTime = MPI_Wtime();
    gatherTime = stopTime-startTime;

    /* Printing data/write to file + time */
    if (rank==0) {
        //printf("Sorted data:\n"); 
        //printMatrix(data,n);
        
        /* Write to file */
        char *output = argv[3];
        FILE *file=fopen(output,"w");
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                fprintf(file,"%6.2f ",data[j+i*n]);
            }
            fprintf(file,"\n");
        }
        fclose(file);
    }
    totalTime = distTime+sortTime+gatherTime;
    if (rank==0) printf("distTime=%f, sortTime=%f, gatherTime=%f, totalTime=%f\n",distTime,sortTime,gatherTime,totalTime);

    /* Free calls */
    if (rank==0) free(data);
    free(row_data);
    free(column_data);
    MPI_Finalize();
    return 0;
}

/* Printing matrix function */
void printMatrix(double *data, int n) {
    int i,j;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            printf("%6.2f ",data[j+i*n]);
        }
        printf("\n");
    }
    printf("\n");
}

/* Quicksort function */
void quicksort(double *data, int left, int right, int n, int sign){
    int pivotIndex, pivotNewIndex;
    if (right > left){
        pivotIndex = left+(right-left)/2;
        pivotNewIndex = partition(data, left, right, pivotIndex, sign);
        quicksort(data,left, pivotNewIndex - 1, n, sign);
        quicksort(data,pivotNewIndex + 1, right, n, sign);
    }
}

/* Partitioning function */
int partition(double *data, int left, int right, int pivotIndex, int sign){
    double pivotValue,temp;
    int storeIndex,i;
    pivotValue = data[pivotIndex];
    temp=data[pivotIndex]; data[pivotIndex]=data[right]; data[right]=temp;
    storeIndex = left;
    for (i=left;i<right;i++)
    if (sign*data[i] <= sign*pivotValue){
        temp=data[i];data[i]=data[storeIndex];data[storeIndex]=temp;
        storeIndex++;
    }
    temp=data[storeIndex];data[storeIndex]=data[right]; data[right]=temp;
    return storeIndex;
}
