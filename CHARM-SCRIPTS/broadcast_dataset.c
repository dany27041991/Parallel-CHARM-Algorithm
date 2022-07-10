#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mpi.h"
#include "common.h"
#include "broadcast_dataset.h"
#include <unistd.h>
#include <math.h>

//Function that allows the file to be read.
void read_data(struct itemset_tids *P, char* namefile, int num_el_P) {
    FILE *fptr;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    if ((fptr = fopen(namefile,"r")) == NULL){
       printf("Error! opening file");
       exit(1);
    }
    int i = 0;
    while ((read = getline(&line, &len, fptr)) != -1) {
        if(i==num_el_P)
            break;
        int counter_tids = 1;
        char *key_ptr = strtok(line, "-");
        P[i].itemsets[0] = 1;
        P[i].itemsets[1] = atoi(key_ptr);
        char *tids_string = strtok(NULL, "-");
        char *tids_ptr = strtok(tids_string, " ");
        while( tids_ptr != NULL ) {
            size_t length = strlen(tids_ptr);
            if((length > 0) && (tids_ptr[length-1] == '\n'))
            {
                tids_ptr[length-1] ='\0';
            }
            P[i].tidsets[counter_tids] = atoi(tids_ptr);
            counter_tids++;
            tids_ptr = strtok(NULL, " ");
        }
        P[i].tidsets[0] = --counter_tids;
        P[i].index = 0;
        if(FLAG_SORT_TIDS == TRUE)
            sort_tids_itemsets(P[i-1].tidsets, P[i-1].tidsets[0]);
        i++;
    }
    fclose(fptr);
    if (line)
        free(line);
}

//Function for broadcasting from master to workers.
void bcast_dataset(struct itemset_tids *P, char* namefile, MPI_Datatype item_type, int rank, int processes, int num_el_P) {
    //If the process has rank 0 (master) I read the data.
    if(rank == 0) {
        if ((processes - 1) > num_el_P) {
            fprintf(stderr, "There are a number of process greater than the number of data present.\n");
            exit(0);
        }
        read_data(P, namefile, num_el_P);
        //I sort the data in order of support.
        sort_by_order_support(P, num_el_P);
    }
    
    //Broadcasts of the data read and inserted in the array from the process with rank "0" (master) to all other processes (workers) 
    //of the communicator MPI_COMM_WORLD.
    int res = MPI_Bcast(P, num_el_P, item_type, 0, MPI_COMM_WORLD);
	if (res != MPI_SUCCESS) {
		fprintf(stderr, "MPI_Bcast failed\n");
		exit(0);
	}
}