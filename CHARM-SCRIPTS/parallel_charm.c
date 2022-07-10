#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mpi.h"
#include "common.h"
#include "broadcast_dataset.h"
#include "closed_itemset.h"

//TO COMPILE: mpicc -O3 parallel_charm.c common.c broadcast_dataset.c closed_itemset.c -o parallel_charm
//TO RUN: mpirun -np 4 parallel_charm dataset.txt 5000 2

/**
 * Structure of a item:
 * - itemset: char[100]
 * - tids: char[100]
 * 
 *
 * How to represent such a structure with an MPI struct:
 *           |            |            |+----- displacement for
 *           |            |            |     block 3: 0
 *           |            |            |    (+ potential padding)
 *           |            |            |            |
 *           +----- displacement for   |            |
 *           |       block 2: 0        |            |
 *           |   (+ potential padding) |            |
 *           |            |            |            |
 *  displacement for      |            |            |
 *    block 1: 0          |            |            |
 * (+ potential padding)  |            |            |
 *           |            |            |            |
 *           V            V            V            V
 *           +------------+------------+------------+
 *           |  itemsets  |  tidsets   |    index   |
 *           +------------+------------+------------+
 *           <-----------> <----------> <----------->
 *               block 1      block 2       block 3
 *            n MPI_INT     n MPI_INT      MPI_INT
 **/

int rank, processes, counter_cfi = 0, NUM_EL_P, MIN_SUP;
double start, end; //Time start-end
char *NAMEFILE;

/*int hash_function(struct itemset_tids itemset) {
    int count = 0;
    int support = itemset.tidsets[0];
    for(int i=1; i<=support; i++) {
        count=count+itemset.tidsets[i];
    }
    return count % FCI_SIZE_ARRAY;
}*/

int hash_function(struct itemset_tids itemset) {
    int support = itemset.tidsets[0];
    char *str = calloc(STRING_HASH_SIZE, sizeof(char));
    int index = 0;
    for(int i=0; i<=support; i++) {
        //printf("%d - %d\n", index, itemset.tidsets[i]);
        index += sprintf(&str[index], "%d", itemset.tidsets[i]);
    }
    //uint64_t hashval = 5381;
    unsigned hashval = 0;
    int k = 0;
    while (str[k] != '\0') {
        //hashval = str[k] + 33*hashval;
        hashval = ((hashval << 5) + hashval) + str[k];
        //hashval = str[k] +(hashval << 6) + (hashval <<16) -hashval;
        k++;
    }
    free(str);
    //printf("%d\n", hashval % FCI_SIZE_ARRAY);
    return (hashval % FCI_SIZE_ARRAY);
}

//CHARM algorithm
int charm(struct itemset_tids *P_ptr, int minsup, struct itemset_tids *C, int flag_starting_sort, int flag_starting_index, int num_el_P, int start_index) {
    int *Xij = (int*) malloc(ITEMSETS_SIZE*sizeof(int));
    int *tij = (int*) malloc(TIDSETS_SIZE*sizeof(int));
    struct itemset_tids *Pi = malloc(num_el_P*sizeof(struct itemset_tids));
    //Check the sort in support order. The first time workers CHARM this check is not performed.
    if(flag_starting_sort == TRUE) {
        sort_by_order_support(P_ptr, num_el_P);
    } 
    int i = start_index;
    //Cycle to the final index identified by the master process. This for loop is the parallelized one.
    while(i < num_el_P) { 
        int j = i+1, num_el_Pi = 0;
        //Internal for loop on all elements following the one in use by the external for.
        while(j < num_el_P) {
            //I create the union set of the itemset relative to the external and internal for using the union_itemset method. 
            //This union does not have duplicate items.
            union_itemset(P_ptr[i].itemsets, P_ptr[j].itemsets, Xij);
            //I create the intersection set of the relative tids relative to the itemset of the external and internal for.
            int dis_i = 0, dis_j = 0, c_i=1, c_j=1, common_element_counter = 1, c1, c2, val1, val2;
            while(P_ptr[i].tidsets[0] >= c_i && P_ptr[j].tidsets[0] >= c_j){
                c1 = P_ptr[i].tidsets[0];
                c2 = P_ptr[j].tidsets[0];
                val1 = P_ptr[i].tidsets[c_i];
                val2 = P_ptr[j].tidsets[c_j];
                if (P_ptr[i].tidsets[c_i] < P_ptr[j].tidsets[c_j]) { 
                    c_i++;
                    dis_i++;
                }else if(P_ptr[j].tidsets[c_j] < P_ptr[i].tidsets[c_i]){
                    c_j++;
                    dis_j++;
                } else {
                    tij[common_element_counter] = P_ptr[i].tidsets[c_i];
                    c_i++;
                    c_j++;
                    common_element_counter++;
                }
                if(c_j > P_ptr[j].tidsets[0] && P_ptr[i].tidsets[0] >= c_i) {
                    dis_i++;
                }
                if(c_i > P_ptr[i].tidsets[0] && P_ptr[j].tidsets[0] >= c_j) {
                    dis_j++;
                }
            }
            tij[0] = --common_element_counter;
            //I enter the if only if the number of tids of the intersection set is greater than MIN_SUPPORT. The total number of tids is stored in index 0 for easier code management.
            if(tij[0] >= minsup) {
                //Property 1
                //I check if the tids relative to the element under analysis of the internal and external for are equal. Use the methods tids_are_equals
                if(dis_i == 0 && dis_j == 0) {
                    //I perform the substitution on P and Pi and finally the removal on P.
                    replace_P(i, P_ptr, Xij);
                    replace_Pi(Pi, Xij, num_el_Pi);
                    remove_itemsets(j, P_ptr, P_ptr[j]);
                    j--;
                } else {
                    //Property 2
                    //I check if the tids relative to the element under analysis of the external for are a subset of the internal one. Use the methods tids_are_subset.
                    if(dis_i == 0 && dis_j > 0) {
                        //I perform the substitution on P and Pi.
                        replace_P(i, P_ptr, Xij);
                        replace_Pi(Pi, Xij, num_el_Pi);
                    } else { //Property 3 
                        //I insert the element into the Pi array.
                        memcpy(Pi[num_el_Pi].itemsets, Xij, (Xij[0]+1) * sizeof(int));
                        memcpy(Pi[num_el_Pi].tidsets, tij, (tij[0]+1) * sizeof(int));
                        num_el_Pi++;
                    }
                }
            }
            j++;
        }

        //If I have elements in Pi I run CHARM recursively.
        if(num_el_Pi > 0) 
            charm(Pi, minsup, C, 1, 1, num_el_Pi, 0);
        //I check if the possible frequent closed itemset generated are a subset of some element present in the array containing the FCIs,
        //we also carry out the verification on the support.
        if(P_ptr[i].tidsets[0] >= minsup) {
            memcpy(C[counter_cfi].itemsets, P_ptr[i].itemsets, (P_ptr[i].itemsets[0]+1) * sizeof(int));
            memcpy(C[counter_cfi].tidsets, P_ptr[i].tidsets, (P_ptr[i].tidsets[0]+1) * sizeof(int));
            C[counter_cfi].index = hash_function(C[counter_cfi]);
            counter_cfi++;
        }
        if(flag_starting_index == TRUE) {
            i++;
        } else {
            i=i+processes;
        }
    }
    free(Xij);
    free(tij);
    free(Pi);
    return counter_cfi;
}

int main (int argc, char* argv[])
{
    NAMEFILE = argv[1];
    NUM_EL_P = atoi(argv[2]);
    MIN_SUP = atoi(argv[3]);

    //Initialize the MPI execution environment
	MPI_Init(&argc, &argv);
    //Determines the rank of the calling process in the communicator
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //Determines the size of the group associated with a communicator
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
    //Create a new type of data.
    MPI_Datatype item_type;
    //This new type of data contains two elements of fixed lenght => two int arrays.
	int lengths[3] = {ITEMSETS_SIZE, TIDSETS_SIZE, 1};
    //MPI_Aint is a portable C data type that can hold memory addresses and it could be larger than the usual int.
    //Used to keep track of the movements between one element and the next.
	MPI_Aint displacements[3];
	struct itemset_tids item;
	MPI_Aint base_address;
    //Get the address of a location in memory.
	MPI_Get_address(&item, &base_address);
	MPI_Get_address(&item.itemsets[0], &displacements[0]);
	MPI_Get_address(&item.tidsets[0], &displacements[1]);
    MPI_Get_address(&item.index, &displacements[2]);
    //Returns the difference between addr1 and addr2.
	displacements[0] = MPI_Aint_diff(displacements[0], base_address);
	displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    //Our structs contains two arrays of integers.
	MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};
    //Create an MPI datatype from a general set of datatypes, displacements, and block sizes.
	MPI_Type_create_struct(3, lengths, displacements, types, &item_type);
    //Commits the datatype to render it effective.
	MPI_Type_commit(&item_type);

    struct itemset_tids *P = calloc(NUM_EL_P, sizeof(struct itemset_tids));
    struct itemset_tids *closedItemset = calloc(CI_SIZE_ARRAY, sizeof(struct itemset_tids));
    struct itemset_tids *frequentClosedItemset = calloc(FCI_SIZE_ARRAY, sizeof(struct itemset_tids));
    //I run the function bcast_dataset which allows the broadcasting of data read from files to all processors after having sorted them by order of support.
    bcast_dataset(P, NAMEFILE, item_type, rank, processes, NUM_EL_P);
    struct itemset_tids *buffer_received = calloc(processes * CI_SIZE_ARRAY, sizeof(struct itemset_tids));
    //Blocks until all processes in the communicator have reached this routine.
    MPI_Barrier(MPI_COMM_WORLD); 
    //Get the start elapsed time on the calling processor.
	start = MPI_Wtime();
    //If the number of processes is 1 I run CHARM sequentially, otherwise it runs in parallel.
    if (processes < 2) {
        int counter_ci = charm(P, MIN_SUP, closedItemset, 0, 1, NUM_EL_P, rank);
        int num_el_FCI = get_frequent_closed_itemsets(frequentClosedItemset, closedItemset, counter_ci);
        printf("The Number of FCI is %d.\n", num_el_FCI);
    } else {
        int num_CI_send = charm(P, MIN_SUP, closedItemset, 0, 0, NUM_EL_P, rank);
        int *num_CI_recv = malloc(processes*sizeof(int));
        MPI_Gather(&num_CI_send, 1, MPI_INT, num_CI_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int tot_CI_el = 0;
        int *displ = malloc(processes*sizeof(int));
        if(rank==0) {
            for (int i = 0; i < processes; i++)
            {
                displ[i] = tot_CI_el;
                tot_CI_el = tot_CI_el + num_CI_recv[i];
            }
        }
        struct itemset_tids *CI_temp = calloc(tot_CI_el, sizeof(struct itemset_tids));
        MPI_Gatherv(closedItemset, num_CI_send, item_type, CI_temp, num_CI_recv, displ, item_type, 0, MPI_COMM_WORLD);
        
        free(num_CI_recv);
        if(rank==0) {
            int num_el_FCI = get_frequent_closed_itemsets(frequentClosedItemset, CI_temp, tot_CI_el);
            printf("The Number of FCI is %d.\n", num_el_FCI);
        }   
        free(CI_temp);
    }

    fflush(stdout);
    end = MPI_Wtime();
    if(rank==0)
        printf("The execution time on %d processes to generate the FCI is %f.\n", processes, end - start);
    
    free(frequentClosedItemset);
    free(closedItemset);
    free(P);
    MPI_Finalize();
	return 0;
}

//Mpi get put accumulate