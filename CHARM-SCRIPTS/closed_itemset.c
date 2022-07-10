#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "common.h"
#include "closed_itemset.h"

//Function that allows the removal of any subsets after receipt of FCI from the different processes.
int get_frequent_closed_itemsets(struct itemset_tids *frequentClosedItemset, struct itemset_tids *closedItemsetProcess, int num_el_CI) {
    int counter = 0;
    if(FLAG_WRITE_FCI==TRUE) {
        int *array_index = malloc(num_el_CI*sizeof(int));
        for (int i = 0; i < num_el_CI; i++) {
            int index = closedItemsetProcess[i].index;
            if(frequentClosedItemset[index].index == 0) {          
                memcpy(frequentClosedItemset[index].itemsets, closedItemsetProcess[i].itemsets, (closedItemsetProcess[i].itemsets[0]+1) * sizeof(int));
                memcpy(frequentClosedItemset[index].tidsets, closedItemsetProcess[i].tidsets, (closedItemsetProcess[i].tidsets[0]+1) * sizeof(int));
                frequentClosedItemset[index].index = index;
                array_index[counter] = index;
                counter++;
            } else {
                if(frequentClosedItemset[index].itemsets[0] < closedItemsetProcess[i].itemsets[0]) {
                    memcpy(frequentClosedItemset[index].itemsets, closedItemsetProcess[i].itemsets, (closedItemsetProcess[i].itemsets[0]+1) * sizeof(int));
                }
            }
        }
        FILE *f = fopen("frequent_closed_itemsets.txt", "w");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }

        for (int i = 0; i < counter; i++) {
            int index = array_index[i];
            fprintf(f, "ITEMSET: ");
            for(int j = 1; j<=frequentClosedItemset[index].itemsets[0]; j++) {
                fprintf(f, "%d ", frequentClosedItemset[index].itemsets[j]);
            }
            fprintf(f, "--> TIDS: ");
            
            for(int k=1; k<=frequentClosedItemset[index].tidsets[0]; k++) {
                fprintf(f, "%d ", frequentClosedItemset[index].tidsets[k]);
            }
            fprintf(f, "\n");
        }
        fclose(f);
        free(array_index);
    } else {
        for (int i = 0; i < num_el_CI; i++) {
            int index = closedItemsetProcess[i].index;
            if(frequentClosedItemset[index].index == 0) {          
                memcpy(frequentClosedItemset[index].itemsets, closedItemsetProcess[i].itemsets, (closedItemsetProcess[i].itemsets[0]+1) * sizeof(int));
                memcpy(frequentClosedItemset[index].tidsets, closedItemsetProcess[i].tidsets, (closedItemsetProcess[i].tidsets[0]+1) * sizeof(int));
                frequentClosedItemset[index].index = index;
                counter++;
            } else {
                if(frequentClosedItemset[index].itemsets[0] < closedItemsetProcess[i].itemsets[0]) {
                    memcpy(frequentClosedItemset[index].itemsets, closedItemsetProcess[i].itemsets, (closedItemsetProcess[i].itemsets[0]+1) * sizeof(int));
                }
            }
        }
    }
    return counter;
}