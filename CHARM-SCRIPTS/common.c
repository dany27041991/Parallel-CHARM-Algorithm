#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "common.h"

//We use quicksort to sort in support order.
void sort_by_order_support(struct itemset_tids *P_ptr, int num_el) {
    int j;
    for(int i=0; i < num_el; i++) {
        struct itemset_tids temp = P_ptr[i]; 
        j=i-1; 
        while(j >= 0 && P_ptr[j].tidsets[0] > temp.tidsets[0]){
            P_ptr[j+1]=P_ptr[j];		
            j--;		
        }
        P_ptr[j+1]=temp;
	}
}

//We use quicksort to sort tids in incremental order of the corresponding integer identifier.
void sort_tids_itemsets(int *tids, int n) {
    int i, key, j;
    for (i = 2; i <= n; i++) {
        key = tids[i];
        j = i - 1;
        while (j >= 1 && tids[j] > key) {
            tids[j + 1] = tids[j];
            j = j - 1;
        }
        tids[j + 1] = key;
    }
}

//Function that allows the removal of the struct being analyzed from the array.
void remove_itemsets(int index, struct itemset_tids *P_ptr, struct itemset_tids Xj) {
    int i = index;
    for(int i=index; P_ptr[i].itemsets[0] != 0; i++) {
        memcpy(P_ptr[i].itemsets, P_ptr[i+1].itemsets, (P_ptr[i+1].itemsets[0]+1) * sizeof(int));
        memcpy(P_ptr[i].tidsets, P_ptr[i+1].tidsets, (P_ptr[i+1].tidsets[0]+1) *  sizeof(int));
    }
}

//Function that allows to obtain the union of the items without considering any duplicates.
void union_itemset(const int *Xi, const int *Xj, int *Xij) {
    if(Xi[0] == 1 && Xj[0] == 1 && Xi[1] != Xj[1]) {
        Xij[0] = 2;
        Xij[1] = Xi[1];
        Xij[2] = Xj[1];
    } else {
        int counter = 0;
        for(int i=1; i <= Xi[0]; i++) {
            Xij[counter+1] = Xi[i];
            counter++;
            Xij[0] = counter;
        }
        for(int k=1; k <= Xj[0]; k++) {
            int flag = TRUE;
            for(int j=1; j <= Xij[0]; j++) {
                if(Xij[j] == Xj[k]) {
                    flag = FALSE;
                    break;
                }
            }
            if(flag == TRUE) {
                Xij[counter+1] = Xj[k];
                counter++;
                Xij[0] = counter;
            }
        }
    }
    if(FLAG_SORT_ITEMSETS==TRUE)
        sort_tids_itemsets(Xij, Xij[0]);
}

//Function that allows the replacement of the i-th element in P.
void replace_P(int index, struct itemset_tids *P_ptr, int *Xij) {
    memcpy(P_ptr[index].itemsets, Xij, (Xij[0]+1) * sizeof(int));
}

//Function that allows the replacement of the i-th element in Pi.
void replace_Pi(struct itemset_tids *P_ptr, int *Xij, int n) {
    for(int q=0; q<n; q++) {
        if(P_ptr[q].itemsets[0] == 1 && Xij[0] == 1 && P_ptr[q].itemsets[1] != Xij[1]) {
            P_ptr[q].itemsets[0] = 2;
            P_ptr[q].itemsets[2] = Xij[1];
        } else {
            for(int k=1; k <= Xij[0]; k++) {
                int flag = TRUE;
                for(int j=1; j <= P_ptr[q].itemsets[0]; j++) {
                    if(P_ptr[q].itemsets[j] == Xij[k]) {
                        flag = FALSE;
                        break;
                    }
                }
                if(flag == TRUE) {
                    int index = P_ptr[q].itemsets[0]+1;
                    P_ptr[q].itemsets[index] = Xij[k];
                    P_ptr[q].itemsets[0] = index;
                }
            }
        }
        if(FLAG_SORT_ITEMSETS==TRUE)
            sort_tids_itemsets(P_ptr[q].itemsets, P_ptr[q].itemsets[0]);
    }
}