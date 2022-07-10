/************************* MACROS **************************/
#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id) + 1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) \
	(BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n) + 1)
#define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))
#define TRUE 1
#define FALSE 0

/************************* CONSTANTS **************************/
#ifndef _COMMON_CONST
#define _COMMON_CONST
static const int CI_SIZE_ARRAY = 500000;
static const int FCI_SIZE_ARRAY = 974775931;
static const int TIDSETS_SIZE = 100;
static const int ITEMSETS_SIZE = 20;
static const int STRING_HASH_SIZE = 1000;
static const int FLAG_SORT_TIDS = TRUE;
static const int FLAG_SORT_ITEMSETS = TRUE;
static const int FLAG_WRITE_FCI = TRUE;
#endif

/************************* STRUCT **************************/
struct itemset_tids {
    int itemsets[ITEMSETS_SIZE];
    int tidsets[TIDSETS_SIZE];
    int index;
};

/************************* COMMON FUNCTIONS **************************/
void sort_tids_itemsets(int *tids, int);
void sort_by_order_support(struct itemset_tids *P_ptr, int num_el);
void remove_itemsets(int index, struct itemset_tids *P_ptr, struct itemset_tids Xj);
void union_itemset(const int *Xi, const int *Xj, int *Xij);
void replace_P(int index, struct itemset_tids *P_ptr, int *Xij);
void replace_Pi(struct itemset_tids *P_ptr, int *Xij, int);
