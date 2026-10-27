#include <stdlib.h>
#include <stdio.h>

void* u_free(void* ptr, const char* file, int line, const char* name) {
    printf("FREEING: %s\n", name);
    free(ptr);
    return NULL;
}

/* OLD buggy version */
#define SFREE_OLD(a) u_free(a,__FILE__,__LINE__,#a)

/* FIXED version */
#define SFREE_FIXED(a) ((a)=u_free(a,__FILE__,__LINE__,#a))

int main() {
    int* ptr1 = malloc(sizeof(int) * 10);
    int* ptr2 = malloc(sizeof(int) * 10);
    
    printf("ptr1 address before free: %p\n", (void*)ptr1);
    
    printf("\n=== Testing OLD SFREE macro ===\n");
    SFREE_OLD(ptr1);
    printf("ptr1 address after free: %p\n", (void*)ptr1);
    
    if (ptr1 != NULL) {
        printf("BUG: ptr1 is not NULL after SFREE!\n");
    }
    
    printf("\n=== Testing FIXED SFREE macro ===\n");
    printf("ptr2 address before free: %p\n", (void*)ptr2);
    SFREE_FIXED(ptr2);
    printf("ptr2 address after free: %p\n", (void*)ptr2);
    
    if (ptr2 == NULL) {
        printf("FIXED: ptr2 is NULL after SFREE!\n");
    }
    
    return 0;
}
