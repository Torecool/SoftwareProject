#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>

#include "vector.h"

/********************** Functions ******************************/
struct vector *VECTOR_allocate(size_t item_size, size_t reserve_items) {
    struct vector *temp_vector = NULL;

    if ((0 == item_size) || (0 == reserve_items)) {
        return NULL;
    }

    temp_vector = malloc(sizeof(struct vector));
    if (NULL == temp_vector) {
        return NULL;
    }

    temp_vector->data = malloc(item_size * reserve_items);
    if (NULL == temp_vector->data) {
        free(temp_vector);
        return NULL;
    }

    temp_vector->item_size = item_size;
    temp_vector->item_capacity = reserve_items;
    temp_vector->item_count = 0;

    return temp_vector;
}


int VECTOR_get_item(struct vector *vector, size_t item_index, void *item_buffer) {
    int status_code = STANDARD_ERROR_CODE;
    
    if (vector->item_count <= item_index) {
        /* Index out of range of initialized items. */
        goto l_cleanup;
    } else {
        memcpy(item_buffer, ((unsigned char *) vector->data) + (item_index * vector->item_size), vector->item_size);
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int VECTOR_set_item(struct vector *vector, size_t item_index, void *item_data) {
    int status_code = STANDARD_ERROR_CODE;
    
    if (vector->item_count <= item_index) {
        /* Index out of range of initialized items. */
        goto l_cleanup;
    } else {
        memcpy(((unsigned char *) vector->data) + (item_index * vector->item_size), item_data, vector->item_size);
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


int VECTOR_add_item(struct vector *vector, void *item_data) {
    int status_code = STANDARD_ERROR_CODE;

    if (vector->item_count == vector->item_capacity) {
        vector->data = realloc(vector->data, vector->item_size * (vector->item_capacity * 2));
        if (NULL == vector->data) {
            goto l_cleanup;
        }

        vector->item_capacity *= 2;
    } 

    memcpy(((unsigned char *) vector->data) + (vector->item_count * vector->item_size), item_data, vector->item_size);
    vector->item_count++;

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    return status_code;
}


void VECTOR_free(struct vector *vector) {
    if (NULL == vector) {
        return;
    }

    free(vector->data);
    free(vector);
}
