#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdlib.h>

/********************** Constants ******************************/
#define STANDARD_ERROR_CODE (-1)
#define STANDARD_SUCCESS_CODE (0)

/********************** Structures *****************************/
/* A dynamically increasing vector for data items of same size. 
 * Note that only adding/setting/clearing items is supported, not popping.
*/
struct vector {
    /* Allocated buffer. */
    void *data;

    /* Size of items in the vector. */
    size_t item_size;

    /* Allocated size in items. */
    size_t item_capacity;
    
    /* Occupied size in items. */
    size_t item_count;
};

/********************** Functions ******************************/
/**
 * Allocate a dynamically resizing vector of same-size items.
 * 
 * @param item_size The size of the items in vector, in bytes.
 * @param reserve_items The initial size to reserve for the vector.
 * @return The allocated vector or NULL if an error occurred. 
 *  The vector must be passed to VECTOR_free() to release memory. 
 */
struct vector *VECTOR_allocate(size_t item_size, size_t reserve_items);

/**
 * Retrieve from the vector the item data at the given index. 
 * 
 * @param vector The vector object.
 * @param item_index The index of the item.
 * @param item_buffer Buffer to which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int VECTOR_get_item(struct vector *vector, size_t item_index, void *item_buffer);

/**
 * Set the item data in the vector at the given index. 
 * 
 * @param vector The vector object.
 * @param item_index The index of the item.
 * @param item_buffer Buffer from which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int VECTOR_set_item(struct vector *vector, size_t item_index, void *item_data);

/**
 * Append an item with given data to the vector.
 * If no capacity is left in the vector, more space will be allocated. 
 * 
 * @param vector The vector object.
 * @param item_buffer Buffer from which the item data will be copied. Must be of size >= vector.item_size.
 * @return The status of the operation. STANDARD_SUCCESS_CODE for success, STANDARD_ERROR_CODE for error.
 */
int VECTOR_add_item(struct vector *vector, void *item_data);

/**
 * Release the vector object and its memory in use.
 * Note that if the items theselves point to allocated memory, their allocations will not be released. 
 * 
 * @param vector The vector object.
 */
void VECTOR_free(struct vector *vector);

#endif /* __VECTOR_H__ */
