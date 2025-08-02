#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>

#include "symnmf.h"
#include "matrix.h"
#include "vector.h"

/********************** Functions **********************/
static int symnmfmodule_numpy_array_to_matrix(
    PyArrayObject *matrix_array, 
    struct matrix **output_matrix
) {
    int status_code = STANDARD_ERROR_CODE;
    double *matrix_buffer = NULL;
    size_t row_count = 0;
    size_t column_count = 0;
    struct matrix *temp_matrix = NULL;
    size_t row_index = 0;
    size_t column_index = 0;
    double item_data = 0;

    /* Check the input is a 2D matrix. */
    if (PyArray_NDIM(matrix_array) != 2) {
        goto l_cleanup;
    }

    /* Extract datapoints and insert into vector. */
    row_count = PyArray_DIM(matrix_array, 0);
    column_count = PyArray_DIM(matrix_array, 1);
    matrix_buffer = (double *) PyArray_DATA(matrix_array);

    temp_matrix = MATRIX_allocate(row_count, column_count);
    if (NULL == temp_matrix) {
        goto l_cleanup;
    }  

    /* Copy cells from first matrix form to the other form. */
    for (row_index = 0; row_index < row_count; row_index++) {
        for (column_index = 0; column_index < column_count; column_index++) {
            item_data = matrix_buffer[row_index * column_count + column_index];

            status_code = MATRIX_set_item(temp_matrix, row_index, column_index, item_data);
            if (STANDARD_SUCCESS_CODE != status_code) {
                goto l_cleanup;
            }
        }
    }

    /* Transfer ownership. */
    if (NULL != output_matrix) {
        *output_matrix = temp_matrix;
        temp_matrix = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_matrix) {
        MATRIX_free(temp_matrix);
        temp_matrix = NULL;
    }

    return status_code;
}


static int symnmfmodule_datapoints_to_vector(
    PyArrayObject *datapoints_array, 
    struct vector **output_vector, 
    size_t *output_data_dimension
) {
    int status_code = STANDARD_ERROR_CODE;
    double *datapoints_buffer = NULL;
    size_t num_datapoints = 0;
    size_t data_dimension = 0;
    size_t datapoint_index = 0;
    struct vector *temp_datapoint_vector = NULL;
    double *temp_coords = NULL;
    struct datapoint_item temp_datapoint = {0};
    
    /* Check the input is a 2D matrix. */
    if (PyArray_NDIM(datapoints_array) != 2) {
        goto l_cleanup;
    }

    /* Extract datapoints and insert into vector. */
    num_datapoints = PyArray_DIM(datapoints_array, 0);
    data_dimension = PyArray_DIM(datapoints_array, 1);
    datapoints_buffer = (double*) PyArray_DATA(datapoints_array);

    temp_datapoint_vector = VECTOR_allocate(sizeof(struct datapoint_item), num_datapoints);
    if (NULL == temp_datapoint_vector) {
        goto l_cleanup;
    }  

    for (datapoint_index = 0; datapoint_index < num_datapoints; datapoint_index++) {
        temp_coords = malloc(sizeof(double) * data_dimension);
        if (NULL == temp_coords) {
            goto l_cleanup;
        }

        memcpy(temp_coords, (datapoints_buffer + data_dimension * datapoint_index), sizeof(double) * data_dimension);
        temp_datapoint.coords = temp_coords;

        status_code = VECTOR_add_item(temp_datapoint_vector, &temp_datapoint);
        if (STANDARD_SUCCESS_CODE != status_code) {
            goto l_cleanup;
        }

        /* Ownership transferred. */
        temp_datapoint.coords = NULL;
        temp_coords = NULL;
    }

    /* Transfer ownership. */
    if (NULL != output_vector) {
        *output_vector = temp_datapoint_vector;
        temp_datapoint_vector = NULL;
    }

    if (NULL != output_data_dimension) {
        *output_data_dimension = data_dimension;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != temp_coords) {
        free(temp_coords);
    }

    if (NULL != temp_datapoint_vector) {
        SYMNMF_free_datapoint_vector(temp_datapoint_vector);
        temp_datapoint_vector = NULL;
    }

    return status_code;
}


static PyObject *symnmfmodule_build_output_matrix(struct matrix *result_matrix) {
    npy_intp matrix_dimensions[2] = {
        (npy_intp) result_matrix->row_count,
        (npy_intp) result_matrix->column_count
    };

    /* Create numpy array wrapping the existing C buffer. */
    PyObject *matrix_array = PyArray_SimpleNewFromData(
        2,                           /* ndim */ 
        matrix_dimensions,           /* shape */
        NPY_DOUBLE,                  /* dtype */
        (void*)result_matrix->data   /* data pointer */
    );
    if (NULL == matrix_array) {
        return NULL;
    }

    /* Notify numpy that is owns the data buffer. */
    PyArray_ENABLEFLAGS((PyArrayObject *) matrix_array, NPY_ARRAY_OWNDATA);

    return matrix_array;
}


/**
 * Python entryoint to execute SymNMF algorithm.
 * 
 * @param self Python caller information.
 * @param args Arguments passed to the Python native call.
 * @return The output Python object, or NULL if an error occurred.
 */
static PyObject* symnmfmodule_symnmf(PyObject *self, PyObject *args) {
    int status_code = STANDARD_ERROR_CODE;
    PyObject *initial_h_matrix_obj = NULL;
    PyObject *w_matrix_obj = NULL;
    PyArrayObject *initial_h_matrix_array = NULL;
    PyArrayObject *w_matrix_array = NULL;
    struct matrix *initial_h_matrix = NULL;
    struct matrix *w_matrix = NULL;
    struct matrix *result_matrix = NULL;
    PyObject *result_obj = NULL;

    if (!PyArg_ParseTuple(args, "OO", &initial_h_matrix_obj, &w_matrix_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Validate inputs. */
    /* Check the objects are numpy arrays. */
    initial_h_matrix_array = (PyArrayObject*) PyArray_FROM_OTF(initial_h_matrix_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (initial_h_matrix_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }
    
    w_matrix_array = (PyArrayObject*) PyArray_FROM_OTF(w_matrix_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (w_matrix_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }
    
    /* Convert input to matrix datastructure. */
    status_code = symnmfmodule_numpy_array_to_matrix(initial_h_matrix_array, &initial_h_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to convert data structures.");
        goto l_cleanup;
    }

    status_code = symnmfmodule_numpy_array_to_matrix(w_matrix_array, &w_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to convert data structures.");
        goto l_cleanup;
    }

    /* Execute algorithm to calculate diagonal degree matrix. */
    status_code = SYMNMF_symnmf_algorithm(initial_h_matrix, w_matrix, &result_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to perform SymNMF algorithm.");
        goto l_cleanup;
    }

    /* Convert the result to numpy format.
     * Note that numpy will own the underlying buffer.
     */
    result_obj = symnmfmodule_build_output_matrix(result_matrix);
    if (NULL == result_obj) {
        PyErr_SetString(PyExc_TypeError, "Failed to build output matrix.");
        goto l_cleanup;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != w_matrix_array) {
        Py_DECREF(w_matrix_array);
        w_matrix_array = NULL;
    }

    if (NULL != initial_h_matrix) {
        MATRIX_free(initial_h_matrix);
        initial_h_matrix = NULL;
    }

    if (NULL != w_matrix) {
        MATRIX_free(w_matrix);
        w_matrix = NULL;
    }

    if (NULL == result_obj) {
        /* Numpy does not yet own the buffer, safe to release. */
        if (NULL != result_matrix) {
            MATRIX_free(result_matrix);
            result_matrix = NULL;
        }
    } else {
        /* Numpy owns the buffer, only release wrapper object. */
        if (NULL != result_matrix) {
            free(result_matrix);
            result_matrix = NULL;
        }
    }

    return result_obj;
}


/**
 * Python entryoint to calculate similarity matrix.
 * 
 * @param self Python caller information.
 * @param args Arguments passed to the Python native call.
 * @return The output Python object, or NULL if an error occurred.
 */
static PyObject* symnmfmodule_sym(PyObject *self, PyObject *args) {
    int status_code = STANDARD_ERROR_CODE;
    PyObject *datapoints_obj = NULL;
    PyArrayObject *datapoints_array = NULL;
    struct vector *datapoints_vector = NULL;
    size_t data_dimension = 0;
    struct matrix *result_matrix = NULL;
    PyObject *result_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &datapoints_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Validate inputs. */
    /* Check the objects are numpy arrays. */
    datapoints_array = (PyArrayObject*) PyArray_FROM_OTF(datapoints_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (datapoints_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Convert datapoints format to vector. */
    status_code = symnmfmodule_datapoints_to_vector(datapoints_array, &datapoints_vector, &data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse datapoints.");
        goto l_cleanup;
    }

    /* Execute algorithm to calculate diagonal degree matrix. */
    status_code = SYMNMF_calc_similarity_matrix(datapoints_vector, data_dimension, &result_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to calculate matrix.");
        goto l_cleanup;
    }

    /* Convert the result to numpy format.
     * Note that numpy will own the underlying buffer.
     */
    result_obj = symnmfmodule_build_output_matrix(result_matrix);
    if (NULL == result_obj) {
        PyErr_SetString(PyExc_TypeError, "Failed to build output matrix.");
        goto l_cleanup;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != datapoints_array) {
        Py_DECREF(datapoints_array);
        datapoints_array = NULL;
    }

    if (NULL != datapoints_vector) {
        SYMNMF_free_datapoint_vector(datapoints_vector);
        datapoints_vector = NULL;
    }

    if (NULL == result_obj) {
        /* Numpy does not yet own the buffer, safe to release. */
        if (NULL != result_matrix) {
            MATRIX_free(result_matrix);
            result_matrix = NULL;
        }
    } else {
        /* Numpy owns the buffer, only release wrapper object. */
        if (NULL != result_matrix) {
            free(result_matrix);
            result_matrix = NULL;
        }
    }

    return result_obj;
}


/**
 * Python entryoint to calculate diagonal degree matrix.
 * 
 * @param self Python caller information.
 * @param args Arguments passed to the Python native call.
 * @return The output Python object, or NULL if an error occurred.
 */
static PyObject* symnmfmodule_ddg(PyObject *self, PyObject *args) {
    int status_code = STANDARD_ERROR_CODE;
    PyObject *datapoints_obj = NULL;
    PyArrayObject *datapoints_array = NULL;
    struct vector *datapoints_vector = NULL;
    size_t data_dimension = 0;
    struct matrix *result_matrix = NULL;
    PyObject *result_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &datapoints_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Validate inputs. */
    /* Check the objects are numpy arrays. */
    datapoints_array = (PyArrayObject*) PyArray_FROM_OTF(datapoints_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (datapoints_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Convert datapoints format to vector. */
    status_code = symnmfmodule_datapoints_to_vector(datapoints_array, &datapoints_vector, &data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse datapoints.");
        goto l_cleanup;
    }

    /* Execute algorithm to calculate diagonal degree matrix. */
    status_code = SYMNMF_calc_diagonal_degree_matrix(datapoints_vector, data_dimension, &result_matrix, NULL);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to calculate matrix.");
        goto l_cleanup;
    }

    /* Convert the result to numpy format.
     * Note that numpy will own the underlying buffer.
     */
    result_obj = symnmfmodule_build_output_matrix(result_matrix);
    if (NULL == result_obj) {
        PyErr_SetString(PyExc_TypeError, "Failed to build output matrix.");
        goto l_cleanup;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != datapoints_array) {
        Py_DECREF(datapoints_array);
        datapoints_array = NULL;
    }

    if (NULL != datapoints_vector) {
        SYMNMF_free_datapoint_vector(datapoints_vector);
        datapoints_vector = NULL;
    }

    if (NULL == result_obj) {
        /* Numpy does not yet own the buffer, safe to release. */
        if (NULL != result_matrix) {
            MATRIX_free(result_matrix);
            result_matrix = NULL;
        }
    } else {
        /* Numpy owns the buffer, only release wrapper object. */
        if (NULL != result_matrix) {
            free(result_matrix);
            result_matrix = NULL;
        }
    }

    return result_obj;
}


/**
 * Python entryoint to calculate normalized similarity matrix.
 * 
 * @param self Python caller information.
 * @param args Arguments passed to the Python native call.
 * @return The output Python object, or NULL if an error occurred.
 */
static PyObject* symnmfmodule_norm(PyObject *self, PyObject *args) {
    int status_code = STANDARD_ERROR_CODE;
    PyObject *datapoints_obj = NULL;
    PyArrayObject *datapoints_array = NULL;
    struct vector *datapoints_vector = NULL;
    size_t data_dimension = 0;
    struct matrix *result_matrix = NULL;
    PyObject *result_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &datapoints_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Validate inputs. */
    /* Check the objects are numpy arrays. */
    datapoints_array = (PyArrayObject*) PyArray_FROM_OTF(datapoints_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (datapoints_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Convert datapoints format to vector. */
    status_code = symnmfmodule_datapoints_to_vector(datapoints_array, &datapoints_vector, &data_dimension);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to parse datapoints.");
        goto l_cleanup;
    }

    /* Execute algorithm to calculate diagonal degree matrix. */
    status_code = SYMNMF_calc_normalized_similarity_matrix(datapoints_vector, data_dimension, &result_matrix);
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_TypeError, "Failed to calculate matrix.");
        goto l_cleanup;
    }

    /* Convert the result to numpy format.
     * Note that numpy will own the underlying buffer.
     */
    result_obj = symnmfmodule_build_output_matrix(result_matrix);
    if (NULL == result_obj) {
        PyErr_SetString(PyExc_TypeError, "Failed to build output matrix.");
        goto l_cleanup;
    }
    
    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (NULL != datapoints_array) {
        Py_DECREF(datapoints_array);
        datapoints_array = NULL;
    }

    if (NULL != datapoints_vector) {
        SYMNMF_free_datapoint_vector(datapoints_vector);
        datapoints_vector = NULL;
    }

    if (NULL == result_obj) {
        /* Numpy does not yet own the buffer, safe to release. */
        if (NULL != result_matrix) {
            MATRIX_free(result_matrix);
            result_matrix = NULL;
        }
    } else {
        /* Numpy owns the buffer, only release wrapper object. */
        if (NULL != result_matrix) {
            free(result_matrix);
            result_matrix = NULL;
        }
    }

    return result_obj;
}


/********************** Structures **********************/
static PyMethodDef symnmf_methods[] = {
    {
      "symnmf",                                      /* the Python method name that will be used */
      (PyCFunction) symnmfmodule_symnmf,             /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,                                  /* flags indicating parameters accepted for this function */
      PyDoc_STR("Perform SymNMF and output H matrix. Expected arguments: (initial_h_matrix: np.ndarray, w_matrix: np.ndarray)")
    },
    {
      "sym",                                         /* the Python method name that will be used */
      (PyCFunction) symnmfmodule_sym,                /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,                                  /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculate and output the similarity matrix. Expected arguments: (datapoints: np.ndarray)")
    },
    {
      "ddg",                                         /* the Python method name that will be used */
      (PyCFunction) symnmfmodule_ddg,                /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,                                  /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculate and output the diagonal degree matrix. Expected arguments: (datapoints: np.ndarray)")
    },
    {
      "norm",                                        /* the Python method name that will be used */
      (PyCFunction) symnmfmodule_norm,               /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,                                  /* flags indicating parameters accepted for this function */
      PyDoc_STR("Calculate and output the normalized similarity matrix. Expected arguments: (datapoints: np.ndarray)")
    },
    {NULL, NULL, 0, NULL}                            /* The last entry must be all NULL as shown to act as a
                                                        sentinel. Python looks for this entry to know that all
                                                        of the functions for the module have been defined. */
};

static struct PyModuleDef symnmf_module = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",              /* name of module */
    NULL,                        /* module documentation, may be NULL */
    -1,                          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    symnmf_methods               /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    
    /* Initialize numpy internal state. */
    import_array();

    m = PyModule_Create(&symnmf_module);
    if (!m) {
        return NULL;
    }
    return m;
}
