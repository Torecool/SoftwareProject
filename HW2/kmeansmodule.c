#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "vector.h"
#include "kmeans.h"


/********************** Functions **********************/
static Py_ssize_t parse_datapoint_dimension(PyObject* datapoints_coords) {
    Py_ssize_t ret_val = STANDARD_ERROR_CODE;
    
    /* Borrowed reference. */
    PyObject* datapoint = PyList_GetItem(datapoints_coords, 0);
    if (!PyTuple_Check(datapoint)) {
        PyErr_SetString(PyExc_TypeError, "Each item in datapoints must be a tuple.");
        goto l_cleanup;
    }

    ret_val = PyTuple_Size(datapoint);
l_cleanup:
    return ret_val;
}


static int add_datapoints_to_vector(struct vector *datapoints_vector, PyObject* datapoints_coords) {
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    double *temp_data = NULL;

    Py_ssize_t num_datapoints = PyList_Size(datapoints_coords);

    /* Load all datapoints. */
    for (Py_ssize_t index = 0; index < num_datapoints; index++) {
        /* Borrowed reference. */
        PyObject* datapoint = PyList_GetItem(datapoints_coords, index); 

        if (!PyTuple_Check(datapoint)) {
            PyErr_SetString(PyExc_TypeError, "Each item in datapoints must be a tuple.");
            goto l_cleanup;
        }

        Py_ssize_t data_dimension = PyTuple_Size(datapoint);

        temp_data = malloc(sizeof(double) * data_dimension);
        if (NULL == temp_data) {
            goto l_cleanup;
        }

        for (Py_ssize_t j = 0; j < data_dimension; ++j) {
            /* Borrowed reference. */
            PyObject* coord = PyTuple_GetItem(datapoint, j); 
            temp_data[j] = PyFloat_AsDouble(coord);
            if (PyErr_Occurred()) {
                /* Conversion failed. */
                goto l_cleanup;
            }
        }

        /* Insert the datapoint into the datastructure. */
        temp_datapoint_item.coords = temp_data;
        status_code = VECTOR_add_item(datapoints_vector, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            PyErr_SetString(PyExc_ValueError, "Failed to add datapoint to vector.");
            goto l_cleanup;
        }
        temp_data = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    free(temp_data);
    return status_code;
}


static PyObject* build_output_centroids(struct vector *clusters, size_t data_dimension) {
    struct cluster_item temp_cluster_item = {0};
    int status_code = STANDARD_ERROR_CODE;
    PyObject* centroid = NULL;
    
    PyObject* centroid_list = PyList_New(clusters->item_count);
    if (NULL == centroid_list) {
        PyErr_SetString(PyExc_ValueError, "Failed to allocate output centroid list.");
        goto l_cleanup;
    }

    for (size_t cluster_index = 0; cluster_index < clusters->item_count; cluster_index++) {
        status_code = VECTOR_get_item(clusters, cluster_index, &temp_cluster_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            PyErr_SetString(PyExc_ValueError, "Failed to read centroid data.");
            goto l_cleanup;
        }

        centroid = PyTuple_New(data_dimension);
        if (NULL == centroid) {
            PyErr_SetString(PyExc_ValueError, "Failed to allocate output centroid.");
            goto l_cleanup;
        }

        for (size_t dimension_index = 0; dimension_index < data_dimension; dimension_index++) {
            
            /* Note: steals reference. */
            PyTuple_SetItem(centroid, dimension_index, Py_BuildValue("f", temp_cluster_item.centroid[dimension_index])); 
        }

        /* Note: steals reference. */
        PyList_SetItem(centroid_list, cluster_index, centroid); 
        centroid = NULL;
    }

    status_code = STANDARD_SUCCESS_CODE;
l_cleanup:
    if (STANDARD_SUCCESS_CODE != status_code && NULL != centroid_list) {
        Py_DECREF(centroid_list);
        centroid_list = NULL;
    }

    return centroid_list;
}


static PyObject* kmeans_fit(PyObject *self, PyObject *args) {
    PyObject *centroids_coords = NULL;
    PyObject *datapoints_coords = NULL;
    size_t num_iterations = 0;
    double convergence_threshold = 0;
    PyObject *output_centroids = NULL;
    int status_code = STANDARD_ERROR_CODE;
    struct datapoint_item temp_datapoint_item = {0};
    struct vector *clusters_vector = NULL;
    struct vector *datapoints_vector = NULL;

    if (!PyArg_ParseTuple(args, "OOnd", &centroids_coords, &datapoints_coords, &num_iterations, &convergence_threshold)) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    /* Validate inputs. */
    if (!PyList_Check(centroids_coords) || !PyList_Check(datapoints_coords) || 
        num_iterations <= MINIMUM_VALID_NUM_ITERATIONS || num_iterations >= MAXIMUM_VALID_NUM_ITERATIONS ||
        convergence_threshold < 0) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    size_t num_centroids = PyList_Size(centroids_coords);
    size_t num_datapoints = PyList_Size(datapoints_coords);

    if (num_centroids >= num_datapoints || num_centroids <= MINIMUM_VALID_K_VALUE) {
        PyErr_SetString(PyExc_TypeError, "Invalid parameters.");
        goto l_cleanup;
    }

    int data_dimension = parse_datapoint_dimension(datapoints_coords);
    if (STANDARD_ERROR_CODE == data_dimension) {
        PyErr_SetString(PyExc_ValueError, "Failed to parse datapoint dimension.");
        goto l_cleanup;
    }

    /* Convert input to vector format. */
    clusters_vector = VECTOR_allocate(sizeof(struct cluster_item), num_centroids);
    if (NULL == clusters_vector) {
        PyErr_SetString(PyExc_EnvironmentError, "Failed to allocate memory.");
        goto l_cleanup;
    }  
    
    datapoints_vector = VECTOR_allocate(sizeof(struct datapoint_item), num_datapoints);
    if (NULL == datapoints_vector) {
        PyErr_SetString(PyExc_EnvironmentError, "Failed to allocate memory.");
        goto l_cleanup;
    }    

    /* First, load centroids into the datapoints vector. */
    status_code = add_datapoints_to_vector(datapoints_vector, centroids_coords);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Create clusters with initial centroids. */
    for (size_t index = 0; index < num_centroids; index++) {
        status_code = VECTOR_get_item(datapoints_vector, index, &temp_datapoint_item);
        if (STANDARD_SUCCESS_CODE != status_code) {
            PyErr_SetString(PyExc_ValueError, "Failed to get initial centroids.");
            goto l_cleanup;
        }

        status_code = KMEANS_create_cluster_item(clusters_vector, &temp_datapoint_item, data_dimension);
        if (STANDARD_SUCCESS_CODE != status_code) {
            PyErr_SetString(PyExc_ValueError, "Failed to create initial centroids.");
            goto l_cleanup;
        }
    }

    /* Load the rest of the datapoints. */
    status_code = add_datapoints_to_vector(datapoints_vector, datapoints_coords);
    if (STANDARD_SUCCESS_CODE != status_code) {
        goto l_cleanup;
    }

    /* Execute clustering algorithm. */
    status_code = KMEANS_perform_kmeans_clustering(
        clusters_vector,
        datapoints_vector, 
        num_centroids, 
        num_iterations, 
        data_dimension,
        convergence_threshold
    );
    if (STANDARD_SUCCESS_CODE != status_code) {
        PyErr_SetString(PyExc_ValueError, "Failed to execute clustering.");
        goto l_cleanup;
    }

    /* Convert result to Python structure. */
    output_centroids = build_output_centroids(clusters_vector, data_dimension);

l_cleanup:
    if (NULL != clusters_vector) {
        KMEANS_free_cluster_vector(clusters_vector);
    }
    if (NULL != datapoints_vector) {
        KMEANS_free_datapoint_vector(datapoints_vector);
    }

    return output_centroids;
}


/********************** Structures **********************/
static PyMethodDef kmeans_methods[] = {
    {"fit",                      /* the Python method name that will be used */
      (PyCFunction) kmeans_fit,  /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,              /* flags indicating parameters accepted for this function */
      PyDoc_STR("TODO")},        /*  The docstring for the function */
    {NULL, NULL, 0, NULL}        /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};

static struct PyModuleDef kmeans_module = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",                /* name of module */
    NULL,                        /* module documentation, may be NULL */
    -1,                          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    kmeans_methods               /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeans_module);
    if (!m) {
        return NULL;
    }
    return m;
}
