#include <Python.h>
#include <numpy/arrayobject.h>
#include "astar.h"
#include <stdio.h>

static char module_docstring[] =
    "This module provides an interface for calculating chi-squared using C.";
static char astar_docstring[] =
    "Calculate the chi-squared of some data given a model.";

static PyObject *astar_astar(PyObject *self, PyObject *args);
static PyObject *astar_astar_bw(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"astar", astar_astar, METH_VARARGS, astar_docstring},
    {"astar_bw", astar_astar_bw, METH_VARARGS, astar_docstring},
    {NULL, NULL, 0, NULL}
};



PyMODINIT_FUNC init_astar(void)
{
    PyObject *m = Py_InitModule3("_astar", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
}

//static PyObject *astar_astar(PyObject *self, PyObject *args)
//{
//    long start_i, start_j, end_i, end_j, size_i, size_j;
//    PyObject *grid_obj, *out_obj;
//
//    /* Parse the input tuple */
//    if (!PyArg_ParseTuple(args, "llllOO", &start_i, &start_j, &end_i, &end_j,
//                             &grid_obj, &out_obj))
//        return NULL;
//
//    /* Interpret the input objects as numpy arrays. */
//    PyObject *grid_array = PyArray_FROM_OTF(grid_obj, NPY_DOUBLE, NPY_IN_ARRAY);
//
//    PyArrayObject *out_array=NULL;
//    out_array = (PyArrayObject*)PyArray_FROM_OTF(out_obj, NPY_DOUBLE, NPY_INOUT_ARRAY);
//
//    /* If that didn't work, throw an exception. */
//    if (grid_array == NULL)  {
//        Py_XDECREF(grid_array);
//        Py_XDECREF(out_array);
//        return NULL;
//    }
//
//    /* How many data points are there? */
//            size_i = PyArray_DIM(grid_array, 0);
//            size_j = PyArray_DIM(grid_array, 1);
//            //size = PyArray_SIZE(grid_array);
//            //type = PyArray_TYPE(grid_array);
//
//    /* Get pointers to the data as C-types. */
//    double *grid    = (double*)PyArray_DATA(grid_array);
//    double *out     = (double*)PyArray_DATA(out_array);
//    
//    //PyObject *out_array_obj=NULL;
//    //PyArrayObject *out_array=NULL;
//    //out_array = (PyArrayObject*)PyArray_FROM_OTF(out_array_obj, NPY_DOUBLE, NPY_INOUT_ARRAY);
//    //double *outgrid    = (double*)PyArray_DATA(out_array);
//
//    /* Call the external C function to compute the chi-squared. */
//    double value = astar(start_i, start_j, end_i, end_j, grid, size_i, size_j, out);
//
//    /* Clean up. */
//    Py_DECREF(grid_array);
//    //Py_DECREF(out_array);
//
//    if (value > 0) {
//        PyErr_SetString(PyExc_RuntimeError,
//                    "Can't find a path.");
//        return NULL;
//    }
//
//    //return out_array_obj;
//    return PyArray_Return(out_array);
//}


static PyObject *astar_astar(PyObject *self, PyObject *args)
{
    long start_i, start_j, end_i, end_j, i_range, j_range;
    PyObject *in_obj, *out_obj; //, *yerr_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "llllOOll", &start_i, &start_j, &end_i, &end_j, 
                                          &in_obj, &out_obj, &i_range, &j_range))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    //PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    PyArrayObject *out_array=NULL;
    out_array = (PyArrayObject*)PyArray_FROM_OTF(out_obj, NPY_DOUBLE, NPY_INOUT_ARRAY);

    /* If that didn't work, throw an exception. */
    //if (x_array == NULL || y_array == NULL || yerr_array == NULL) {
    if (in_array == NULL || out_array == NULL ) {
        Py_XDECREF(in_array);
        Py_XDECREF(out_array);
        //Py_XDECREF(yerr_array);
        return NULL;
    }


    /* How many data points are there? */
    //int N = (int)PyArray_DIM(in_array, 0);

    /* Get pointers to the data as C-types. */
    double *in    = (double*)PyArray_DATA(in_array);
    double *out    = (double*)PyArray_DATA(out_array);
    //double *yerr = (double*)PyArray_DATA(yerr_array);

    /* Call the external C function to compute the chi-squared. */
    //double value = chi2(start_i, start_j, in, out, N);
    double value = astar(start_i, start_j, end_i, end_j, in, out, i_range, j_range);

    /* Clean up. */
    Py_DECREF(in_array);
    //Py_DECREF(y_array);
    //Py_DECREF(yerr_array);

    //if (value < 0.0) {
    //    PyErr_SetString(PyExc_RuntimeError,
    //                "Chi-squared returned an impossible value.");
    //    return NULL;
    //}

    if (value > 0) {
        PyErr_SetString(PyExc_RuntimeError,
                    "Can't find a path.");
        return NULL;
    }


    return PyArray_Return(out_array);
}


static PyObject *astar_astar_bw(PyObject *self, PyObject *args)
{
    long start_i, start_j, end_i, end_j, i_range, j_range;
    PyObject *in_obj, *out_obj; //, *yerr_obj;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "llllOOll", &start_i, &start_j, &end_i, &end_j, 
                                          &in_obj, &out_obj, &i_range, &j_range))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *in_array = PyArray_FROM_OTF(in_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    //PyObject *y_array = PyArray_FROM_OTF(y_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    PyArrayObject *out_array=NULL;
    out_array = (PyArrayObject*)PyArray_FROM_OTF(out_obj, NPY_DOUBLE, NPY_INOUT_ARRAY);

    /* If that didn't work, throw an exception. */
    //if (x_array == NULL || y_array == NULL || yerr_array == NULL) {
    if (in_array == NULL || out_array == NULL ) {
        Py_XDECREF(in_array);
        Py_XDECREF(out_array);
        //Py_XDECREF(yerr_array);
        return NULL;
    }


    /* How many data points are there? */
    //int N = (int)PyArray_DIM(in_array, 0);

    /* Get pointers to the data as C-types. */
    double *in    = (double*)PyArray_DATA(in_array);
    double *out    = (double*)PyArray_DATA(out_array);
    //double *yerr = (double*)PyArray_DATA(yerr_array);

    /* Call the external C function to compute the chi-squared. */
    //double value = chi2(start_i, start_j, in, out, N);
    double value = astar_bw(start_i, start_j, end_i, end_j, in, out, i_range, j_range);

    /* Clean up. */
    Py_DECREF(in_array);
    //Py_DECREF(y_array);
    //Py_DECREF(yerr_array);

    //if (value < 0.0) {
    //    PyErr_SetString(PyExc_RuntimeError,
    //                "Chi-squared returned an impossible value.");
    //    return NULL;
    //}

    if (value > 0) {
        PyErr_SetString(PyExc_RuntimeError,
                    "Can't find a path.");
        return NULL;
    }


    return PyArray_Return(out_array);
}


