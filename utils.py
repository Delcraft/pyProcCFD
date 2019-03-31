#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 19:23:06 2019

@author: delcraft
"""
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read_vtu(filename):
    """
    Reads .vtu file and returns vtkUnstructuredGrid object
    Parameters:
        filename: str
            path to a .vtu file
    Returns:
        out: vtkCommonDataModelPython.vtkUnstructuredGrid
            result of GetOutput() attrubute of vtkXMLUnstructuredGridReader()
    """
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(filename)
    r.Update()
    return r.GetOutput()

def get_node_field_names(grid):
    """
    Obtain the list of node-stored field names
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        out: list of str
            list of node-stored field names
    """
    pdata = grid.GetPointData()
    num_fields = pdata.GetNumberOfArrays()
    field_names = []
    for i in range(num_fields):
        field_names.append(pdata.GetArrayName(i))
    return field_names

def get_node_coords(grid):
    """
    Obtain the coordinate of the nodes in mesh
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        (x,y,z): tuple, np.float32
            tuple of the node coordinates
    """
    pts = grid.GetPoints().GetData()
    pts_array = vtk_to_numpy(pts)
    x,y,z = (pts_array[:,0],pts_array[:,1],pts_array[:,2])
    return x,y,z

def get_node_field(grid,fieldname):
    """
    Obtain the node-stored field by its name
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
        fieldname: str
            name of the field to retrieve
    Returns:
        out: numpy.ndarray
            numpy array with the node values of the field
    """
    pdata = grid.GetPointData()
    try:
        fpos = get_node_field_names(grid).index(fieldname)
    except:
        raise Exception('The field {0:s} is not stored in the nodes of the mesh.'.format(fieldname))
    field = vtk_to_numpy(pdata.GetArray(fpos))
    return field