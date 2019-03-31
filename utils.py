#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Generic functions for data manipulation on vtkUnstructuredGrid object
@author: Delcraft
"""
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk

def readVtu(filename):
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

def getNodeFieldNames(grid):
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

def getCellFieldNames(grid):
    """
    Obtain the list of node-stored field names
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        out: list of str
            list of node-stored field names
    """
    cdata = grid.GetCellData()
    num_fields = cdata.GetNumberOfArrays()
    field_names = []
    for i in range(num_fields):
        field_names.append(cdata.GetArrayName(i))
    return field_names

def getNodeCoords(grid):
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

def getNodeField(grid,fieldname):
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
        fpos = getNodeFieldNames(grid).index(fieldname)
    except:
        raise Exception('The field {0:s} is not stored in the nodes of the mesh.'.format(fieldname))
    field = vtk_to_numpy(pdata.GetArray(fpos))
    return field

def getCellField(grid,fieldname):
    """
    Obtain the cell-stored field by its name
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
        fieldname: str
            name of the field to retrieve
    Returns:
        out: numpy.ndarray
            numpy array with the node values of the field
    """
    pdata = grid.GetCellData()
    try:
        fpos = getNodeFieldNames(grid).index(fieldname)
    except:
        raise Exception('The field {0:s} is not stored in the cells of the mesh.'.format(fieldname))
    field = vtk_to_numpy(pdata.GetArray(fpos))
    return field

def removeNodeField(grid,fieldNames):
    """
    Remove node field from the grid
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
        fieldNames: list of str
            node field names to remove
    Returns:
        out: None
    """
    pdata = grid.GetPointData()
    for fieldname in fieldNames:
        pdata.RemoveArray(fieldname)
    
def removeAllNodeFields(grid):
    """
    Remove all node fields from the grid
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        out: None
    """
    pdata = grid.GetPointData()
    fieldnames = getNodeFieldNames(grid)
    for field in fieldnames:
        pdata.RemoveArray(field)
        
def getGridSize(grid):
    """
    Get number of points and cells in the given grid
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        (cells,points) : tuple of int
            Tuple with the number of cells and points
    """
    return int(grid.GetNumberOfCells()),int(grid.GetNumberOfPoints())

def nodeFieldsToCells(grid):
    """
    Project node-stored fields on cells of the given grid
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        out: None
    """
    filt = vtk.vtkPointDataToCellData()
    filt.SetInputData(grid)
    filt.PassPointDataOn()
    filt.Update()
    data = filt.GetOutput()
    cdata = data.GetCellData()
    gridcdata = grid.GetCellData()
    for cellArrayIdx in range(cdata.GetNumberOfArrays()):
        gridcdata.AddArray(cdata.GetArray(cellArrayIdx))

def cellFieldsToNodes(grid):
    """
    Project node-stored fields on cells of the given grid
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
    Returns:
        out: None
    """
    filt = vtk.vtkCellDataToPointData()
    filt.SetInputData(grid)
    filt.PassCellDataOn()
    filt.Update()
    data = filt.GetOutput()
    pdata = data.GetPointData()
    gridpdata = grid.GetPointData()
    for nodeArrayIdx in range(pdata.GetNumberOfArrays()):
        gridpdata.AddArray(pdata.GetArray(nodeArrayIdx))
        
def computeGradient(grid,fieldNames,datatype='cell'):
    """
    Compute gradients of the cell- or node-based fields and add them
    as 1D arrays to the current mesh. The name of the new arrays will be
    'grad_[var_name]_[x/y/z]'. For each field, three new arays will be
    provided: gradients in X,Y and Z directions.
    Parameters:
        grid: vtkUnstructuredGrid
            grid vtk object
        fieldNames: list of str
            list of field names for which to compute the gradients
        datatype: str, default: 'cell'
            basis on which to compute the gradients
    Returns:
        out: None
    """
    grad = vtk.vtkGradientFilter()
    grad.SetInputData(grid)
    grad.SetComputeGradient(1)
    if datatype == 'cell':
        field_ass = vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS
    elif datatype == 'node':
        field_ass = vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS
    for fieldname in fieldNames:
        gradname = 'grad_{0:s}'.format(fieldname)
        print('---> Compute gradient for {0:s} at {1:s}s'.format(fieldname,datatype))
        grad.SetResultArrayName(gradname)
        grad.SetInputArrayToProcess(0,0,0,field_ass,fieldname)
        grad.Update()
        data = grad.GetOutput()
        if datatype == 'cell':
            cellfields = getCellFieldNames(data)
            idx = cellfields.index(gradname)
            arr = data.GetCellData().GetArray(idx)
            arr_np = vtk_to_numpy(arr)
            del arr
            for jdx,jname in enumerate(['_x','_y','_z']):
                print('Assign ' + gradname+jname + ' to mesh.')
                arr = numpy_to_vtk(arr_np[:,jdx],array_type=vtk.VTK_FLOAT)
                arr.SetName(gradname+jname)
                grid.GetCellData().AddArray(arr)
        elif datatype == 'node':
            nodefields = getNodeFieldNames(data)
            idx = nodefields.index(gradname)
            arr = data.GetPointData().GetArray(idx)
            arr_np = vtk_to_numpy(arr)
            del arr
            for jdx,jname in enumerate(['_x','_y','_z']):
                print('Assign ' + gradname+jname + ' to mesh.')
                arr = numpy_to_vtk(arr_np[:,jdx],array_type=vtk.VTK_FLOAT)
                arr.SetName(gradname+jname)
                grid.GetPointData().AddArray(arr)
            
        
        
        