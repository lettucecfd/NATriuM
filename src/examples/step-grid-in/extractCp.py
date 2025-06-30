import numpy as np
import vtk
import matplotlib.pyplot as plt
from scipy.spatial import KDTree


# Example Usage:
jobid = "9682930"
pvtu_file_path = jobid + "/output/t_0.200000.pvtu"  # Replace with the actual path to your .pvtu file

"""
Reads a .pvtu file and extracts values along the x-axis.

If no points lie exactly on the x-axis (y=0, z=0), it finds the point with the lowest x-coordinate
and uses that point's data for the x-axis value.

Args:
    pvtu_file (str): Path to the .pvtu file.

Returns:
    list: A list of values extracted along the x-axis.  Returns an empty list if the file
            cannot be read or if no data is found.
"""

reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName(pvtu_file_path)
reader.Update()

data = reader.GetOutput()

if data is None:
    print(f"Error: Could not read data from {pvtu_file_path}")
    exit()

points = data.GetPoints()
point_data = data.GetPointData().GetArray("rho")  # points.GetData()

xList = []
yList = []
dataList = []
xFoilList = []
yFoilList = []
dataFoilList = []

# If no points on the x-axis, find the point with the lowest x-coordinate
for i in range(points.GetNumberOfPoints()):
    xi, yi, _ = points.GetPoint(i)
    if (abs(yi) < 1e-6) and (-1 < xi < 1.5):# or 0 < x < 1:
        xList.append(xi)
        yList.append(yi)
        dataList.append(point_data.GetValue(i))
    if (0 < xi < 1) and (0 < yi < 0.1):
        xFoilList.append(xi)
        yFoilList.append(yi)
        dataFoilList.append(point_data.GetValue(i))

valuesAxis = np.array([xList,yList,dataList])
uniqueCoords, uniqueCounts = np.unique(np.array([xList,yList]), axis=1, return_counts=True)
dataAxis = []
for i in range(len(uniqueCoords[0,:])):
    x, y = uniqueCoords[:,i]
    valuesAtI = valuesAxis[:,np.isclose(x, valuesAxis[0,:]) * np.isclose(y, valuesAxis[1,:])]
    dataAxis.append([x,y,np.mean(valuesAtI[2,:])])
dataAxis = np.array(dataAxis)

# get lowest y-coordinate points
dataAroundFoil = np.unique(np.array([xFoilList, yFoilList, dataFoilList]), axis=1)
dataFoil = []
n = 1000
dx = 1/n
for xi in np.linspace(0,1,n):
    dataAtX = dataAroundFoil[:,np.isclose(xi, dataAroundFoil[0,:], atol=dx)]
    # dataAtX = dataAtX[:,dataAtX[1]>np.min(dataAtX[1])+1e-10]
    dataFoil.append(dataAtX[:,np.argmin(dataAtX[1])])
dataFoil = np.array(dataFoil)
    

results = np.concat((dataAxis, dataFoil))
x1 = results[:,0]
rhoLU = results[:,2]
Cp = (rhoLU-.9)/(0.5*1.4*1.5*1.5)*2
results = np.hstack((results,np.expand_dims(Cp,1)))
np.savetxt(jobid + "/Cp.txt", results)
plt.scatter(results[:,0], results[:,3])
plt.show()