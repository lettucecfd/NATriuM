import numpy as np
import vtk
import matplotlib.pyplot as plt
import os


basepath = "/mnt/c/Users/phili/Desktop/"
jobids = ["9682930", "dw872713"]
vtupaths = ["9682930/output", "NACA/dw872713-natrium/step-grid-in-old/Re10000-Ma1.5-reflevel0-time1695800330"]
imgpath = basepath + "images/"
if not os.path.exists(imgpath):
    os.mkdir(imgpath)
xFoil = []
yFoil = []
coordsfile = "/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/varyRefinement/naca0012_res100.txt"
with open(coordsfile, "r") as f:
    string = f.read().splitlines()[1:]
    for line in string:
      if line not in ['', ' ']:
        line = line.split(' ')
        line = [l for l in line if l != '']
        xFoil.append(float(line[0]))
        yFoil.append(float(line[1]))
foilCoords = np.array([xFoil,yFoil])
foilCoords = foilCoords[:,foilCoords[1,:] >= 0]
foilCoords = foilCoords[:,:-2]

for jobid, vtupath in zip(jobids, vtupaths):
    pvtu_file_path = basepath + vtupath + "/t_0.200000.pvtu"

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
        if jobid == "dw872713":
            xi += 0.5
        if ((-1 < xi < 0) or (1 < xi < 1.5)) and (abs(yi) < 1e-6):# or 0 < x < 1:
            xList.append(xi)
            yList.append(yi)
            dataList.append(point_data.GetValue(i))
        # if (0 < xi < 1) and (0 < yi < 0.1):
        #     xFoilList.append(xi)
        #     yFoilList.append(yi)
        #     dataFoilList.append(point_data.GetValue(i))

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(data)
    locator.BuildLocator()

    dataFoilList = []
    for coord in foilCoords.T:
        xi, yi = coord
        zi = 0
        point_id_list = vtk.vtkIdList()
        point_id = locator.FindClosestPoint(xi, yi, zi)
        dataFoilList.append(point_data.GetValue(point_id))
    dataFoil = np.hstack((foilCoords.T,np.expand_dims(np.array(dataFoilList),1)))

    valuesAxis = np.array([xList,yList,dataList])
    uniqueCoords, uniqueCounts = np.unique(np.array([xList,yList]), axis=1, return_counts=True)
    dataAxis = []
    for i in range(len(uniqueCoords[0,:])):
        x, y = uniqueCoords[:,i]
        valuesAtI = valuesAxis[:,np.isclose(x, valuesAxis[0,:]) * np.isclose(y, valuesAxis[1,:])]
        dataAxis.append([x,y,np.mean(valuesAtI[2,:])])
    dataAxis = np.array(dataAxis)

    # # get lowest y-coordinate points
    # dataAroundFoil = np.unique(np.array([xFoilList, yFoilList, dataFoilList]), axis=1)
    # dataFoil = []
    # n = 1000
    # dx = 1/n
    # for xi in np.linspace(0,1,n):
    #     dataAtX = dataAroundFoil[:,np.isclose(xi, dataAroundFoil[0,:], atol=dx)]
    #     # dataAtX = dataAtX[:,dataAtX[1]>np.min(dataAtX[1])+1e-10]
    #     dataFoil.append(dataAtX[:,np.argmin(dataAtX[1])])
    # dataFoil = np.array(dataFoil)        

    results = np.vstack((dataAxis, dataFoil))
    x1 = results[:,0]
    rhoLU = results[:,2]
    rho0 = 1
    Cp = (rhoLU-rho0)/(0.5*1.4*1.5*1.5)
    results = np.hstack((results,np.expand_dims(Cp,1)))
    results = results[np.argsort(results[:,0])]
    np.savetxt(imgpath + "Cp_" + jobid + ".txt", results)

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,1])
    # ax.set_ylim((-.5,.5))
    # ax.set_xlim((-1, 1.5))
    fig.savefig(imgpath + "coords_" + jobid + ".png")

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,3])
    # ax.set_ylim((1.5,-0.1))
    # ax.set_xlim((-1, 1.5))
    fig.savefig(imgpath + "Cp_" + jobid + ".png")

    print(f"Finished jobid {jobid}")

    