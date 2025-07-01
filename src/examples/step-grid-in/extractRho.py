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
    point_rho = data.GetPointData().GetArray("rho")
    point_ux = data.GetPointData().GetArray("ux")
    point_uy = data.GetPointData().GetArray("uy")
    point_T = data.GetPointData().GetArray("T")

    xListAxis = []
    yListAxis = []
    rhoListAxis = []
    uxListAxis = []
    uyListAxis = []
    TListAxis = []
    xListFoil = []
    yListFoil = []
    rhoListFoil = []
    uxListFoil = []
    uyListFoil = []
    TListFoil = []

    # If no points on the x-axis, find the point with the lowest x-coordinate
    for i in range(points.GetNumberOfPoints()):
        xi, yi, _ = points.GetPoint(i)
        if jobid == "dw872713":
            xi += 0.5
        if ((-1 <= xi < 0) or (1 < xi <= 1.5)) and (abs(yi) < 1e-8):# or 0 < x < 1:
            xListAxis.append(xi)
            yListAxis.append(yi)
            rhoListAxis.append(point_rho.GetValue(i))
            uxListAxis.append(point_ux.GetValue(i))
            uyListAxis.append(point_uy.GetValue(i))
            TListAxis.append(point_T.GetValue(i))
        # if (0 < xi < 1) and (0 < yi < 0.1):
        #     xListFoil.append(xi)
        #     yListFoil.append(yi)
        #     rhoListFoil.append(point_rho.GetValue(i))

    locator = vtk.vtkPointLocator()
    locator.SetDataSet(data)
    locator.BuildLocator()

    coordsList = foilCoords.T
    for coord in coordsList:
        xi, yi = coord
        if jobid == "dw872713":
            xi -= 0.5
        zi = 0
        point_id_list = vtk.vtkIdList()
        point_id = locator.FindClosestPoint(xi, yi, zi)
        rhoListFoil.append(point_rho.GetValue(point_id))
        uxListFoil.append(point_ux.GetValue(point_id))
        uyListFoil.append(point_uy.GetValue(point_id))
        TListFoil.append(point_T.GetValue(point_id))
    dataFoil = np.hstack((coordsList,
                          np.expand_dims(np.array(rhoListFoil),1),
                          np.expand_dims(np.array(uxListFoil),1),
                          np.expand_dims(np.array(uyListFoil),1),
                          np.expand_dims(np.array(TListFoil),1)))

    valuesAxis = np.array([xListAxis,yListAxis,rhoListAxis,uxListAxis,uyListAxis,TListAxis])
    uniqueCoords, uniqueCounts = np.unique(np.array([xListAxis,yListAxis]), axis=1, return_counts=True)
    dataAxis = []
    for i in range(len(uniqueCoords[0,:])):
        x, y = uniqueCoords[:,i]
        valuesAtI = valuesAxis[:,np.isclose(x, valuesAxis[0,:]) * np.isclose(y, valuesAxis[1,:])]
        dataAxis.append([x,y,np.mean(valuesAtI[2,:]),np.mean(valuesAtI[3,:]),np.mean(valuesAtI[4,:]),np.mean(valuesAtI[5,:])])
    dataAxis = np.array(dataAxis)

    # # get lowest y-coordinate points
    # dataAroundFoil = np.unique(np.array([xListFoil, yListFoil, rhoListFoil]), axis=1)
    # rhoFoil = []
    # n = 1000
    # dx = 1/n
    # for xi in np.linspace(0,1,n):
    #     dataAtX = dataAroundFoil[:,np.isclose(xi, dataAroundFoil[0,:], atol=dx)]
    #     # dataAtX = dataAtX[:,dataAtX[1]>np.min(dataAtX[1])+1e-10]
    #     rhoFoil.append(dataAtX[:,np.argmin(dataAtX[1])])
    # rhoFoil = np.array(rhoFoil)        

    results = np.vstack((dataAxis, dataFoil))
    results = results[np.argsort(results[:,0])]
    np.savetxt(imgpath + "rhoUxUy_" + jobid + ".txt", results)

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,1])
    # ax.set_ylim((-.5,.5))
    # ax.set_xlim((-1, 1.5))
    fig.savefig(imgpath + "coords_" + jobid + ".png")

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,2])
    fig.savefig(imgpath + "rho_" + jobid + ".png")

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,3])
    fig.savefig(imgpath + "ux_" + jobid + ".png")

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,4])
    fig.savefig(imgpath + "uy_" + jobid + ".png")

    fig, ax = plt.subplots()
    ax.plot(results[:,0], results[:,5])
    fig.savefig(imgpath + "T_" + jobid + ".png")

    print(f"Finished jobid {jobid}")

