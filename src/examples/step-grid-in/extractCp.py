import numpy as np
import vtk
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import itertools

rc('font', **{'size': 9, 'family': 'serif', 'serif': ['Charis SIL']})
plt.rcParams['text.usetex'] = True
plt.rcParams["savefig.dpi"] = 600
plt.rcParams['markers.fillstyle'] = 'none'
plt.rcParams['figure.constrained_layout.use'] = True

sllbm_marker = 'x'
sllbm_color = 'red'
sllbm_size = 1.5*plt.rcParams['lines.markersize']

Tref = 1  # TODO: may differ
rho0LU = 1
RLU = 1

basepath = "/mnt/c/Users/phili/Desktop/"

# jobids = ["dw872713"]
# jobpaths = [basepath+"NACA/dw872713-natrium/step-grid-in-old/"]
# vtupaths = [basepath+"NACA/dw872713-natrium/step-grid-in-old/Re10000-Ma1.5-reflevel0-time1695800330/"]

# jobids = ["9694158", "9694159", "9694162"]
# jobids = ["9700609", "9700727", "9693848", "9684155"]
# coordsfile = "/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/varyRefinement/naca0012_res100.txt"
    # jobids = ["9687726"]
    # coordsfile = "/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/varyRefinement/naca0012_res60.txt"
jobids = ["9836194_final_final"]#["9796896"]
jobpaths = [basepath + jobid + "/" for jobid in jobids]
vtupaths = [jobpath + "output/" for jobpath in jobpaths]
imgpaths = [jobpath + "images/" for jobpath in jobpaths]
coordsfile = "/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/varyRefinement/naca0012_res100_combined.txt"

for imgpath in imgpaths:
    if not os.path.exists(imgpath):
        os.mkdir(imgpath)
xFoil = []
yFoil = []
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

xListAxis = [xi for xi in np.linspace(-1,0,200)]
[xListAxis.append(xi) for xi in np.linspace(1,1.5,100)]
yListAxis = [0 for _ in xListAxis]
axisCoords = np.array([xListAxis,yListAxis])

ref = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/reference.txt',delimiter=';',skiprows=0)
frap = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/frapolli.txt',delimiter=';',skiprows=0)
# hafezOld = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/reference_hafez_from_frapolli.txt',delimiter=';',skiprows=0)
hafez = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets.csv',delimiter=';',skiprows=1)
hafez = hafez[hafez[:,0].argsort()]
# hafez = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets(1).csv',delimiter=',',skiprows=2)

for jobid, jobpath, vtupath, imgpath in zip(jobids, jobpaths, vtupaths, imgpaths):
    iTlist = [iT.removeprefix("t_0.").removesuffix(".pvtu") for iT in os.listdir(vtupath) if iT.endswith(".pvtu")]
    for iT in iTlist:#["100000","200000"]:
    # for iT in ["94000"]:
        ref_m = itertools.cycle(('o', 'v', '^', 's', 'p', 'h', 'D'))
        ref_c = itertools.cycle(('royalblue', 'green', 'grey', 'black', 'cyan', 'magenta'))
    
        if jobid != "dw872713":
            file = open(jobpath + "/slurm_natrium_naca.out", "r")
        else:
            file = open(jobpath + "/natrium.log", "r")
        lines = file.read().splitlines()
        cs = float([line for line in lines if "::::Sound speed:" in line][0].removeprefix("::::Sound speed:              "))
        gamma = float([line for line in lines if "::::Heat capacity ratio:" in line][0].removeprefix("::::Heat capacity ratio:      "))
        Ma = float([line for line in lines if "::::Mach number:" in line][0].removeprefix("::::Mach number:              "))/np.sqrt(gamma)
        file.close()

        if not os.path.exists(imgpath + "Cp_" + jobid + "_iT" + iT + ".txt"):
            pvtu_file_path = vtupath + "t_0." + iT + ".pvtu"

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
            x1 = np.max(foilCoords[0])
            y1 = np.min(foilCoords[1])
            x2 = 6
            y2 = 0

            locator = vtk.vtkPointLocator()
            locator.SetDataSet(data)
            locator.BuildLocator()

            ### X-AXIS SOLUTION 1
            ### get lowest y-coordinate points
            # dataAroundFoil = np.unique(np.array([xListFoil, yListFoil, rhoListFoil]), axis=1)
            # rhoFoil = []
            # n = 1000
            # dx = 1/n
            # for xi in np.linspace(0,1,n):
            #     dataAtX = dataAroundFoil[:,np.isclose(xi, dataAroundFoil[0,:], atol=dx)]
            #     # dataAtX = dataAtX[:,dataAtX[1]>np.min(dataAtX[1])+1e-10]
            #     rhoFoil.append(dataAtX[:,np.argmin(dataAtX[1])])
            # rhoFoil = np.array(rhoFoil)

            ### X-AXIS SOLUTION 2
            # getting points on x-axis
            # # If no points on the x-axis, find the point with the lowest x-coordinate
            # for i in range(points.GetNumberOfPoints()):
            #     xi, yi, _ = points.GetPoint(i)
            #     if jobid == "dw872713":
            #         xi += 0.5
            #     outletCentreI = (y2-y1)/(x2-x1)*(xi-x1)+y1
            #     if ((-1 <= xi < 0) and (abs(yi) < 1e-10)) or ((1 < xi <= 1.5) and (abs(yi-outletCentreI) < 1e-5)):# or 0 < x < 1:
            #         xListAxis.append(xi)
            #         yListAxis.append(yi)
            #         rhoListAxis.append(point_rho.GetValue(i))
            #         uxListAxis.append(point_ux.GetValue(i))
            #         uyListAxis.append(point_uy.GetValue(i))
            #         TListAxis.append(point_T.GetValue(i))

            # valuesAxis = np.array([xListAxis,yListAxis,rhoListAxis,uxListAxis,uyListAxis,TListAxis])
            # uniqueCoords, uniqueCounts = np.unique(np.array([xListAxis,yListAxis]), axis=1, return_counts=True)
            # dataAxis = []
            # for i in range(len(uniqueCoords[0,:])):
            #     x, y = uniqueCoords[:,i]
            #     valuesAtI = valuesAxis[:,np.isclose(x, valuesAxis[0,:]) * np.isclose(y, valuesAxis[1,:])]
            #     dataAxis.append([x,y,np.mean(valuesAtI[2,:]),np.mean(valuesAtI[3,:]),np.mean(valuesAtI[4,:]),np.mean(valuesAtI[5,:])])
            # dataAxis = np.array(dataAxis)

            ### X-AXIS SOLUTION 3
            # getting points closest to samples on x-axis
            coordsList = axisCoords.T
            for coord in coordsList:
                xi, yi = coord
                if jobid == "dw872713":
                    xi -= 0.5
                zi = 0
                point_id_list = vtk.vtkIdList()
                point_id = locator.FindClosestPoint(xi, yi, zi)
                rhoListAxis.append(point_rho.GetValue(point_id))
                uxListAxis.append(point_ux.GetValue(point_id))
                uyListAxis.append(point_uy.GetValue(point_id))
                TListAxis.append(point_T.GetValue(point_id))
            dataAxis = np.hstack((coordsList,
                                np.expand_dims(np.array(rhoListAxis),1),
                                np.expand_dims(np.array(uxListAxis),1),
                                np.expand_dims(np.array(uyListAxis),1),
                                np.expand_dims(np.array(TListAxis),1)))

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

            results = np.vstack((dataAxis, dataFoil))
            results = results[np.argsort(results[:,0])]
            
            rhoLU = results[:,2]
            uxLU = results[:,3]
            uyLU = results[:,4]
            TLU = results[:,5]
            
            UmagLU = np.sqrt(uxLU*uxLU+uyLU*uyLU)
            MaLocal = UmagLU*1.5*1.5*1.5*gamma*Tref/3

            pLU = rhoLU*RLU*TLU
            p0LU = rho0LU*RLU*Tref
            Cp = (pLU-p0LU)/(0.5*gamma*p0LU*Ma*Ma)

            results = np.hstack((results,np.expand_dims(Cp,1)))
            np.savetxt(imgpath + "Cp_" + jobid + "_iT" + iT + ".txt", results)

        else:
            results = np.loadtxt(imgpath + "Cp_" + jobid + "_iT" + iT + ".txt")
        
        uxLU = results[:,3]
        uyLU = results[:,4]
        Cp = results[:,6]
        UmagLU = np.sqrt(uxLU*uxLU+uyLU*uyLU)
        MaLocal = UmagLU*1.5*1.5*1.5*gamma*Tref/3

        fig, ax = plt.subplots()
        ax.plot(results[:,0], results[:,1])
        # ax.set_ylim((-.5,.5))
        # ax.set_xlim((-1, 1.5))
        fig.savefig(imgpath + "coords_" + iT + ".png")

        fig, ax = plt.subplots()
        ax.plot(results[:,0], results[:,2])
        fig.savefig(imgpath + "rho_" + iT + ".png")

        fig, ax = plt.subplots()
        ax.plot(results[:,0], results[:,3])
        fig.savefig(imgpath + "ux_" + iT + ".png")

        fig, ax = plt.subplots()
        ax.plot(results[:,0], results[:,4])
        fig.savefig(imgpath + "uy_" + iT + ".png")

        fig, ax = plt.subplots()
        ax.plot(results[:,0], results[:,5])
        fig.savefig(imgpath + "T_" + iT + ".png")

        fig, ax = plt.subplots()
        ax.plot(results[:,0], MaLocal)
        fig.savefig(imgpath + "MaLocal_" + jobid + ".png")

        fig, ax = plt.subplots(figsize=[5.8, 2.3])
        ax.plot(ref[:,0],ref[:,1],
                linewidth=0, color=next(ref_c), marker=next(ref_m),
                label='Latt et al.')
        ax.plot(frap[:,0],frap[:,1],
                linewidth=0, color=next(ref_c), marker=next(ref_m),
                label='Frapolli et al.')
        # ax.plot(hafez[:,0],hafez[:,1],label=r'Hafez & Wahba',color='black')
        ax.set_ylim(1.8,-0.3)
        ax.set_xlim(-1,1.5)
        ax.set_xlabel(r'$x/C$')
        ax.set_ylabel(r'$C_P$')
        ax.plot(results[:,0], Cp,    color='red',    label='SLLBM')
        ax.legend(frameon=False, loc='lower right')
        fig.savefig(imgpath + "Cp_" + jobid + "_iT" + iT + ".png")
        fig.savefig(imgpath + "Cp_" + jobid + "_iT" + iT + ".pdf")

        plt.close("all")
    print(f"Finished jobid {jobid}")

