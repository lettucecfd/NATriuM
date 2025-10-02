import numpy as np
import vtk
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import itertools
import matplotlib.image as mpimg

rc('font', **{'size': 11, 'family': 'sans-serif', 'sans-serif': ['Myriad Pro', 'Arial', 'Tahoma']})
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

# hafezOld = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/reference_hafez_from_frapolli.txt',delimiter=';',skiprows=0)
hafez = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets.csv',delimiter=';',skiprows=1)
hafez = hafez[hafez[:,0].argsort()]
latt = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/reference.txt',delimiter=';',skiprows=0)
latt = latt[latt[:,0].argsort()]
frap = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/frapolli.txt',delimiter=';',skiprows=0)
frap = frap[frap[:,0].argsort()]
dorschner = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_dorschner.csv',delimiter=',',skiprows=2)
dorschner = dorschner[dorschner[:,0].argsort()]
saadat = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_saadat.csv',delimiter=',',skiprows=2)
saadat = saadat[saadat[:,0].argsort()]
tran = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_tran2.csv',delimiter=',',skiprows=2)
tran = tran[tran[:,0].argsort()]
hafez_from_dorschner = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_hafez_from_dorschner2.csv',delimiter=',',skiprows=2)
hafez_from_dorschner = hafez_from_dorschner[hafez_from_dorschner[:,0].argsort()]
thyagarajan = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_thyagarajan.csv',delimiter=',',skiprows=2)
thyagarajan = thyagarajan[thyagarajan[:,0].argsort()]
noh = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/wpd_datasets_noh.csv',delimiter=',',skiprows=2)
noh = noh[noh[:,0].argsort()]


for jobid, jobpath, vtupath, imgpath in zip(jobids, jobpaths, vtupaths, imgpaths):
    iTlist = [iT.removeprefix("t_0.").removesuffix(".pvtu") for iT in os.listdir(vtupath) if iT.endswith(".pvtu")]
    for iT in iTlist:#["100000","200000"]:
    # for iT in ["94000"]:
        ref_m = itertools.cycle(('o', 'v', '^', 's', 'p', 'h', 'D'))
        ref_c = itertools.cycle(('royalblue', 'green', 'cyan', 'magenta', 'grey', 'black'))
    
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

            ### X-AXIS
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

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], results[:,1])
        # fig.savefig(imgpath + "coords_" + iT + ".png")

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], results[:,2])
        # fig.savefig(imgpath + "rho_" + iT + ".png")

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], results[:,3])
        # fig.savefig(imgpath + "ux_" + iT + ".png")

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], results[:,4])
        # fig.savefig(imgpath + "uy_" + iT + ".png")

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], results[:,5])
        # fig.savefig(imgpath + "T_" + iT + ".png")

        # fig, ax = plt.subplots()
        # ax.plot(results[:,0], MaLocal)
        # fig.savefig(imgpath + "MaLocal_" + jobid + ".png")

        # img = mpimg.imread(imgpath + 'CpField_9836194_final_final_iT50000.png')
        fig, ax = plt.subplots(figsize=[7, 3.5])
        ax.plot(hafez[::2,0],hafez[::2,1],
                color='black', linestyle='-', linewidth=1.5,
                label=r'Hafez \& Wahba [1]')
        clbm = ('green', .7)
        mlbm = '.'#next(ref_m)
        mslbm = 2
        lw = 1.2
        ax.plot(latt[:,0],latt[:,1],
                # color=next(ref_c), marker=next(ref_m), linewidth=0, markersize=2,
                # color=clbm, marker=mlbm, linewidth=0, markersize=mslbm,
                color=clbm, linewidth=lw,
                label='LBM results [2-5,7,8,10,12]'
                )
        ax.plot(frap[:,0],frap[:,1],color=clbm, linewidth=lw)
        ax.plot(dorschner[:,0],dorschner[:,1],color=clbm,linewidth=lw)
        ax.plot(saadat[:,0],saadat[:,1],color=clbm,linewidth=lw)
        ax.plot(tran[:,0],tran[:,1],color=clbm,linewidth=lw)
        ax.plot(thyagarajan[:,0],thyagarajan[:,1],color=clbm,linewidth=lw)
        ax.plot(noh[:,0],noh[:,1],color=clbm,linewidth=lw)
        ax.plot(results[:,0], Cp[:],color=clbm, linewidth=lw)
        ax.plot(hafez_from_dorschner[:,0],hafez_from_dorschner[:,1],
                # color=next(ref_c), linestyle='-.',
                # color=next(ref_c), marker=next(ref_m), linewidth=0,
                color='red', linewidth=1.5,
                label='[1] as displayed in [2-5]')
        ax.set_ylim(-0.25,1.8)
        ax.set_xlim(-1,1.5)
        # ax.set_xticks([-.75, -0.5, -.25, 0, .25, 0.5, .75, 1])
        ax.set_xticks([-1, -0.5, 0, .5, 1, 1.5])
        ax.set_yticks([0, .5, 1, 1.5])
        ax.set_xlabel(r'$\mathbf{x/C}$', fontdict={'size': 14})
        ax.set_ylabel(r'$\mathbf{C_P}$', fontdict={'size': 14})
        ax.legend(frameon=False, loc='upper left')
        # ax.imshow(img, extent=(-1.7, 3.8, -1.5, 1.5))
        # plt.show()
        fig.savefig(imgpath + "Cp_" + jobid + "_iT" + iT + "_high.png", transparent=True, dpi=300)
        fig.savefig(imgpath + "Cp_" + jobid + "_iT" + iT + "_high.pdf", transparent=True)

        plt.close("all")
    print(f"Finished jobid {jobid}")

