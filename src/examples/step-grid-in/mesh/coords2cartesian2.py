import numpy as np
import matplotlib.pyplot as plt
nx = 400

for filename in ["nonUni_closed/NACA0012_0deg.geo"]:
  with open(filename, 'r') as file:
    string = file.read()
    string = string.split("Point(1)")
    string = string[1].split("Line(1)")
    string = string[0].splitlines()
    x = []
    y = []
    for line in string:
      line = line.split("{")
      coords = line[1].split(",")
      x.append(coords[0].strip(' '))
      y.append(coords[1].strip(' '))

    x = np.array([float(xi) for xi in x])
    y = np.array([float(yi) for yi in y])
    xTop = np.flip(x[y>=0])
    yTop = np.flip(y[y>=0])
    xBot = x[y<0][1:]
    yBot = y[y<0][1:]
    xLinTop = np.round(np.linspace(start=1, stop=0, num=nx+1, endpoint=True), 4)
    yIntTop = np.interp(xLinTop, xTop, yTop, left=0, right=0)
    xLinBot = np.round(np.linspace(start=0, stop=1, num=nx+1, endpoint=True), 4)
    yIntBot = np.flip(np.interp(xLinTop, xBot, yBot, left=0, right=0))
    
    dy = 1/nx
    yLinTop = []
    for yInt in yIntTop:
      y = 0
      while y < yInt:
        y += dy
      yLinTop.append(np.round(y,5))
    yLinBot = []
    for yInt in yIntBot:
      y = 0
      while y > yInt:
        y -= dy
      yLinBot.append(np.round(y,5))
    yLinBot = np.array(yLinBot)
    # xLinBot = xLinBot[yLinBot != 0]
    # yLinBot = yLinBot[yLinBot != 0]

    # yLinTop = np.round(yIntTop*nx,0)/nx
    # yLinBot = (np.round(yIntBot*nx,0)-1)/nx

    xStrings = []
    yStrings = []
    for i in range(len(xLinTop)-1):
      xStrings.append(xLinTop[i])
      yStrings.append(yLinTop[i])
      if yLinTop[i] != yLinTop[i+1]:
        xStrings.append(xLinTop[i])
        yStrings.append(yLinTop[i+1])
    # xStrings.append(xLinTop[-1])
    # yStrings.append(yLinTop[-1])
    for i in range(len(xLinBot)-2):
      xStrings.append(xLinBot[i])
      yStrings.append(yLinBot[i])
      if yLinBot[i] != yLinBot[i+1]:
        xStrings.append(xLinBot[i])
        yStrings.append(yLinBot[i+1])
    xStrings.append(xLinBot[-2])
    yStrings.append(yLinBot[-3])
    xStrings.append(xLinBot[-2])
    yStrings.append(yLinBot[-2])
    xStrings.append(xLinBot[-2])
    yStrings.append(yLinBot[-1])
    xStrings = np.array(xStrings)
    yStrings = np.array(yStrings)

    pidTop = np.argwhere(yStrings == max(yStrings)).max()
    pidFront = np.argwhere(xStrings == min(xStrings)).min()
    pidBot = np.argwhere(yStrings == min(yStrings)).min()
    
    print("x={" + ''.join([str(xi) + ", " for xi in xStrings[:-2]]) + str(xStrings[-1]) + "};")
    print("y={" + ''.join([str(yi) + ", " for yi in yStrings[:-2]]) + str(yStrings[-1]) + "};")
    print("point_id_top = " + str(pidTop+1)  + ";\npoint_id_front = " + str(pidFront+1) + ";\npoint_id_bot = " + str(pidBot+1) + ";")
    # print("y={" + str(yIntTop[0]) + ", " + ''.join([str(yi) + ", " + str(yi) + ", " for yi in yIntTop[2:]]) + ''.join([str(yi) + ", " for yi in yIntBot[2:-2]]) + str(yIntBot[-2]) + "};")
    plt.plot(xLinTop, yLinTop, 'r-')
    plt.plot(xLinBot, yLinBot, 'b-')
    # print(xStrings.shape, yStrings.shape)
    plt.scatter(xStrings, yStrings)
    plt.show()