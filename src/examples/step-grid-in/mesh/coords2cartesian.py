import numpy as np
import matplotlib.pyplot as plt
nx = 100

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
      x.append(coords[0].strip(' '))#float(coords[0]))
      y.append(coords[1].strip(' '))#float(coords[1]))

    # writing all entries except second and last (so that trailing edge is not quite so fine)
    #print("x={" + x[0] + ", " + ''.join([xi + ", " for xi in x[2:-2]]) + x[-2] + "};")
    #print("y={" + y[0] + ", " + ''.join([yi + ", " for yi in y[2:-2]]) + y[-2] + "};")

    x = np.array([float(xi) for xi in x])
    y = np.array([float(yi) for yi in y])
    xTop = np.flip(x[y>=0])
    yTop = np.flip(y[y>=0])
    xBot = x[y<0][1:]
    yBot = y[y<0][1:]
    xLinTop = np.round(np.linspace(start=1, stop=0, num=nx+1, endpoint=True), 4)
    yIntTop = np.interp(xLinTop, xTop, yTop, left=0, right=0)
    xLinBot = np.round(np.linspace(start=0, stop=1, num=nx+1, endpoint=True), 4)
    # print(xBot)
    # print(yBot)
    # print(xLinBot)
    yIntBot = np.flip(np.interp(xLinTop, xBot, yBot, left=0, right=0))
    # print(yIntBot)
    xStrings = []
    yStrings = []
    for i in range(len(xLinTop)-1):
      xStrings.append(str(xLinTop[i]))
      xStrings.append(str(xLinTop[i]))
      yStrings.append(str(yIntTop[i]))
      yStrings.append(str(yIntTop[i+1]))
    xStrings.append(str(xLinTop[-1]))
    yStrings.append(str(yIntTop[-1]))
    for i in range(len(xLinBot)-1):
      xStrings.append(str(xLinBot[i]))
      xStrings.append(str(xLinBot[i]))
      yStrings.append(str(yIntBot[i]))
      yStrings.append(str(yIntBot[i+1]))
    xStrings.append(str(xLinBot[-1]))
    yStrings.append(str(yIntBot[-1]))
    
    print("x={" + ''.join([xi + ", " for xi in xStrings[:-2]]) + xStrings[-1] + "};")
    print("y={" + ''.join([yi + ", " for yi in yStrings[:-2]]) + yStrings[-1] + "};")
    # print("y={" + str(yIntTop[0]) + ", " + ''.join([str(yi) + ", " + str(yi) + ", " for yi in yIntTop[2:]]) + ''.join([str(yi) + ", " for yi in yIntBot[2:-2]]) + str(yIntBot[-2]) + "};")
    # plt.plot(xLinTop, yIntTop, 'r-')
    # plt.plot(xLinBot, yIntBot, 'b-')
    # plt.show()