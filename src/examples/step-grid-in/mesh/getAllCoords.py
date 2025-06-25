for filename in ["uni_closed/NACA0012_0deg.geo"]:
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

    print("x={" + ''.join([xi + ", " for xi in x[:-1]]) + x[-1] + "};")
    print("y={" + ''.join([yi + ", " for yi in y[:-1]]) + y[-1] + "};")