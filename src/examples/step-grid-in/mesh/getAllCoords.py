for filename in ["nonUni_closed/NACA0012_4deg.geo"]:
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
    print("x={" + x[0] + ", " + ''.join([xi + ", " for xi in x[2:-2]]) + x[-2] + "};")
    print("y={" + y[0] + ", " + ''.join([yi + ", " for yi in y[2:-2]]) + y[-2] + "};")