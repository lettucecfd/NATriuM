import numpy as np

with open("/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/newNonUni/NACA0012_0deg.msh", 'r') as f:
  string = f.read()
  string = string.split("$Elements")[1]
  string = string.split("$EndElements")[0]
  #print(string)
  #a = np.array([np.array([int(n) for n in s.split(' ') if n != '']) for s in string.splitlines() if s != ''])
  a = []
  for s in string.splitlines():
    line = [int(n) for n in s.split(' ') if n != '']
    if len(line) == 5:
      a.append(line)
  a = np.array(a)
  # a = np.array([np.array([int(n) for n in s.split(' ') if n != '']) for s in string.splitlines() if len(s.split(' ')[s.split(' ') != '']) == 4])
  print(a)

