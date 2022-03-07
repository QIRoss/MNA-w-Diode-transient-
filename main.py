import sys
import numpy as np
import math
import matplotlib.pyplot as plt

def solve(Gn, I ):
  e2 = np.linalg.solve(Gn,I)
  return e2

#open netlist text file given as argv[1]
f = open(sys.argv[1])
#split file in a list of strings using 
textFile = f.read().split('\n')
#remove empty lines from list
textFile = list( filter(lambda x:x !='', textFile) )
#remove comments from list
netlist = list( filter(lambda x:x[0] !='*', textFile))

simulationTime = 3
deltaT         = 0.0005
epsilon        = 0.0001

iterations     = int(simulationTime/deltaT)

vc = np.zeros(iterations)
ic = np.zeros(iterations)

#count total number of nodes
nodeCount = 0
#count total number of currents
currentCount = 0

prevActiveCurrentList = 0
c = 0
#map each string from list into an array of netlist parameters as strings
for i in range(len(netlist)):
  netlist[i] = netlist[i].split()

#transform every parameter into integer and floats and count number of nodes
for i in range(len(netlist)):
  if(netlist[i][0][0] == 'D'):
    aux = netlist[i]
    aux[1] = np.uintc(aux[1])
    aux[2] = np.uintc(aux[2])
    aux[3] = np.double(aux[3])
    aux[4] = np.double(aux[4])
  
  elif(netlist[i][0][0] == 'R'):
    aux = netlist[i]
    aux[1] = np.uintc(aux[1])
    aux[2] = np.uintc(aux[2])
    aux[3] = np.double(aux[3])
    #swap nodes if a > b
    if(aux[1] > aux[2]):
      aux[1], aux[2] = aux[2], aux[1]
    #count node value
    if(aux[2] > nodeCount):
      nodeCount = aux[2]

  elif(netlist[i][0][0] == 'V'):
    aux = netlist[i]
    aux[1] = np.uintc(aux[1])
    aux[2] = np.uintc(aux[2])
    aux[4] = np.double(aux[4])
    aux[5] = np.double(aux[5])
    aux[6] = np.double(aux[6])
    aux[7] = np.double(aux[7])
    #count node value
    if(aux[2] > nodeCount):
      nodeCount = aux[2]
    currentCount += 1

  elif(netlist[i][0][0] == 'C'):
    aux = netlist[i]
    aux[1] = np.uintc(aux[1])
    aux[2] = np.uintc(aux[2])
    aux[3] = np.double(aux[3])
    aux[4] = np.double(aux[4])
    c = netlist[i]
    vc[0] = aux[4]
    #count node value
    if(aux[2] > nodeCount):
      nodeCount = aux[2]
    

Gn = np.zeros((nodeCount+1+currentCount, nodeCount+1+currentCount))
I = np.zeros(nodeCount+1+currentCount)

results = np.zeros((iterations, nodeCount + currentCount))


for actualIteration in range(1,iterations):
  previousIteration = 0
  for nrIteration in range(30):
    Gn = np.zeros((nodeCount+1+currentCount, nodeCount+1+currentCount))
    I = np.zeros(nodeCount+1+currentCount)
    t = actualIteration * deltaT
    actualCurrent = 1
    for i in range(len(netlist)):
      aux = netlist[i]
      #insert resistor stamp
      if(aux[0][0] == 'D'):
        a = aux[1]
        b = aux[2]
        Is = aux[3]
        nVt = aux[4]
        vn = 0.1
        if nrIteration == 0:
          vn = np.random.rand()
        else:
          vn = previousIteration[a-1] - previousIteration[b-1]
          #print(previousIteration,previousIteration[a],previousIteration[b],a,b,vn)
        exponential = np.exp(vn/nVt)
        G0 = Is * exponential /nVt
        I0 = Is * (exponential - 1) - G0 * vn
        #add diode current source stamp
        I[a] -= I0
        I[b] += I0
        #add diode resistance source stamp
        Gn[a][a] += G0
        Gn[a][b] -= G0
        Gn[b][a] -= G0
        Gn[b][b] += G0
          
      elif(aux[0][0] == 'R'):
        a = aux[1]
        b = aux[2]
        conductance = np.double(1/aux[3])
        Gn[a][a] += conductance
        Gn[a][b] -= conductance
        Gn[b][a] -= conductance
        Gn[b][b] += conductance

      #insert V stamp
      elif(aux[0][0] == 'V'):
        a = aux[1]
        b = aux[2]
        omega = 2 * math.pi * aux[6]
        phase = aux[7] * math.pi / 180
        value = aux[4] + aux[5] * math.cos(omega * t + phase + math.pi/2)
        I[nodeCount + actualCurrent] -= value
        Gn[nodeCount + actualCurrent ][a] -= 1
        Gn[nodeCount + actualCurrent ][b] += 1
        Gn[a][nodeCount + actualCurrent]  += 1
        Gn[b][nodeCount + actualCurrent]  -= 1
        actualCurrent += 1
        
      #insert capacitor stamp
      elif(aux[0][0] == 'C'):
        a = aux[1]
        b = aux[2]
        conductance = 2 * aux[3]/deltaT
        it0 = ic[actualIteration]
        Gn[a][a] += conductance
        Gn[a][b] -= conductance
        Gn[b][a] -= conductance
        Gn[b][b] += conductance
        vt0 = vc[actualIteration - 1]
        it0 = ic[actualIteration - 1]
        trapezoidalMethodCurrent = conductance * vt0 + it0
        I[a] += trapezoidalMethodCurrent
        I[b] -= trapezoidalMethodCurrent
        

    #remove ground line/column
    Gn = Gn[1:,1:]
    I = I[1:]
    
    e = solve(Gn,I)
    #breaking newton raphson iterations:
    if nrIteration > 1:
      if abs(e[0] - previousIteration[0]) < epsilon and abs(e[1] - previousIteration[1]) < epsilon:
        #print('breaking at:',nrIteration)
        break
    previousIteration = e
  if c != 0:
    vc[actualIteration] = e[1]
    ic[actualIteration] = e[1]*2*c[3]/deltaT-((2*c[3]/deltaT)*vc[actualIteration-1]+ic[actualIteration -1])
  results[actualIteration] = e

e1 = [i[0] for i in results]
e2 = [i[1] for i in results]

timeline = [deltaT * i for i in range(iterations)]

plt.plot(timeline, e1, color='red')
plt.plot(timeline, e2, color='blue')
plt.savefig('Result_' + sys.argv[1][0:8] + '.png')