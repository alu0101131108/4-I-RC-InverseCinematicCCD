#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Robótica Computacional - Curso 2014/2015
# Grado en Ingeniería Informática (Cuarto)
# Práctica: Resolución de la cinemática inversa mediante CCD
#           (Cyclic Coordinate Descent).

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt
import colorsys as cs

# ******************************************************************************
# Declaración de funciones

def muestra_origenes(O,final=0):
  # Muestra los orígenes de coordenadas para cada articulación
  print('Origenes de coordenadas:')
  for i in range(len(O)):
    print('(O'+str(i)+')0\t= '+str([round(j,3) for j in O[i]]))
  if final:
    print('E.Final = '+str([round(j,3) for j in final]))

def muestra_robot(O,obj):
  # Muestra el robot graficamente
  plt.figure(1)
  plt.xlim(-L,L)
  plt.ylim(-L,L)
  T = [np.array(o).T.tolist() for o in O]
  for i in range(len(T)):
    plt.plot(T[i][0], T[i][1], '-o', color=cs.hsv_to_rgb(i/float(len(T)),1,1))
  plt.plot(obj[0], obj[1], '*')
  plt.show()
  input()
  plt.clf()

def matriz_T(d,th,a,al):
  # Calcula la matriz T (ángulos de entrada en grados)
  
  return [[cos(th), -sin(th)*cos(al),  sin(th)*sin(al), a*cos(th)]
         ,[sin(th),  cos(th)*cos(al), -sin(al)*cos(th), a*sin(th)]
         ,[      0,          sin(al),          cos(al),         d]
         ,[      0,                0,                0,         1]
         ]

def cin_dir(th,a):
  #Sea 'th' el vector de thetas
  #Sea 'a'  el vector de longitudes
  T = np.identity(4)
  o = [[0,0]]
  for i in range(len(th)):
    T = np.dot(T,matriz_T(0,th[i],a[i],0))
    tmp=np.dot(T,[0,0,0,1])
    o.append([tmp[0],tmp[1]])
  return o

# ******************************************************************************
# Cálculo de la cinemática inversa de forma iterativa por el método CCD

# valores articulares arbitrarios para la cinemática directa inicial
th=[0.,0.,0.]
a =[5.,1.,5.]

# rango de valores máximo para cada articulación de rotación
range_th =  [[-0.785398, 0.785398],
           [-0.785398, 0.785398],
           [-0.785398, 0.785398]]

# rango de valores máximo para cada articulación prismática
range_a =  [[0, 10.],
          [0, 10.],
          [0, 10.]]

# indicar si las articulaciones son de rotación (0) o prismaticas (1)
articulation_type = [0, 1, 0]

L = sum(a) # variable para representación gráfica
EPSILON = .01 # error
plt.ion() # modo interactivo

# introducción del punto para la cinemática inversa
if len(sys.argv) != 3:
  sys.exit("python " + sys.argv[0] + " x y")
objetivo=[float(i) for i in sys.argv[1:]]

O=list(range(len(th)+1)) # Reservamos estructura en memoria
O[0]=cin_dir(th,a) # Calculamos la posicion inicial
print ("- Posicion inicial:")
muestra_origenes(O[0])

dist = float("inf")
prev = 0.
iteracion = 1
while (dist > EPSILON and abs(prev-dist) > EPSILON/100.):
  prev = dist
  # Para cada combinación de articulaciones:
  for i in range(len(th)):
    # cálculo de la cinemática inversa:
    print("- Objetivo", objetivo)
    print("- Posiciones", O[i])
    print("- Angulos", th)
    print("- Longitudes", a)

    targetIndex = len(th) - i - 1
    last = O[i][-1]
    target = O[i][targetIndex]

    # Articulacion rotativa.
    if articulation_type[targetIndex] == 0:
      # Target to Last and Target to Objective vectors.
      TL = np.subtract(last, target)
      TO = np.subtract(objetivo, target)

      # Angle needed to correct thetha.
      alfa = atan2(TO[1], TO[0]) - atan2(TL[1], TL[0])
      
      # Apply correction.
      th[targetIndex] += alfa

      # Check constraints.
      th[targetIndex] = th[targetIndex] % (pi * 2)
      if th[targetIndex] < range_th[targetIndex][0]:
        th[targetIndex] = range_th[targetIndex][0]
      elif th[targetIndex] > range_th[targetIndex][1]:
        th[targetIndex] = range_th[targetIndex][1]


      # angulo = th[targetIndex] % (pi * 2)
      # if angulo < range_th[targetIndex][0]:
      #   angulo = range_th[targetIndex][0]
      # elif angulo > range_th[targetIndex][1]:
      #   angulo = range_th[targetIndex][1]

    # Articulacion prismática.
    elif articulation_type[targetIndex] == 1:
      # Angle associated to the prismatic articulation direction.
      omega = 0 
      for j in range(targetIndex):
        if articulation_type[j] == 0:
          omega += th[j]

      # Calculate distance (d) to correct by projecting LastToObjetive vector (LO)
      # into the prismatic articulation direction vector (PAD).
      PAD = [cos(omega), sin(omega)]  
      LO = np.subtract(objetivo, last)
      d = np.dot(PAD, LO)

      # Apply correction
      a[targetIndex] += d

      # Check constraints.
      if a[targetIndex] < range_a[targetIndex][0]:
        a[targetIndex] = range_a[targetIndex][0]
      elif a[targetIndex] > range_a[targetIndex][1]:
        a[targetIndex] = range_a[targetIndex][1]

    O[i+1] = cin_dir(th,a)

  dist = np.linalg.norm(np.subtract(objetivo,O[-1][-1]))
  print ("\n- Iteracion " + str(iteracion) + ':')
  muestra_origenes(O[-1])
  muestra_robot(O,objetivo)
  print ("Distancia al objetivo = " + str(round(dist,5)))
  iteracion+=1
  O[0]=O[-1]

if dist <= EPSILON:
  print ("\n" + str(iteracion) + " iteraciones para converger.")
else:
  print ("\nNo hay convergencia tras " + str(iteracion) + " iteraciones.")
print ("- Umbral de convergencia epsilon: " + str(EPSILON))
print ("- Distancia al objetivo:          " + str(round(dist,5)))
print ("- Valores finales de las articulaciones:")
for i in range(len(th)):
  print ("  theta" + str(i+1) + " = " + str(round(th[i],3)))
for i in range(len(th)):
  print ("  L" + str(i+1) + "     = " + str(round(a[i],3)))