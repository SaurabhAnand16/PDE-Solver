import matplotlib.pyplot as plt
import numpy as np
from csv import reader


x_values = np.array([])
y_values = np.array([])
u_values = np.array([])

fpath = ""  # Enter relative filepath of the csv file here.

with open('fname','r') as f:
    csv_reader = reader(f)
    
    for row in csv_reader:
        if row[0] != 'x':
         x_values = np.append(x_values,float(row[0]))
        if row[1] != 'y':
         y_values = np.append(y_values,float(row[1]))
        if row[2] != 'z':
         u_values = np.append(u_values,float(row[2]))   
         

X = x_values.reshape((51,51))
       
Y = y_values.reshape((51,51))

Z = u_values.reshape((51,51))
print(Z)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

#ax.plot_surface(X,Y,Z)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap="Spectral", linewidth=0, antialiased=True)
ax.set(xlabel = 'X', ylabel='Y',zlabel='U(x,y)', title = 'Surface of Function U vs x & y')
plt.show()
