#####Code for question #3a #####

import numpy as np

#latitude = y
#longitude = x

sf = np.loadtxt("LongIslandData/longdata.txt")
#data format array([latitude,longitude,height])

nx, ny = (100, 100) #size of domain/grid
grid=np.zeros((nx,ny)) #create grid of zeros
grid=grid.tolist() #convert grid to list
x=np.linspace(min(sf[:,0]),max(sf[:,0]),nx) #array of longitudes from min to max of data values
y=np.linspace(min(sf[:,1]),max(sf[:,1]),ny) #array of latitiudes from min to max of data values
dx=0.029780860000002463/2.0 #half the difference between two longitude values
dy=0.0096799899999950867/2.0 #half the difference between two longitude values

frame=10.0 #consider every tenth data points (i.e. use only 10th of data)

#the following loop checks to see if a given latitude/longitude is within a certain range of array y/x
#if it is, then it appends the height value to that grid point
#if that grid point is empty, it fills it with an empty list

for i in range(0,len(sf)):
    if i%frame==0.0:
        for j in range(0,len(x)):
            if np.abs(np.abs(sf[i][0])-np.abs(x[j]))<dx:
                for k in range(0,len(y)):
                    if np.abs(np.abs(sf[i][1])-np.abs(y[j]))<dy:
                        
                        if grid[j][k]==0.0:
                            grid[j][k]=[]
                        else:
                            grid[j][k].append(sf[i][2])
        
grid=np.asarray(grid) #the grid is converted back to array form

#Each grid point has several height values, the average is now taken and
#the average height filled into the grid point
for m in range(0,len(grid)):
    for n in range(0,len(grid[0])):
        grid[m][n]=(np.sum(grid[m][n])/len(grid[m][n]))
        
    print m
    
#the 100 by 100 grid is now filled with height values and saved
#np.savetxt('longisland1.txt', grid)