import os
import sys 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel
from numpy.random import rand


path_in = '../OUTPUT_FILES/'

nrecs_x = 72
nrecs_y = 1
nrecs = nrecs_x * nrecs_y

tsteps = 2001
dt = 0.0003

#****************************************************************************
# input
#****************************************************************************
data = np.zeros ((nrecs_x, tsteps))
dat_type = 'semv'
comp = 'FXZ'
sn = 1
prec = 'single'

num_traces = [1, 36, 72]
#print(num_traces)

index = 0
for rec_x in range(0,nrecs_x):
#for rec_x in num_traces[:]: 
    file_name_in = path_in + 'DB.X' + str(rec_x+1) + '.' + comp + '.' + dat_type
    xz = np.fromfile(file_name_in, dtype=np.float32, count=-1, sep='')
    data[index,:] = xz
    index += 1
    #print(index, rec_x)

norm = np.linalg.norm(data)
data = data / norm

x1 = 0
x2 = nrecs_x
t1 = 0
t2 = tsteps


##****************************************************************************
# Plot pcolor/pcolormesh
##****************************************************************************
# number of plots (rows and columns)
Nr = len(num_traces)
Nc = 1

fig = figure(figsize=(8, 6))
cmap = cm.seismic

# title
figtitle = 'Vz'
t = fig.text(0.5, 0.95, figtitle,
             horizontalalignment='center',
             fontproperties=FontProperties(size=18))

# create axis
t_axis = np.arange(t1*dt,t2*dt,dt)

#plt.plot(t_axis, data[1,:])
#plt.xlabel('Time (s)')
#plt.ylabel('Amplitude')
#plt.grid()
#show()
#sys.exit("Error message")
#print(len(t_axis))

w = 0.8
h = 0.25
ax = []
images = []
vmin = 1e40
vmax = -1e40

for i in range(Nr):
    for j in range(Nc):
        pos = [0.15 + j*1.1*w, 0.1 + i*1.1*h, w, h]
        
        a = fig.add_axes(pos)
        if i > 0:
            a.set_xticklabels([])
        
        dd = ravel(data)
        print(num_traces[i])
        images.append(plt.plot(t_axis, data[(num_traces[i])-1,:]))
        plt.grid()
        plt.ylim((-2e-2,2e-2))        

        ## labels
        if i == 0:
            plt.xlabel('Time (s)', fontsize=16)
        plt.ylabel('Amplitude', fontsize=16)
        ax.append(a)

#show()
#sys.exit("Error message")

# Set the first image as the master, with all the others
# observing it for changes in cmap or norm.

class ImageFollower:
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())



# save fig
path_out = '../fig'
if not os.path.exists(path_out):
    os.makedirs(path_out)

plt.savefig(path_out + '/traces_' + dat_type + '_s' + str(sn) + '_rec_slice' + '.png',figsize=(8,4),dpi=300)

show()
