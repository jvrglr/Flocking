
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import moviepy.editor as mpy
from moviepy.editor import VideoClip
from math import sqrt
from moviepy.video.io.bindings import mplfig_to_npimage
np.random.seed(1994)
#import pdb; pdb.set_trace() #debugger

def initialize_positions(N):
    """
    generate N random 2D positions (x,y) in squared canvas of [-lat/2:lat/2]x[-lat/2:lat/2]
    """
    mat=np.zeros((N,2))
    for i in range(N):
        lat=400.0
        mat[i][0]=np.random.rand()*lat-lat/2.0
        mat[i][1]=np.random.rand()*lat-lat/2.0
    return mat

def initialize_velocity(N):
    """
    generate N random 2D positions (vx,vy) in [-1,1]X[-1,1] (uniformly distrinuted)
    """
    mat=np.zeros((N,2))
    for i in range(N):
        v_cota=100.0
        mat[i][0]=np.random.rand()*v_cota-v_cota/2.0
        mat[i][1]=np.random.rand()*v_cota-v_cota/2.0
    return mat

def cohesion(state,N,i):
    """
    The speed of this movement is parametrized by the speed_c parameter (higher speed as speed_c is lower).
    """
    r_cm=np.zeros(2)
    v= np.zeros(2)
    ri=np.array([state[i][1],state[i][2]])
    contador=0
    radius=50
    #CONTROLZ HASTA AQUÍ
    #-------------compute center of mass, each void has the same mass -------------------------------------
    for j in range(N):
        if j!=i:
            contador=contador+1
            rj=np.array([state[j][1],state[j][2]])
            distance=sqrt((ri[0]-rj[0])**2+(ri[1]-rj[1])**2)
            if distance < radius:
                r_cm=r_cm+rj
    r_cm=r_cm/float(contador)
    #-------------------------------------------------------------------------------
    speed_c=300
    speed=sqrt(r_cm[0]**2+r_cm[1]**2)
    if speed<=0.0:
        v=np.zeros(2)
    else:
        v=(r_cm-ri)/speed_c
    return v

def aligment(state,N,i):
    """
    align speed of boid-i with the averaged "perceived" velocity
    """
    v= np.zeros(2)
    ri=np.array([state[i][1],state[i][2]])
    radius=50
    contador=0
    for j in range(N):
        if j!= i:
            contador=contador+1
            rj=np.array([state[j][1],state[j][2]])
            distance=sqrt((ri[0]-rj[0])**2+(ri[1]-rj[1])**2)
            if distance < radius:
                vj=np.array([state[j][3],state[j][4]])
                v=v+vj
    v=v/float(contador)
    vi=np.array([state[i][3],state[i][4]])
    speed=sqrt(v[0]**2+v[1]**2)
    if speed<=0.0:
        v=np.zeros(2)
    else:
        v=(v-vi)/5.0
    return v


def exclusion(state,N,i):
    """
    Boids within a distance parametrized by radius are repelled.
    """
    ri=np.array([state[i][1],state[i][2]])
    c=np.zeros(2)
    radius=15.0
    for j in range(N):
        rj=np.array([state[j][1],state[j][2]])
        distance=sqrt((ri[0]-rj[0])**2+(ri[1]-rj[1])**2)
        if distance < radius:
            c=c+(ri-rj)/200.0
    return c

def cruise_v(state,i):
    """
    Module of velocity is constant
    """
    vi=np.array([state[i][3],state[i][4]])
    v_cruise=2.0
    speed=sqrt(vi[0]**2+vi[1]**2)
    v=vi*(v_cruise-speed)/speed
    return v

def updating_position (state,N,i):
    """
    Create New position for the boid i
    """
    #I HAVE ADDED /100
    x=state[i][1]+state[i][3]/1.0
    y=state[i][2]+state[i][4]/1.0
    return x,y

def updating_velocity (state,N,i):
    """
    Create New position for the boid i
    """
    #I HAVE ADDED /100
    x=state[i][3]+state[i][5]
    y=state[i][4]+state[i][6]
    return x,y


# Creates blank canvas
fig = plt.figure()
axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#limits of the graph
L=400
axes.set_xlim([-L, L])
axes.set_ylim([-L, L])
axes.axes.get_xaxis().set_ticks([]) #ERASE XTICS ON X_AXE
axes.axes.get_yaxis().set_ticks([])

#duration = 2


N=60 #NUMBER OF VOIDS

state=np.zeros((N,7))
#-------------------------------------------------------------
# state [i][0] = i name of the boid
# state [i][1],state [i][2] = 2D position of the boid i
# state [i][3],state [i][4] = 2D velocity of the boid i
# state [i][5],state [i][6] = 2D acceleration of the boid i
#-------------------------------------------------------------

#------------------INITIAL CONDITION-------------------------------------
aux= initialize_positions(N)
aux2= initialize_velocity(N)
for i in range(N):
    pass
    state[i][0]=i
    state[i][1]=aux[i][0]
    state[i][2]=aux[i][1]
    state[i][3]=aux2[i][0]
    state[i][4]=aux2[i][1]
    state[i][5]=0.0
    state[i][6]=0.0




#-------------------------------------------------------------------------


def make_frame(t):
    axes.clear()

    #wind=np.array([0.1,0.1])
    wind=np.zeros(2)
    saving=np.zeros((N,2))
# ------------------------UPDATING ACCELERATION---------------------------------------
    for i in range(N):
        aux=cohesion(state,N,i)
        #aux=np.zeros(2)
        aux2=exclusion(state,N,i)
        #aux2=np.zeros(2)
        aux3=aligment(state,N,i)
        #aux3=np.zeros(2)
        saving[i][0]=aux[0]+aux2[0]+aux3[0]+0.0*random()
        saving[i][1]=aux[1]+aux2[1]+aux3[1]+0.0*random()
    #UPDATING
    for i in range(N):
        state[i][5]=saving[i][0]
        state[i][6]=saving[i][1]
#---------------------------------------------------------------------------------

    #Updating VELOCITIES
    for i in range(N):
        aux=updating_velocity(state,N,i)
        state[i][3]=aux[0]+wind[0]
        state[i][4]=aux[1]+wind[1]
        aux4=cruise_v(state,i)
        state[i][3]=state[i][3]+aux4[0]
        state[i][4]=state[i][4]+aux4[1]


    #Updating positions
    for i in range(N):
        aux=updating_position(state,N,i)
        state[i][1]=aux[0]
        if ((state[i][1] >= float(L)) or (state[i][1] <= float(-L))): #BOUNDARY CONDITIONS
            state[i][3]=-state[i][3]
            state[i][1]=state[i][1]+2.0*state[i][3]
        state[i][2]=aux[1]
        if ((state[i][2] >= float(L)) or (state[i][2] <= float(-L))):
            state[i][4]=-state[i][4]
            state[i][2]=state[i][2]+2.0*state[i][4]

    #SAVE DATE ON FILE
    #if controls==True:
    #data=pd.DataFrame(state,columns="#N x1 x2 v1 v2".split())
    #data.to_csv('data.dat', header=True, index=False, sep='\t', mode='a')


    #CREATE ARROWS
    for i in range (N):
        pass
        speed=sqrt(state[i][3]**2+state[i][4]**2)
        plt.arrow(state[i][1],state[i][2],state[i][3]/speed,state[i][4]/speed,head_width=10)
    axes.axes.get_xaxis().set_ticks([]) #ERASE XTICS ON X_AXE
    axes.axes.get_yaxis().set_ticks([])
    #fig.savefig("filename.png", dpi=200) # plot canvas
    #plt.show()
    return mplfig_to_npimage(fig)

animation = VideoClip(make_frame, duration=14) # duration= duration of video in seconds
animation.write_videofile("my_animation.mp4", fps=24) # export as video, fps= frames per second

