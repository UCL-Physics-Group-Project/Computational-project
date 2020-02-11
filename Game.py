#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 21:51:52 2020

@author: theoreignier
"""

###################   Importing data   #########################################################################################

#Importing the numpy library
import numpy as np

#Importing functions from the vpython library
from vpython import sphere, color, rate, canvas, vector, curve, label, box , cross, mag, random

####################   ANGRY BIRDS   ############################################################################################################################################

### Setting up initial data ###

g = 9.81 
mass_bird = 0.1
radius_bird = 0.5
t_applied = 0
p = 200
contact_time = 0.01 
dt = 0.01
d_angle = 0.01
x0 = 0.0
y0 = 0.0
n = 1
r = 1


### Defining the function of the trajectory y(x) ###

def trajectory_bird(x_pos):
    ''''''
    trajectory_y = -g*(x_pos-x0)**2 / (2*v0**2*np.cos(theta)**2) + np.tan(theta)*(x_pos-x0) + y0
    
    return trajectory_y


### Main code of the game ### 

while n < 10 :

    ## Setting up the scene ##
    scene = canvas(width=640, height=480, center=vector(8,0,0), range=8)
    scene.camera.pos = vector(4,2,12)
    scene.forward = vector(8,0,0) - scene.camera.pos
    scene.camera.axis = vector(8,0,0) - scene.camera.pos
    ground = box(pos=vector(8,0,0), length=16, height=0.1, width=7, color=color.green)
    target = box(pos=vector((2*random()+1)*5,1,0), axis=vector(1,0,0), length=0.5, height=2, width=0.5)
    label_initial = label(pos=vector(0,0,0), xoffset= 75, yoffset= 150)
    
    ## Initial data with respect to the scene set up ##
    mass_target = p*(target.length*target.width*target.height)
    hit_tolerance = 0.2 + target.length/2
    x_impact = mag(target.pos-vector(0,1,0)) - hit_tolerance
    t_restoring = mass_target*g*target.length/2
    
    print ('-------------------- LEVEL', n, '--------------------')
    
    while t_applied <= t_restoring :
    
        print('------ Lauch', r, '------')
        
        bird = sphere(pos = vector(x0,y0,0), radius = 0.2, color=color.red, make_trail=True)
        theta_deg = float(input("Input the launch angle in degrees: "))
        v0 = float(input("Input the speed in m/s: "))
        theta = np.radians(theta_deg)
        label_initial.text='Last launch initial values: \n speed={0:0.1f} m/s \n angle={1:0.1f} deg'.format(v0, theta_deg)
    
        if 2 + 0.2 > trajectory_bird(x_impact) > 0 :
            
            x = x0
            y = y0
            t = 0
            
            while x <= x_impact :
                rate(40)
                t += dt
                x = x0 + v0*t*np.cos(theta) 
                y = y0 + v0*t*np.sin(theta) - 0.5*g*t**2 
                bird.pos = vector(x,y,0)
            
            c_rot = target.pos + vector(target.length/2,-target.height/2,0)
            bird_momentum = vector(mass_bird*v0*np.cos(theta) , mass_bird*v0*np.sin(theta) - mass_bird*g*t , 0)
            vector_da = vector(x_impact,trajectory_bird(x_impact),0) - c_rot
            t_applied_vector = cross(bird_momentum/contact_time,vector_da)
            t_applied = mag(t_applied_vector)
            
            if t_applied > t_restoring :
                
                s = target.pos
                w = bird.pos
                alpha = 0
                betha = np.pi - np.arccos((target.length/2)/mag(s-c_rot))
                gamma = np.pi - np.arccos((hit_tolerance + target.length/2) / mag(w-c_rot))
                
                x = x_impact
                y = trajectory_bird(x_impact)
                t = 0
                
                while alpha >= -np.pi/2:
                    rate(100)
                    betha -= d_angle
                    alpha -= d_angle
                    gamma -= d_angle
                    o = vector(np.cos(alpha),np.sin(alpha),0)
                    u = c_rot + vector(np.cos(betha),np.sin(betha),0)*mag(s-c_rot)
                    z = c_rot + vector(np.cos(gamma),np.sin(gamma),0)*mag(w-c_rot)
                    target.pos = u
                    target.axis = o
                    bird.pos = z
                    target.length = 0.5
                    
            else :
                x = x_impact
                y = trajectory_bird(x_impact)
                t = 0
                while y > 0 :
                    rate(40)
                    t += dt
                    y = -g*t**2 + trajectory_bird(x_impact)
                    bird.pos = vector(x,y,0)
                    
        else :
                                  
            x = x0
            y = y0
            t = 0
            
            while x < ground.length and y >= 0:
                rate(40)
                t += dt
                x = x0 + v0*t*np.cos(theta)
                y = y0 + v0*t*np.sin(theta) - 0.5*g*t**2 
                bird.pos = vector(x,y,0)
                    
        r += 1
        
    print('Congratulations! You finished level', n, '.')
    t_applied = 0
    r = 1
    n += 1

print('Amazing! You finished the game!')

##################################################################################################################################################################################