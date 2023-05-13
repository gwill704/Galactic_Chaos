#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import random as rd

class Chaos:
    def __init__(self, vx0, vy0, vz0, dt, b, c, rc):
        self.vx0 = vx0
        self.vy0 = vy0
        self.vz0 = vz0
        self.dt = dt
        self.b = b
        self.c = c
        self.rc = rc
        self.v0 = self.v(self.vx0, self.vy0, self.vz0)
    
    def v(self, vx, vy, vz):
        return np.sqrt(vx**2+vy**2+vz**2)

    def phi(self, x, y, z):
        return (self.v0**2/2)*np.log(np.abs((x**2+(y/self.b)**2+(z/self.c)**2)/self.rc**2))

    def a(self, x, y, z):
        dphix = (self.v0**2)*(x/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        dphiy = (self.v0**2)*((y/self.b)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2)) 
        dphiz = (self.v0**2)*((z/self.c)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        grad = [-dphix, -dphiy, -dphiz]
        return grad

    def verlet(self, x0, y0, z0, k):
        x = []
        y = []
        z = []
        E = []
        Lz = []
        t = []
        vx = []
        vy = []
        vz = []
        ace = []
        poinx = []
        poinvx = []

        #if (1/2*self.v(self.vx0, self.vy0, self.vz0)**2 + self.phi(x0, y0, z0)) < 0:
        x.append(x0 + self.vx0*self.dt + 0.5*self.a(x0, y0, z0)[0]*self.dt**2)
        y.append(y0 + self.vy0*self.dt + 0.5*self.a(x0, y0, z0)[1]*self.dt**2)
        z.append(z0 + self.vz0*self.dt + 0.5*self.a(x0, y0, z0)[2]*self.dt**2)
        print(x[0], y[0], z[0])


        vx.append(self.vx0 + 0.5*self.a(x0, y0, z0)[0]*self.dt)
        vy.append(self.vy0 + 0.5*self.a(x0, y0, z0)[1]*self.dt)
        vz.append(self.vz0 + 0.5*self.a(x0, y0, z0)[2]*self.dt)

        E.append(1/2*self.v(self.vx0, self.vy0, self.vz0)**2 + self.phi(x0, y0, z0))
        
        Lz.append(x[0]*vy[0] - y[0]*vx[0])

        t.append(0)

        for i in range(50000):

           # if (1/2*self.v(vx[i], vy[i], vz[i])**2 + self.phi(x[i], y[i], z[i])) < 0:
            x.append(x[i] + vx[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[0]*self.dt**2)
            y.append(y[i] + vy[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[1]*self.dt**2)
            z.append(z[i] + vz[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[2]*self.dt**2)

            vx.append(vx[i] + 0.5*self.a(x[i], y[i], z[i])[0]*self.dt)
            vy.append(vy[i] + 0.5*self.a(x[i], y[i], z[i])[1]*self.dt)
            vz.append(vz[i] + 0.5*self.a(x[i], y[i], z[i])[2]*self.dt)

            E.append(1/2*self.v(vx[i], vy[i], vz[i])**2 + self.phi(x[i], y[i], z[i]))

            Lz.append(x[i]*vy[i] - y[i]*vx[i])
            
            t.append(self.dt*i)

            if (y[i-1]<0.) and (y[i]>0.) and (z[i] > 0.):# and Lz[i] == Lz[0]:
                poinx.append(x[i])
                poinvx.append(vx[i])

        print('\nLz_inicial = ',Lz[0])
        print('\nE_inicial = ',E[0])
                
        print('\nLz_mid = ',Lz[25000])
        print('\nE_mid = ',E[25000])
        
        print('\nLz_ final = ',Lz[50000])
        print('\nE_final = ',E[50000])
            
        pos = [x,y,z]
        vel = [vx,vy,vz]
        energy = [E, t, Lz]
        poin = [poinx, poinvx]
        return pos, vel, energy, poin






