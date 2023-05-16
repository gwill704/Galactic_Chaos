#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import random as rd
import os

class Chaos:
    def __init__(self, x0, y0, z0, dt, b, c, rc, E, omega):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dt = dt
        self.b = b
        self.c = c
        self.rc = rc
        self.E = E
        self.omega = omega
    def v0(self,x, y, z):
        return np.sqrt((2*self.E)/(1 + np.log(np.abs((self.x0**2+(self.y0/self.b)**2+(self.z0/self.c)**2)/(self.rc**2)))))
    
    def v0_w(self,x, y, z):
        return np.sqrt((2*self.E + self.x0**2 * self.omega**2)/(1 + np.log(np.abs((self.x0**2+(self.y0/self.b)**2+(self.z0/self.c)**2)/(self.rc**2)))))
    
    def v(self, vx, vy, vz):
        return np.sqrt(vx**2+vy**2+vz**2)

    def phi(self, x, y, z):
        return (self.v0(self.x0,self.y0,self.z0)**2/2)*np.log(np.abs((x**2+(y/self.b)**2+(z/self.c)**2)/self.rc**2))
    
    def phi_w(self, x, y, z):
        return (self.v0_w(self.x0,self.y0,self.z0)**2/2)*np.log(np.abs((x**2+(y/self.b)**2+(z/self.c)**2)/self.rc**2))

    def a(self, x, y, z):
        dphix = (self.v0(self.x0,self.y0,self.z0)**2)*(x/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        dphiy = (self.v0(self.x0,self.y0,self.z0)**2)*((y/self.b)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2)) 
        dphiz = (self.v0(self.x0,self.y0,self.z0)**2)*((z/self.c)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        grad = [-dphix, -dphiy, -dphiz]
        return grad
    
    def a_w(self, x, y, z):
        dphix = (self.v0_w(self.x0,self.y0,self.z0)**2)*(x/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        dphiy = (self.v0_w(self.x0,self.y0,self.z0)**2)*((y/self.b)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2)) 
        dphiz = (self.v0_w(self.x0,self.y0,self.z0)**2)*((z/self.c)/np.abs(x**2+(y/self.b)**2+(z/self.c)**2))
        grad = [-dphix, -dphiy, -dphiz]
        return grad

    def verlet(self, k):
        
        x = []
        y = []
        z = []
        En = []
        Lzm = []
        t = []
        vx = []
        vy = []
        vz = []
        ace = []
        poinx = []
        poinvx = []
        
        # Creation of .txt files to use gnuplot
        orbit = 'orbits b = '+str(self.b)+' c = '+str(self.c)+' E = '+str(self.E)
        os.makedirs(orbit, exist_ok = True)
        g = open(orbit+'/poin'+str(k)+'.txt', 'w')
        f = open(orbit+'/datos'+str(k) +'.txt', 'w')
        f.write(str(self.x0) + str('\t') + str(self.y0) + str('\t') + str(self.z0) + str('\n'))
        
        
        # Initial velocities and initial Lz
        
        Lz = rd.uniform(0,np.sqrt(self.v0(self.x0, self.y0, self.z0)*self.x0**2))
        vx0 = 0
        vy0 = Lz / self.x0
        vz0 = np.sqrt((-vy0**2 + self.v0(self.x0,self.y0,self.z0)**2))
        
        
        x.append(self.x0 + vx0*self.dt + 0.5*self.a(self.x0, self.y0, self.z0)[0]*self.dt**2)
        y.append(self.y0 + vy0*self.dt + 0.5*self.a(self.x0, self.y0, self.z0)[1]*self.dt**2)
        z.append(self.z0 + vz0*self.dt + 0.5*self.a(self.x0, self.y0, self.z0)[2]*self.dt**2)
        f.write(str(x[0]) + str('\t') + str(y[0]) + str('\t') + str(z[0]) + str('\n'))

        vx.append(vx0 + 0.5*self.a(x[0], y[0], z[0])[0]*self.dt)
        vy.append(vy0 + 0.5*self.a(x[0], y[0], z[0])[1]*self.dt)
        vz.append(vz0 + 0.5*self.a(x[0], y[0], z[0])[2]*self.dt)
        
        Lzm.append(x[0]*vy[0] - y[0]*vx[0])
        
        En.append(1/2*self.v(vx[0], vy[0], vz[0])**2 + self.phi(x[0], y[0], z[0]))

        t.append(0)

        for i in range(5000):

            x.append(x[i] + vx[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[0]*self.dt**2)
            y.append(y[i] + vy[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[1]*self.dt**2)
            z.append(z[i] + vz[i]*self.dt + 0.5*self.a(x[i], y[i], z[i])[2]*self.dt**2)
            f.write(str(x[i+1]) + str('\t') + str(y[i+1]) + str('\t') + str(z[i+1]) + str('\n'))

            vx.append(vx[i] + 0.5*(self.a(x[i], y[i], z[i])[0])*self.dt)
            vy.append(vy[i] + 0.5*(self.a(x[i], y[i], z[i])[1])*self.dt)
            vz.append(vz[i] + 0.5*(self.a(x[i], y[i], z[i])[2])*self.dt)
            '''vy.append((Lz[i]-y[i+1]*vx[i+1])/x[i+1])
            vz.append(np.sqrt((2*(self.E-self.phi(x[i+1], y[i+1], z[i+1])))-vx[i+1]**2-vy[i+1]**2))'''

            En.append(1/2*self.v(vx[i+1], vy[i+1], vz[i+1])**2 + self.phi(x[i+1], y[i+1], z[i+1]))

            Lzm.append(x[i+1]*vy[i+1] - y[i+1]*vx[i+1])
            
            t.append(self.dt*i)

            if (y[i-1]<0.) and (y[i]>0.) and (z[i] > 0.):
                poinx.append(x[i])
                poinvx.append(vx[i])
                g.write(str(x[i]) + str('\t') +str(vx[i]) + str('\n'))
        
        f.close()
        g.close()
        
        # Remove ''' if you want display of the iterations
        '''print('PARTICULE:' , k)

        print('\nLz_inicial = ',Lzm[0])
        print('\nE_inicial = ',En[0])
                
        print('\nLz_mid = ',Lzm[25000])
        print('\nE_mid = ',En[25000])
        
        print('\nLz_ final = ',Lzm[50000])
        print('\nE_final = ',En[50000])'''
            
        pos = [x,y,z]
        vel = [vx,vy,vz]
        energy = [En, t, Lzm]
        poin = [poinx, poinvx]
        return pos, vel, energy, poin
    
    def angularspeed(self, k):
        x = []
        y = []
        z = []
        En = []
        Lzm = []
        t = []
        vx = []
        vy = []
        vz = []
        ace = []
        poinx = []
        poinvx = []
        
        # Creation of .txt files to use gnuplot
        orbit_w = 'orbits_with_angular_speed b = '+str(self.b)+' c = '+str(self.c)+' E = '+str(self.E)
        os.makedirs(orbit_w, exist_ok = True)
        g = open(orbit_w+'/poin'+str(k)+'.txt', 'w')
        f = open(orbit_w+'/datos'+str(k) +'.txt', 'w')
        f.write(str(self.x0) + str('\t') + str(self.y0) + str('\t') + str(self.z0) + str('\n'))
        
        
        # Initial velocities and initial Lz
        
        Lz = rd.uniform(0,np.sqrt((2*(self.E-self.phi(self.x0,self.y0,self.z0)**2))))
        print('Lz = ', Lz)
        vx0 = 0.
        vy0 = Lz / self.x0
        vz0 = np.sqrt((-vy0**2 + ((2*self.E + self.x0**2*self.omega**2)/(1+np.log(np.abs(self.x0**2/self.rc**2))))))
        
        
        x.append(self.x0 + vx0*self.dt + 0.5*self.a_w(self.x0, self.y0, self.z0)[0]*self.dt**2)
        y.append(self.y0 + vy0*self.dt + 0.5*self.a_w(self.x0, self.y0, self.z0)[1]*self.dt**2)
        z.append(self.z0 + vz0*self.dt + 0.5*self.a_w(self.x0, self.y0, self.z0)[2]*self.dt**2)
        f.write(str(x[0]) + str('\t') + str(y[0]) + str('\t') + str(z[0]) + str('\n'))

        vx.append(vx0 + 0.5*self.a_w(x[0], y[0], z[0])[0]*self.dt)
        vy.append(vy0 + 0.5*self.a_w(x[0], y[0], z[0])[1]*self.dt)
        vz.append(vz0 + 0.5*self.a_w(x[0], y[0], z[0])[2]*self.dt)
        
        Lzm.append(x[0]*vy[0] - y[0]*vx[0])
        
        En.append(1/2*self.v(vx[0], vy[0], vz[0])**2 + self.phi_w(x[0], y[0], z[0])-1/2*(x[0]**2+y[0]**2+z[0]**2)*self.omega**2)
        

        t.append(0)

        for i in range(50000):

            x.append(x[i] + vx[i]*self.dt + 0.5*self.a_w(x[i], y[i], z[i])[0]*self.dt**2)
            y.append(y[i] + vy[i]*self.dt + 0.5*self.a_w(x[i], y[i], z[i])[1]*self.dt**2)
            z.append(z[i] + vz[i]*self.dt + 0.5*self.a_w(x[i], y[i], z[i])[2]*self.dt**2)
            f.write(str(x[i+1]) + str('\t') + str(y[i+1]) + str('\t') + str(z[i+1]) + str('\n'))

            vx.append(vx[i] + 0.5*(self.a_w(x[i], y[i], z[i])[0])*self.dt)
            vy.append(vy[i] + 0.5*(self.a_w(x[i], y[i], z[i])[1])*self.dt)
            vz.append(vz[i] + 0.5*(self.a_w(x[i], y[i], z[i])[2])*self.dt)

            En.append(1/2*self.v(vx[i+1], vy[i+1], vz[i+1])**2 + self.phi_w(x[i+1], y[i+1], z[i+1])-1/2*((x[i+1]**2+y[i+1]**2+z[i+1]**2)*self.omega**2))

            Lzm.append(x[i+1]*vy[i+1] - y[i+1]*vx[i+1])
            
            t.append(self.dt*i)

            if (y[i-1]<0.) and (y[i]>0.) and (z[i] > 0.):
                poinx.append(x[i])
                poinvx.append(vx[i])
                g.write(str(x[i]) + str('\t') +str(vx[i]) + str('\n'))
        
        f.close()
        
        # Remove ''' if you want display of the iterations
        '''
        print('PARTICULE:' , k)

        print('\nLz_inicial = ',Lzm[0])
        print('\nE_inicial = ',En[0])
                
        print('\nLz_mid = ',Lzm[25000])
        print('\nE_mid = ',En[25000])
        
        print('\nLz_ final = ',Lzm[50000])
        print('\nE_final = ',En[50000])
            '''
        pos = [x,y,z]
        vel = [vx,vy,vz]
        energy = [En, t, Lzm]
        poin = [poinx, poinvx]
        return pos, vel, energy, poin
        
    