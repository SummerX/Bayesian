#!/usr/bin/python
import numpy as np

def J3HNHa(phi=0,psi=0):
    results=7.09*(np.cos(np.radians(phi-60)))**2-1.42*np.cos(np.radians(phi-60))+1.55
    return results

def J3HNC(phi=0,psi=0):
    results=4.29*(np.cos(np.radians(phi+180)))**2-1.01*np.cos(np.radians(phi+180))
    return results

def J3HaC(phi=0,psi=0):
    results=3.72*(np.cos(np.radians(phi+120)))**2-2.18*np.cos(np.radians(phi+120))+1.28
    return results

def J3CC(phi=0,psi=0):
    results=1.36*(np.cos(np.radians(phi)))**2-0.93*np.cos(np.radians(phi))+0.60
    return results

def J3HNCb(phi=0,psi=0):
    results=3.06*(np.cos(np.radians(phi+60)))**2-0.74*np.cos(np.radians(phi+60))+0.13
    return results

def J1NCa(phi=0,psi=0):
    results=1.70*(np.cos(np.radians(psi)))**2-0.98*np.cos(np.radians(psi))+9.51
    return results

def J2NCa(phi=0,psi=0):
    results=-0.66*(np.cos(np.radians(psi)))**2-1.52*np.cos(np.radians(psi))+7.85
    return results

def J3HNCa(phi=0,psi=0):
    results=-0.23*np.cos(np.radians(phi))-0.20*np.cos(np.radians(psi))+0.07*np.sin(np.radians(phi))+0.08*np.sin(np.radians(psi))+0.07*np.cos(np.radians(phi))*np.cos(np.radians(psi))+0.12*np.cos(np.radians(phi))*np.sin(np.radians(psi))-0.08*np.sin(np.radians(phi))*np.cos(np.radians(psi))-0.14*np.sin(np.radians(phi))*np.sin(np.radians(psi))+0.54
    return results