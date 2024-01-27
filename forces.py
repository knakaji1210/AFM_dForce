# Import of several modules

import math
import time
import numpy as np
from scipy import integrate
from scipy.signal import sawtooth
from numpy import sqrt,pi,cos,sin,tan,exp


def punch(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    if (z + zc) < a0:                
        I = -z - zc + a0        
        force = 2 * Rt * E_eff * I    
    else:
        force = 0

    return force

# ================================= More forces ====================
def sneddon(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    tan = np.tan(angle)   
    if (z + zc) < a0:                
        I = -z - zc + a0
        force = 2 * E_eff * np.tan(angle) * I**2 / (np.pi)        
        
        
    else:
        force = 0

    return force

def kv_sphere_BEC(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):    
    
    if (z + zc) < a0:
        I = -z - zc + a0
        R = Rt
        
        h = sahe
        """
        eta_G = eta/(2*(1+nu))
        
        force = 16/9*sqrt(R*I)*(-9/2*eta_G*v+E_sample*I)  \
                + 1.133*16/9*1/h*R*I*(-6*eta_G*v+E_sample*I) \
                + 1.497*16/9*1/h**2*sqrt(R*I)**3*(-15/2*eta_G*v+E_sample*I) \
                + 1.469*16/9*1/h**3*(R*I)**2*(-9*eta_G*v+I*E_sample) \
                + 0.755*16/9*1/h**4*sqrt(R*I)**5*(-21/2*eta_G*v+I*E_sample) 
        """
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        force = 4*I**0.5*sqrt(R)*(E_sample*I - 1.5*v*eta)/(3*(1 - nu**2)) \
            -8*I**1.0*R*alpha0*(E_sample*I - 2.0*v*eta)/(3*pi*h*(1 - nu**2))\
            +I**1.5*(E_sample*I - 2.5*v*eta)*(5.33333333333333*R**(3/2)*alpha0**2/pi**2 + 8*R**1.5*alpha0**2/(9*pi**2))/(h**2*(1 - nu**2)) \
            +I**2.0*(E_sample*I - 3.0*v*eta)*(-15.8024691358025*R**2*alpha0**3/pi**3 - 2.84444444444444*R**2*beta0/pi)/(h**3*(1 - nu**2))     \
            +I**2.5*(E_sample*I - 3.5*v*eta)*(39.9012345679012*R**2.5*alpha0**4/pi**4 + 17.0666666666667*R**2.5*alpha0*beta0/pi**2)/(h**4*(1 - nu**2))
     
    else:
        force = 0

    return force

def kv_cone_BEC(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    tan = np.tan(angle)    
    h = sahe
    if (z + zc) < a0:                
        I = -z - zc + a0
        
        """
        eta_G = eta/(2*(1+nu))
        force = 8*tan/(3*pi)*I*(-6*eta_G*v+E_sample*I) \
                    +0.721*8*tan**2/(3*pi*h)*I**2*(-9*eta_G*v+E_sample*I) \
                    +0.650*8*tan**3/(3*pi*h**2)*I**3*(-12*eta_G*v+E_sample*I) \
                    +0.491*8*tan**4/(3*pi*h**3)*I**4*(-15*eta_G*v+E_sample*I) \
                    +0.225*8*tan**5/(3*pi*h**4)*I**5*(-18*eta_G*v+E_sample*I) 
        """
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        force = I**5*(E_sample*I - 6*v*eta)*(1216*alpha0**4*tan**5/pi**9 + 448*alpha0*beta0*tan**5/pi**7)/(h**4*(1 - nu**2)) \
            + I**4*(E_sample*I - 5*v*eta)*(-224*alpha0**3*tan**4/pi**7 - 32*beta0*tan**4/pi**5)/(h**3*(1 - nu**2)) \
            + 40*I**3*alpha0**2*tan**3*(E_sample*I - 4*v*eta)/(pi**5*h**2*(1 - nu**2)) \
            - 8*I**2*alpha0*tan**2*(E_sample*I - 3*v*eta)/(pi**3*h*(1 - nu**2)) \
            + 2*I*tan*(E_sample*I - 2*v*eta)/(pi*(1 - nu**2))

    else:
        force = 0

    return force

def kv_punch_BEC(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    a = Rt
    h = sahe
    
    # This is assuming nu = 0.5
    if (z + zc) < a0:                
        I = -z - zc + a0
        
        """
        if nu == 0.5 and bonded:
            force = 8*a/3*(-eta*v+E_sample*I) \
            + 1.133*8/3*a**2/h*(-eta*v+E_sample*I) \
            + 1.283*8/3*a*(a/h)**2*(-eta*v+E_sample*I)  \
            + 0.598*8/3*a*(a/h)**3*(-eta*v+E_sample*I)  \
            - 0.291*8/3*a*(a/h)**4*(-eta*v+E_sample*I)  
        """
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        force = 2*Rt*(E_sample*I - eta * v)/(1-nu**2) * (1  \
                            - 2*alpha0 * a/(pi*h)  \
                            + (2*alpha0/pi * a/h)**2    \
                            -8*a**3*(3*alpha0**3+pi**2+beta0)/(3*pi**3)* (a/h)**3  
                            -16*alpha0*(3*alpha0**3+2*pi**2*beta0)/(3*pi**4)*(a/h)**4)
            
            
    else:
        force = 0

    return force

def kv_bec_nanowire(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    if (z + zc) < a0:
        
        h = sahe
        I = -z - zc + a0
        a = Rt         
        
        alpha0 = - 1.7795
        
        
        eta_G = eta/3
            
        # Critical indentation between sphere and punch
        I_critical = 3*(-h*pi*(4*Rt*alpha0 - 3*h*pi) - sqrt(3)*sqrt(-h**3*pi**3*(8*Rt*alpha0 - 3*h*pi)))/(8*Rt*alpha0**2)  
        
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
        
        def F_sphere(I,R,v):           
                
            force = 4*I**0.5*sqrt(R)*(E_sample*I - 1.5*v*eta)/(3*(1 - nu**2)) \
                -8*I**1.0*R*alpha0*(E_sample*I - 2.0*v*eta)/(3*pi*h*(1 - nu**2))\
                +I**1.5*(E_sample*I - 2.5*v*eta)*(5.33333333333333*R**(3/2)*alpha0**2/pi**2 + 8*R**1.5*alpha0**2/(9*pi**2))/(h**2*(1 - nu**2)) \
                +I**2.0*(E_sample*I - 3.0*v*eta)*(-15.8024691358025*R**2*alpha0**3/pi**3 - 2.84444444444444*R**2*beta0/pi)/(h**3*(1 - nu**2))     \
                +I**2.5*(E_sample*I - 3.5*v*eta)*(39.9012345679012*R**2.5*alpha0**4/pi**4 + 17.0666666666667*R**2.5*alpha0*beta0/pi**2)/(h**4*(1 - nu**2))
                    
            return force
        
        def F_punch(I,a,v):
                
            force = 2*Rt*(E_sample*I - eta * v)/(1-nu**2) * (1  \
                                - 2*alpha0 * a/(pi*h)  \
                                + (2*alpha0/pi * a/h)**2    \
                                -8*a**3*(3*alpha0**3+pi**2+beta0)/(3*pi**3)* (a/h)**3  
                                -16*alpha0*(3*alpha0**3+2*pi**2*beta0)/(3*pi**4)*(a/h)**4)
            return force
        
        if I<I_critical:
            # Only the a spherical area is in contact,the tip behaves like a sphere
            # a = np.sqrt(Rt*I) - 2*alpha0*Rt*I/(3*h*pi)
            force = F_sphere(I,Rt,v)           
            
        else:
            # All the sphere in contact, the tip is now a punch
            # a = Rt
            force = F_sphere(I_critical,Rt  ,v=0) + F_punch(I-I_critical,Rt, v)
            
                
        
    else:
        force = 0

    return force





# ================================= Use custom forces ====================

def custom1(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    # for the magnetic results
    if (z + zc) < a0:
        force = 0
    else:
        alpha = 1E-45
        force = -alpha / (z + zc)**4
    return force


def custom2(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    if (z + zc) < a0:
        force = 0
    else:
        eta = V1
        force = visco(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
                      landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4)
    return force


# ================================= Long range forces ====================

def vdw(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    #g_n = Rt/6    
    #n = 2
    if (z + zc) < a0:
        force = -H * Rt / (6 * a0**2)
    else:
        force = -H * Rt / (6 * (z + zc) * (z + zc))
    return force

def vdw_cone(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    # Hartmann,1990, van der Waals interactions between sharp probes and Aat sample surfaces
    #g_n = np.tan(angle)**2/6
    #n = 1    
    if (z + zc) < a0:
        force = -H * np.tan(angle)**2/(6*a0)
    else:
        force = -H * np.tan(angle)**2/(6*(z+zc))
    return force

def vdw_punch(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    # Hartmann,1990, van der Waals interactions between sharp probes and Aat sample surfaces
    #g_n = Rt**2/6
    #n = 3
    if (z + zc) < a0:
        force = -H * Rt**2/ (6 * (a0)**3 )
    else:
        force = -H * Rt**2/ (6 * (z + zc)**3 )
    return force


def dlvo(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    epsilon = 8.854187 * 10**-12
    if (z + zc) < a0:
        force = 4 * math.pi * Rt / \
            (epsilon * epsilonr) * sigmat * sigmas * landadeb * np.exp(-(a0) / landadeb) - H * Rt / (6 * a0 * a0)
    else:
        force = -H * Rt / (6 * (z + zc)**2) + 4 * math.pi * Rt / (
            epsilon * epsilonr) * sigmat * sigmas * landadeb * np.exp(-(z + zc) / landadeb) -H * Rt / (6 * (z + zc) * (z + zc))

    return force


# ================================= Short range forces ===================


def hertz(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    if (z + zc) < a0:
        force = 4.0 / 3.0 * E_eff * math.sqrt(Rt) * math.pow(-(z + zc - a0), 3.0 / 2.0)
    else:
        force = 0

    return force


def dmt(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    if (z + zc) < a0:
        force = 4.0 / 3.0 * E_eff * \
            math.sqrt(Rt) * math.pow(-(z + zc - a0), 3.0 / 2.0) - H * Rt / (6 * a0 * a0)
    else:
        force = -H * Rt / (6 * (z + zc) * (z + zc))

    return force


def tatara(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    if ((z + zc) < a0) and nc > 0:
        force = alfa * math.pow(-(z + zc - a0),
                                1.5) / math.pow(2,
                                                1.5) + 3 * alfa**2 * (z + zc - a0)**2 / (8 * nc) + 15 * alfa**3 / (math.pow(2,
                                                                                                                            11 / 2.0) * nc**2) * math.pow(-(z + zc - a0),
                                                                                                                                                          2.5) - H * Rt / (6 * a0 * a0)
    else:
        force = 0

    return force

# Bottom Effect correction for the cone
def becc(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    if (z + zc) < a0:        
        h = sahe
        I = -z - zc + a0
        
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        tan = np.tan(angle) # Tangent of the angle
        cot = 1 / np.tan(angle) # Cotangent of the angle
        a = 2*tan*I/pi * (1-2*alpha0*tan*I/(h*pi**2))
        force = -(a*E_eff*(-4*I+a*pi*cot))/(2) \
            + (a*E_eff*(-4*I+a*pi*cot)*2*alpha0*a)/(2*pi*h) \
            - (a*E_eff*(-4*I+a*pi*cot))/(2) * (2*alpha0*a/(pi*h))**2 \
            + ( a**4*E_eff*(-16*I*(3*alpha0**3+pi**2*beta0)+a*pi*(12*alpha0**3+5*pi**2*beta0)*cot) ) / (3*pi**3*h**3) \
            -(2*a**5*E_eff*alpha0*(-16*I*(3*alpha0**3+2*pi**2*beta0)+3*a*pi*(4*alpha0**3+3*pi**2*beta0)*cot))/(3*pi**4*h**4) 
        
    else:
        force = 0

    return force


def behc(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    if (z + zc) < a0:
        
        h = sahe
        I = -z - zc + a0
        R = Rt 
        
        
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        
        
        a = sqrt(R*I)* (1-2*alpha0*sqrt(R*I)/(3*h*pi))
 
        force = 2*E_eff*(3*a*I*R-a**3)/(3*R) \
            + 2*E_eff*(a**3-3*a*I*R)*2*alpha0*a/(3*R*pi*h**1.0) \
            -2*E_eff*(a**3-3*a*I*R)/(3*R)*((2*alpha0*a)/(h*pi))**2.0 \
            +8*E_eff* ( 2/15*a**6*(15*alpha0**3+7*pi**2*beta0) - 2*a**4*I*(3*alpha0**3+pi**2*beta0)*R )/(3*pi**3*R)/h**3.0 \
            -32*a**5*E_eff*alpha0*(a**2*(5*alpha0**3+4*pi**2*beta0)-5*I*(3*alpha0**3+2*pi**2*beta0)*R)/(15*pi**4*R)/h**4.0

      

    else:
        force = 0

    return force

def bepc(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    if (z + zc) < a0:
        
        h = sahe
        I = -z - zc + a0
        a = Rt         
        
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        force = 2*a*I*E_eff * (1  \
                                - 2*alpha0*a/(pi*h) \
                                + (2*alpha0/pi)**2*(a/h)**2 \
                                -8*a**3*(3*alpha0**3+pi**2*beta0)/(3*pi**3)*(a/h)**3 \
                                -16*alpha0*(3*alpha0**3+2*pi**2*beta0)/(3*pi**4) * (a/h)**4 \
                                 )
                
    else:
        force = 0

    return force

def nanowire(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    I = -z - zc + a0
    a = Rt       
    I_critical = Rt
    
    def F_sphere(I,R):
        a = sqrt(R*I)      
        return 2*E_eff*(3*a*I*R-a**3)/(3*R)             
    
    def F_punch(I,a):
        return 2*a*I*E_eff  
    
    # All the sphere in contact, the tip is now a punch 
    if I>I_critical:                         
        force = F_sphere(I_critical,Rt ) + F_punch(I-I_critical,Rt)   
    elif I<I_critical and I>0:
        # Only the spherical area is in contact,the tip behaves like a sphere            
        force = F_sphere(I,Rt)
    else:
        force = 0
    
        
    return force
   
        
   


def bec_nanowire(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    if (z + zc) < a0:
        
        h = sahe
        I = -z - zc + a0
        a = Rt         
        
        if bonded:
            # bonded case
            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
        else:
            # Unbonded case
            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
            
        # Critical indentation between sphere and punch
        I_critical = 3*(-h*pi*(4*Rt*alpha0 - 3*h*pi) - sqrt(3)*sqrt(-h**3*pi**3*(8*Rt*alpha0 - 3*h*pi)))/(8*Rt*alpha0**2)     
        
        def F_sphere(I,R,alpha0,beta0):
            a = sqrt(R*I)* (1-2*alpha0*sqrt(R*I)/(3*h*pi))
     
            force = 2*E_eff*(3*a*I*R-a**3)/(3*R) \
                + 2*E_eff*(a**3-3*a*I*R)*2*alpha0*a/(3*R*pi*h**1.0) \
                -2*E_eff*(a**3-3*a*I*R)/(3*R)*((2*alpha0*a)/(h*pi))**2.0 \
                +8*E_eff* ( 2/15*a**6*(15*alpha0**3+7*pi**2*beta0) - 2*a**4*I*(3*alpha0**3+pi**2*beta0)*R )/(3*pi**3*R)/h**3.0 \
                -32*a**5*E_eff*alpha0*(a**2*(5*alpha0**3+4*pi**2*beta0)-5*I*(3*alpha0**3+2*pi**2*beta0)*Rt)/(15*pi**4*R)/h**4.0
            return force
        
        def F_punch(I,a,alpha0,beta0):
            force = 2*a*I*E_eff * (1  \
                                    - 2*alpha0*a/(pi*h) \
                                    + (2*alpha0/pi)**2*(a/h)**2 \
                                    -8*a**3*(3*alpha0**3+pi**2*beta0)/(3*pi**3)*(a/h)**3 \
                                    -16*alpha0*(3*alpha0**3+2*pi**2*beta0)/(3*pi**4) * (a/h)**4 \
                                     )
            return force
        
        if I<I_critical:
            # Only the a spherical area is in contact,the tip behaves like a sphere
            a = np.sqrt(Rt*I) - 2*alpha0*Rt*I/(3*h*pi)
            force = F_sphere(I,Rt,alpha0,beta0)           
            
        else:
            # All the sphere in contact, the tip is now a punch
            a = Rt
            force = F_sphere(I_critical,Rt,alpha0,beta0) + F_punch(I-I_critical,Rt,alpha0,beta0)
            
                
        
    else:
        force = 0

    return force





# ================================= Viscosity forces =====================

def visco(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    if (z + zc) < a0:
        I = -z - zc + a0        
        
        #force = 4.0 / 3.0 * E_eff * math.sqrt(Rt) * math.pow(-(z + zc - a0), 3.0 / 2.0) - eta * sqrt(Rt * (-z - zc + a0)) * v
        force = 4.0/3.0 * np.sqrt(Rt*I)/(1-nu**2)*(E_sample*I - 3.0/2.0*eta*v)
    else:
        force = 0
    return force

def visco_punch(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    if (z + zc) < a0:
        I = -z - zc + a0
        
        #force = 2 * Rt * E_eff * I - eta * Rt * v
        force = 2 * Rt/(1-nu**2)*(E_sample*I - eta * v)
    else:
        force = 0

    
    return force

def visco_cone(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    if (z + zc) < a0:
        I = (-z - zc + a0)
        force = 2*np.tan(angle)*I/(np.pi*(1-nu**2)) * ( E_sample*I - 2*eta*v)
    else:
        force = 0

    return force

def visco_nanowire(z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    I = -z - zc + a0    
    I_critical = Rt
    
    def F_sphere(I,R,v):
                     
        return 4.0/3.0 * np.sqrt(Rt*I)/(1-nu**2)*(E_sample*I - 3.0/2.0*eta*v)
    def F_punch(I,Rt,v):
        return  2 * Rt/(1-nu**2)*(E_sample*I - eta * v)
    
    # All the sphere in contact, the tip is now a punch 
    if I>I_critical:                            
        force = F_sphere(I_critical,Rt,v=0) + F_punch(I-I_critical,Rt,v)   
    elif I<I_critical and I>0:
        #a = sqrt(Rt*I)
        # Only the spherical area is in contact,the tip behaves like a sphere            
        force = F_sphere(I,Rt,v)
    else: # I<0
        force = 0
    
        
    return force


def Lennard_Jones_Potential(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    
    d = z + zc
    
    if d<=0:
        return np.nan
    else:
        return 4*epsilon_lj * ( 12 * sigma_lj**12/d**13 - 6 * sigma_lj**6/d**7)


# -----------------------------------jkr-------------------------------------------------


forcepri = np.arange(-1, 200, 0.1)
identipri = 3 * np.power(forcepri + 2 + 2 * np.sqrt(forcepri + 1), 2.0 / 3.0) - \
    4 * np.power(forcepri + 2 + 2 * np.sqrt(forcepri + 1), 1 / 6.0)


def jkr(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    ident_a = np.power(Rt * H**2 * 9 / (3 * E**2 * 64 * 36 * a0**4), 1 / 3.0)
    force_a = 9 * Rt * H / (8 * 6 * a0**2)
    if (z + zc) < a0:
        force = np.interp(
            (-z - zc + ident_a) / ident_a,
            identipri,
            forcepri) * force_a
    elif (z + zc) < ident_a and v > 0:
        force = np.interp(
            (-z - zc + ident_a) / ident_a,
            identipri,
            forcepri) * force_a
    else:
        force = -H * Rt / (6 * (z + zc) * (z + zc))

    return force


# -----------------------------------magnetic-------------------------------------------------


def magnetic_dipole(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    mu_0 = 4*np.pi*1E-7
        
    if (z + zc) < a0:
        force = 0
    else:
        # 1 emu → 10−3 A·m2 = 10−3 J/T
        force = 3/(2*np.pi)* mu_0 * (magnetic_m_tip*1E-14*1E-3) * (magnetic_m_sample*1E-14*1E-3)/(z + zc)**4
        
    return force

def magnetic_exponential(z,v,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,sigmat,
        landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):
    mu_0 = 4*np.pi*1E-7
    force = (magnetic_m_tip*1E-14*1E-3) *(magnetic_k_sample*1E9) * (magnetic_B0_sample*1E-3) * np.exp(-(z+zc)*(magnetic_k_sample*1E9))
    
    return force

# == Bryant L. Doss, Kiarash Rahmani Eliato, Keng-hui Lind  and  Robert Ros ===

def DELR(z,zc,Rt,a0,E_eff,E_sample,E2,sahe,nu):
    # https://doi.org/10.1039/C8SM02121J
    
    
    h = sahe
    def A(x):                
        return np.min([1,0.72-0.34*x+0.51*x**2])
    
    def B(x):        
        return 0.85*x+3.36*x**2
    if (z + zc) < a0:        
        # This method assumes nu_1 = nu_2 = 0.5
        E1 = E_eff*(1-nu**2)        
        I = (-z - zc + a0)
        F0 = E1/(1-nu**2) *4.0/3.0*np.sqrt(Rt*I**3)
        a0 = np.sqrt(Rt*I)
        x = a0/h
        force = F0*( B(x) + 1 ) /( B(x)*(E1/E2)**A(x)+1)
    else:
        force = 0
    
    return force

# ================================= Driven forces ========================


def driven(F0, f0, t):
    if f0.size == 1:
        force = F0 * math.cos(f0 * 2 * math.pi * t)
    else:
        force = F0[1] * math.cos(f0[1] * 2 * math.pi * t) + \
            F0[2] * math.cos(f0[2] * 2 * math.pi * t)

    return force







def sumatory(option_mask,z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,H2,angle,sahe,bonded,nu,m,D,F0,
        zc,eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4):

    
    if H2 > 0:
        H1 = H
        H = (v > 0) * H2 + (v < 0) * H1  
        
    args = (z,v,t,f0,k,Q,Rt,a0,E_eff,E_sample,E2,H,angle,sahe,bonded,nu,m,D,F0,zc,eta,
            sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, 
            magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,
            epsilon_lj,sigma_lj,
            V1,V2,V3,V4)
    
    force = driven(F0, f0, t) 
    
    if option_mask[0]:
        force += vdw(*args)
    if option_mask[1]:
        force += dlvo(*args)
    if option_mask[2]:
        force += magnetic_dipole(*args)
    if option_mask[3]:
        force += magnetic_exponential(*args)
    if option_mask[4]:
        force += hertz(*args)
    if option_mask[5]:
        force += sneddon(*args)
    if option_mask[6]:
        force += punch(*args)
    if option_mask[7]:
        force += dmt(*args)
    if option_mask[8]:
        force += jkr(*args)
    if option_mask[9]:
        force += tatara(*args)
    if option_mask[10]:
        force += behc(*args)
    if option_mask[11]:
        force += becc(*args)
    if option_mask[12]:
        force += bepc(*args)
    if option_mask[13]:
        force += visco(*args)
    if option_mask[14]:
        force += visco_cone(*args)
    if option_mask[15]:
        force += visco_punch(*args)
    if option_mask[16]:
        force += kv_sphere_BEC(*args)
    if option_mask[17]:
        force += kv_cone_BEC(*args)
    if option_mask[18]:
        force += kv_punch_BEC(*args)
    if option_mask[19]:
        force += Lennard_Jones_Potential(*args)
    if option_mask[20]:
        force += custom1(*args)
    if option_mask[21]:
        force += custom2(*args)
    if option_mask[22]:
        force += bec_nanowire(*args)
    if option_mask[23]:
        force += nanowire(*args)    
    if option_mask[24]:
        force += DELR(z,zc,Rt,a0,E_eff,E_sample,E2,sahe,nu)  
    if option_mask[25]:
        force +=visco_nanowire(*args)
    if option_mask[26]:
        force +=kv_bec_nanowire(*args)
    if option_mask[27]:
        force +=vdw_cone(*args)
    if option_mask[28]:
        force +=vdw_punch(*args)
        

    

    return force


    
    
