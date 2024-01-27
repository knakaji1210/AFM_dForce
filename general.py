#!/usr/bin/python
# =======================dFORCE==============================
# The code written below is part of the dFORCE project and has
# been developed at the Forcetool laboratory in Madrid. This
# software will enable students and researchers from the
# AFM sphere to easily get an insight into the AFM theory
# through numerical simulations.
# ======================CREDITS==============================
#
# Prof. Ricardo Garcia
# PhD Victor G Gisbert


#
# We thank you in advance for sending your feedback and/or
# suggestions to:
#             dforce@forceforfuture.es


# ======================dForce LICENSE==============================

# Copyright (C) 2023  # Victor G. Gisbert, Ricardo Garcia

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
# @@@@@@@@@@@@@@@@@@@@@2014@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# What we observe is not nature itself, but nature exposed to
# our method of questioning.
#                                  W. Heisenberg
# =========================Updates=================================
#
# 19-06-2018 Victor Garcia Gisbert
#  -Python update for the 2.7 no longers support indexation with non-integers. A values like nper_sp, npp_sp, and
#       nperfin_sp are set to intergers and loops that used np.linspace(0,ns-1,ns) are substituted by range(ns)
#  -Error when using Hysteresis solved, the program no longer gives an error when this force is used
#
# 06-03-2023 Victor Garcia Gisbert
# - Migrated to python 3
# - Switched gui to tkiner
# - Added support for magnetic forces
# - Fixed forces that use bottom effect correction
# - Fixed adhesion
# - Added the punch geometry, the cone geomtry and the nanowire with BEC and KV
# - Added Lennard Jhones
# - Added option to autocalculate simulation parameteres.
# ===========================================================
#
# ========================Code lines========================
# IMPORT MODULES
import os
#import shutil
#import gtk
#import gtk.glade
#import pygtk  # temp
from scipy.integrate import odeint
from pylab import plot,xlabel,ylabel,figure
import numpy as np
import math
import forces
import analysis
import f_curve
#import time
#import tkinter as tk

from scipy import integrate
import importlib


ROOT = os.path.realpath(os.path.dirname(__file__))


# ========FUNCTIONS=================
# |||||||||||||||||||||||||||||||||||||||||||||||


# PUNTUAL MASS
def deriv(y,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,H2,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,
          sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4,option_mask):
    uprime = y[1]
    wprime = - k * y[0] / m - D * y[1] / m + forces.sumatory(
        option_mask,y[0],y[1],t,f0,k,Q,Rt,a0,E,E_sample,E2,H,H2,angle,sahe,bonded,nu,m,D,F0,zc,
        eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4) / m
    yprime = np.array([uprime, wprime])

    return yprime


def runner(b, date, iteration, self):
    
    
    
    importlib.reload(forces)
    """
    forces.custom_reset(self.toplevel1.custom1_sp_var,self.toplevel1.custom2_sp_var,
        self.toplevel1.custom3_sp_var,self.toplevel1.custom4_sp_var)
    """
    print("Single mode")

    # Transform the readed data in the rigth units

    # Frequency of maximum response in amplitude tking in account the Q factor
    # (the thermal tunning gives us that). The driving frequency will be this.
    f0 = b[0] * 1000
    k = b[2]
    Q = b[3]
    
    if self.fcorrection_q1_var.get():
        # Frequency of resonance if there weren't viscosity (in the case Q -->
        # infinity)
        f_ori = f0 / math.sqrt(1 - 1 / (2 * Q**2))
    else:
        f_ori = f0
        
    A01 = b[4]
    Rt = b[5]
    a0 = b[6]
    E = b[7]
    E2 = self.toplevel1.saym2_sp_var * 1E6
    E_sample = self.toplevel1.saym_sp_var * 1E6
    H = b[8]
    H2 = self.toplevel1.ham2_sp_var * 1E-20
    angle = self.toplevel1.angle_sp_var * np.pi / 180.0 # Angle to radians
    sahe = self.toplevel1.sahe_sp_var * 10**-9 # Height of the sample for BEC
    bonded = self.toplevel1.bec_bonded_var # Sample bonded for BEC
    nu = self.toplevel1.musam_sp_var
    eta = b[9]
    sigmas = b[10]
    sigmat = b[11]
    landadeb = b[12]
    epsilonr = b[13]
    alfa = b[14]
    nc = b[15]
    m = k / (f_ori * f_ori * 4 * np.pi * np.pi)
    D = math.sqrt(m * k) / Q

    if self.copy_fre1_var.get() == 0:
        fd = float(self.fd1_sp.get()) * 1000
        fd = np.array(fd)
    else:
        fd = f0
    
    F0 = m * A01 * math.sqrt((f_ori**2 - fd**2)**2 +(fd * f_ori / Q)**2) * 4 * np.pi * np.pi

    magnetic_m_tip = self.toplevel1.magnetic_m_tip
    magnetic_m_sample = self.toplevel1.magnetic_m_sample
    magnetic_B0_sample = self.toplevel1.magnetic_B0_sample
    magnetic_k_sample = self.toplevel1.magnetic_k_sample
    
    
    sigma_lj = self.toplevel1.sigma_lj_var * 1E-9
    epsilon_lj = self.toplevel1.epsilon_lj_var * 1E-9* 1E-9
    

    V1 = self.toplevel1.custom1_sp_var
    V2 = self.toplevel1.custom2_sp_var
    V3 = self.toplevel1.custom3_sp_var
    V4 = self.toplevel1.custom4_sp_var
    


    if self.copy_fre1_var.get():  
        # The value of fo is, from here, the value of the driven frequency, 
        # not the resonant frequency
        f0 = fd

    

    tolerance = b[23]

    start = 0
    end = b[16] / f0
    numsteps = int(b[16] * b[17])
    numsteps_ss = int(b[17] * b[18])
    n_cut = int(numsteps - numsteps_ss)
    t = np.linspace(start, end, numsteps)

    y01 = np.array([0.0000, A01 * 2 * np.pi * f0])  # Initial values

    # ns = math.floor((b[19] - b[20])/b[21])   #number of steps for zc

    #zc = np.linspace(b[20], b[19], ns)
    zc_max = b[19]
    zc_min = b[20]
    zc_step = b[21]
    zc = np.arange(zc_max,zc_min, -zc_step)
    
    
    if self.simulate_retrace_var.get():
        print("simulating retrace")        
        zc = np.concatenate([zc,zc[::-1]])
        
    ns = zc.size

    z_output = np.zeros([ns, 10])
    z_output[:, 0] = zc * 10**9
    t_output = np.zeros([int(numsteps_ss), 8])
    t_output[:, 0] = t[n_cut:].copy() * f0

    # Which step in zc is closer to the choosen zc to do ana anlysis
    closer = np.argmin(np.fabs(b[22] * np.ones(len(zc)) - zc))
    
  


    self.abort = 0

    #  Here it starts the loop for the differents values of zc

    for counter in range(ns):

       
        
        self.progressbar1['value'] = int(100*float(counter+1) / float(ns))        
        percentage = int(100*float(counter+1) / float(ns)) 
        self.style.configure('text.Horizontal.TProgressbar', 
                        text='{:g} %'.format(percentage))  # update label
        
        self.toplevel1.update_idletasks()
        self.toplevel1.update()
        
        if self.abort == 1:
            return

       
        option_mask = self.get_option_mask()

        print(("zc step: %i" % counter))

        a = t[1] - t[0]
        #forces.sls_reset()
        args = (f0, k, Q, Rt, a0, E, E_sample, E2, H, H2, angle, sahe,bonded, nu, m, D, 
                F0,zc[counter], eta, sigmas, sigmat, landadeb, epsilonr, alfa,
                nc, magnetic_m_tip, magnetic_m_sample,
                magnetic_B0_sample,magnetic_k_sample,
                epsilon_lj,sigma_lj,
                V1, V2, V3, V4, option_mask)
        y, info = odeint(deriv, y01, t, args=args, full_output=1, 
                                   rtol=10**-tolerance, atol=10**-tolerance)
        
        stationary = y[n_cut:, 0].copy()
        #v = analysis.velocity(stationary,a)
        v = y[n_cut:, 1].copy()
        f_ts = analysis.force_ts(option_mask,stationary,v,t[n_cut:],f0,k,Q,Rt,
                                 a0,E,E_sample, E2,H,H2,angle,sahe,bonded,nu,m,D,F0,zc[counter],
                                 eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, 
                                 magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4)

        amplitude = analysis.amplitude(stationary)
        z_output[counter, 1] = amplitude * 10**9
        phase = analysis.phase_shift(
            stationary, analysis.force_dr(F0, f0, t[n_cut:]))
        z_output[counter, 2] = phase
        minimum_dis = analysis.minimum_distance(stationary, zc[counter], a0)
        z_output[counter, 3] = minimum_dis * 10**9
        f_max = analysis.max_force(f_ts)
        z_output[counter, 4] = f_max * 10**12
        dis_en = analysis.dissipated_energy(v, f_ts)
        z_output[counter, 5] = a * dis_en * 10**20 / (b[18])
        vir = analysis.virial(stationary, f_ts)
        z_output[counter, 6] = vir * 10**20 / (b[18] * b[17])
        meand = analysis.mean_deflection(stationary)
        z_output[counter, 7] = meand * 10**9
        t_con = analysis.contact_time(stationary, zc[counter], a0)
        z_output[counter, 8] = t_con
        dis_pow = analysis.dissipated_power(v, f_ts)
        z_output[counter, 9] = dis_pow * 10**16

        
        if counter == closer:
            t_output[:, 1] = y[n_cut:, 0] * 10**9
            t_output[:, 2] = v * 10**6
            t_output[:, 3] = f_ts * 10**12
            identation = analysis.identation(stationary, zc[counter], a0)
            t_output[:, 4] = identation * 10**9
            fampt, fphat, freq_bas = analysis.fourier(stationary, a)
            t_output[:, 5] = fampt
            t_output[:, 6] = fphat
            t_output[:, 7] = freq_bas

            if A01 == 0:  # If the amplitude is 0 a default curve is calculated with an Amplitude of 1 nm
                w1, w2, fre_amplitude = f_curve.curve(
                    option_mask, 1, v, t[n_cut:], f0, k, Q, Rt, a0, E,E_sample, E2, H, H2, angle, sahe,bonded, nu, m, D, F0, zc[counter], eta, sigmas, sigmat, landadeb, epsilonr, alfa, nc, magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1, V2, V3, V4)
            else:
                w1, w2, fre_amplitude = f_curve.curve(
                    option_mask, A01, v, t[n_cut:], f0, k, Q, Rt, a0, E,E_sample, E2, H, H2, angle, sahe,bonded, nu, m, D, F0, zc[counter], eta, sigmas, sigmat, landadeb, epsilonr, alfa, nc, magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1, V2, V3, V4)

            w1 = np.array(w1)
            w2 = np.array(w2)
            fre_amplitude = np.array(fre_amplitude)

            if w1.ndim > 1:
                # this is only to avoid an error, the results are still wrong
                w1 = w1[:, 0]
                w2 = w2[:, 0]

            w_output = np.zeros([len(fre_amplitude), 3])
            w_output[:, 0] = w1
            w_output[:, 1] = w2
            w_output[:, 2] = fre_amplitude

        y01 = np.array([amplitude * np.cos(-phase * np.pi / 180),
                        amplitude * 2 * np.pi * f0 * np.cos((-phase + 90) * np.pi / 180)])
        y01 = y[-1]


    
    name = "dForceproject"
    name = name + date
    name = name + "/tdom"
    name = name + iteration
    name = name + ".csv"
    header = "time zdeflection velocity force indentation famp fphat freq_bas"
    delimiter = ","
    comments=""
    header = header.replace(" ",delimiter)
    np.savetxt(name, t_output,header = header,comments = comments,delimiter = delimiter)

    
    name = "dForceproject"
    name = name + date
    name = name + "/zcdom"
    name = name + iteration
    name = name + ".csv"
    header = "zc Amplitude Phase Distance_min Force_max Energy_diss Virial1 Deflection Contact_time Power_diss"
    header = header.replace(" ",delimiter)
    np.savetxt(name, z_output,header = header,comments = comments,delimiter = delimiter)

    name = "dForceproject"
    name = name + date
    name = name + "/wdom"
    name = name + iteration
    name = name + ".csv"
    header = "w1 w2 fre_amplitude"
    header = header.replace(" ",delimiter)
    np.savetxt(name, w_output,header = header,comments = comments,delimiter = delimiter)
    
    self.style.configure('text.Horizontal.TProgressbar',text='Finished')
    self.toplevel1.update_idletasks()
    self.toplevel1.update()
    
    print("Finished simulation")

    return t_output


# ==================== Euler-Bernoulii ================
def derivEB(y,t,f0,k,Q,Rt,a0,E,E_sample, E2,H,H2,angle,sahe,bonded,nu,m,D,F0,zc,eta,sigmas,
            sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4,option_mask):
    f = [y[1],
         - k[1] * y[0] / m[1] - D[1] * y[1] / m[1] + 
         forces.sumatory(option_mask,y[0] + y[2],y[1] + y[3],t,f0,k[1],Q[1],
                         Rt,a0,E,E_sample, E2,H,H2,angle,sahe,bonded,nu,m[1],D[1],F0,zc,eta,
                         sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4) / m[1],
         y[3],
         - k[2] * y[2] / m[2] - D[2] * y[3] / m[2] + forces.sumatory(option_mask,
                         y[2] + y[0],y[3] + y[1],t,f0,k[2],Q[2],Rt,a0,E,E_sample,E2,H,
                         H2,angle,sahe,bonded,nu,m[2],D[2],F0,zc,eta,sigmas,sigmat,
                         landadeb,epsilonr,alfa,nc,magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4) / m[2]]

    return f


def runnerEB(b, date, iteration,self):

    
    importlib.reload(forces)
    print("....EB....")

    # Transform the readed data in the rigth units

    f0 = np.empty(3)
    f_ori = np.empty(3)
    k = np.empty(3)
    Q = np.empty(3)
    A0 = np.empty(3)
    m = np.empty(3)
    D = np.empty(3)
    F0 = np.empty(3)
    fd = np.empty(3)
    
    #f01,f02,kc,kc2,q,q2,a0,a02,rad,a00,ebarra,ham,eta,sigmas,sigmat,landadeb,epsilon,alfa,nc,nper,npp,nperfin,zcmax,zcmin,dzc,fixzc,tolerance = b
    """
    b
    0 f01
    1 f02
    2 kc
    3 kc2
    4 q
    5 q2
    6 a0
    7 a02
    8 rad
    9 a00
    10 ebarra
    11 ham
    12 eta
    13 sigmas
    14 sigmat
    15 landadeb
    16 epsilon
    17 alfa
    18 nc
    19 nper
    20 npp
    21 nperfin
    22 zcmax
    23 zcmin
    24 dzc
    25 fixzc
    26 tolerance
    """

    f0[1] = b[0] * 1000
    f0[2] = b[1] * 1000
    k[1] = b[2]
    k[2] = b[3]
    Q[1] = b[4]
    Q[2] = b[5]
    

        
    if self.fcorrection_q1_var.get():
        # Frequency of resonance if there weren't viscosity (in the case Q -->
        # infinity)
        f_ori[1] = f0[1] / math.sqrt(1 - 1 / (2 * Q[1]**2))
    else:
        f_ori[1] = f0[1]
    if self.fcorrection_q2_var.get():
        # Frequency of resonance if there weren't viscosity (in the case Q -->
        # infinity)
        f_ori[2] = f0[2] / math.sqrt(1 - 1 / (2 * Q[2]**2))
    else:
        f_ori[2] = f0[2]
    
    A0[1] = A01 = b[6]
    A0[2] = A02 = b[7]
    Rt = b[8]
    a0 = b[9]
    E = b[10]
    E_sample = self.toplevel1.saym_sp_var * 1E6
    E2 = self.toplevel1.saym2_sp_var * 1E6
    H = b[11]
    H2 = self.toplevel1.ham2_sp_var * 1E-20
    angle = self.toplevel1.angle_sp_var * np.pi / 180.0
    sahe = self.toplevel1.sahe_sp_var * 1E-9
    bonded = self.toplevel1.bec_bonded_var # Sample bonded for BEC
    nu = self.toplevel1.musam_sp_var # poisson coefficient of the sample
    
    sigma_lj = self.toplevel1.sigma_lj_var * 1E-9
    epsilon_lj = self.toplevel1.epsilon_lj_var * 1E-9* 1E-9
        
    eta = b[12]
    sigmas = b[13]
    sigmat = b[14]
    landadeb = b[15]
    epsilonr = b[16]
    alfa = b[17]
    nc = b[18]
    
    
    m[1] = k[1] / (f_ori[1] * f_ori[1] * 4 * np.pi * np.pi)
    m[2] = k[2] / (f_ori[2] * f_ori[2] * 4 * np.pi * np.pi)
    D[1] = math.sqrt(m[1] * k[1]) / Q[1]
    D[2] = math.sqrt(m[2] * k[2]) / Q[2]

    if self.copy_fre1_var.get() == 0:
        fd[1] = float(self.fd1_sp.get()) * 1000
    else:
        fd[1] = f0[1]
        
    if self.copy_fre2_var.get() == 0:
        fd[2] = float(self.fd2_sp.get()) * 1000
    else:
        fd[2] = f0[2]

    F0[1] = m[1] * A0[1] * math.sqrt((f_ori[1]**2 - fd[1]**2) **
                                     2 + (fd[1] * f_ori[1] / Q[1])**2) * 4 * np.pi * np.pi
    F0[2] = m[2] * A0[2] * math.sqrt((f_ori[2]**2 - fd[2]**2) **
                                     2 + (fd[2] * f_ori[2] / Q[2])**2) * 4 * np.pi * np.pi
    
    
    magnetic_m_tip = self.toplevel1.magnetic_m_tip
    magnetic_m_sample = self.toplevel1.magnetic_m_sample
    magnetic_B0_sample = self.toplevel1.magnetic_B0_sample
    magnetic_k_sample = self.toplevel1.magnetic_k_sample
    
    V1 = self.toplevel1.custom1_sp_var
    V2 = self.toplevel1.custom2_sp_var
    V3 = self.toplevel1.custom3_sp_var
    V4 = self.toplevel1.custom4_sp_var
    #nu = 4  # dummy variable

    tolerance = b[26]

    if not self.copy_fre1_var.get() :  # The value of fo is, from here, the value of the driven frequency, not the resonant frequency
        f0[1] = fd[1]
    if not self.copy_fre2_var.get() :  # The value of fo is, from here, the value of the driven frequency, not the resonant frequency
        f0[2] = fd[2]

    nc = b[18] # what is this?
    nper = b[19]
    npp = b[20]
    nperfin = b[21]
    
    assert nper == float( self.nper_sp.get()) 
    assert npp == float( self.npp_sp.get()) 
    assert nperfin == float( self.nperfin_sp.get()) 
    assert nper == float( self.nper_sp.get()) 
    
  
    start = 0
    end = nper / f0[1]
    numsteps = int(nper * npp)
    numsteps_ss = int(nperfin * npp)
    n_cut = int(numsteps - numsteps_ss)
    t = np.linspace(start, end, numsteps)
    
    
    

    y01 = np.array([0.0000, A0[1] * 2 * np.pi * f0[1], 0.0000,
                   A0[2] * 2 * np.pi * f0[2]])  # Initial values

    

    #zc = np.linspace(b[20], b[19], ns)
    zc_max = b[22]
    zc_min = b[23]
    zc_step = b[24] # is this actually what we want?
    zc = np.arange(zc_max,zc_min, -zc_step)
    if self.simulate_retrace_var.get():
        print("simulating retrace")
        zc = np.concatenate([zc,zc[::-1]])
    ns = zc.size
    

    z_output = np.zeros([int(ns), 16])
    z_output[:, 0] = zc * 10**9
    t_output = np.zeros([numsteps_ss, 16])
    t_output[:, 0] = t[n_cut:].copy() * f0[1]

    # Which step in zc is closer to the choosen zc to do ana anlysis
    closer = np.argmin(np.fabs(b[25] * np.ones(len(zc)) - zc))

    self.abort = 0
    
    option_mask = self.get_option_mask()


    #  Here it starts the loop for the differents values of zc

    for counter in range(ns):
        if self.abort == 1:
            return
        
        percentage = int(100*float(counter+1) / float(ns))
        
        self.progressbar1['value'] = percentage       
        percentage = int(100*float(counter+1) / float(ns)) 
        self.style.configure('text.Horizontal.TProgressbar', 
                        text='{:g} %'.format(percentage))  # update label        
        self.toplevel1.update_idletasks()
        self.toplevel1.update()
  

        print(("zc step: %i" % counter))
    
    
        a = t[1] - t[0]
        #forces.sls_reset()
        args = (f0,k,Q,Rt,a0,E,E_sample, E2,H,H2,angle,sahe,bonded,nu,m,D,F0,zc[counter],
            eta,sigmas,sigmat,landadeb,epsilonr,alfa,nc,magnetic_m_tip, 
            magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1,V2,V3,V4,
            option_mask)
    
        y = odeint(derivEB,y01,t,args=args,rtol=10**-tolerance,
            atol=10**-tolerance)

        stationary1 = y[n_cut:, 0].copy()
        stationary2 = y[n_cut:, 2].copy()
        stationary_sum = stationary1 + stationary2

        v1 = analysis.velocity(stationary1, a)
        v2 = analysis.velocity(stationary2, a)
        v_sum = analysis.velocity(stationary_sum, a)

        f_ts = analysis.force_ts(option_mask,
                                 stationary_sum,
                                 v_sum,
                                 t[n_cut:],
                                 f0,
                                 k[1],
                                 Q[1],
                                 Rt,
                                 a0,
                                 E,
                                 E_sample,
                                 E2,
                                 H,
                                 H2,
                                 angle,
                                 sahe,bonded,
                                 nu,
                                 m[1],
                                 D[1],
                                 F0,
                                 zc[counter],
                                 eta,
                                 sigmas,
                                 sigmat,
                                 landadeb,
                                 epsilonr,
                                 alfa,
                                 nc,
                                 magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,
                                 epsilon_lj,sigma_lj,
                                 V1,V2,V3,V4)

        amplitude1 = analysis.amplitude(stationary1)
        z_output[counter, 1] = amplitude1 * 10**9
        phase1 = analysis.phase_shift(
            stationary1, analysis.force_dr(F0[1], f0[1], t[n_cut:]))
        z_output[counter, 2] = phase1
        amplitude2 = analysis.amplitude(stationary2)
        z_output[counter, 3] = amplitude2 * 10**9
        phase2 = analysis.phase_shift(
            stationary2, analysis.force_dr(F0[2], f0[2], t[n_cut:]))
        z_output[counter, 4] = phase2
        minimum_dis = analysis.minimum_distance(
            stationary_sum, zc[counter], a0)
        z_output[counter, 5] = minimum_dis * 10**9
        maximum_dis = analysis.maximum_distance(stationary_sum, zc[counter])
        z_output[counter, 6] = maximum_dis * 10**9
        f_max = analysis.max_force(f_ts)
        z_output[counter, 7] = f_max * 10**12
        edisp1 = analysis.dissipated_energy(v1, f_ts)
        z_output[counter, 8] = a * b[20] * \
            edisp1 * 10**20 / (b[21] * f0[1] / f0[2])
        edisp2 = analysis.dissipated_energy(v2, f_ts)
        z_output[counter, 9] = a * b[20] * edisp2 * 10**20 / (b[21] * b[20])
        viri1 = analysis.virial(stationary1, f_ts)
        z_output[counter, 10] = viri1 * 10**20 / (b[21] * f0[1] / f0[2])
        viri2 = analysis.virial(stationary2, f_ts)
        z_output[counter, 11] = viri2 * 10**20 / (b[21] * b[20])
        meand = analysis.mean_deflection(stationary_sum)
        z_output[counter, 12] = meand * 10**9
        t_con = analysis.contact_time(stationary_sum, zc[counter], a0)
        z_output[counter, 13] = t_con
        pow1 = analysis.dissipated_power(v1, f_ts)
        z_output[counter, 14] = pow1 * 10**16
        pow2 = analysis.dissipated_power(v2, f_ts)
        z_output[counter, 15] = pow2 * 10**16

        if counter == closer:

            t_output[:, 1] = stationary_sum.copy() * 10**9
            t_output[:, 2] = v_sum.copy() * 10**6
            t_output[:, 3] = stationary1.copy() * 10**9
            t_output[:, 4] = v1.copy() * 10**6
            t_output[:, 5] = stationary2.copy() * 10**9
            t_output[:, 6] = v2.copy() * 10**6
            t_output[:, 7] = f_ts.copy() * 10**12
            identation = analysis.identation(stationary_sum, zc[counter], a0)
            t_output[:, 8] = identation.copy() * 10**9
            famp1, fpha1, freq_bas = analysis.fourier(stationary1, a)
            t_output[:, 9] = famp1.copy()
            t_output[:, 10] = fpha1.copy()
            famp2, fpha2, freq_bas = analysis.fourier(stationary2, a)
            t_output[:, 11] = famp2.copy()
            t_output[:, 12] = fpha2.copy()
            fampt, fphat, freq_bas = analysis.fourier(stationary_sum, a)
            t_output[:, 13] = fampt.copy()
            t_output[:, 14] = fphat.copy()
            t_output[:, 15] = freq_bas.copy()
            
            if A01 == 0:  # If the amplitude is 0 a default curve is calculated with an Amplitude of 1 nm
            
                w1, w2, fre_amplitude = f_curve.curve(
                    option_mask, 1, v1, t[n_cut:], f0, k, Q, Rt, a0, E,E_sample, E2, H, H2, angle, sahe,bonded, nu, m, D, F0, zc[counter], eta, sigmas, sigmat, landadeb, epsilonr, alfa, nc, magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1, V2, V3, V4)
            else:
                w1, w2, fre_amplitude = f_curve.curve(
                    option_mask, A01, v1, t[n_cut:], f0, k, Q, Rt, a0, E,E_sample, E2, H, H2, angle, sahe,bonded, nu, m, D, F0, zc[counter], eta, sigmas, sigmat, landadeb, epsilonr, alfa, nc, magnetic_m_tip, magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,epsilon_lj,sigma_lj,V1, V2, V3, V4)

            w1 = np.array(w1)
            w2 = np.array(w2)
            fre_amplitude = np.array(fre_amplitude)

            if w1.ndim > 1:
                # this is only to avoid an error, the results are still wrong
                w1 = w1[:, 0]
                w2 = w2[:, 0]

            w_output = np.zeros([len(fre_amplitude), 3])
            w_output[:, 0] = w1
            w_output[:, 1] = w2
            w_output[:, 2] = fre_amplitude
            

        y01 = np.array([amplitude1 * np.cos(-phase1 * np.pi / 180),
                        amplitude1 * 2 * np.pi * f0[1] * np.cos((-phase1 + 90) * np.pi / 180),
                        amplitude2 * np.cos(-phase2 * np.pi / 180),
                        amplitude2 * 2 * np.pi * f0[2] * np.cos((-phase2 + 90) * np.pi / 180)])
        y01 = y[-1]

    name = "dForceproject"
    name = name + date
    name = name + "/tdom"
    name = name + iteration
    name = name + ".csv"
    header = "time zdeflection velocity zdeflection1 velocity1 zdeflection2 velocity2 force indentation famp1 fpha1 famp2 pha2 fampt fphat freq_bas"
    delimiter = ","
    comments=""
    header = header.replace(" ",delimiter)
    np.savetxt(name, t_output,header = header,comments = comments,delimiter = delimiter)
    

    name = "dForceproject"
    name = name + date
    name = name + "/zcdom"
    name = name + iteration
    name = name + ".csv"    
    header = "zc Amplitude1 Phase1 Amplitude2 Phase2 Distance_min1 Distance_min2 Force_max Energy_diss1 Eenergy_diss2 Virial1 Virial2 Deflection Contact_time Power_diss1 Power_diss2"
    header = header.replace(" ",delimiter)
    np.savetxt(name, z_output,header = header,comments = comments,delimiter = delimiter)

    w_output = np.zeros([len(fre_amplitude), 3])
    w_output[:, 0] = w1
    w_output[:, 1] = w2
    w_output[:, 2] = fre_amplitude
    
    name = "dForceproject"
    name = name + date
    name = name + "/wdom"
    name = name + iteration
    name = name + ".csv"
    header = "w1 w2 fre_amplitude"
    header = header.replace(" ",delimiter)
    np.savetxt(name, w_output,header = header,comments = comments,delimiter = delimiter)
    
 
    self.style.configure('text.Horizontal.TProgressbar', text='Finished')  
    
    print("Finished simulation")
  

    return t_output



