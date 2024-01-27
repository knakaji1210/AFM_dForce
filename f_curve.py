import forces
import numpy as np
from pylab import *
import warnings

def curve(option_mask, A0, v, t, f0, k, Q, Rt, a0, E,E_sample, E2, H, H2, angle, sahe, 
          bonded, nu, m, D, F0, zc, eta, sigmas, sigmat, landadeb, epsilonr, 
          alfa, nc, magnetic_m_tip, magnetic_m_sample, magnetic_B0_sample, 
          magnetic_k_sample,epsilon_lj,sigma_lj, V1, V2, V3, V4):
    
    # TODO revise for dissipative forces
    
    #0 self.dlvo_check_var.get(),
    #1 self.vdw_check_var.get(),
    #2 self.hertz_check_var.get(),
    #3 self.dmt_check_var.get(),
    #4 self.tatara_check_var.get(),
    #5 self.viscosity_check_var.get(),
    option_mask[5] = False  # It works just for non-disipative forces 
    #6 self.becc_check_var.get(),
    #7 self.behc_check_var.get(),
    #8 self.sls_check_var.get(),
    option_mask[8] = False  # It works just for non-disipative forces
    #9 self.jkr_check_var.get(),           
    option_mask[9] = False  # It works just for non-disipative forces 
    #10 self.dipole_check_var.get(),
    #11 self.exp_check_var.get(),
    #12 self.sneddon_check_var.get(),            
    #13 self.kv_sphere_check_var.get(),
    option_mask[13] = False  # It works just for non-disipative forces
    #14 self.kv_cone_check_var.get(),            
    option_mask[14] = False  # It works just for non-disipative forces
    #15 self.custom1_check_var.get(),
    #16 self.custom2_check_var.get()
    
    
    
    
    H2 = 0  # It works just for non-disipative forces
    f0 = 1.0
    f0 = np.array(f0)

    def ies():
        f_ave = 0
        f_vir = 0

        for t,vi in zip(time,v):
            z = A * np.cos(2 * np.pi * t)
            f_ts = forces.sumatory(option_mask,z,vi,t,f0,k,Q,Rt,a0,E,E_sample,E2,H,H2,
                        angle,sahe,bonded,nu,m,D,0,zc,eta,sigmas,sigmat,
                        landadeb,epsilonr,alfa,nc,magnetic_m_tip, 
                        magnetic_m_sample,magnetic_B0_sample,magnetic_k_sample,
                        epsilon_lj,sigma_lj,
                        V1, V2, V3, V4)
            f_ave = f_ave + f_ts
            f_vir = f_vir + f_ts * (z)
        f_ave = f_ave / time.size
        f_vir = f_vir / time.size
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            i = 2 / (k * A**2) * (f_ave**2 / k - f_vir)
        return i

    delta_t = 0.001
    time = np.arange(0, 1 + delta_t, delta_t)
    delta_A = 0.01 * A0
    amplitude = np.arange(0.1 * A0, 1.2 * A0, delta_A)
    w1 = []
    w2 = []
    for A in amplitude:
        i = ies()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            omega1 = f0 * np.sqrt(i + 1 - 1 / (2 * Q**2) *
                                  (1 + np.sqrt(1 + 4 * Q**2 * (A0**2 / A**2 - 1 - i))))
            omega2 = f0 * np.sqrt(i + 1 - 1 / (2 * Q**2) *
                                  (1 - np.sqrt(1 + 4 * Q**2 * (A0**2 / A**2 - 1 - i))))
        w1.append(omega1)
        w2.append(omega2)
    amplitude = amplitude * 1e9  # In nanometers    
    return w1, w2, amplitude
