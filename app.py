print("Loading dForce 2.0")
#!/usr/bin/python3
import tkinter as tk
import tkinter.ttk as ttk
import os
import numpy as np

import time
import matplotlib
from matplotlib.pyplot import plot,figure,xlabel,ylabel,savefig,grid,title,show
from general import runner, runnerEB
from guidforceapp import GuiDforceApp as MainApp
from force_parameters_menu import ForceParametersMenuApp
from credits_app import credits_app
from examples_app import examples_app

import matplotlib.pyplot as plt

from tkinter import NORMAL,DISABLED
from custom_widgets import CheckbuttonWithFormula
   
class dforce_APP(MainApp):
    
    def __init__(self,master=None):
        super().__init__(master)
        
        
        
        self.toplevel1.wm_title("dForce 2.0")
        
        self.img = tk.PhotoImage(file="icons/icon.png")
        self.toplevel1.wm_iconphoto(False, self.img)
                
        #============ Some style changes====================================
        style = ttk.Style()
        # This will change the foreground color of ALL tabs to red.
        style.configure('TNotebook.Tab', foreground='green',font='-weight bold -size 9')
        
        # TODO fix this?
        """
        style.configure('TLabelFrame', background='SystemWindow')
        style.configure('Labelframe', background='SystemWindow')
        style.configure('TFrame', background='SystemWindow')
        style.configure('Frame', background='SystemWindow')
        """
        # Tabs of the forces        
        
        self.label33.config(font = '-weight bold -size 10')
        self.label48.config(font = '-weight bold -size 10')
        self.label73.config(font = '-weight bold -size 10')
        self.label23.config(font = '-weight bold -size 10')
        self.label74.config(font = '-weight bold -size 10')
        
        self.label35.config(font = '-weight bold -size 9')
        self.label49.config(font = '-weight bold -size 9')
        self.label40.config(font = '-weight bold -size 9')
        # ============ Checkbox for zc need some fancy configuration
       
        self.ampst_cb.grid_forget()
        self.ampst_cb = CheckbuttonWithFormula(self.frame5,"Amplitude",self.ampst_cb_var)
        self.ampst_cb.grid(column=0, row=0, sticky="ew")

        self.fasest_cb.grid_forget()
        self.fasest_cb = CheckbuttonWithFormula(self.frame5,"Phase shift",self.fasest_cb_var)
        self.fasest_cb.grid(column=0, row=1, sticky="ew")

        self.etst_cb.grid_forget()
        self.etst_cb = CheckbuttonWithFormula(self.frame5,"Dissipated energy",self.etst_cb_var)
        self.etst_cb.grid(column=0, row=2, sticky="ew")

        self.ptst_cb.grid_forget()
        self.ptst_cb = CheckbuttonWithFormula(self.frame5,"Dissipated power",self.ptst_cb_var)
        self.ptst_cb.grid(column=0, row=3, sticky="ew")
        
        self.virialt_cb.grid_forget()
        self.virialt_cb = CheckbuttonWithFormula(self.frame5,"Virial",self.virialt_cb_var)
        self.virialt_cb.grid(column=0, row=4, sticky="ew")

        self.mdef_cb.grid_forget()
        self.mdef_cb = CheckbuttonWithFormula(self.frame5,"Mean deflection",self.mdef_cb_var)
        self.mdef_cb.grid(column=0, row=5, sticky="ew")

        self.mdis_cb.grid_forget()
        self.mdis_cb = CheckbuttonWithFormula(self.frame5,"Minimum distance",self.mdis_cb_var)
        self.mdis_cb.grid(column=0, row=6, sticky="ew")
        
        self.maforce_cb.grid_forget()
        self.maforce_cb = CheckbuttonWithFormula(self.frame5,"Maximum force",self.maforce_cb_var)
        self.maforce_cb.grid(column=0, row=7, sticky="ew")
        
        self.ctime_cb.grid_forget()
        self.ctime_cb = CheckbuttonWithFormula(self.frame5,"Contact time",self.ctime_cb_var)
        self.ctime_cb.grid(column=0, row=8, sticky="ew")
        
        # 2nd mode and others
        self.amps2_cb.grid_forget()
        self.amps2_cb = CheckbuttonWithFormula(self.frame5,"2nd Amplitude",self.amps2_cb_var)
        self.amps2_cb.grid(column=1, row=0, sticky="ew")

        self.fases2_cb.grid_forget()
        self.fases2_cb = CheckbuttonWithFormula(self.frame5,"2nd Phase shift",self.fases2_cb_var)
        self.fases2_cb.grid(column=1, row=1, sticky="ew")

        self.ets2_cb.grid_forget()
        self.ets2_cb = CheckbuttonWithFormula(self.frame5,"2nd Dissipated energy",self.ets2_cb_var)
        self.ets2_cb.grid(column=1, row=2, sticky="ew")

        self.pts2_cb.grid_forget()
        self.pts2_cb = CheckbuttonWithFormula(self.frame5,"2nd Dissipated power",self.pts2_cb_var)
        self.pts2_cb.grid(column=1, row=3, sticky="ew")

        self.virial2_cb.grid_forget()
        self.virial2_cb = CheckbuttonWithFormula(self.frame5,"2nd Virial",self.virial2_cb_var)
        self.virial2_cb.grid(column=1, row=4, sticky="ew")
        
        self.indent_cb.grid_forget()
        self.indent_cb = CheckbuttonWithFormula(self.frame5,"Maximum indentation",self.indent_cb_var)
        self.indent_cb.grid(column=1, row=5, sticky="ew")
        
        self.lateralresolution_cb.grid_forget()
        
        self.lateralresolution_cb= CheckbuttonWithFormula(self.frame5,
                                                          "Lateral resolution/Contact radius",
                                                          self.lateralresolution_cb_var,width = 35)
        self.lateralresolution_cb.grid(column=1, row=6, sticky="ew")
        
        
        #============ Values of the forces====================================
        # Mechanical
    
        self.toplevel1.sr_sp_var = 0.00
        self.toplevel1.saym_sp_var = 100.00
        self.toplevel1.saym2_sp_var = 100.0
        self.toplevel1.musam_sp_var = 0.3
        # Electrical
        self.toplevel1.epsilon_sp_var = 1.00
        self.toplevel1.sigmas_sp_var = 0.00
        self.toplevel1.sigmat_sp_var = 0.00
        self.toplevel1.landadeb_sp_var = 1.00
        # Other
        self.toplevel1.ham_sp_var = 1.00
        self.toplevel1.ham2_sp_var = 0.00
        self.toplevel1.a00_sp_var = 0.169
        self.toplevel1.angle_sp_var = 10.00
        self.toplevel1.sahe_sp_var = 10.00
        self.toplevel1.bec_bonded_var = True
        # Viscosity
        self.toplevel1.mvis_sp_var = 0.00
        
        # Lennad Jhones
        self.toplevel1.sigma_lj_var = 1.00
        self.toplevel1.epsilon_lj_var = 1.00
        
        # User customized
        self.toplevel1.custom1_sp_var = 0.00
        self.toplevel1.custom2_sp_var = 0.00
        self.toplevel1.custom3_sp_var = 0.00
        self.toplevel1.custom4_sp_var = 0.00
        
        # Magnetic
        self.toplevel1.magnetic_m_tip = 2.5 # SSS-MFMR
        self.toplevel1.magnetic_m_sample = 0.01 # 
        self.toplevel1.magnetic_B0_sample = 1.00 # 1 mT
        self.toplevel1.magnetic_k_sample = 0.1 # nota kappa = 2*np.pi/lambda
        
        # Parameters used by the runner
        self.flagin = 0
        self.flagpm = 1
        self.flagfile = 0
        self.flagpnam = 0
        self.flagaux = 0
        self.compareflag = 0     # Flag for compare two different simulations
        self.date = time.strftime("20%y.%m.%d-%H.%M")
        self.abort = 0
        
        # Preloading some parameters for the simulation
        self.copy_fre1_var.set(1) # Copy the frequency 
        self.copy_fre2_var.set(1) # Copy the frequency 
        
        self.tdomload = None
        self.tdomload1 = None
        self.zdomload = None
        self.zdomload1 = None
        self.wdomload = None
        self.wdomload1 = None
        
        self.point_mass_activated()
        self.single_mode_activated()
        
        self.piat = None
        self.button_plot.config(state=DISABLED)
        
      
    def gtk_main_quit(self):                
        self.toplevel1.destroy()
        
    def on_calculate_simulation_parameters(self):
        
        f01 = float( self.f01_sp.get() ) # Natural resonant frequency
        
        Q1 = float( self.q1_sp.get() )
        
        f02 = float( self.f02_sp.get() )         
        Q2 = float( self.q2_sp.get() )    
        
        A01 = float( self.a01_sp.get() ) 
        A02 = float( self.a02_sp.get() ) 
        zmax = (A01+A02)*1.25
        zmin = (A01+A02)*0.25
        dz = (zmax-zmin)/20.0 
        z_fixed = (zmin + zmax)/2
        
        
        w01 = f01 * np.pi * 2
        w02 = f02 * np.pi * 2
        tau_1 = 2*Q1/w01
        tau_2 = 2*Q2/w02
        
        # if single mode
        if self.multifrequency_mode.get():
            tau = tau_1
        # Multifrequency mode
        else:            
            tau = tau_1 + tau_2
        # %10 error
        t_calculation = -tau *  np.log(0.1)
        # %1 error
        t_calculation = -tau *  np.log(0.01)
        n_calculation = int(t_calculation*f01) 
        n_calculation = n_calculation + 1 # one period more, just in case
        #n_steady = n_calculation//10 
        n_steady = n_calculation//20 
        n_steady = n_steady + 1 # one period more, just in case
        n_total = n_calculation + n_steady
        
        
                
        self.nper_sp # Number of periods
        self.npp_sp # Numbr of points per period        
        self.nperfin_sp # Interval of periods to calculate the stady state
        
        self.nper_sp.set(str( n_total) ) 
        self.nperfin_sp.set(str( n_steady))
        
        if zmax > zmin and zmin >0:
            self.zcmax_sp.set('{0:.2f}'.format(zmax))
            self.zcmin_sp.set('{0:.2f}'.format(zmin))
            self.deltazc_sp.set('{0:.2f}'.format(dz))
            self.fixedzc_sp.set('{0:.2f}'.format(z_fixed))
        
        
        
    def on_run_ex_clicked(self):
        
        
        """
        self.abort").set_visible(True)
        self.savingd.show()

        while gtk.events_pending():
            gtk.main_iteration(False)  # Refresh the GUI
        """
        self.button_abort.config(state=NORMAL)
        self.button_run.config(state=DISABLED)
        self.button_plot.config(state=DISABLED)
        self.button_about.config(state=DISABLED)
        self.button_quit.config(state=DISABLED)
        self.button_examples.config(state=DISABLED)
        self.toplevel1.update_idletasks()
        self.toplevel1.update()
        
        
        self.flagin += 1
        
        self.flagfile = 0   # There will not print the load files
        
        #========================== CANTILEVER    =============================
        f01 = float( self.f01_sp.get() ) # Natural resonant frequency
        fd1 = float( self.fd1_sp.get() )  # Driving frequency
        q = float( self.q1_sp.get() )
        kc = float( self.kc1_sp.get() )
        a0 = float( self.a01_sp.get() ) * 1E-9
        
        #========================== SAMPLE   ==================================
        
        emuestraa = self.toplevel1.saym_sp_var
        emuestra = emuestraa * 1E6  # in Mpa now
        epuntaa = float( self.ymt_sp.get() )   
        epunta = epuntaa * 1E9  # in Gpa now            
        rada = float( self.tr_sp.get() )   
        rad = rada * 1e-9
        rs = self.toplevel1.sr_sp_var  # Surface radius
        rs = rs * 1e-9
        a00a = self.toplevel1.a00_sp_var # a00 interatomic distance
        a00 = a00a * 1e-9
        hama = self.toplevel1.ham_sp_var # Hamaker constant
        ham = hama * 1e-20
        eta = self.toplevel1.mvis_sp_var    # Viscosity coeficient
        sigmasa = self.toplevel1.sigmas_sp_var  # Surface charge density
        sigmas = sigmasa * 10**-3
        sigmata = self.toplevel1.sigmat_sp_var  # Tip charge density
        sigmat = sigmata * 10**-3
        landadeba = self.toplevel1.landadeb_sp_var  # landa debye (nm)
        landadeb = landadeba * 10**-9
        epsilon = self.toplevel1.epsilon_sp_var  # epsilon
        
        mutip = float(self.mutip_sp.get())  # poisson tip        
        
        musample = self.toplevel1.musam_sp_var  # poisson sample
        
        #sigma_lj = self.toplevel1.sigma_lj_var
        #epsilon_lj = self.toplevel1.epsilon_lj_var
        
     
        
        nper = int(self.nper_sp.get())
        npp = int(self.npp_sp.get())
        #naux = int(nper * npp) 				# need C
        nperfin = int(self.nperfin_sp.get())
        dzca = float( self.deltazc_sp.get() )
        dzc = dzca * 1.e-9
        zcmina = float( self.zcmin_sp.get())
        zcmin = zcmina * 1.e-9
        zcmaxa = float( self.zcmax_sp.get()) 
        zcmax = zcmaxa * 1.e-9
        fixzca = float( self.fixedzc_sp.get())
        tolerance = float( self.tolerance_sp.get())
        #==========================BIMODAL CASE    ============================
        f02 = float( self.f02_sp.get() ) 
        #float( self.fd2_sp.get() )  # Driving frequency
        q2 = float( self.q2_sp.get() )    
        kc2 = float( self.kc2_sp.get() )        
        a02 = float( self.a02_sp.get() ) * 1E-9 # AMPLITUDE OF THE SECOND MODE
              
        
        ###FIRST
        if fixzca == 0.0:
            fixzc = zcmax
        else:
            fixzc = fixzca * 1.e-9
        if emuestra > 0 and epunta > 0:		# NEW
            ebarra = (1. / ((1. - musample**2) /
                      emuestra + (1. - mutip**2) / epunta))
        else:
            ebarra = 0.
        
        # -----------------tatara----------------------
        if rs > 0 and self.tatara_check_var.get():
            radeff = (1.0 / (1.0 / rs + 1.0 / rad))
            alfa = 4.0 / 3.0 * ebarra * np.sqrt(radeff)                    
            nc = 4 * np.pi * rs * rad * emuestra * epunta / \
                (6 + mutip + musample - 2 * mutip**2 - 2 * musample**2)
            rad = radeff
        else:
            alfa = 0
            nc = 0
        
        # =========================Files handling w Python ============
        # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        if self.flagpnam == 1 and self.flagin == 1:
            self.palo1 = "dForceproject" + self.date
            try:
                os.mkdir(self.palo1)
            except:
                print("Folder already exist, files will be overwriten")
                
        elif self.flagpnam == 0 and self.flagin == 1:
            self.palo1 = "dForceproject" + self.date
            try:
                os.mkdir(self.palo1)
            except:
                print("Folder already exist, files will be overwriten")
        elif self.flagin > 1:
            print("Project folder has been already created")
        
        date = self.date
        iteration = str(self.flagin)
  
        
        # if single mode
        #if self.multifrequency_mode.get():
        if self.mass_model.get() == "point_mass":
            
            # Point Mass, right?
            self.flagpm == 1            

             # =========================== Saving Sim Inputs ===============

             
            # TODO add more parameters
            self.f01 = f01
            self.fd1 = fd1
            self.kc = kc
            self.q = q
            self.a0 = a0
            self.rad = rad
            self.a00 = a00
            self.ebarra = ebarra
            self.ham = ham
            self.eta = eta
            self.sigmas = sigmas
            self.sigmat = sigmat
            self.landadeb = landadeb
            self.epsilon = epsilon
            self.alfa = alfa
            self.nc = nc
            self.nper = nper
            self.npp = npp
            self.nperfin = nperfin
            self.zcmax = zcmax
            self.zcmin = zcmin
            self.dzc = dzc
            self.fixzc = fixzc
            self.tolerance = tolerance
            
            parameters = [f01,fd1,kc,q,a0,rad,a00,ebarra,ham,eta,sigmas,
                 sigmat,landadeb,epsilon,alfa,nc,nper,npp,nperfin,zcmax,
                 zcmin,dzc,fixzc,tolerance]  # Converted to an unique array
            parameters = np.array(parameters)
            
            # Save the parameters too an ouput file
           
            
            with open("dForceproject" +date +"/inputs" +iteration +".txt","w+") as file1:
                file1.write(f"""res freq 1 (kHz) : { f01 :.2f}
drive freq 1 (kHz) : { fd1 :.2f}
k1 (N/m) : { kc :.2f}
Q1 (adim) : { q :.2f}
A01 (nm) : { a0*1E9 :.2f}
Tip Radius (nm) : { rad *1E9:.2f}
a0 (nm) : { a00 * 1E9 :.2f}
Young's Modulus (MPa) : { emuestraa :.2f}
Hamaker Const(10^-20 J) : { ham *1E20:.2f}
Viscosity (Pa s) : { eta :.2f}
sigma surface (adim) : { sigmas :.2f}
sigma tip(adim) : { sigmat :.2f},
Sample thickness (nm): {self.toplevel1.sahe_sp_var:.2f},
Debye length : { landadeb :.2f}
Epsilon : { epsilon :.2f}
alfa : { alfa :.2f}
nc : { nc :.2f}
N periods osc : { nper }
N point period : { npp }
N per to the ss : { nperfin}
zc max (nm) : { zcmax * 1E9 :.2f}
zc min (nm) : { zcmin * 1E9 :.2f}
delta zc (nm) : { dzc * 1E9:.2f}
zc fixed td(nm) : { fixzc * 1E9:.2f}
tolerance : { tolerance :.2f}
""")

            
            t_output = runner(parameters, date, iteration, self)

         


# =========================Files handling w Python =======================
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            self.piat = os.getcwd()
            self.naind = "inputs" + str(self.flagin) + ".txt"
            self.nazcd = "zcdom" + str(self.flagin) + ".csv"
            self.natd = "tdom" + str(self.flagin) + ".csv"
            self.nawd = "wdom" + str(self.flagin) + ".csv"
           
# ========================Obtaining outputs zc domain ====================
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        # Bimodal mode
        else:    
            print("simulating bimodal")
            self.flagpm = 0
           
            if (self.flagfile == 0):
                # =========================== Saving Sim Inputs ===========
                # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                
                parameters = [f01,f02,kc,kc2,q,q2,a0,a02,rad,a00,ebarra,ham,
                    eta,sigmas,sigmat,landadeb,epsilon,alfa,nc,nper,npp,nperfin,
                    zcmax,zcmin,dzc,fixzc,tolerance]  # Converted to an unique array
                parameters = np.array(parameters)
                
                with open("dForceproject" + date + "/inputs" + iteration + ".txt","w+")  as file1:
                    file1.write(f"""res freq 1 (kHz) : { f01 :.2f}
res freq 2 (kHz) : { f02 :.2f}
k1 (N/m) : { kc :.2f}
k2 (N/m) : { kc2 :.2f}
Q1 (adim) : { q :.2f}
Q2 (adim) : { q2 :.2f}
A01 (nm) : { a0 *1E9 :.2f}
A02 (nm) : { a02 * 1E9 :.2f}
Tip Radius (nm) : { rad :.2f}
a0 (nm) : { a00*1E9 :.2f}
Young's Modulus (MPa) : { emuestraa :.2f}
Hamaker Const(10^-20 J) : { ham *1E20:.2f}
Viscosity (Pa s) : { eta :.2f}
sigma surface (adim) : { sigmas :.2f}
sigma tip(adim) : { sigmat :.2f},
Sample thickness (nm): {self.toplevel1.sahe_sp_var:.2f},
LJdepth (nN) : { landadeb :.2f}
LJlength (nm) : { epsilon :.2f}
Tip magn (A m2) : { alfa :.2f}
Surf magn (A m2) : { nc :.2f}
N periods osc : { nper }
N point period : { npp }
N per to the ss : { nperfin}
zc max (nm) : { zcmax *1E9:.2f}
zc min (nm) : { zcmin *1E9:.2f}
delta zc (nm) : { dzc * 1E9 :.2f}
zc fixed td (nm) : { fixzc *1E9:.2f}
""")
                

                date = self.date
                iteration = str(self.flagin)
                
                t_output = runnerEB(parameters, date, iteration,self)

# =========================Files handling w Python =======================
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                self.piat = os.getcwd()
                self.naind = "inputs" + str(self.flagin) + ".txt"
                self.nazcd = "zcdom" + str(self.flagin) + ".csv"
                self.natd = "tdom" + str(self.flagin) + ".csv"
                self.nawd = "wdom" + str(self.flagin) + ".csv"
              
            else:
                # self.flagin=0
                print("Mist hyper Mist")
        
        self.button_abort.config(state=DISABLED)
        self.button_run.config(state=NORMAL)
        self.button_plot.config(state=NORMAL)
        self.button_about.config(state=NORMAL)
        self.button_quit.config(state=NORMAL)
        self.button_examples.config(state=NORMAL)
        self.toplevel1.update_idletasks()
        self.toplevel1.update()
        
           
    # TODO: not implemented
    def save_input(self):
        return f"""A01 (nm): {self.a01_sp.get()}
A02(nm): {self.a02_sp.get()}

kc1_sp: {self.kc1_sp.get()}
kc2_sp: {self.kc2_sp.get()}
f01_sp: {self.f01_sp.get()}
f02_sp: {self.f02_sp.get()}
q1_sp: {self.q1_sp.get()}
q2_sp: {self.q2_sp.get()}
tolerance_sp: {self.tolerance_sp.get()}
tr_sp: {self.tr_sp.get()}

copy_fre1: {self.copy_fre1_var.get()}
copy_fre2: {self.copy_fre2_var.get()}

amp_x1: {self.amp_x1_var.get()}
amp_x2: {self.amp_x2_var.get()}
amp_x3: {self.amp_x3_var.get()}
amp_y1: {self.amp_y1_var.get()}
amp_y2: {self.amp_y2_var.get()}
amp_y3: {self.amp_y3_var.get()}
amps2_cb: {self.amps2_cb_var.get()}
ampst_cb: {self.ampst_cb_var.get()}


bias: {self.bias_var.get()}
compareflag: {self.compareflag}

ctime_cb: {self.ctime_cb_var.get()}

deltazc_sp: {self.deltazc_sp.get()}


dmin_x1: {self.dmin_x1_var.get()}
dmin_x2: {self.dmin_x2_var.get()}
dmin_x3: {self.dmin_x3_var.get()}
dmin_y1: {self.dmin_y1_var.get()}
dmin_y2: {self.dmin_y2_var.get()}
dmin_y3: {self.dmin_y3_var.get()}

ets2_cb: {self.ets2_cb_var.get()}
ets_x1: {self.ets_x1_var.get()}
ets_x2: {self.ets_x2_var.get()}
ets_x3: {self.ets_x3_var.get()}
ets_y1: {self.ets_y1_var.get()}
ets_y2: {self.ets_y2_var.get()}
ets_y3: {self.ets_y3_var.get()}
etst_cb: {self.etst_cb_var.get()}

famp1: {self.famp1_var.get()}
famp2: {self.famp2_var.get()}
fampt: {self.fampt_var.get()}
fases2_cb: {self.fases2_cb_var.get()}
fasest_cb: {self.fasest_cb_var.get()}
fcorrection_q1: {self.fcorrection_q1_var.get()}
fcorrection_q2: {self.fcorrection_q2_var.get()}
fd1_sp: {self.fd1_sp.get()}
fd2_sp: {self.fd2_sp.get()}
fd_cb: {self.fd_cb_var.get()}
fixedzc_sp: {self.fixedzc_sp.get()}
flagaux: {self.flagaux}
flagfile: {self.flagfile}
flagin: {self.flagin}
flagpm: {self.flagpm}
flagpnam: {self.flagpnam}
fmax_x1: {self.fmax_x1_var.get()}
fmax_x2: {self.fmax_x2_var.get()}
fmax_x3: {self.fmax_x3_var.get()}
fmax_y1: {self.fmax_y1_var.get()}
fmax_y2: {self.fmax_y2_var.get()}
fmax_y3: {self.fmax_y3_var.get()}
fpha1: {self.fpha1_var.get()}
fpha2: {self.fpha2_var.get()}
fphat: {self.fphat_var.get()}
ft_cb: {self.ft_cb_var.get()}

indent_cb: {self.indent_cb_var.get()}
it_cb: {self.it_cb_var.get()}





mass_model: {self.mass_model.get()}
maforce_cb: {self.maforce_cb_var.get()}


mdef_cb: {self.mdef_cb_var.get()}
mdef_x1: {self.mdef_x1_var.get()}
mdef_x2: {self.mdef_x2_var.get()}
mdef_x3: {self.mdef_x3_var.get()}
mdef_y1: {self.mdef_y1_var.get()}
mdef_y2: {self.mdef_y2_var.get()}
mdef_y3: {self.mdef_y3_var.get()}
mdis_cb: {self.mdis_cb_var.get()}
multifrequency_mode: {self.multifrequency_mode.get()}
mutip_sp: {self.mutip_sp.get()}
nper_sp: {self.nper_sp.get()}
nperfin_sp: {self.nperfin_sp.get()}
npp_sp: {self.npp_sp.get()}
phase_x1: {self.phase_x1_var.get()}
phase_x2: {self.phase_x2_var.get()}
phase_x3: {self.phase_x3_var.get()}
phase_y1: {self.phase_y1_var.get()}
phase_y2: {self.phase_y2_var.get()}
phase_y3: {self.phase_y3_var.get()}
phi1a1_cb: {self.phi1a1_cb_var.get()}
phi1a2_cb: {self.phi1a2_cb_var.get()}
phi2a1_cb: {self.phi2a1_cb_var.get()}
phi2a2_cb: {self.phi2a2_cb_var.get()}
ps_cb: {self.ps_cb_var.get()}
pts2_cb: {self.pts2_cb_var.get()}
pts_x1: {self.pts_x1_var.get()}
pts_x2: {self.pts_x2_var.get()}
pts_x3: {self.pts_x3_var.get()}
pts_y1: {self.pts_y1_var.get()}
pts_y2: {self.pts_y2_var.get()}
pts_y3: {self.pts_y3_var.get()}
ptst_cb: {self.ptst_cb_var.get()}

rescur: {self.rescur_var.get()}
simulate_retrace: {self.simulate_retrace_var.get()}


tcd_x1: {self.tcd_x1_var.get()}
tcd_x2: {self.tcd_x2_var.get()}
tcd_x3: {self.tcd_x3_var.get()}
tcd_y1: {self.tcd_y1_var.get()}
tcd_y2: {self.tcd_y2_var.get()}
tcd_y3: {self.tcd_y3_var.get()}


virial2_cb: {self.virial2_cb_var.get()}
virialt_cb: {self.virialt_cb_var.get()}



vt_cb: {self.vt_cb_var.get()}
vts_x1: {self.vts_x1_var.get()}
vts_x2: {self.vts_x2_var.get()}
vts_x3: {self.vts_x3_var.get()}
vts_y1: {self.vts_y1_var.get()}
vts_y2: {self.vts_y2_var.get()}
vts_y3: {self.vts_y3_var.get()}
ymt_sp: {self.ymt_sp.get()}
z1_cb: {self.z1_cb_var.get()}
z2_cb: {self.z2_cb_var.get()}
z3_cb: {self.z3_cb_var.get()}
z4_cb: {self.z4_cb_var.get()}
zc_x1: {self.zc_x1_var.get()}
zc_x2: {self.zc_x2_var.get()}
zc_x3: {self.zc_x3_var.get()}
zc_y1: {self.zc_y1_var.get()}
zc_y2: {self.zc_y2_var.get()}
zc_y3: {self.zc_y3_var.get()}
zcmax_sp: {self.zcmax_sp.get()}
zcmin_sp: {self.zcmin_sp.get()}
zt_cb: {self.zt_cb_var.get()}

# Forces
hertz_check: {self.hertz_check_var.get()}
adhy_check: {self.adhy_check_var.get()}
magnetic_dipole_check: {self.magnetic_dipole_check_var.get()}
magnetic_exponential_check: {self.magnetic_exponential_check_var.get()}
becc_check: {self.becc_check_var.get()}
behc_check: {self.behc_check_var.get()}
bepc_check: {self.bepc_check_var.get()}
custom1_check: {self.custom1_check_var.get()}
custom2_check: {self.custom2_check_var.get()}
dlvo_check: {self.dlvo_check_var.get()}
dmt_check: {self.dmt_check_var.get()}
jkr_check: {self.jkr_check_var.get()}
kv_cone_check: {self.kv_cone_check_var.get()}
kv_punch_check: {self.kv_punch_check_var.get()}
kv_sphere_check: {self.kv_sphere_check_var.get()}
lennardJones_check: {self.lennardJones_check_var.get()}
punch_check: {self.punch_check_var.get()}
sneddon_check: {self.sneddon_check_var.get()}
tatara_check: {self.tatara_check_var.get()}
vdw_check: {self.vdw_check_var.get()}
viscosity_check: {self.viscosity_check_var.get()}
viscosity_cone_check: {self.viscosity_cone_check_var.get()}
viscosity_punch_check: {self.viscosity_punch_check_var.get()}

#Sample parameters
sr_sp: {self.toplevel1.sr_sp_var}
saym_sp : {self.toplevel1.saym_sp_var }
saym2_sp: {self.toplevel1.saym2_sp_var}
musam_sp: {self.toplevel1.musam_sp_var}
# Electrical
epsilon_sp: {self.toplevel1.epsilon_sp_var}
sigmas_sp: {self.toplevel1.sigmas_sp_var}
sigmat_sp: {self.toplevel1.sigmat_sp_var}
landadeb_sp: {self.toplevel1.landadeb_sp_var}
# Other
ham_sp: {self.toplevel1.ham_sp_var}
ham2_sp: {self.toplevel1.ham2_sp_var}
a00_sp: {self.toplevel1.a00_sp_var}
angle_sp: {self.toplevel1.angle_sp_var}
sahe_sp: {self.toplevel1.sahe_sp_var}
bec_bonded: {self.toplevel1.bec_bonded_var}
# Viscosity
mvis_sp: {self.toplevel1.mvis_sp_var}

# Lennad Jhones
sigma_lj: {self.toplevel1.sigma_lj_var}
epsilon_lj: {self.toplevel1.epsilon_lj_var}

# User customized
custom1_sp: {self.toplevel1.custom1_sp_var}
custom2_sp: {self.toplevel1.custom2_sp_var}
custom3_sp : {self.toplevel1.custom3_sp_var }
custom4_sp: {self.toplevel1.custom4_sp_var}

# Magnetic
magnetic_m_tip: {self.toplevel1.magnetic_m_tip}
magnetic_m_sample: {self.toplevel1.magnetic_m_sample}
magnetic_B0_sample: {self.toplevel1.magnetic_B0_sample}
magnetic_k_sample: {self.toplevel1.magnetic_k_sample}
        
        """
        
    def on_plot_ex_clicked(self):       
        
        palette = ['b.--', 'r.--']
        # http://matplotlib.org/users/customizing.html
        matplotlib.rcParams.update(
            {'font.size': 16, 'figure.figsize': [12, 10]})       
        
        if self.piat == None:
            print("Nothing to plot, run a simulation first")
            return
      

        for comp_counter in range(self.compareflag + 1):
            if self.flagpm == 1:
                # ========Reading files for point-mass outputs===========
                # ||||||||||||||||||||||||||||||||||||||||||                
          
                if self.flagfile == 0:                 # Is there any file load?
                
                    # if zero, there has been a simulation runned by the user
                    os.chdir(self.piat + "/" + self.palo1)
                    tdom = np.genfromtxt(self.natd,skip_header = 1,delimiter = ",")
                    zcdom = np.genfromtxt(self.nazcd,skip_header= 1,delimiter = ",")
                    wdom = np.genfromtxt(self.nawd,skip_header= 1,delimiter = ",")
                else:
                    # I think this is the case were the user loaded something
                    # This is not used anymore
                    directionstdom = [self.tdomload, self.tdomload1]
                    directionszdom = [self.zdomload, self.zdomload1]
                    directionswdom = [self.wdomload, self.wdomload1]
                    tdom = np.genfromtxt(directionstdom[comp_counter],skip_header = 1,delimiter = ",")
                    zcdom = np.genfromtxt(directionszdom[comp_counter],skip_header = 1,delimiter = ",")
                    wdom = np.genfromtxt(directionswdom[comp_counter],skip_header = 1,delimiter = ",")
                tt = tdom[:, 0]
                zt = tdom[:, 1]
                vt = tdom[:, 2]
                forcet = tdom[:, 3]
                ident = tdom[:, 4]
                famptex = tdom[:, 5]
                fphatex = tdom[:, 6]
                freq_basex = tdom[:, 7]
                
                
                zc = zcdom[:, 0]
                amp1d = zcdom[:, 1]
                fase1d = zcdom[:, 2]
                dmind = zcdom[:, 3]
                #dmaxd = zcdom[:, 4]  # kill
                fmaxd = zcdom[:, 4]
                ets1d = zcdom[:, 5]
                vts1d = zcdom[:, 6]
                #fmedd = zcdom[:, 8]  # kill
                defld = zcdom[:, 7]
                tcd = zcdom[:, 8]  # new
                pts1d = zcdom[:, 9]  # new
                w1 = wdom[:, 0]
                w2 = wdom[:, 1]
                fre_amplitude = wdom[:, 2]

            else:
                os.chdir(self.piat + "/" + self.palo1)  # new trying
                # ========Reading files for Euler Bernoulli outputs===========
                # ||||||||||||||||||||||||||||||||||||||||||
                tdom = np.genfromtxt(self.natd,skip_header = 1,delimiter = ",")
                wdom = np.genfromtxt(self.nawd,skip_header= 1,delimiter = ",")
                zcdom = np.genfromtxt(self.nazcd,skip_header = 1,delimiter = ",")
                
                tt = tdom[:, 0]
                zt = tdom[:, 1]
                vt = tdom[:, 2]
                z1 = tdom[:, 3]
                v1 = tdom[:, 4]
                z2 = tdom[:, 5]
                v2 = tdom[:, 6]                
                forcet = tdom[:, 7]
                ident = tdom[:, 8]
                famp1ex = tdom[:, 9]
                fpha1ex = tdom[:, 10]
                famp2ex = tdom[:, 11]
                fpha2ex = tdom[:, 12]
                famptex = tdom[:, 13]
                fphatex = tdom[:, 14]
                freq_basex = tdom[:, 15]
                # Bimodal data
                print("Creating vectors from files for Bimodal ARRAYS ZC distance")
                # ========Material1=================
                
                
                
                zc = zcdom[:, 0]
                amp1d = zcdom[:, 1]
                fase1d = zcdom[:, 2]
                amp2d = zcdom[:, 3]
                fase2d = zcdom[:, 4]
                dmind = zcdom[:, 5]
                dmaxd = zcdom[:, 6]
                fmaxd = zcdom[:, 7]
                ets1d = zcdom[:, 8]
                ets2d = zcdom[:, 9]
                vts1d = zcdom[:, 10]
                vts2d = zcdom[:, 11]
                
                defld = zcdom[:, 12]
                tcd = zcdom[:, 13]
                pts1d = zcdom[:, 14]  # new
                pts2d = zcdom[:, 15]  # new
                
                w1 = wdom[:, 0]
                w2 = wdom[:, 1]
                fre_amplitude = wdom[:, 2]
                
              
            # ==================================================
            # End of getting arrays from simulated DATA
            # ======================PLOT========================
            
            if self.flagin != 0 or self.flagfile == 1:

                # ========Building the big Switch===========
                # ||||||||||||||||||||||||||||||||||||||||||
                
                def amp1(self):                    
                    # ========Amp 1 vs. zc=================
                    figure(1)
                    grid(False)
                    # title('$A_1$($z_c$)')
                    ylabel('$A_1 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, amp1d, palette[comp_counter])
                    savefig('amp1Vzc' + str(self.flagin) + '.png')

                def pha1(self):
                    # ========Phase 1 vs. zc=================
                    figure(2)
                    grid(False)
                    # title('$\phi_1$($z_c$)')
                    ylabel('$\\phi_1 \\; (deg)$',weight='bold',size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, fase1d, palette[comp_counter])
                    savefig('phi1Vzc' + str(self.flagin) + '.png')

                def den1(self):
                    # ========Diss Energy 1st vs. zc=================
                    figure(3)
                    grid(False)
                    #title('First disp. energy vs. Average distance')
                    ylabel('$E_{ts1}  \\; (10^{-20} J)$',weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, ets1d, palette[comp_counter])
                    savefig('ets1Vzc' + str(self.flagin) + '.png')

                def dpo1(self):
                    # ========Diss Power 1st vs. zc=================
                    figure(4)
                    grid(False)
                    #title('First disp. power vs. Average distance')
                    ylabel(
                        '$P_{ts1} \\; (10^{-20} W)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, pts1d, palette[comp_counter])
                    savefig('pts1Vzc' + str(self.flagin) + '.png')

                def vir1(self):
                    # ========Virial 1st vs. zc=================
                    figure(5)
                    grid(False)
                    #title('Virial first vs. Average distance')
                    ylabel(
                        '$V_{ts1} \\; (10^{-20} J )$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, vts1d, palette[comp_counter])
                    savefig('vts1Vzc' + str(self.flagin) + '.png')

                def def1(self):
                    # ========Deflection vs. zc=================
                    """
                    figure(6)
                    grid(False)
                    
                    
                    #title('Deflection vs. Average distance')
                    ylabel(
                        '$Deflection \\; (nm)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, defld, palette[comp_counter])
                    """
                    
                    
                    fig,ax = plt.subplots()
                    ax.grid(False)
                    ax.set_ylabel(
                        '$Deflection \\; (nm)$',
                        weight='bold',
                        size='x-large')
                    ax.set_xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    ax.plot(zc, defld, palette[comp_counter])
                    
                    
                    """
                    # TODO
                    ax2 = ax.twinx()
                    #limits = ax2.get_ylim()
                    #ax2.set_ylim(limits[0]*self.kc,limits[1]*self.kc)
                    ax2.set_ylabel(
                        '$Average Force \\; (nN)$',
                        weight='bold',
                        size='x-large')
                    """
                    
      
                def dmin1(self):
                    # ========Dmin vs. zc=================
                    figure(7)
                    grid(False)
                    #title('Minimum Distance vs. Average distance')
                    ylabel('$d_{min} \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, dmind, palette[comp_counter])
                    savefig('dminVzc' + str(self.flagin) + '.png')

                def maf(self):
                    # ========Fmax vs. zc=================
                    figure(8)
                    grid(False)
                    #title('Maximum Force vs. Average distance')
                    ylabel('$F_{max} \\; (pN)$', weight='bold', size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, fmaxd, palette[comp_counter])
                    savefig('fmaxVzc' + str(self.flagin) + '.png')

                def tco(self):
                    # ========Contact time vs. zc=================
                    figure(9)
                    grid(False)
                    #title('Contact Time (tc/T) vs. Average distance (nm)')
                    ylabel(
                        '$Contact \\;Time \\; (t/T)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, tcd, palette[comp_counter])
                    savefig('ctimeVzc' + str(self.flagin) + '.png')

                def amp2(self):
                    # ========Amp 2 vs. zc=================
                    figure(10)
                    grid(False)
                    # title('$A_2$($z_c$)')
                    ylabel('$A_2 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, amp2d, palette[comp_counter])
                    savefig('amp2Vzc' + str(self.flagin) + '.png')

                def pha2(self):
                    # ========Phase 2 vs. zc=================
                    figure(11)
                    grid(False)
                    # title('$\phi_2$($z_c$)')
                    ylabel(
                        '$\\phi_2 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, fase2d, palette[comp_counter])
                    savefig('phi2Vzc' + str(self.flagin) + '.png')

                def den2(self):
                    # ========Diss Energy 2nd vs. zc=================
                    figure(12)
                    grid(False)
                    #title('Second disp. energy vs. Average distance')
                    ylabel(
                        '$E_{ts2}  \\; (10^{-20} J)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, ets2d, palette[comp_counter])
                    savefig('ets2Vzc' + str(self.flagin) + '.png')

                def dpo2(self):
                    # ========Diss Power 1st vs. zc=================
                    figure(13)
                    grid(False)
                    #title('Second disp. power vs. Average distance')
                    ylabel(
                        '$P_{ts2} \\; (10^{-20} W)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, pts2d, palette[comp_counter])
                    savefig('pts2Vzc' + str(self.flagin) + '.png')

                def vir2(self):
                    # ========Virial 2nd vs. zc=================
                    figure(14)
                    grid(False)
                    #title('Virial second vs. Average distance ')
                    ylabel(
                        '$V_{ts2} \\; ( 10^{-20} J)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    plot(zc, vts2d, palette[comp_counter])
                    savefig('vts2Vzc' + str(self.flagin) + '.png')
                    
                def indentation_graph(self):
                    # ========Indentation vs. zc=================
                    figure(35)
                    grid(False)
                    #title('Virial second vs. Average distance ')
                    ylabel('$I_{max} \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    I_max = - dmind
                    I_max[I_max<=0] = 0
                    plot(zc, I_max, palette[comp_counter])
                    savefig('ImaxVzc' + str(self.flagin) + '.png')
                    # ----------------------------------------------------------------
                    # ===================Time Domain===========================
                    # ----------------------------------------------------------------
                    
                def lateralresolution_graph(self):
                    # ========lateral resolution vs. zc=================
                    
                    I_max = - dmind 
                    I_max[I_max<=0] = 0
                                         
                    if self.punch_check_var.get() or self.viscosity_punch_check_var.get() or self.bepc_check_var.get() or self.kv_punch_check_var.get():
                        # Punch case
                        R_tip = float(self.tr_sp.get())
                        contact_radius = R_tip * (I_max>0)
                    
                    # Spherical case
                    elif self.hertz_check_var.get() or self.viscosity_check_var.get():                        
                        R_tip = float(self.tr_sp.get()) 
                        contact_radius = np.sqrt(R_tip*I_max)
                        
                    # Other sperical-like cases
                    elif self.jkr_check_var.get() or self.tatara_check_var.get() or self.DELR_check_var.get():                        
                        R_tip = float(self.tr_sp.get()) 
                        contact_radius = np.sqrt(R_tip*I_max)
                        
                    # Spherical tip with bec
                    elif self.behc_check_var.get() or self.kv_sphere_check_var.get():
                        
                        R_tip = float(self.tr_sp.get()) 
                        sahe = self.toplevel1.sahe_sp_var # Height of the sample for BEC
                        bonded = self.toplevel1.bec_bonded_var # Sample bonded for BEC
                        nu = self.toplevel1.musam_sp_var
                        
                        if bonded:
                            # bonded case
                            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
                            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
                        else:
                            # Unbonded case
                            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
                            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
                        
                        contact_radius = np.sqrt(R_tip*I_max) * (1-2 * alpha0 * np.sqrt(R_tip*I_max) / (3 * np.pi * sahe))
                        
                    # Conical case
                    elif self.sneddon_check_var.get() or self.viscosity_cone_check_var.get():
                        
                        angle_tip = self.toplevel1.angle_sp_var * np.pi / 180
                        contact_radius = 2 * np.tan(angle_tip) / np.pi * I_max    
                    
                    # Conical tip with bec
                    elif self.becc_check_var.get() or self.kv_cone_check_var.get():
                        
                        angle_tip = self.toplevel1.angle_sp_var * np.pi / 180
                        sahe = self.toplevel1.sahe_sp_var # Height of the sample for BEC
                        bonded = self.toplevel1.bec_bonded_var # Sample bonded for BEC
                        nu = self.toplevel1.musam_sp_var
                        if bonded:
                            # bonded case
                            alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
                            beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
                        else:
                            # Unbonded case
                            alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
                            beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
                        # Conical case
                         
                        angle_tip = self.toplevel1.angle_sp_var * np.pi / 180
                        contact_radius = 2 * np.tan(angle_tip) / np.pi * I_max * (1- 2 *alpha0*np.tan(angle_tip) * I_max / (np.pi**2*sahe) )
                        
                        
                    # Nanowire case
                    elif self.nanowire_check_var.get() or self.viscosity_nanowire_check_var.get():                        
                        R_tip = float(self.tr_sp.get()) 
                        Ic = R_tip  
                        def calculate_contact_radius_nanowire(I_max):
                            if I_max<Ic:
                                contact_radius = np.sqrt(R_tip*I_max)
                            else:
                                contact_radius = np.sqrt(R_tip*Ic)
                            return contact_radius
                        contact_radius = np.array( [calculate_contact_radius_nanowire(i_max) for i_max in I_max] )
                        
                    # Nanowire with bec
                    elif self.benc_check_var.get() or self.kv_nanowire_check_var.get():
                        
                        def calculate_contact_radius_nanowire_bec(I_max):
                            R_tip = float(self.tr_sp.get()) 
                            sahe = self.toplevel1.sahe_sp_var # Height of the sample for BEC
                            bonded = self.toplevel1.bec_bonded_var # Sample bonded for BEC
                            nu = self.toplevel1.musam_sp_var
                            
                            if bonded:
                                # bonded case
                                alpha0 = -(1.2876-1.4678*nu+1.34442*nu**2)/(1-nu)
                                beta0 = (0.6387-1.0277*nu+1.516*nu**2)/(1-nu)
                            else:
                                # Unbonded case
                                alpha0 = -0.374 * (3 - 2 * nu) / (1 - nu)
                                beta0 = 0.056 * (5 - 2 * nu) / (1 - nu)
                            
                            h = sahe
                            Rt = R_tip
                            pi = np.pi
                            sqrt = np.sqrt
                            I_critical = 3*(-h*pi*(4*Rt*alpha0 - 3*h*pi) - sqrt(3)*sqrt(-h**3*pi**3*(8*Rt*alpha0 - 3*h*pi)))/(8*Rt*alpha0**2)  
                            
                            if I_max<I_critical:
                                contact_radius = np.sqrt(R_tip*I_max) * (1-2 * alpha0 * np.sqrt(R_tip*I_max) / (3 * np.pi * sahe))
                            else:
                                contact_radius = np.sqrt(R_tip*I_critical) * (1-2 * alpha0 * np.sqrt(R_tip*I_critical) / (3 * np.pi * sahe))
                            return contact_radius
                        contact_radius = np.array( [calculate_contact_radius_nanowire_bec(i_max) for i_max in I_max] )
                            
                        
                        
                    else:
                        contact_radius = I_max * np.nan
                        
                    lateral_resolution = 2 * contact_radius
                    
                    #plt.figure(37)
                    fig,ax = plt.subplots()
                    ax.grid(False)
                    #ax = plt.gca()
                    ax2 = ax.twinx()
                    
                    ax.set_ylabel('$L_{r} \\; (nm)$', weight='bold', size='x-large')
                    ax.set_xlabel('$z_c\\; (nm)$', weight='bold', size='x-large')
                    ax2.set_ylabel('$a \\; (nm)$', weight='bold', size='x-large')                    
                    
                    ax.plot(zc, lateral_resolution, palette[comp_counter])
                    
                    lims = ax.get_ylim()
                    ax2.set_ylim(0.5 * lims[0],0.5*lims[1])
                    
                    
                    fig.savefig('LateralResolutionVzc' + str(self.flagin) + '.png')
                    # ----------------------------------------------------------------
                    # ===================Time Domain===========================
                    # ----------------------------------------------------------------

                def ztt(self):
                    # ========Instantaneous position vs. time=================
                    figure(15)
                    grid(False)
                    # title("$z_T$($t/T$)")
                    ylabel('$z_T \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, zt, palette[comp_counter])
                    savefig('totalposTime' + str(self.flagin) + '.png')

                def vtt(self):
                    # ========Instantaneous velocity vs. time=================
                    figure(16)
                    grid(False)
                    # title("$v_T$($t/T$)")
                    ylabel(
                        '$v_T \\; (nm/\\mu s)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, vt, palette[comp_counter])
                    savefig('totalveloTime' + str(self.flagin) + '.png')
                    
                def phase_space_graph(self):
                    # ==Instantaneous position vs. Instantaneous velocity======
                    figure(36)
                    grid(False)
                    # title("$v_T$($t/T$)")
                    ylabel(
                        '$v_T \\; (nm/\\mu s)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$z_T \\; (nm)$', weight='bold', size='x-large')
                    plot(zt, vt, palette[comp_counter])
                    savefig('phase_space' + str(self.flagin) + '.png')

                def ftt(self):
                    # ========Instantaneous force vs. time=================
                    figure(17)
                    grid(False)
                    # title('$f_T$($t/T$)')
                    ylabel('$f_T \\; (pN)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, forcet, palette[comp_counter])
                    savefig('totalforceTime' + str(self.flagin) + '.png')

                def zt1(self):
                    # ========Instantaneous position 1 vs. time================
                    figure(18)
                    grid(False)
                    # title("$z_1$($t/T$)")
                    ylabel('$z_1 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, z1, palette[comp_counter])
                    savefig('FirstposTime' + str(self.flagin) + '.png')

                def zt2(self):
                    # ========Instantaneous position 2 vs. time================
                    figure(19)
                    grid(False)
                    # title("$z_2$($t/T$)")
                    ylabel('$z_2 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, z2, palette[comp_counter])
                    savefig('SecondposTime' + str(self.flagin) + '.png')

                def zt3(self):
                    # ========Instantaneous position 3 vs. time================
                    import pdb
                    pdb.set_trace()
                    figure(20)
                    grid(False)
                    # title("$z_3$($t/T$)")
                    ylabel('$z_3 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, zt, palette[comp_counter])
                    savefig('ThirdposTime' + str(self.flagin) + '.png')

                def zt4(self):
                    # ========Instantaneous position 4 vs. time================
                    figure(21)
                    grid(False)
                    # title("$z_4$($t/T$)")
                    ylabel('$z_4 \\; (nm)$', weight='bold', size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, zt, palette[comp_counter])
                    savefig('FourthposTime' + str(self.flagin) + '.png')

                def p1a1(self):
                    # ========Phase 1 vs. Amp 21=================
                    figure(22)
                    grid(False)
                    # title("$\phi_1$($A_1$)")
                    ylabel(
                        '$\\phi_1 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$A_1\\; (nm)$', weight='bold', size='x-large')
                    plot(amp1d, fase1d, palette[comp_counter])
                    savefig('phi1Vamp1' + str(self.flagin) + '.png')

                def p2a2(self):
                    # ========Phase 2 vs. Amp 2=================
                    figure(23)
                    grid(False)
                    # title("$\phi_2$($A_2$)")
                    ylabel(
                        '$\\phi_2 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$A_2\\; (nm)$', weight='bold', size='x-large')
                    plot(amp2d, fase2d, palette[comp_counter])
                    savefig('phi2Vamp2' + str(self.flagin) + '.png')

                def p2a1(self):
                    # ========Phase 2 vs. Amp 1=================
                    figure(24)
                    grid(False)
                    # title("$\phi_2$($A_1$)")
                    ylabel(
                        '$\\phi_2 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$A_1\\; (nm)$', weight='bold', size='x-large')
                    plot(amp1d, fase2d, palette[comp_counter])
                    savefig('phi2Vamp1' + str(self.flagin) + '.png')

                def p1a2(self):
                    # ========Phase 1 vs. Amp 2=================
                    figure(25)
                    grid(False)
                    # title("$\phi_1$($A_2$)")
                    ylabel(
                        '$\\phi_1 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$A_2\\; (nm)$', weight='bold', size='x-large')
                    plot(amp2d, fase1d, palette[comp_counter])
                    savefig('phi1Vamp2' + str(self.flagin) + '.png')

                def famp1graph(self):
                    # ========Fourier amplitude 1=================
                    figure(26)
                    grid(False)
                    #title("Amplitude 1")
                    ylabel(
                        '$A_1\\; (logarithmic  units)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         np.log(famp1ex[:n]),
                         palette[comp_counter])
                    savefig('famp1' + str(self.flagin) + '.png')

                def fpha1graph(self):
                    # ========Fourier phase 1=================
                    figure(27)
                    grid(False)
                    #title("Phase 1")
                    ylabel(
                        '$\\phi_1 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         fpha1ex[:n],
                         palette[comp_counter])
                    savefig('fpha1' + str(self.flagin) + '.png')

                def famp2graph(self):
                    # ========Fourier amplitude 2=================
                    figure(28)
                    grid(False)
                    #title("Amplitude 2")
                    ylabel(
                        '$A_2\\; (logarithmic  units)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         np.log(famp2ex[:n]),
                         palette[comp_counter])
                    savefig('famp2' + str(self.flagin) + '.png')

                def fpha2graph(self):
                    # ========Fourier phase 2=================
                    figure(29)
                    grid(False)
                    #title("Phase 2")
                    ylabel(
                        '$\\phi_2 \\; (deg)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         fpha2ex[:n],
                         palette[comp_counter])
                    savefig('fpha2' + str(self.flagin) + '.png')

                def famptgraph(self):
                    # ========Fourier amplitude total=================
                    figure(30)
                    grid(False)
                    #title("Amplitude total")
                    ylabel(
                        '$A \\; (logarithmic \\, units)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         np.log(famptex[:n]),
                         palette[comp_counter])
                    savefig('fampT' + str(self.flagin) + '.png')

                def fphatgraph(self):
                    # ========Fourier phase total=================
                    figure(31)
                    grid(False)
                    #title("Phase total")
                    ylabel('$\\phi \\; (deg)$', weight='bold', size='x-large')
                    xlabel('$Freq  \\; (Hz)$', weight='bold', size='x-large')
                    n = int(freq_basex.size / 2)
                    plot(freq_basex[:n],
                         fphatex[:n],
                         palette[comp_counter])
                    savefig('fphaT' + str(self.flagin) + '.png')

                def iden(self):
                    # ========Fourier phase total=================
                    figure(32)
                    grid(False)
                    # title("Identation")
                    ylabel(
                        '$Indentation \\; (nm)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$(t/T)$', weight='bold', size='x-large')
                    plot(tt, ident, palette[comp_counter])
                    savefig('iden' + str(self.flagin) + '.png')

                def res_cur(self):
                    # ========Resonance curve=================
                    figure(33)
                    grid(False)
                    # title("Identation")
                    ylabel(
                        '$Amplitude \\; (nm)$',
                        weight='bold',
                        size='x-large')
                    xlabel('$f/f_{0}$', weight='bold', size='x-large')
                    plot(w1, fre_amplitude, palette[comp_counter])
                    plot(w2, fre_amplitude, palette[comp_counter])
                    savefig('res' + str(self.flagin) + '.png')

                def f_d(self):
                    # ========Fourier phase total=================
                    
                    figure(34)
                    grid(False)
                    # title("Identation")
                    ylabel('$Force \\; (pN)$', weight='bold', size='x-large')
                    xlabel(
                        '$Distance \\; (nm)$',
                        weight='bold',
                        size='x-large')
                    z_temp = zt + float(self.fixedzc_sp.get())                    
                    #plot(z_temp,forcet,palette[comp_counter])
                    plot(z_temp,forcet,"b-")
                    savefig('f_d' + str(self.flagin) + '.png')

                # ========Custom figures===========

                def multi_draw(opt_mixed_x,opt_mixed_y,curves_mixed,names_mixed,counter):

                    # deleting axix before move them
                    def make_patch_spines_invisible(ax):
                        ax.set_frame_on(True)
                        ax.patch.set_visible(False)
                        for sp in list(ax.spines.values()):
                            sp.set_visible(False)

                    fig = figure(100 + counter)
                    grid(False)
                    title("Custom graph")
                    ax = plt.gca()

                    fig.subplots_adjust(right=0.8, left=0.2)

                    x_axis = np.arange(len(curves_mixed))
                    counter = 0
                    for ii in range(len(curves_mixed)):
                        if opt_mixed_x[ii] == 1:
                            x_axis = curves_mixed[ii]
                            ax.set_xlabel(names_mixed[ii])

                    for ii in range(len(curves_mixed)):
                        if opt_mixed_y[ii] == 1:
                            counter = counter + 1
                            y_axis = curves_mixed[ii]
                            if counter == 1:
                                lns1 = ax.plot(
                                    x_axis, y_axis, 'k-', label=names_mixed[ii])
                                ax.set_ylabel(names_mixed[ii])
                                ax.yaxis.label.set_color('k')
                                ax.tick_params(axis='y', colors='k')
                            elif counter == 2:
                                ax2 = ax.twinx()
                                lns2 = ax2.plot(
                                    x_axis, y_axis, 'r-', label=names_mixed[ii])
                                ax2.set_ylabel(names_mixed[ii])
                                ax2.yaxis.label.set_color('r')
                                ax2.tick_params(axis='y', colors='r')
                                ax2.spines['right'].set_color('r')
                            elif counter == 3:
                                ax3 = ax.twinx()
                                lns2 = ax3.plot(
                                    x_axis, y_axis, 'g-', label=names_mixed[ii])
                                ax3.spines["right"].set_position(("axes", 1.2))
                                make_patch_spines_invisible(ax3)
                                ax3.spines["right"].set_visible(True)
                                ax3.set_ylabel(names_mixed[ii])
                                ax3.yaxis.label.set_color('g')
                                ax3.tick_params(axis='y', colors='g')
                                ax3.spines['right'].set_color('g')
                            elif counter == 4:
                                ax4 = ax.twinx()
                                lns2 = ax4.plot(
                                    x_axis, y_axis, 'b-', label=names_mixed[ii])
                                ax4.spines["right"].set_position(
                                    ("axes", -0.3))
                                make_patch_spines_invisible(ax4)
                                ax4.spines["right"].set_visible(True)
                                ax4.set_ylabel(names_mixed[ii])
                                ax4.yaxis.label.set_color('b')
                                ax4.tick_params(axis='y', colors='b')
                                ax4.spines['right'].set_color('b')
                    fig.tight_layout()
                    # END OF multi_draw function
                    
                
                # Create a Tuple for the command self.glade_ui.get_widget
                if self.ampst_cb_var.get():
                    amp1(self)
                if self.fasest_cb_var.get():
                    pha1(self)
                if self.etst_cb_var.get():
                    den1(self)
                if self.ptst_cb_var.get():
                    dpo1(self)
                if self.virialt_cb_var.get():
                    vir1(self)
                if self.mdef_cb_var.get():
                    def1(self)
                if self.mdis_cb_var.get():
                    dmin1(self)
                if self.maforce_cb_var.get():
                    maf(self)
                if self.ctime_cb_var.get():
                    tco(self)
                if self.amps2_cb_var.get():
                    amp2(self)
                if self.fases2_cb_var.get():
                    pha2(self)
                if self.ets2_cb_var.get():
                    den2(self)
                if self.pts2_cb_var.get():
                    dpo2(self)
                if self.virial2_cb_var.get():
                    vir2(self)
                if self.indent_cb_var.get():
                    indentation_graph(self)
                if self.zt_cb_var.get():
                    ztt(self)
                if self.vt_cb_var.get():
                    vtt(self)
                if self.ft_cb_var.get():
                    ftt(self)
                if self.z1_cb_var.get():
                    zt1(self)
                if self.z2_cb_var.get():
                    zt2(self)
                if self.z3_cb_var.get():
                    zt3(self)
                if self.z4_cb_var.get():
                    zt4(self)
                if self.phi1a1_cb_var.get():
                    p1a1(self)
                if self.phi2a2_cb_var.get():
                    p2a2(self)
                if self.phi2a1_cb_var.get():
                    p2a1(self)
                if self.phi1a2_cb_var.get():
                    p1a2(self)
                if self.famp1_var.get():
                    famp1graph(self)
                if self.fpha1_var.get():
                    fpha1graph(self)
                if self.famp2_var.get():
                    famp2graph(self)
                if self.fpha2_var.get():
                    fpha2graph(self)
                if self.fampt_var.get():
                    famptgraph(self)
                if self.fphat_var.get():
                    fphatgraph(self)
                if self.it_cb_var.get():
                    iden(self)
                if self.rescur_var.get():
                    res_cur(self)
                if self.fd_cb_var.get():
                    f_d(self)
                if self.ps_cb_var.get():
                    phase_space_graph(self)
                if self.lateralresolution_cb_var.get():
                    lateralresolution_graph(self)
                    
                curves_to_plot = [amp1,pha1,den1,dpo1,vir1,def1,dmin1,maf,tco,amp2,
                    pha2,den2,dpo2,vir2,indentation_graph,
                    ztt,vtt,ftt,zt1,zt2,zt3,zt4,p1a1,p2a2,
                    p2a1,p1a2,famp1graph,fpha1graph,famp2graph,fpha2graph,
                    famptgraph,fphatgraph,iden,res_cur,f_d,
                    phase_space_graph,
                    lateralresolution_graph]
                

                
                # ===========Draw custom figures=================
                
                if self.flagpm == 1:
                    
                    # calculation of the indentation from the minimum distance:
                    imax = -dmind
                    imax[imax<0] = 0
                    names_mixed = [
                        "$A \\; (nm)$",
                        "$\\phi\\; (deg)$",
                        "zc (nm)",
                        "$P \\; ( 10^{-16} J)$",
                        "$E \\; ( 10^{-20} J)$",
                        "Minimum distance (nm)",
                        "Maximum force (pN)",
                        "Contact time (t/T)",
                        "$V \\; ( 10^{-20} J)$",
                        "$I_{max} (nm)$"]
                    
                    curves_mixed = [amp1d, fase1d, zc, pts1d, ets1d, dmind, fmaxd, tcd, vts1d,imax]
                    ii = 1
                    opt_mixed_x = [
                        self.amp_x1_var.get(),
                        self.phase_x1_var.get() ,
                        self.zc_x1_var.get() ,
                        self.pts_x1_var.get() ,
                        self.ets_x1_var.get() ,
                        self.dmin_x1_var.get() ,
                        self.fmax_x1_var.get() ,
                        self.tcd_x1_var.get() ,
                        self.vts_x1_var.get(),
                        self.imax_x1_var.get()]
                    opt_mixed_y = [
                        self.amp_y1_var.get() ,
                        self.phase_y1_var.get() ,
                        self.zc_y1_var.get() ,
                        self.pts_y1_var.get() ,
                        self.ets_y1_var.get() ,
                        self.dmin_y1_var.get() ,
                        self.fmax_y1_var.get() ,
                        self.tcd_y1_var.get() ,
                        self.vts_y1_var.get(),
                        self.imax_y1_var.get()]                    
                    
                    if sum(opt_mixed_x) > 0:
                        multi_draw(opt_mixed_x, opt_mixed_y, curves_mixed, names_mixed, ii)
                        
                    ii = 2
                    opt_mixed_x = [
                        self.amp_x2_var.get(),
                        self.phase_x2_var.get() ,
                        self.zc_x2_var.get() ,
                        self.pts_x2_var.get() ,
                        self.ets_x2_var.get() ,
                        self.dmin_x2_var.get() ,
                        self.fmax_x2_var.get() ,
                        self.tcd_x2_var.get() ,
                        self.vts_x2_var.get(),
                        self.imax_x2_var.get()]
                    opt_mixed_y = [
                        self.amp_y2_var.get() ,
                        self.phase_y2_var.get() ,
                        self.zc_y2_var.get() ,
                        self.pts_y2_var.get() ,
                        self.ets_y2_var.get() ,
                        self.dmin_y2_var.get() ,
                        self.fmax_y2_var.get() ,
                        self.tcd_y2_var.get() ,
                        self.vts_y2_var.get(),
                        self.imax_y2_var.get()]
                    
                    
                    if sum(opt_mixed_x) > 0:
                        multi_draw(opt_mixed_x, opt_mixed_y, curves_mixed, names_mixed, ii)
                        
                    ii = 3
                    opt_mixed_x = [
                        self.amp_x3_var.get(),
                        self.phase_x3_var.get() ,
                        self.zc_x3_var.get() ,
                        self.pts_x3_var.get() ,
                        self.ets_x3_var.get() ,
                        self.dmin_x3_var.get() ,
                        self.fmax_x3_var.get() ,
                        self.tcd_x3_var.get() ,
                        self.vts_x3_var.get(),
                        self.imax_x3_var.get()]
                    opt_mixed_y = [
                        self.amp_y3_var.get() ,
                        self.phase_y3_var.get() ,
                        self.zc_y3_var.get() ,
                        self.pts_y3_var.get() ,
                        self.ets_y3_var.get() ,
                        self.dmin_y3_var.get() ,
                        self.fmax_y3_var.get() ,
                        self.tcd_y3_var.get() ,
                        self.vts_y3_var.get() , 
                        self.imax_y3_var.get()]
                    
                    
                    if sum(opt_mixed_x) > 0:
                        multi_draw(opt_mixed_x, opt_mixed_y, curves_mixed, names_mixed, ii)


                else:

                    opt_mixed_x = [
                       self.amp_x1_var.get(),
                        self.phase_x1_var.get() ,
                        self.zc_x1_var.get() ,
                        self.pts_x1_var.get() ,
                        self.ets_x1_var.get() ,
                        self.dmin_x1_var.get() ,
                        self.fmax_x1_var.get() ,
                        self.tcd_x1_var.get() ,
                        self.vts_x1_var.get(),
                        self.amp_x2_var.get(),
                        self.phase_x2_var.get() ,
                        self.zc_x2_var.get() ,
                        self.pts_x2_var.get() ,
                        self.ets_x2_var.get() ,
                        self.dmin_x2_var.get() ,
                        self.fmax_x2_var.get() ,
                        self.tcd_x2_var.get() ,
                        self.vts_x2_var.get()]
                    opt_mixed_y = [
                        self.amp_y1_var.get() ,
                        self.phase_y1_var.get() ,
                        self.zc_y1_var.get() ,
                        self.pts_y1_var.get() ,
                        self.ets_y1_var.get() ,
                        self.dmin_y1_var.get() ,
                        self.fmax_y1_var.get() ,
                        self.tcd_y1_var.get() ,
                        self.vts_y1_var.get(),                        
                        self.amp_y2_var.get() ,
                        self.phase_y2_var.get() ,
                        self.zc_y2_var.get() ,
                        self.pts_y2_var.get() ,
                        self.ets_y2_var.get() ,
                        self.dmin_y2_var.get() ,
                        self.fmax_y2_var.get() ,
                        self.tcd_y2_var.get() ,
                        self.vts_y2_var.get()]
                    curves_mixed = [
                        amp1d,
                        fase1d,
                        zc,
                        pts1d,
                        ets1d,
                        dmind,
                        fmaxd,
                        tcd,
                        vts1d,
                        amp2d,
                        fase2d,
                        zc,
                        pts2d,
                        ets2d,
                        dmind,
                        fmaxd,
                        tcd,
                        vts2d]
                    names_mixed = [
                        "$A_{1} \\; (nm)$",
                        "$\\phi_1 \\; (deg)$",
                        "zc (nm)",
                        "$$P_{1} \\; ( 10^{-16} W)$",
                        "$E_{1} \\; ( 10^{-20} J)$",
                        "dmin (nm)",
                        "maximum force (pN)",
                        "contact time (t/T)",
                        "$V_{1} \\; ( 10^{-20} J)$",
                        "A_2 (nm)",
                        "$\\phi_2 \\; (deg)$",
                        " (nm)",
                        "$P_{2} \\; ( 10^{-16} J)$",
                        "$E_{2} \\; ( 10^{-20} J)$",
                        "dmin (nm)",
                        "Maximum force (pN)",
                        "Contact time (t/T)",
                        "$V_{2} \\; ( 10^{-20} J)$"]
                    if sum(opt_mixed_x) > 0:
                        multi_draw(
                            opt_mixed_x,
                            opt_mixed_y,
                            curves_mixed,
                            names_mixed,
                            4)
                
                # ============Bias figure=====================
                if self.bias_var.get():
                    
                    if self.flagpm == 1:
                        opt_mixed_x = [1, 0, 0, 0, 0]
                        opt_mixed_y = [
                            0,
                            self.zt_cb_var.get(),
                            self.vt_cb_var.get(),
                            self.ft_cb_var.get(),
                            self.it_cb_var.get()]
                        curves_mixed = [tt, zt, vt, forcet, ident]
                        names_mixed = [
                            "time (t/T)",
                            "z (nm)",
                            "$v_T \\; (nm/\\mu s)$",
                            "forc (pN)",
                            "indentation (nm)"]
                        
                        multi_draw(opt_mixed_x, 
                                   opt_mixed_y, 
                                   curves_mixed, 
                                   names_mixed, 2)
                    else:
                        opt_mixed_x = [1, 0, 0, 0, 0, 0, 0]
                        opt_mixed_y = [
                            0,
                            self.zt_cb_var.get(),
                            self.vt_cb_var.get(),
                            self.ft_cb_var.get(),
                            self.it_cb_var.get(),
                            self.z1_cb_var.get(),
                            self.z2_cb_var.get()]
                        curves_mixed = [tt, zt, vt, forcet, ident, z1, z2]
                        names_mixed = [
                            "time (t/T)",
                            "z (nm)",
                            "$v_T \\; (nm/\\mu s)$",
                            "force (pN)",
                            "indentation (nm)",
                            "z1 (nm)",
                            "z2 (nm)"]
                        if sum(opt_mixed_x) > 0:
                            multi_draw(
                                opt_mixed_x, opt_mixed_y, curves_mixed, names_mixed, 3)

                # ||||||||||||||||||||||||||||||||||||||||||
                # =========Big switch ready to use==============
                
                
                
                # self.flagin+=1
                
                if self.flagfile == 0:
                    os.chdir(self.piat)
            else:
                print("You have nothing to visualize")  # -----Message
        show()		# intented for plotting the bounch of curves...
    
    def on_show_value_clicked(self):                
        window = ForceParametersMenuApp(self,
                            title = "All parameters",
                            image_path="value table",
                            sr_show = True,
                            saym_show = True,
                            saym2_show = True,
                            musam_show = True,
                            epsilon_show = True,
                            sigmas_show = True,
                            sigmat_show = True,
                            landadeb_show = True,
                            ham_show = True,
                            ham2_show = True,
                            a00_show = True,
                            angle_show = True,
                            sahe_show = True,
                            sample_bonded_show = True,
                            mvis_show = True,
                            custom_show = True,
                            rtip_show = True,
                            ymt_show = True,
                            magnetic_exp_show = True,
                            magnetic_dipole_show = True,
                            lennardjones_show = True,
                            text_width = 40)
        window.toplevel.grab_set()
        
        
    def on_about_clicked(self):
        window = credits_app(self.toplevel1)        
        window.credits_menu.grab_set()
        
    def on_examples_ex_clicked(self):
        window = examples_app(self)        
        window.examples_menu.grab_set()
        
        
        
    
            
            
    def single_mode_activated(self):
        self.pm_cb["state"] = 'normal' # point-mass activated
        self.eb_cb["state"] = 'normal' # Euler-bernuilly activated
        
        
        #Options of excitation bimodal are deactivated
        self.a02_sp["state"] = 'disabled'
        self.fd2_sp["state"] = 'disabled'
        self.copy_fre2["state"] = 'disabled'
        
        
        
        self.amps2_cb.disable()
        self.fases2_cb.disable()
        self.ets2_cb.disable()
        self.pts2_cb.disable()
        self.virial2_cb.disable()
                        
        #self.z1_cb["state"] = 'disabled'
        #self.z2_cb["state"] = 'disabled'
        
        
    def bimodal_mode_activated(self):
        self.mass_model.set("euler")
        
        self.pm_cb["state"] = 'disabled' # point-mass activated
        self.eb_cb["state"] = 'normal' # Euler-bernuilly activated
        
        self.a02_sp["state"] = 'normal'
        self.fd2_sp["state"] = 'normal'
        self.copy_fre2["state"] = 'normal'
        
        self.z1_cb["state"] = 'normal'
        self.z2_cb["state"] = 'normal'
        
     
        self.amps2_cb.enable()
        self.fases2_cb.enable()
        self.ets2_cb.enable()
        self.pts2_cb.enable()
        self.virial2_cb.enable()
        self.z2_cb["state"] = 'normal'
        self.euler_activated()
        
        
    def point_mass_activated(self):
        #Options for graphs that use bimodal are deactivated
        self.amps2_cb.disable()
        self.fases2_cb.disable()
        self.ets2_cb.disable()
        self.pts2_cb.disable()
        self.virial2_cb.disable()
        self.z1_cb["state"] = 'disabled'
        self.z2_cb["state"] = 'disabled'
        self.phi2a2_cb["state"] = 'disabled'
        self.phi2a1_cb["state"] = 'disabled'
        self.phi1a2_cb["state"] = 'disabled'
        self.famp1["state"] = 'disabled'
        self.famp2["state"] = 'disabled'
        self.fpha1["state"] = 'disabled'
        self.fpha2["state"] = 'disabled'
        
        self.kc2_sp["state"] = 'disabled'
        self.f02_sp["state"] = 'disabled'
        self.fcorrection_q2["state"] = 'disabled'
        self.q2_sp["state"] = 'disabled'
        
        
    def euler_activated(self):
        self.bmafm_cb["state"] = 'normal'
        
        #Options for graphs that use bimodal are activated
        #self.amps2_cb.enable()
        #self.fases2_cb.enable()
        #self.ets2_cb.enable()
        #self.pts2_cb.enable()
        #self.virial2_cb.enable()
        self.z1_cb["state"] = 'normal'
        self.z2_cb["state"] = 'normal'
        self.phi2a2_cb["state"] = 'normal'
        self.phi2a1_cb["state"] = 'normal'
        self.phi1a2_cb["state"] = 'normal'
        self.famp1["state"] = 'normal'
        self.famp2["state"] = 'normal'
        self.fpha1["state"] = 'normal'
        self.fpha2["state"] = 'normal'
            
        self.kc2_sp["state"] = 'normal'
        self.f02_sp["state"] = 'normal'
        self.fcorrection_q2["state"] = 'normal'
        self.q2_sp["state"] = 'normal'
        
    def on_abort_clicked(self):
        print("aborted")
        self.abort = 1
        
    
    def on_activate_all_checks(self):
        self.vdw_check.config(state=NORMAL)
        self.vdw_cone_check.config(state=NORMAL)
        self.vdw_punch_check.config(state=NORMAL)
        self.dlvo_check.config(state=NORMAL)
        self.magnetic_dipole_check.config(state=NORMAL)
        self.magnetic_exponential_check.config(state=NORMAL)
        self.lennardJones_check.config(state=NORMAL)
        self.hertz_check.config(state=NORMAL)
        self.sneddon_check.config(state=NORMAL)
        self.punch_check.config(state=NORMAL)
        self.dmt_check.config(state=NORMAL)
        self.jkr_check.config(state=NORMAL)
        self.tatara_check.config(state=NORMAL)
        self.behc_check.config(state=NORMAL)
        self.becc_check.config(state=NORMAL)
        self.bepc_check.config(state=NORMAL)
        self.viscosity_check.config(state=NORMAL)
        self.viscosity_cone_check.config(state=NORMAL)
        self.viscosity_punch_check.config(state=NORMAL)
        self.kv_sphere_check.config(state=NORMAL)
        self.kv_cone_check.config(state=NORMAL)
        self.kv_punch_check.config(state=NORMAL)
        self.adhy_check.config(state=NORMAL)
        self.custom1_check.config(state=NORMAL)
        self.custom2_check.config(state=NORMAL)
        self.benc_check.config(state=NORMAL)
        self.nanowire_check.config(state=NORMAL)
        self.DELR_check.config(state=NORMAL)
        self.viscosity_nanowire_check.config(state=NORMAL)
        self.kv_nanowire_check.config(state=NORMAL)
        
        # If all are empty then all are set to normal
        if not np.any(self.get_option_mask()):
            return
        # If the user is using customized forces all options are available
        elif self.custom1_check_var.get():
            return
        elif self.custom2_check_var.get():
            return
        
        self.vdw_check.config(state=DISABLED)
        self.vdw_cone_check.config(state=DISABLED)
        self.vdw_punch_check.config(state=DISABLED)        
        self.dlvo_check.config(state=DISABLED)
        self.magnetic_dipole_check.config(state=DISABLED)
        self.magnetic_exponential_check.config(state=DISABLED)
        self.lennardJones_check.config(state=DISABLED)
        self.hertz_check.config(state=DISABLED)
        self.sneddon_check.config(state=DISABLED)
        self.punch_check.config(state=DISABLED)
        self.dmt_check.config(state=DISABLED)
        self.jkr_check.config(state=DISABLED)
        self.tatara_check.config(state=DISABLED)
        self.behc_check.config(state=DISABLED)
        self.becc_check.config(state=DISABLED)
        self.bepc_check.config(state=DISABLED)
        self.viscosity_check.config(state=DISABLED)
        self.viscosity_cone_check.config(state=DISABLED)
        self.viscosity_punch_check.config(state=DISABLED)
        self.kv_sphere_check.config(state=DISABLED)
        self.kv_cone_check.config(state=DISABLED)
        self.kv_punch_check.config(state=DISABLED)
        self.adhy_check.config(state=DISABLED)
        self.benc_check.config(state=DISABLED)
        self.nanowire_check.config(state=DISABLED)
        self.DELR_check.config(state=DISABLED)
        self.viscosity_nanowire_check.config(state=DISABLED)
        self.kv_nanowire_check.config(state=DISABLED)
        
        self.custom1_check.config(state=NORMAL)
        self.custom2_check.config(state=NORMAL)
        self.magnetic_dipole_check.config(state=NORMAL)
        self.magnetic_exponential_check.config(state=NORMAL)        
        #self.lennardJones_check.config(state=NORMAL)
        
        if self.lennardJones_check_var.get():        
            self.magnetic_dipole_check.config(state=DISABLED)
            self.magnetic_exponential_check.config(state=DISABLED)   
            self.lennardJones_check.config(state=NORMAL)
            
        if self.dmt_check_var.get():
            #self.vdw_check.config(state=NORMAL)
            #self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            self.dmt_check.config(state=NORMAL)
            
        
           
        if self.adhy_check_var.get():
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            self.dmt_check.config(state=NORMAL)      
            self.hertz_check.config(state=NORMAL)
            self.behc_check.config(state=NORMAL)
            self.viscosity_check.config(state=NORMAL)
            self.kv_sphere_check.config(state=NORMAL)
            
        if self.vdw_check_var.get() :            
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=DISABLED)
            self.adhy_check.config(state=NORMAL)
            #self.dmt_check.config(state=NORMAL)
            
            self.hertz_check.config(state=NORMAL)
            self.behc_check.config(state=NORMAL)
            self.viscosity_check.config(state=NORMAL)
            self.kv_sphere_check.config(state=NORMAL)
        
        if self.vdw_cone_check_var.get() :            
            self.vdw_cone_check.config(state=NORMAL)
            self.dlvo_check.config(state=DISABLED)
            self.adhy_check.config(state=NORMAL)
            #self.dmt_check.config(state=NORMAL)            
            self.sneddon_check.config(state=NORMAL)
            self.becc_check.config(state=NORMAL)
            self.viscosity_cone_check.config(state=NORMAL)
            self.kv_cone_check.config(state=NORMAL)
            
        if self.vdw_punch_check_var.get() :            
            self.vdw_punch_check.config(state=NORMAL)
            self.dlvo_check.config(state=DISABLED)
            self.adhy_check.config(state=NORMAL)
            #self.dmt_check.config(state=NORMAL)
            
            self.punch_check.config(state=NORMAL)
            self.bepc_check.config(state=NORMAL)
            self.viscosity_punch_check.config(state=NORMAL)
            self.kv_punch_check.config(state=NORMAL)
            
        if self.dlvo_check_var.get():            
            self.dlvo_check.config(state=DISABLED)
            self.adhy_check.config(state=NORMAL)            
            self.hertz_check.config(state=NORMAL)
            self.behc_check.config(state=NORMAL)
            self.viscosity_check.config(state=NORMAL)
            self.kv_sphere_check.config(state=NORMAL)
            
            self.vdw_check.config(state=DISABLED)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            self.dmt_check.config(state=NORMAL)
            
            
            
        if self.hertz_check_var.get():
            self.hertz_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            
            self.behc_check.config(state=DISABLED)
            self.viscosity_check.config(state=DISABLED)
            self.kv_sphere_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.DELR_check_var.get():
            self.DELR_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            
        if self.sneddon_check_var.get():
            self.sneddon_check.config(state=NORMAL)
            self.vdw_cone_check.config(state=NORMAL)            
            self.becc_check.config(state=DISABLED)
            self.viscosity_cone_check.config(state=DISABLED)
            self.kv_cone_check.config(state=DISABLED)
            
        if self.punch_check_var.get():            
            self.vdw_punch_check.config(state=NORMAL)            
            self.punch_check.config(state=NORMAL)
            self.bepc_check.config(state=DISABLED)
            self.viscosity_punch_check.config(state=DISABLED)
            self.kv_punch_check.config(state=DISABLED)
            
        if self.jkr_check_var.get():
            self.jkr_check.config(state=NORMAL)
            
        if self.tatara_check_var.get():
            self.tatara_check.config(state=NORMAL)
            
        if self.behc_check_var.get():
            self.behc_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            
            self.hertz_check.config(state=DISABLED)            
            self.viscosity_check.config(state=DISABLED)
            self.kv_sphere_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.becc_check_var.get():
            self.becc_check.config(state=NORMAL)
            self.vdw_cone_check.config(state=NORMAL)
            self.sneddon_check.config(state=DISABLED)            
            self.viscosity_cone_check.config(state=DISABLED)
            self.kv_cone_check.config(state=DISABLED)
            
        if self.bepc_check_var.get():            
            self.vdw_punch_check.config(state=NORMAL)
            self.punch_check.config(state=DISABLED)
            self.bepc_check.config(state=NORMAL)
            self.viscosity_punch_check.config(state=DISABLED)
            self.kv_punch_check.config(state=DISABLED)
            
            
        if self.viscosity_check_var.get():
            self.viscosity_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            self.hertz_check.config(state=DISABLED)
            self.behc_check.config(state=DISABLED)            
            self.kv_sphere_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.viscosity_nanowire_check_var.get():
            self.viscosity_nanowire_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            self.hertz_check.config(state=DISABLED)
            self.behc_check.config(state=DISABLED)            
            self.kv_sphere_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.viscosity_cone_check_var.get():
            self.viscosity_cone_check.config(state=NORMAL)
            self.vdw_cone_check.config(state=NORMAL)
            self.sneddon_check.config(state=DISABLED)
            self.becc_check.config(state=DISABLED)            
            self.kv_cone_check.config(state=DISABLED)
            
            
        if self.viscosity_punch_check_var.get():
            
            self.vdw_punch_check.config(state=NORMAL)
            self.punch_check.config(state=DISABLED)
            self.bepc_check.config(state=DISABLED)
            self.viscosity_punch_check.config(state=NORMAL)
            self.kv_punch_check.config(state=DISABLED)
            
            
        if self.kv_sphere_check_var.get():
            self.kv_sphere_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            
            self.hertz_check.config(state=DISABLED)
            self.behc_check.config(state=DISABLED)
            self.viscosity_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.kv_nanowire_check_var.get():
            self.kv_nanowire_check.config(state=NORMAL)
            self.vdw_check.config(state=NORMAL)
            self.dlvo_check.config(state=NORMAL)
            self.adhy_check.config(state=NORMAL)
            
            self.hertz_check.config(state=DISABLED)
            self.behc_check.config(state=DISABLED)
            self.viscosity_check.config(state=DISABLED)
            self.dmt_check.config(state=DISABLED)
            
        if self.kv_cone_check_var.get():
            self.kv_cone_check.config(state=NORMAL)
            self.vdw_cone_check.config(state=NORMAL)
            self.sneddon_check.config(state=DISABLED)
            self.becc_check.config(state=DISABLED)
            self.viscosity_cone_check.config(state=DISABLED)
            
            
        if self.kv_punch_check_var.get():            
            self.vdw_punch_check.config(state=NORMAL)            
            self.punch_check.config(state=DISABLED)
            self.bepc_check.config(state=DISABLED)
            self.viscosity_punch_check.config(state=DISABLED)
            self.kv_punch_check.config(state=NORMAL)
            
            
        if self.nanowire_check_var.get():
            self.nanowire_check.config(state=NORMAL)            
        if self.benc_check_var.get():
            self.benc_check.config(state=NORMAL)
        
        
    def on_vdw_activated(self):
        self.on_activate_all_checks()
        if self.vdw_check_var.get():
            
            window = ForceParametersMenuApp(self,
                                            title = "van der Wals",
                                            image_path="vdw.png",
                                            sr_show = False,
                                            saym_show = False,
                                            saym2_show = False,
                                            musam_show = False,
                                            epsilon_show = False,
                                            sigmas_show = False,
                                            sigmat_show = False,
                                            landadeb_show = False,
                                            ham_show = True,
                                            ham2_show = False,
                                            a00_show = False,
                                            angle_show = False,
                                            sahe_show = False,
                                            mvis_show = False,
                                            rtip_show = True)
            window.toplevel.grab_set()
            
    def on_vdw_cone_activated(self):
        self.on_activate_all_checks()
        if self.vdw_cone_check_var.get():
            
            window = ForceParametersMenuApp(self,
                                            title = "Cone van der Wals",
                                            image_path="vdw_cone.png",
                                            sr_show = False,
                                            saym_show = False,
                                            saym2_show = False,
                                            musam_show = False,
                                            epsilon_show = False,
                                            sigmas_show = False,
                                            sigmat_show = False,
                                            landadeb_show = False,
                                            ham_show = True,
                                            ham2_show = False,
                                            a00_show = False,
                                            angle_show = True,
                                            sahe_show = False,
                                            mvis_show = False,
                                            rtip_show = False)
            window.toplevel.grab_set()
            
            
            
    def on_vdw_punch_activated(self):
        self.on_activate_all_checks()
        if self.vdw_punch_check_var.get():
            
            window = ForceParametersMenuApp(self,
                                            title = "Flat van der Wals",
                                            image_path="vdw_flat.png",
                                            sr_show = False,
                                            saym_show = False,
                                            saym2_show = False,
                                            musam_show = False,
                                            epsilon_show = False,
                                            sigmas_show = False,
                                            sigmat_show = False,
                                            landadeb_show = False,
                                            ham_show = True,
                                            ham2_show = False,
                                            a00_show = False,
                                            angle_show = False,
                                            sahe_show = False,
                                            mvis_show = False,
                                            rtip_show = True)
            window.toplevel.grab_set()

    def on_dlvo_activated(self):
        self.on_activate_all_checks()
        if self.dlvo_check_var.get():
               window = ForceParametersMenuApp(self,
                                    title = "DLVO",
                                    image_path="dlvo.png",
                                    sr_show = False,
                                    saym_show = False,
                                    saym2_show = False,
                                    musam_show = False,
                                    epsilon_show = True,
                                    sigmas_show = True,
                                    sigmat_show = True,
                                    landadeb_show = True,
                                    ham_show = True,
                                    ham2_show = False,
                                    a00_show = False,
                                    angle_show = False,
                                    sahe_show = False,
                                    mvis_show = False,
                                    rtip_show = True)
               window.toplevel.grab_set()
           

    def on_hertz_activated(self):
        self.on_activate_all_checks()
        if self.hertz_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Hertz",
                                            image_path="hertz.png",                                            
                                            saym_show = True,                                            
                                            musam_show = True,
                                            rtip_show = True,
                                            ymt_show = True)
            window.toplevel.grab_set()
            
    def on_punch_activated(self):
        self.on_activate_all_checks()
        if self.punch_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Flat cylinder",
                                            image_path="punch.png",
                                            saym_show = True,
                                            musam_show = True,
                                            rtip_show = True,
                                            ymt_show = True)
            window.toplevel.grab_set()
            
    def on_sneddon_activated(self): 
        
        self.on_activate_all_checks()
        if self.sneddon_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Sneddon cone",
                                            image_path="sneddon.png",                                            
                                            saym_show = True,                                            
                                            musam_show = True,
                                            angle_show = True,
                                            rtip_show = True,
                                            ymt_show = True)
            window.toplevel.grab_set()
            
    def on_nanowire_activated(self): 
        
        self.on_activate_all_checks()
        if self.nanowire_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Nanowire",
                                            image_path="benc.png",                                            
                                            saym_show = True,                                            
                                            musam_show = True,
                                            rtip_show = True,
                                            ymt_show = True)
            window.toplevel.grab_set()
            
            

    def on_dmt_activated(self):
        self.on_activate_all_checks()
        if self.dmt_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "DMT",
                                image_path="dmt.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = True,
                                ham2_show = False,
                                a00_show = True,
                                angle_show = False,
                                sahe_show = False,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True)
            window.toplevel.grab_set()

    def on_jkr_activated(self):
        self.on_activate_all_checks()
        if self.jkr_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "JKR",
                                image_path="jkr.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = True,
                                ham2_show = False,
                                a00_show = True,
                                angle_show = False,
                                sahe_show = False,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True)
            window.toplevel.grab_set()

    def on_becc_activated(self):
        self.on_activate_all_checks()
        if self.becc_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "bec Cone",
                                image_path="becc.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = False,
                                ham2_show = False,
                                a00_show = False,
                                angle_show = True,
                                sahe_show = True,
                                sample_bonded_show = True,
                                mvis_show = False,
                                rtip_show = False,
                                ymt_show = True,
                                text_width = 40)
            window.toplevel.grab_set()

    def on_behc_activated(self):
        self.on_activate_all_checks()
        if self.behc_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "bec Sphere",
                                image_path="behc.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = False,
                                ham2_show = False,
                                a00_show = False,
                                angle_show = False,
                                sahe_show = True,
                                sample_bonded_show = True,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True,
                                text_width = 40)
            window.toplevel.grab_set()

    def on_bepc_activated(self):
        self.on_activate_all_checks()
        if self.bepc_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "bec Flat cylinder",
                                image_path="bepc.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = False,
                                ham2_show = False,
                                a00_show = False,
                                angle_show = False,
                                sahe_show = True,
                                sample_bonded_show = True,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True,
                                text_width = 40)
            window.toplevel.grab_set()
            
    def on_benc_activated(self):
        self.on_activate_all_checks()
        if self.benc_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "bec Nanowire",
                                image_path="benc.png",
                                sr_show = False,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = False,
                                ham2_show = False,
                                a00_show = False,
                                angle_show = False,
                                sahe_show = True,
                                sample_bonded_show = True,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True,
                                text_width = 40)
            window.toplevel.grab_set()
            
    def on_tatara_activated(self):
        self.on_activate_all_checks()
        if self.tatara_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Tatara",
                                image_path="tatara.png",
                                sr_show = True,
                                saym_show = True,
                                saym2_show = False,
                                musam_show = True,
                                epsilon_show = False,
                                sigmas_show = False,
                                sigmat_show = False,
                                landadeb_show = False,
                                ham_show = False,
                                ham2_show = False,
                                a00_show = False,
                                angle_show = False,
                                sahe_show = False,
                                mvis_show = False,
                                rtip_show = True,
                                ymt_show = True)
            window.toplevel.grab_set()
            
    def on_DELR_activated(self):
        self.on_activate_all_checks()
        if self.DELR_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Two layered material", # "Doss-Eliato-Lin-Ros",
                                image_path="Doss-Eliato-Lin-Ros.png",      
                                saym_show = True,
                                saym2_show = True,
                                musam_show = True,
                                sahe_show = True,                                                            
                                rtip_show = True,
                                ymt_show = False,
                                
                                text_width = 40)
            window.toplevel.grab_set()
    def on_viscosity_activated(self):
        self.on_activate_all_checks()
        if self.viscosity_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Viscosity",
                                image_path="voigt_sphere.png",
                                mvis_show = True,
                                rtip_show = True,
                                saym_show = True, 
                                musam_show = True,                                                                                                
                                ymt_show = True)
            window.toplevel.grab_set()
    
    def on_viscosity_cone_activated(self):
        self.on_activate_all_checks()
        if self.viscosity_cone_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Viscosity",
                                image_path="voigt_cone.png",
                                mvis_show = True,
                                angle_show = True,
                                musam_show = True,                                
                                ymt_show = True,
                                saym_show = True                                                                
                                )
            window.toplevel.grab_set()
            
    def on_viscosity_punch_activated(self):
        self.on_activate_all_checks()
        if self.viscosity_punch_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Viscosity",
                                image_path="voigt_punch.png",
                                mvis_show = True,
                                rtip_show = True,
                                musam_show = True,                                
                                ymt_show = True,
                                saym_show = True,                                            
                                )
            window.toplevel.grab_set()

    def on_adhesion_activated(self):        
        self.on_activate_all_checks()
        if self.adhy_check_var.get(): 
            window = ForceParametersMenuApp(self,
                        title = "Adhesion hysteresis",
                        image_path="hyster.png",
                        sr_show = False,
                        saym_show = False,
                        saym2_show = False,
                        musam_show = False,
                        epsilon_show = False,
                        sigmas_show = False,
                        sigmat_show = False,
                        landadeb_show = False,
                        ham_show = True,
                        ham2_show = True,
                        a00_show = False,
                        angle_show = False,
                        sahe_show = False,
                        mvis_show = False,
                        rtip_show = False,
                        ymt_show = False)
            window.toplevel.grab_set()   
 
            
    def on_dipole_activated(self):
        self.on_activate_all_checks()
        if self.magnetic_dipole_check_var.get():
            window = ForceParametersMenuApp(self,
                        title = "Magnetic Dipole-Dipole",
                        image_path="magnetic_dipole-dipole.png",
                        magnetic_dipole_show = True)
            window.toplevel.grab_set()

    def on_exp_activated(self):
        self.on_activate_all_checks()
        if self.magnetic_exponential_check_var.get():
            window = ForceParametersMenuApp(self,
                        title = "Magnetic Exponential Decay",
                        image_path="magnetic-exponential.png",
                        magnetic_exp_show = True)
            window.toplevel.grab_set()
    
    def on_lennardJones_activated(self):
        self.on_activate_all_checks()
        if self.lennardJones_check_var.get():
            window = ForceParametersMenuApp(self,
                        title = "Lennard–Jones potential",
                        image_path="Lennard_Jones.png",
                        lennardjones_show = True)
            window.toplevel.grab_set()

    def on_custom1_activated(self):
        self.on_activate_all_checks()
        if self.custom1_check_var.get():
            window = ForceParametersMenuApp(self,
                        title = "Customize force",
                        image_path=None,
                        custom_show = True)
            window.toplevel.grab_set()

    def on_custom2_activated(self):
        self.on_activate_all_checks()
        if self.custom2_check_var.get():
            window = ForceParametersMenuApp(self,
                        title = "Customize force",
                        image_path=None,
                        custom_show = True)
            window.toplevel.grab_set()
    

    # THIS ONES ARE viscosity AND Bottom Effect correction
    def on_kv_sphere_activated(self):
        self.on_activate_all_checks()
        if self.kv_sphere_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Kelvin-Voigt sphere with bec",
                                            image_path="kv_sphere_bec.png",
                                            saym_show = True,
                                            musam_show = True,
                                            sahe_show = True,
                                            mvis_show = True,
                                            rtip_show = True,  
                                            sample_bonded_show = True,
                                            text_width = 40)
            window.toplevel.grab_set()


    def on_kv_cone_activated(self):
        self.on_activate_all_checks()
        if self.kv_cone_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Kelvin-Voigt cone with bec",
                                            image_path="kv_cone_bec.png",                                            
                                            saym_show = True,
                                            musam_show = True,
                                            sahe_show = True,
                                            mvis_show = True,
                                            angle_show = True,
                                            sample_bonded_show = True,                                            
                                            text_width = 40)
            window.toplevel.grab_set()
        
    def on_kv_punch_activated(self):
        self.on_activate_all_checks()
        if self.kv_punch_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Kelvin-Voigt flat cylinder with bec",
                                            image_path="kv_punch_bec.png",
                                            saym_show = True,
                                            musam_show = True,
                                            sahe_show = True,
                                            mvis_show = True,
                                            rtip_show = True,
                                            sample_bonded_show = True,
                                            text_width = 40)
            window.toplevel.grab_set()
            
    def on_viscosity_nanowire_activated(self):
        self.on_activate_all_checks()
        if self.viscosity_nanowire_check_var.get():
            window = ForceParametersMenuApp(self,
                                title = "Viscosity nanowire",
                                image_path="benc.png",
                                mvis_show = True,
                                rtip_show = True,
                                saym_show = True, 
                                musam_show = True,                                                                                                
                                ymt_show = True)
            window.toplevel.grab_set()
    
    def on_kv_nanowire_activated(self):
        self.on_activate_all_checks()
        if self.kv_nanowire_check_var.get():
            window = ForceParametersMenuApp(self,
                                            title = "Kelvin-Voigt nanowire with bec",
                                            image_path="benc.png",
                                            saym_show = True,
                                            musam_show = True,
                                            sahe_show = True,
                                            mvis_show = True,
                                            rtip_show = True,
                                            sample_bonded_show = True,                                            
                                            text_width = 40)
            window.toplevel.grab_set()

   
    
    def get_option_mask(self):
        self.option_mask = [            
                self.vdw_check_var.get(),
                self.dlvo_check_var.get(),
                self.magnetic_dipole_check_var.get(),
                self.magnetic_exponential_check_var.get(),
                
                self.hertz_check_var.get(),
                self.sneddon_check_var.get(),
                self.punch_check_var.get(),
                
                self.dmt_check_var.get(),
                self.jkr_check_var.get(),           
                self.tatara_check_var.get(),
                
                self.behc_check_var.get(),
                self.becc_check_var.get(),
                self.bepc_check_var.get(),
                
                self.viscosity_check_var.get(),
                self.viscosity_cone_check_var.get(),
                self.viscosity_punch_check_var.get(),
                
                self.kv_sphere_check_var.get(),
                self.kv_cone_check_var.get(),            
                self.kv_punch_check_var.get(),    
                self.lennardJones_check_var.get(),
                
                self.custom1_check_var.get(),
                self.custom2_check_var.get(),
                self.benc_check_var.get(),
                self.nanowire_check_var.get(),
                self.DELR_check_var.get(),
                self.viscosity_nanowire_check_var.get(),
                self.kv_nanowire_check_var.get(),
                self.vdw_cone_check_var.get(),
                self.vdw_punch_check_var.get()]
        return self.option_mask
    
if __name__ == "__main__":    
    app = dforce_APP()   
    app.run()
