# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 12:16:32 2022

@author: Victor
"""
    
#!/usr/bin/python3
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import font

class ForceParametersMenuApp:
    def __init__(self, master,title,image_path,
                 sr_show = False,
                 saym_show = False,
                 saym2_show = False,
                 musam_show = False,
                 epsilon_show = False,
                 sigmas_show = False,
                 sigmat_show = False,
                 landadeb_show = False,
                 ham_show = False,
                 ham2_show = False,
                 a00_show = False,
                 angle_show = False,
                 sahe_show = False,
                 sample_bonded_show = False,
                 mvis_show = False,
                 custom_show = False,
                 rtip_show = False,
                 ymt_show = False,
                 magnetic_exp_show = False,
                 magnetic_dipole_show = False,
                 lennardjones_show = False,
                 lock_poisson_coeffient = False,
                 text_width = 30):
        
        self.sr_show = sr_show
        self.saym_show = saym_show
        self.saym2_show = saym2_show
        self.musam_show = musam_show
        self.epsilon_show = epsilon_show
        self.sigmas_show = sigmas_show
        self.sigmat_show = sigmat_show
        self.landadeb_show = landadeb_show
        self.ham_show = ham_show
        self.ham2_show = ham2_show
        self.a00_show = a00_show
        self.angle_show = angle_show
        self.sahe_show = sahe_show
        self.sample_bonded_show = sample_bonded_show
        self.mvis_show = mvis_show
        self.custom_show = custom_show
        self.rtip_show = rtip_show
        self.ymt_show = ymt_show
        
        
        self.magnetic_exp_show = magnetic_exp_show
        self.magnetic_dipole_show = magnetic_dipole_show
        self.lennardjones_show = lennardjones_show

        # build ui
        self.master = master
        self.toplevel = tk.Toplevel(self.master.toplevel1)
        self.toplevel.wm_title(title + " force")        
        self.frame = ttk.Frame(self.toplevel)        
        
        
        
    
        # Label for the name of the force
        self.label_title = ttk.Label(self.frame)
        self.label_title.configure(text=title)
        self.label_title.grid(column=0, row=0)
        
        
        # Formula of the force is displayed along with the text
        self.img = tk.PhotoImage(file="icons/icon.png")
        self.toplevel.wm_iconphoto(False, self.img)        
                
        self.label_formula = ttk.Label(self.frame)
        if image_path =="value table":
            self.label_formula.configure(text="Value Table")
        elif not image_path == None:
            self.img_hyster = tk.PhotoImage(file="icons/"+image_path)
            self.label_formula.configure(image=self.img_hyster, text="label_image")
        else:
            # When user selects custom force we display a text as there is not an image
            self.label_formula.configure(text="To configure your own force, open and edit the file 'force.py'\n Within the file you shall find instructions about how to do it")
        
        self.label_formula.grid(column=1, row=0,sticky = "we")
        
        
        
        #======================CONFIGURATION OPTIONS FOR LABELS================
        text_entry_conf = {"width":text_width,"height":1.5,
                                    "borderwidth":0,
                                    "font":"TkDefaultFont",
                                    "background" : self.toplevel.cget("background")}
        
 
        
        subscript_conf = {"tagName":"subscript", "offset":-4,"font":("TkDefaultFont",7) }
        superscript_conf = {"tagName":"superscript", "offset":+4,"font":("TkDefaultFont",7) }
        
     
        
        #============================ENTRY LABELS==============================        
        
        if self.sr_show:
            #self.label_sr_sp = ttk.Label(self.frame)
            #self.label_sr_sp.configure(text="R_sample (nm)")

            self.label_sr_sp = tk.Text(self.frame,**text_entry_conf)
            
            self.label_sr_sp.tag_configure(**subscript_conf)
            self.label_sr_sp.insert("insert", "R", "", "sample", "subscript")
            self.label_sr_sp.configure(state="disabled")
            self.label_sr_sp.grid(column=0, row=1)            
            self.label_sr_sp.bindtags((str(self.label_sr_sp), str(master), "all"))
            
            self.sr_sp = ttk.Spinbox(self.frame)
            self.sr_sp.configure(from_=0, increment=0.01, justify="right", to=1E9)
            _text_ = str(self.master.toplevel1.sr_sp_var)
            self.sr_sp.delete("0", "end")
            self.sr_sp.insert("0", _text_)
            self.sr_sp.grid(column=1, row=1, sticky="e")
            
            
        if self.saym_show:
            #self.label_saym_sp = ttk.Label(self.frame)
            #self.label_saym_sp.configure(text="E_sample (MPa)")
            
            self.label_saym_sp = tk.Text(self.frame,**text_entry_conf)
            
            
            self.label_saym_sp.tag_configure(**subscript_conf)
            self.label_saym_sp.insert("insert", "E", "", "sample", "subscript"," (MPa)","")
            self.label_saym_sp.configure(state="disabled")
            
            self.label_saym_sp.grid(column=0, row=2)
            self.saym_sp = ttk.Spinbox(self.frame)
            self.saym_sp.configure(from_=0, increment=0.01, justify="right", to=1E9)
            _text_ = str(self.master.toplevel1.saym_sp_var)
            self.saym_sp.delete("0", "end")
            self.saym_sp.insert("0", _text_)
            self.saym_sp.grid(column=1, row=2, sticky="e")
            
        if self.saym2_show:
            #self.label_saym2_sp = ttk.Label(self.frame)
            #self.label_saym2_sp.configure(text="E_{2 sample} (MPa)")
            
            self.label_saym2_sp = tk.Text(self.frame,**text_entry_conf)
            
            self.label_saym2_sp.tag_configure(**subscript_conf)
            self.label_saym2_sp.insert("insert", "E", "", "2 sample", "subscript"," (MPa)","")
            self.label_saym2_sp.configure(state="disabled")
            
            
            self.label_saym2_sp.grid(column=0, row=3)
            self.saym2_sp = ttk.Spinbox(self.frame)
            self.saym2_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.saym2_sp_var)
            self.saym2_sp.delete("0", "end")
            self.saym2_sp.insert("0", _text_)
            self.saym2_sp.grid(column=1, row=3, sticky="e")
        
        if self.musam_show:
            
            #self.label_musam_sp = ttk.Label(self.frame)
            #self.label_musam_sp.configure(text="v_sample (adimensional)")
             
            self.label_musam_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_musam_sp.tag_configure(**subscript_conf)
            self.label_musam_sp.insert("insert", "ν", "", "sample", "subscript"," (adimensional)","")
            self.label_musam_sp.configure(state="disabled")
             
            self.label_musam_sp.grid(column=0, row=4)
            self.musam_sp = ttk.Spinbox(self.frame)
            self.musam_sp.configure(from_=0, increment=0.01, justify="right", to=0.5)
             
             
            if lock_poisson_coeffient:
                _text_ =  str(0.5)
                self.musam_sp.delete("0", "end")
                self.musam_sp.insert("0", _text_)
                self.musam_sp.grid(column=1, row=4, sticky="e")
                self.musam_sp.config(state=tk.DISABLED)
            else:
                _text_ =  str(self.master.toplevel1.musam_sp_var)
                self.musam_sp.delete("0", "end")
                self.musam_sp.insert("0", _text_)
                self.musam_sp.grid(column=1, row=4, sticky="e")
        
        if self.epsilon_show:
            #self.label_epsilon_sp = ttk.Label(self.frame)
            #self.label_epsilon_sp.configure(text="epsilon_r (adimensional)")
            
            self.label_epsilon_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_epsilon_sp.tag_configure(**subscript_conf)
            self.label_epsilon_sp.insert("insert", "ϵ", "", "r", "subscript"," (adimensional)","")
            self.label_epsilon_sp.configure(state="disabled")
            
            
            self.label_epsilon_sp.grid(column=0, row=7)
            self.epsilon_sp = ttk.Spinbox(self.frame)
            self.epsilon_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.epsilon_sp_var)
            self.epsilon_sp.delete("0", "end")
            self.epsilon_sp.insert("0", _text_)
            self.epsilon_sp.grid(column=1, row=7, sticky="e")
        
        if self.sigmas_show:
            #self.label_sigmas_sp = ttk.Label(self.frame)
            #self.label_sigmas_sp.configure(text="sigma_sample (mC/m^2)")
            
            self.label_sigmas_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_sigmas_sp.tag_configure(**subscript_conf)
            self.label_sigmas_sp.tag_configure(**superscript_conf)
            self.label_sigmas_sp.insert("insert", "σ", "", "sample", "subscript"," (mC/m","","2","superscript",")","")
            self.label_sigmas_sp.configure(state="disabled")
            
            self.label_sigmas_sp.grid(column=0, row=8)
            self.sigmas_sp = ttk.Spinbox(self.frame)
            self.sigmas_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.sigmas_sp_var)
            self.sigmas_sp.delete("0", "end")
            self.sigmas_sp.insert("0", _text_)
            self.sigmas_sp.grid(column=1, row=8, sticky="e")
            
        if self.sigmat_show:
            #self.label_sigmat_sp = ttk.Label(self.frame)
            #self.label_sigmat_sp.configure(text="sigma_tip (mC/m^2)")
            self.label_sigmat_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_sigmat_sp.tag_configure(**subscript_conf)
            self.label_sigmat_sp.tag_configure(**superscript_conf)
            self.label_sigmat_sp.insert("insert", "σ", "", "tip", "subscript"," (mC/m","","2","superscript",")","")
            self.label_sigmat_sp.configure(state="disabled")
            
            self.label_sigmat_sp.grid(column=0, row=9)
            self.sigmat_sp = ttk.Spinbox(self.frame)
            self.sigmat_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.sigmat_sp_var)
            self.sigmat_sp.delete("0", "end")
            self.sigmat_sp.insert("0", _text_)
            self.sigmat_sp.grid(column=1, row=9, sticky="e")
        
        if self.landadeb_show:
            #self.label_landadeb_sp = ttk.Label(self.frame)
            #self.label_landadeb_sp.configure(text="1/kappa_D = \\labmda_D (nm)")
            
            self.label_landadeb_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_landadeb_sp.tag_configure(**subscript_conf)
            self.label_landadeb_sp.tag_configure(**superscript_conf)
            self.label_landadeb_sp.insert("insert", "1/κ", "", "D", "subscript","=λ","","D","subscript"," (nm)","")
            self.label_landadeb_sp.configure(state="disabled")
            
            self.label_landadeb_sp.grid(column=0, row=10)
            self.landadeb_sp = ttk.Spinbox(self.frame)
            self.landadeb_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.landadeb_sp_var)
            self.landadeb_sp.delete("0", "end")
            self.landadeb_sp.insert("0", _text_)
            self.landadeb_sp.grid(column=1, row=10, sticky="e")
            
        if self.ham_show:
            #self.label_ham_sp = ttk.Label(self.frame)
            #self.label_ham_sp.configure(text="Hamaker constant (10^-20 J)")
            
            self.label_ham_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_ham_sp.tag_configure(**subscript_conf)
            self.label_ham_sp.tag_configure(**superscript_conf)
            self.label_ham_sp.insert("insert", "Hamaker Constant (10", "", "-20", "superscript"," J)","")
            self.label_ham_sp.configure(state="disabled")
            
            
            self.label_ham_sp.grid(column=0, row=13)
            self.ham_sp = ttk.Spinbox(self.frame)
            self.ham_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.ham_sp_var)
            self.ham_sp.delete("0", "end")
            self.ham_sp.insert("0", _text_)
            self.ham_sp.grid(column=1, row=13, sticky="e")
            
        if self.ham2_show:
            #self.label_ham2_sp = ttk.Label(self.frame)
            #self.label_ham2_sp.configure(text="Retractive Hamaker (10^-20 J)")
            
            self.label_ham2_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_ham2_sp.tag_configure(**subscript_conf)
            self.label_ham2_sp.tag_configure(**superscript_conf)
            self.label_ham2_sp.insert("insert", "Retractive Hamaker (10", "", "-20", "superscript"," J)","")
            self.label_ham2_sp.configure(state="disabled")
            
            
            self.label_ham2_sp.grid(column=0, row=14)
            self.ham2_sp = ttk.Spinbox(self.frame)
            self.ham2_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.ham2_sp_var)
            self.ham2_sp.delete("0", "end")
            self.ham2_sp.insert("0", _text_)
            self.ham2_sp.grid(column=1, row=14, sticky="e")
        
        if self.a00_show:
            #self.label_a00_sp = ttk.Label(self.frame)
            #self.label_a00_sp.configure(text="a_0 (nm)")
            
            self.label_a00_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_a00_sp.tag_configure(**subscript_conf)
            self.label_a00_sp.tag_configure(**superscript_conf)
            self.label_a00_sp.insert("insert", "a", "", "0", "subscript"," (nm)","")
            self.label_a00_sp.configure(state="disabled")
            
            
            self.label_a00_sp.grid(column=0, row=15)
            self.a00_sp = ttk.Spinbox(self.frame)
            self.a00_sp.configure(
                from_=0, increment=0.001, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.a00_sp_var)
            self.a00_sp.delete("0", "end")
            self.a00_sp.insert("0", _text_)
            self.a00_sp.grid(column=1, row=15, sticky="e")
        
        if self.angle_show:
            #self.label_angle_sp = ttk.Label(self.frame)
            #self.label_angle_sp.configure(text="Tip cone angle (degress)")
            
            self.label_angle_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_angle_sp.tag_configure(**subscript_conf)
            self.label_angle_sp.tag_configure(**superscript_conf)
            self.label_angle_sp.insert("insert","Tip cone angle (degrees)","")
            self.label_angle_sp.configure(state="disabled")
            
            
            self.label_angle_sp.grid(column=0, row=16)
            self.angle_sp = ttk.Spinbox(self.frame)
            self.angle_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.angle_sp_var)
            self.angle_sp.delete("0", "end")
            self.angle_sp.insert("0", _text_)
            self.angle_sp.grid(column=1, row=16, sticky="e")
        
        if self.sahe_show:
            #self.label_sahe_sp = ttk.Label(self.frame)
            #self.label_sahe_sp.configure(text="Sample height (nm)")
            
            self.label_sahe_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_sahe_sp.tag_configure(**subscript_conf)
            self.label_sahe_sp.tag_configure(**superscript_conf)
            self.label_sahe_sp.insert("insert", "Layer thickness (nm)", "")
            self.label_sahe_sp.configure(state="disabled")            
            self.label_sahe_sp.grid(column=0, row=17)
            
            
            self.sahe_sp = ttk.Spinbox(self.frame)
            self.sahe_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.sahe_sp_var)
            self.sahe_sp.delete("0", "end")
            self.sahe_sp.insert("0", _text_)
            self.sahe_sp.grid(column=1, row=17, sticky="e")
            
            # Bonded, yes/no
            #self.label_bonded_sp = ttk.Label(self.frame)
            #self.label_bonded_sp.configure(text="Sample bonded")
        if self.sample_bonded_show:
            
            self.frame2 = ttk.Frame(self.frame) 
            self.label_bonded_sp = tk.Text(self.frame,**text_entry_conf)                        
            self.label_bonded_sp.insert("insert", "Boundary conditions: sample bonded")
            self.label_bonded_sp.configure(state="disabled")            
            
            self.label_bonded_sp.grid(column=0, row=18)
            
            self.bonded_sp_value = tk.BooleanVar(self.frame,value = self.master.toplevel1.bec_bonded_var)
            """
            self.bonded_sp = ttk.Checkbutton(self.frame,
                    text="",
                variable=self.bonded_sp_value)
            
            self.bonded_sp.grid(column=1, row=18, sticky="e")
            """
            self.bonded_sp1 = ttk.Radiobutton(self.frame2,text="Bonded",value=True,variable=self.bonded_sp_value)
            self.bonded_sp1.grid(column=0, row=0)
            self.bonded_sp1 = ttk.Radiobutton(self.frame2,text="Frictionless",value=False,variable=self.bonded_sp_value)
            self.bonded_sp1.grid(column=1, row=0,sticky="e")
            
            self.frame2.grid(column=1,row = 18,sticky = "e")
            
            
            
            

        
            
        if self.mvis_show:
            #self.label_mvis_sp = ttk.Label(self.frame)
            #self.label_mvis_sp.configure(text="eta (Pa \\dot s)")
            
            self.label_mvis_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_mvis_sp.tag_configure(**subscript_conf)
            self.label_mvis_sp.tag_configure(**superscript_conf)
            self.label_mvis_sp.insert("insert", "η", "","e","subscript"," (Pa•s)","")
            self.label_mvis_sp.configure(state="disabled")
            
            
            self.label_mvis_sp.grid(column=0, row=20)
            self.mvis_sp = ttk.Spinbox(self.frame)
            self.mvis_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            self.mvis_sp.configure(validate="none")
            _text_ =  str(self.master.toplevel1.mvis_sp_var)
            self.mvis_sp.delete("0", "end")
            self.mvis_sp.insert("0", _text_)
            self.mvis_sp.grid(column=1, row=20, sticky="e")
            
        if self.custom_show:
            #============================CUSTOM 1==============================
            #self.label_custom1_sp = ttk.Label(self.frame)
            #self.label_custom1_sp.configure(text="V_1")
            self.label_custom1_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_custom1_sp.tag_configure(**subscript_conf)
            self.label_custom1_sp.tag_configure(**superscript_conf)
            self.label_custom1_sp.insert("insert", "V", "","1","subscript")
            self.label_custom1_sp.configure(state="disabled")
            
            self.label_custom1_sp.grid(column=0, row=23)
            self.custom1_sp = ttk.Spinbox(self.frame)
            self.custom1_sp.configure(
                from_=-1E9, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.custom1_sp_var)
            self.custom1_sp.delete("0", "end")
            self.custom1_sp.insert("0", _text_)
            self.custom1_sp.grid(column=1, row=23, sticky="e")
            
            #=============================CUSTOM 2"============================
            #self.label_custom2_sp = ttk.Label(self.frame)
            #self.label_custom2_sp.configure(text="V_2")
            self.label_custom2_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_custom2_sp.tag_configure(**subscript_conf)
            self.label_custom2_sp.tag_configure(**superscript_conf)
            self.label_custom2_sp.insert("insert", "V", "","2","subscript")
            self.label_custom2_sp.configure(state="disabled")
            self.label_custom2_sp.grid(column=0, row=24)
            self.custom2_sp = ttk.Spinbox(self.frame)
            self.custom2_sp.configure(
                from_=-1E9, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.custom2_sp_var)
            self.custom2_sp.delete("0", "end")
            self.custom2_sp.insert("0", _text_)
            self.custom2_sp.grid(column=1, row=24, sticky="e")
                        
            #=============================CUSTOM 2"============================
            #self.label_custom3_sp = ttk.Label(self.frame)
            #self.label_custom3_sp.configure(text="V_3")
            self.label_custom3_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_custom3_sp.tag_configure(**subscript_conf)
            self.label_custom3_sp.tag_configure(**superscript_conf)
            self.label_custom3_sp.insert("insert", "V", "","3","subscript")
            self.label_custom3_sp.configure(state="disabled")
            
            self.label_custom3_sp.grid(column=0, row=25)
            self.custom3_sp = ttk.Spinbox(self.frame)
            self.custom3_sp.configure(
                from_=-1E9, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.custom3_sp_var)
            self.custom3_sp.delete("0", "end")
            self.custom3_sp.insert("0", _text_)
            self.custom3_sp.grid(column=1, row=25, sticky="e")
            
            #=============================CUSTOM 4"============================
            #self.label_custom4_sp = ttk.Label(self.frame)
            #self.label_custom4_sp.configure(text="V_4")
            self.label_custom4_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_custom4_sp.tag_configure(**subscript_conf)
            self.label_custom4_sp.tag_configure(**superscript_conf)
            self.label_custom4_sp.insert("insert", "V", "","4","subscript")
            self.label_custom4_sp.configure(state="disabled")
            
            self.label_custom4_sp.grid(column=0, row=26)
            self.custom4_sp = ttk.Spinbox(self.frame)
            self.custom4_sp.configure(
                from_=-1E9, increment=0.01, justify="right", to=1E9
            )
            self.custom4_sp.configure(validate="focusout")
            _text_ =  str(self.master.toplevel1.custom4_sp_var)
            self.custom4_sp.delete("0", "end")
            self.custom4_sp.insert("0", _text_)
            self.custom4_sp.grid(column=1, row=26, sticky="e")
            
        if self.rtip_show:
            #self.label1 = ttk.Label(self.frame)
            #self.label1.configure(text="Radius tip (nm)")
       
            self.label1 = tk.Text(self.frame,**text_entry_conf)            
            self.label1.tag_configure(**subscript_conf)
            self.label1.tag_configure(**superscript_conf)
            self.label1.insert("insert", "Tip Radius (nm)","")
            
            self.label1.configure(state="disabled")
            
            self.label1.grid(column=0, row=27)
            self.rtip_sp = ttk.Spinbox(self.frame)
            self.rtip_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.rtip_sp.configure(validate="focusout")            
            _text_ = str(self.master.tr_sp.get())            
            self.rtip_sp.delete("0", "end")
            self.rtip_sp.insert("0", _text_)
            self.rtip_sp.grid(column=1, row=27, sticky="e")
        
        if self.ymt_show:
            #self.label2 = ttk.Label(self.frame)
            #self.label2.configure(text="Tip Young's modulus")
            
            self.label2 = tk.Text(self.frame,**text_entry_conf)            
            self.label2.tag_configure(**subscript_conf)
            self.label2.tag_configure(**superscript_conf)
            self.label2.insert("insert", "Tip Young's Modulus (GPa)","")
            
            self.label2.grid(column=0, row=28)
            self.ymt_sp = ttk.Spinbox(self.frame)
            self.ymt_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.ymt_sp.configure(validate="focusout")
            _text_ = str(self.master.ymt_sp.get())            
            self.ymt_sp.delete("0", "end")
            self.ymt_sp.insert("0", _text_)
            self.ymt_sp.grid(column=1, row=28, sticky="e")
            
            #self.label3 = ttk.Label(self.frame)
            #self.label3.configure(text="Tip Poisson coefficient")
            self.label3 = tk.Text(self.frame,**text_entry_conf)            
            self.label3.tag_configure(**subscript_conf)
            self.label3.tag_configure(**superscript_conf)
            self.label3.insert("insert", "Tip Poisson coefficient","")
            
            self.label3.grid(column=0, row=29)
            
            self.mutip_sp = ttk.Spinbox(self.frame)
            self.mutip_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.mutip_sp.configure(validate="focusout")
            _text_ = str(self.master.mutip_sp.get())            
            self.mutip_sp.delete("0", "end")
            self.mutip_sp.insert("0", _text_)
            self.mutip_sp.grid(column=1, row=29, sticky="e")
        
        if self.magnetic_exp_show or self.magnetic_dipole_show:
            #=====================MAGNETIC tip dipole=======================
            #self.label4 = ttk.Label(self.frame)
            #self.label4.configure(text="m_tip (10^−14 emu)")
            
            self.label4 = tk.Text(self.frame,**text_entry_conf)            
            self.label4.tag_configure(**subscript_conf)
            self.label4.tag_configure(**superscript_conf)
            self.label4.insert("insert", "m","","tip","subscript"," (10","","−14","superscript"," emu)","")
            
            self.label4.grid(column=0, row=30)
            
            self.magnetic_m_tip_sp = ttk.Spinbox(self.frame)
            self.magnetic_m_tip_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.magnetic_m_tip_sp.configure(validate="focusout")
            _text_ =  str(self.master.toplevel1.magnetic_m_tip)
            self.magnetic_m_tip_sp.delete("0", "end")
            self.magnetic_m_tip_sp.insert("0", _text_)
            self.magnetic_m_tip_sp.grid(column=1, row=30, sticky="e")
        
        if self.magnetic_dipole_show:
            
            #=====================MAGNETIC sample dipole=======================
            #self.label5 = ttk.Label(self.frame)
            #self.label5.configure(text="m_sample (10^−14 emu)")
            self.label5 = tk.Text(self.frame,**text_entry_conf)            
            self.label5.tag_configure(**subscript_conf)
            self.label5.tag_configure(**superscript_conf)
            self.label5.insert("insert", "m","","sample","subscript"," (10","","−14","superscript"," emu)","")
            
            self.label5.grid(column=0, row=31)
            
            self.magnetic_m_sample_sp = ttk.Spinbox(self.frame)
            self.magnetic_m_sample_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.magnetic_m_sample_sp.configure(validate="focusout")
            _text_ =  str(self.master.toplevel1.magnetic_m_sample)
            self.magnetic_m_sample_sp.delete("0", "end")
            self.magnetic_m_sample_sp.insert("0", _text_)
            self.magnetic_m_sample_sp.grid(column=1, row=31, sticky="e")
            
        if self.magnetic_exp_show:
            
            #=====================MAGNETIC B0    ==============================
            #self.label6 = ttk.Label(self.frame)
            #self.label6.configure(text="B_0 sample (mT)")
            self.label6 = tk.Text(self.frame,**text_entry_conf)            
            self.label6.tag_configure(**subscript_conf)
            self.label6.tag_configure(**superscript_conf)
            self.label6.insert("insert", "B","","0","subscript"," sample (mT)","")
            
            self.label6.grid(column=0, row=32)
            
            self.magnetic_B0_sample_sp = ttk.Spinbox(self.frame)
            self.magnetic_B0_sample_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.magnetic_B0_sample_sp.configure(validate="focusout")
            _text_ =  str(self.master.toplevel1.magnetic_B0_sample)
            self.magnetic_B0_sample_sp.delete("0", "end")
            self.magnetic_B0_sample_sp.insert("0", _text_)
            self.magnetic_B0_sample_sp.grid(column=1, row=32, sticky="e")
            
            #=====================MAGNETIC KAPPA==============================
            #self.label7 = ttk.Label(self.frame)
            #self.label7.configure(text="kappa (nm^-1)")
            self.label7 = tk.Text(self.frame,**text_entry_conf)            
            self.label7.tag_configure(**subscript_conf)
            self.label7.tag_configure(**superscript_conf)
            self.label7.insert("insert", "κ",""," (nm","","-1","superscript",")","")
            
            self.label7.grid(column=0, row=33)
            
            self.magnetic_k_sample_sp = ttk.Spinbox(self.frame)
            self.magnetic_k_sample_sp.configure(
                from_=-10000000000000, increment=0.01, justify="right", to=10000000000000
            )
            self.magnetic_k_sample_sp.configure(validate="focusout")
            _text_ =  str(self.master.toplevel1.magnetic_k_sample)
            self.magnetic_k_sample_sp.delete("0", "end")
            self.magnetic_k_sample_sp.insert("0", _text_)
            self.magnetic_k_sample_sp.grid(column=1, row=33, sticky="e")
  
        if self.lennardjones_show:
            #=====================Lennard-Jones======================= part 1
            self.label_sigma_lj_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_sigma_lj_sp.tag_configure(**subscript_conf)
            self.label_sigma_lj_sp.tag_configure(**superscript_conf)
            self.label_sigma_lj_sp.insert("insert", "σ (nm)", "")
            self.label_sigma_lj_sp.configure(state="disabled")            
            self.label_sigma_lj_sp.grid(column=0, row=34)
            
            
            self.sigma_lj_sp = ttk.Spinbox(self.frame)
            self.sigma_lj_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.sigma_lj_var)
        
            self.sigma_lj_sp.delete("0", "end")
            self.sigma_lj_sp.insert("0", _text_)
            self.sigma_lj_sp.grid(column=1, row=34, sticky="e")
            
            #=====================Lennard-Jones======================= part 2
            self.label_epsilon_lj_sp = tk.Text(self.frame,**text_entry_conf)            
            self.label_epsilon_lj_sp.tag_configure(**subscript_conf)
            self.label_epsilon_lj_sp.tag_configure(**superscript_conf)
            self.label_epsilon_lj_sp.insert("insert", "ϵ  (nN⋅nm)","")
            self.label_epsilon_lj_sp.configure(state="disabled")            
            self.label_epsilon_lj_sp.grid(column=0, row=35)
            
            
            self.epsilon_lj_sp = ttk.Spinbox(self.frame)
            self.epsilon_lj_sp.configure(
                from_=0, increment=0.01, justify="right", to=1E9
            )
            _text_ =  str(self.master.toplevel1.epsilon_lj_var)
            self.epsilon_lj_sp.delete("0", "end")
            self.epsilon_lj_sp.insert("0", _text_)
            self.epsilon_lj_sp.grid(column=1, row=35, sticky="e")
           
            
            
            
            
        # Other frame configurations and button application
        self.frame.configure(height=200, width=200)
        self.frame.pack(side="top")
        #self.frame.rowconfigure("all", pad=10, weight=1)
        self.frame.columnconfigure("all", pad=10, weight=1)
        self.button10 = ttk.Button(self.toplevel)
        self.button10.configure(text="OK")
        self.button10.pack(side="top")
        self.button10.configure(command=self.value_table_exit)
        # If in focus when pressing enter the button will work
        self.button10.bind("<Return>", ((lambda event: self.value_table_exit())))
        # Focus on the button when opening the window to make sure exiting is fast
        self.button10.focus()
        self.toplevel.configure(height=200, width=200)

        # Main widget
        self.mainwindow = self.toplevel

    def run(self):
        # NEVER USED
        self.mainwindow.mainloop()

    def value_table_exit(self):
        if self.sr_show:            
            self.master.toplevel1.sr_sp_var = float( self.sr_sp.get() )
        if self.saym_show:            
            self.master.toplevel1.saym_sp_var = float( self.saym_sp.get() )
        if self.saym2_show:            
            self.master.toplevel1.saym2_sp_var = float( self.saym2_sp.get() )
        if self.musam_show:            
            self.master.toplevel1.musam_sp_var = float( self.musam_sp.get() )
        if self.epsilon_show:            
            self.master.toplevel1.epsilon_sp_var = float( self.epsilon_sp.get() )
        if self.sigmas_show:            
            self.master.toplevel1.sigmas_sp_var = float( self.sigmas_sp.get() )
        if self.sigmat_show:            
            self.master.toplevel1.sigmat_sp_var = float( self.sigmat_sp.get() )
        if self.landadeb_show:            
            self.master.toplevel1.landadeb_sp_var = float( self.landadeb_sp.get() )
        if self.ham_show:            
            self.master.toplevel1.ham_sp_var = float( self.ham_sp.get() )
        if self.ham2_show:            
            self.master.toplevel1.ham2_sp_var = float( self.ham2_sp.get() )
        if self.a00_show:            
            self.master.toplevel1.a00_sp_var = float( self.a00_sp.get() )
        if self.angle_show:            
            self.master.toplevel1.angle_sp_var = float( self.angle_sp.get() )
        if self.sahe_show:            
            self.master.toplevel1.sahe_sp_var = float( self.sahe_sp.get() )
        if self.sample_bonded_show:
            self.master.toplevel1.bec_bonded_var = self.bonded_sp_value.get()            
        if self.mvis_show:            
            self.master.toplevel1.mvis_sp_var = float( self.mvis_sp.get() )  
            
        if self.magnetic_dipole_show:
            self.master.toplevel1.magnetic_m_tip = float(self.magnetic_m_tip_sp.get())
            self.master.toplevel1.magnetic_m_sample = float(self.magnetic_m_sample_sp.get())
        if self.magnetic_exp_show:
            self.master.toplevel1.magnetic_m_tip = float(self.magnetic_m_tip_sp.get())
            self.master.toplevel1.magnetic_B0_sample = float(self.magnetic_B0_sample_sp.get())
            self.master.toplevel1.magnetic_k_sample = float(self.magnetic_k_sample_sp.get())
            
        if self.custom_show:            
            self.master.toplevel1.custom1_sp_var = float( self.custom1_sp.get() )        
            self.master.toplevel1.custom2_sp_var = float( self.custom2_sp.get() )        
            self.master.toplevel1.custom3_sp_var = float( self.custom3_sp.get() )        
            self.master.toplevel1.custom4_sp_var = float( self.custom4_sp.get() )
            
        if self.rtip_show:            
            self.master.tr_sp.set(self.rtip_sp.get())
        if self.ymt_show:
            self.master.ymt_sp.set(self.ymt_sp.get())
            self.master.mutip_sp.set(self.mutip_sp.get())
            
        if self.lennardjones_show:
            self.master.toplevel1.sigma_lj_var = float(self.sigma_lj_sp.get())
            self.master.toplevel1.epsilon_lj_var = float(self.epsilon_lj_sp.get())
            
        self.toplevel.destroy()



