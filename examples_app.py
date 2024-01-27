
import tkinter as tk
import tkinter.ttk as ttk
import webbrowser

text1 = """Monomodal
A_{01} = 15 mn
E_{sample} = 100
R_{tip} = 5 nm
k_{tip} = 0.9 N/m
Q_{tip} = 200
f_{tip} = 100 kHz"""

text2 = """Monomodal
A_{01} = 25 mn
E_{sample} = 1 GPa
Hamaker = 300 zJ 
R_{tip} = 5 nm
k_{tip} = 0.9 N/m
Q_{tip} = 200
f_{tip} = 100 kHz
Including retraction """

text3 = """Bimodal 
A_{01} = 100 nm
A_{02} = 0.25 nm
E_{sample} = 100
R_{tip} = 5 nm
k_{1} = 3.0 N/m
Q_{1} = 220
f_{1} = 82 kHz
k_{2} = 152 N/m
f_{2} = 518 kHz
Q_{2} = 1000"""

text4 = """Monomodal
A_{01} = 10 mn
E_{sample} = 100 MPa
η_{sample} = 2 Pa•s
h_{sample}= 10 nm
R_{tip} = 5 nm
k_{tip} = 0.9 N/m
Q_{tip} = 200
f_{tip} = 100 kHz"""


text5 = """Monomodal
A_{01} = 41 mn
E_{sample} = 100
R_{tip} = 50 nm
k_{tip} = 0.9 N/m
Q_{tip} = 200
f_{tip} = 100 kHz
Hamaker =  100 zJ 
B_{sample} = 10 
k_{sample} = 0.05
m_{tip} = 5.0"""
text6 =  """Monomodal
A_{01} = 41 mn
R_{tip} = 50 nm
k_{tip} = 0.9 N/m
Q_{tip} = 200
f_{tip} = 100 kHz
Hamaker =  100 zJ
m_{tip} = 5.0
m_{sample} = 0.0003"""


class examples_app():
    
    def refine(self,text):
        text = text.replace("_{","/split/split")
        text = text.replace("}","/splitsubscript/split ")
        text = text.split("/split")
        text = tuple(text)
        return text

    def __init__(self, master):     
        
        
        
        self.master = master           
        self.examples_menu = tk.Toplevel(master.toplevel1)
        self.examples_menu.wm_title("Examples")
        
        
        
        
        self.img = tk.PhotoImage(file="icons/icon.png")
        self.examples_menu.wm_iconphoto(False, self.img)
        self.frame = ttk.Frame(self.examples_menu)
        self.label1 = ttk.Label(self.frame)        
        self.img = tk.PhotoImage(file="icons/dForceLogo4-1.png")
        
        
        # Textboxes with details of the simulations
        # All this details are to show subscripts in the textbox
        # https://stackoverflow.com/questions/62230669/how-to-make-a-text-subscript-or-superscript-in-text-widget-of-python-tkinter
        
        for n,text in enumerate([text1,text2,text3,text4,text5,text6]):
            self.label_sr_sp = tk.Text(self.examples_menu,width=20,height=14,
                                        borderwidth=0,
                                        font=("TkDefaultFont",15),
                                        background = self.master.toplevel1.cget("background"))        
            self.label_sr_sp.tag_configure(tagName="subscript", offset=-4,font=("TkDefaultFont",7) )
            self.label_sr_sp.insert("insert", *self.refine(text))
            self.label_sr_sp.configure(state="disabled")
            self.label_sr_sp.grid(column=n, row=1)            
            self.label_sr_sp.bindtags((str(self.label_sr_sp), str(master), "all"))
            
            

        self.img_0 = tk.PhotoImage(file="icons/demo.png")
        self.label_display_0 = tk.Label(self.examples_menu,image = self.img_0)#,text = "Hertz",compound=tk.BOTTOM,font=("TkDefaultFont",15))   
        self.label_display_0.grid(column=0, row=0,columnspan=6)

        width = 20
        height=3
        pady=5
        # Buttons to run the simulations
        self.button1 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button1.configure(text="Hertz")        
        self.button1.grid(column=0, row=2, pady=pady)
        self.button1.configure(command=self.on_button_1_pressed)
        
        self.button2 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button2.configure(text="DMT")
        self.button2.grid(column=1, row=2, pady=pady)
        self.button2.configure(command=self.on_button_2_pressed)
        
        self.button3 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button3.configure(text="Bimodal")
        self.button3.grid(column=2, row=2, pady=pady)
        self.button3.configure(command=self.on_button_3_pressed)
        
        self.button4 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button4.configure(text="Kelvin-Voigt + bec")
        self.button4.grid(column=3, row=2, pady=pady)
        self.button4.configure(command=self.on_button_4_pressed)
        
        self.button5 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button5.configure(text="Magnetic transfer function")
        self.button5.grid(column=4, row=2, pady=pady)
        self.button5.configure(command=self.on_button_5_pressed)
        
        self.button6 = tk.Button(self.examples_menu,height=height,width = width,font=("TkDefaultFont",15))
        self.button6.configure(text="Magnetic dipole-dipole")        
        self.button6.configure(command=self.on_button_6_pressed)
        self.button6.grid(column=5, row=2, pady=pady)
        
   
        
        
        

    
        self.mainwindow = self.examples_menu

    def run(self):
        self.mainwindow.mainloop()
        
    def on_button_1_pressed(self):
        self.reset()
        
        self.master.multifrequency_mode.set(True) 
        self.master.mass_model.set('point_mass') #
    
    
        self.master.a01_sp.set("15.00")
        self.master.copy_fre1_var.set(1)
        self.master.nperfin_sp.set(10)
    
        self.master.hertz_check_var.set(1) # Hertz model
        self.master.ampst_cb_var.set(1) # Plot amplitude           
        self.master.fasest_cb_var.set(1) # Plot phase        
        self.master.zt_cb_var.set(1)
        self.master.ft_cb_var.set(1)
        
        self.master.on_activate_all_checks()
        self.master.on_calculate_simulation_parameters()
                    
        self.examples_menu.destroy()
        self.master.on_run_ex_clicked()
        self.master.on_plot_ex_clicked()
         
    def on_button_2_pressed(self):
         self.reset()
         
         self.master.multifrequency_mode.set(True) 
         self.master.mass_model.set('point_mass') #
     
     
         self.master.a01_sp.set("20.00")
         self.master.copy_fre1_var.set(1)
         self.master.nperfin_sp.set(50)
         self.master.deltazc_sp.set("1.0")
         self.master.dmt_check_var.set(1) # DMT model
         self.master.vdw_check_var.set(1) # VdW model
         self.master.toplevel1.ham_sp_var = 3.0
         self.master.toplevel1.saym_sp_var = 1000.0    
         
         self.master.ampst_cb_var.set(1) # Plot amplitude           
         self.master.fasest_cb_var.set(1) # Plot phase        
         self.master.zt_cb_var.set(1)
         self.master.ft_cb_var.set(1)
         self.master.simulate_retrace_var.set(1)
         
         self.master.on_activate_all_checks()
         self.master.on_calculate_simulation_parameters()
         
         self.examples_menu.destroy()
         self.master.on_run_ex_clicked()
         self.master.on_plot_ex_clicked()
    
        
        
    def on_button_3_pressed(self):
        
        self.reset()
         
        
        self.master.mass_model.set('euler') #
        self.master.multifrequency_mode.set(False) 
        self.master.euler_activated()
        self.master.bimodal_mode_activated()
        
        
        #self.master.a01_sp.set("10.00")
        #self.master.a02_sp.set("0.10")
        
        self.master.a01_sp.set("100.00")
        self.master.a02_sp.set("0.25")
        
        self.master.kc1_sp.set("40") # k1
        self.master.kc2_sp.set("1600") # k2
        
        self.master.copy_fre1_var.set(1)
        self.master.copy_fre2_var.set(1)
                   
        self.master.hertz_check_var.set(1) # Hertz model
        self.master.ampst_cb_var.set(1) # Plot amplitude
            
        self.master.ampst_cb_var.set(1)
        self.master.fasest_cb_var.set(1)
            
        #BIMODAL
        self.master.amps2_cb_var.set(1)
        self.master.fases2_cb_var.set(1)
        
        # As in nature protocols
        self.master.a01_sp.set("100.00")
        self.master.a02_sp.set("0.25")        
        self.master.kc1_sp.set("3.00") # k1
        
        self.master.f01_sp.set("82.00") # f1
        self.master.q1_sp.set("220.00") # Q1
        self.master.kc2_sp.set("152.00") # k2
        self.master.f02_sp.set("518.00") # f2
        self.master.q2_sp.set("1000") # Q1
        
                 
        self.master.on_activate_all_checks()
        self.master.on_calculate_simulation_parameters()
        
        
        self.examples_menu.destroy()
        self.master.on_run_ex_clicked()
        self.master.on_plot_ex_clicked()
    
        
    def on_button_4_pressed(self):
         
         self.reset()        
         self.master.multifrequency_mode.set(True) 
         self.master.mass_model.set('point_mass') #
         
         self.master.ampst_cb_var.set(1)
         self.master.fasest_cb_var.set(1)
         self.master.maforce_cb_var.set(1)
         self.master.zt_cb_var.set(1)
         self.master.ft_cb_var.set(1)
         
         
         
         self.master.tr_sp.set("5.00") # tip radius
         self.master.a01_sp.set("10.00") # A01
         self.master.kc1_sp.set("0.3") # k1
         self.master.f01_sp.set("350") # f1
         self.master.q1_sp.set("2.00") # Q1
         self.master.zcmin_sp.set("-1.00")
         self.master.zcmax_sp.set("3.00")
         self.master.fixedzc_sp.set("1.00")
         self.master.deltazc_sp.set("0.1")
         self.master.npp_sp.set("512")
         self.master.copy_fre1_var.set(1)
         self.master.nperfin_sp.set(10)
         
         self.master.on_calculate_simulation_parameters()      
         
         self.master.kv_sphere_check_var.set(1)
         self.master.toplevel1.musam_sp_var = 0.5 # poission coefficient of the sample
         self.master.toplevel1.saym_sp_var = 100 # sample young's modulus in MPa
         self.master.toplevel1.sahe_sp_var = 10 # height of the sample
         self.master.toplevel1.mvis_sp_var = 2.00 # Viscosity
         
         self.master.ampst_cb_var.set(1)
         self.master.maforce_cb_var.set(1)
         self.master.etst_cb_var.set(1)
         
         self.master.on_activate_all_checks()
         self.master.on_calculate_simulation_parameters()
         
         self.examples_menu.destroy()
         self.master.on_run_ex_clicked()
         self.master.on_plot_ex_clicked()
         
    def on_button_5_pressed(self):
        self.reset()
        
        self.master.multifrequency_mode.set(True) # This is confusing
        self.master.mass_model.set('point_mass') #
        
        self.master.magnetic_exponential_check_var.set(1)
        self.master.magnetic_dipole_check_var.set(0)
        # Magnetic
        
        self.master.toplevel1.magnetic_m_tip = 5.0 # PPP-MFMR
        self.master.toplevel1.magnetic_m_sample =  0.0003 # 
        
        self.master.toplevel1.magnetic_k_sample = 0.05 # nota kself.mastera = 2*np.pi/lambda
        
        
        
        self.master.vdw_check_var.set(1)
        self.master.toplevel1.ham_sp_var = 10.0
        
        self.master.tr_sp.set("50.00") # tip radius
        self.master.a01_sp.set("41.00") # A01
        #self.master.a02_sp.set("10.0") # A02
        
        self.master.kc1_sp.set("3.4") # k1
        self.master.f01_sp.set("79") # f1
        self.master.q1_sp.set("160") # Q1
        
        self.master.kc2_sp.set("172") # k2
        self.master.f02_sp.set("350") # f2
        self.master.q2_sp.set("502") # Q2
        
        self.master.on_calculate_simulation_parameters()
        self.master.zcmax_sp.set("70.00")
        self.master.fixedzc_sp.set("45.00")
        self.master.zcmin_sp.set("40.00")
        self.master.zcmax_sp.set("60.00")
        self.master.zcmin_sp.set("41.00")
        
        
        self.master.deltazc_sp.set("0.5")
        self.master.npp_sp.set("512")
        self.master.copy_fre1_var.set(1)
        self.master.nperfin_sp.set("128")
        
        self.master.copy_fre1_var.set(1)
        self.master.copy_fre2_var.set(1)
        
        self.master.ampst_cb_var.set(0) # Plot amplitude
        self.master.fasest_cb_var.set(0)
        
        self.master.toplevel1.magnetic_B0_sample = 10.0
        
        self.master.fasest_cb_var.set(1)
        
        self.master.on_activate_all_checks()
        self.master.on_calculate_simulation_parameters()
        self.master.nperfin_sp.set("128")        
        self.examples_menu.destroy()
        self.master.on_run_ex_clicked()
        self.master.on_plot_ex_clicked()
    
    def on_button_6_pressed(self):
        self.reset()
        
        self.master.multifrequency_mode.set(True) 
        self.master.mass_model.set('point_mass') #
        
        self.master.magnetic_exponential_check_var.set(0)
        self.master.magnetic_dipole_check_var.set(1)
        # Magnetic
        self.master.toplevel1.magnetic_m_tip = 5.0 # PPP-MFMR
        self.master.toplevel1.magnetic_m_sample = 0.0003
        
        self.master.toplevel1.magnetic_k_sample = 0.05 # nota kself.mastera = 2*np.pi/lambda
        self.master.vdw_check_var.set(1)
        self.master.toplevel1.ham_sp_var =  10.0
        
        
        self.master.tr_sp.set("50.00") # tip radius
        self.master.a01_sp.set("41.00") # A01
        self.master.a02_sp.set("10.0") # A02
        
        self.master.kc1_sp.set("3.4") # k1
        self.master.f01_sp.set("79") # f1
        self.master.q1_sp.set("160") # Q1
        
        self.master.kc2_sp.set("172") # k2
        self.master.f02_sp.set("350") # f2
        self.master.q2_sp.set("502") # Q2
        
        self.master.on_calculate_simulation_parameters()
        self.master.fixedzc_sp.set("45.00")
        self.master.zcmax_sp.set("50.00")
        self.master.zcmin_sp.set("40.00")
        
        
        self.master.zcmin_sp.set("42.00")
        self.master.zcmax_sp.set("60.00")
        self.master.deltazc_sp.set("0.25")
        
        self.master.npp_sp.set("1024")
        
        self.master.copy_fre1_var.set(1)
        
        
        
        self.master.copy_fre1_var.set(1)
        self.master.copy_fre2_var.set(1)
        
        self.master.ampst_cb_var.set(0) # Plot amplitude
        self.master.fasest_cb_var.set(0)
        
        self.master.toplevel1.magnetic_m_sample = 0.0003
        
        self.master.fasest_cb_var.set(1)
        
        
        
        self.master.on_activate_all_checks()
        self.master.on_calculate_simulation_parameters()        
        self.master.nperfin_sp.set("128")
        self.examples_menu.destroy()
        self.master.on_run_ex_clicked()
        self.master.on_plot_ex_clicked()
        
   
        
    def reset(self):
        
        self.master.flagpm = 1
        self.master.multifrequency_mode.set(True) 
        self.master.mass_model.set('point_mass') #
        self.master.point_mass_activated()
        self.master.single_mode_activated()
        self.master.multifrequency_mode.set(True) 
        self.master.mass_model.set('point_mass') #
        self.master.point_mass_activated()
        self.master.single_mode_activated()
        
        
        self.master.a01_sp.set('0.00')
        self.master.a02_sp.set('0.00')
        self.master.adhy_check_var.set(0)
        self.master.amp_x1_var.set(False)
        self.master.amp_x2_var.set(False)
        self.master.amp_x3_var.set(False)
        self.master.amp_y1_var.set(False)
        self.master.amp_y2_var.set(False)
        self.master.amp_y3_var.set(False)
        self.master.amps2_cb_var.set(0)
        self.master.ampst_cb_var.set(0)
        self.master.nanowire_check_var.set(0)
        self.master.becc_check_var.set(0)
        self.master.behc_check_var.set(0)
        self.master.bepc_check_var.set(0)
        self.master.benc_check_var.set(0)        
        self.master.DELR_check_var.set(0)        
        self.master.bias_var.set(0)
        self.master.copy_fre1_var.set(1)
        self.master.copy_fre2_var.set(1)
        self.master.ctime_cb_var.set(0)
        self.master.custom1_check_var.set(0)
        self.master.custom2_check_var.set(0)
        self.master.deltazc_sp.set('2.5')
        self.master.magnetic_dipole_check_var.set(0)
        self.master.dlvo_check_var.set(0)
        self.master.dmin_x1_var.set(False)
        self.master.dmin_x2_var.set(False)
        self.master.dmin_x3_var.set(False)
        self.master.dmin_y1_var.set(False)
        self.master.dmin_y2_var.set(False)
        self.master.dmin_y3_var.set(False)
        self.master.dmt_check_var.set(0)
        self.master.ets2_cb_var.set(0)
        self.master.ets_x1_var.set(False)
        self.master.ets_x2_var.set(False)
        self.master.ets_x3_var.set(False)
        self.master.ets_y1_var.set(False)
        self.master.ets_y2_var.set(False)
        self.master.ets_y3_var.set(False)
        self.master.etst_cb_var.set(0)
        self.master.magnetic_exponential_check_var.set(0)
        self.master.f01_sp.set('100.00')
        self.master.f02_sp.set('600')
        self.master.famp1_var.set(False)
        self.master.famp2_var.set(False)
        self.master.fampt_var.set(False)
        self.master.fases2_cb_var.set(0)
        self.master.fasest_cb_var.set(0)
        self.master.fcorrection_q1_var.set(0)
        self.master.fcorrection_q2_var.set(0)
        self.master.fd1_sp.set('0.00')
        self.master.fd2_sp.set('0.00')
        self.master.fd_cb_var.set(0)
        self.master.fixedzc_sp.set('10.00')
        self.master.fmax_x1_var.set(False)
        self.master.fmax_x2_var.set(False)
        self.master.fmax_x3_var.set(False)
        self.master.fmax_y1_var.set(False)
        self.master.fmax_y2_var.set(False)
        self.master.fmax_y3_var.set(False)
        self.master.fpha1_var.set(False)
        self.master.fpha2_var.set(False)
        self.master.fphat_var.set(False)
        self.master.ft_cb_var.set(0)
        self.master.hertz_check_var.set(0)
        self.master.indent_cb_var.set(0)
        self.master.lateralresolution_cb_var.set(0)
        self.master.it_cb_var.set(0)
        self.master.jkr_check_var.set(0)
        self.master.kc1_sp.set('0.90')
        self.master.kc2_sp.set('35.20')
        self.master.kv_cone_check_var.set(0)
        self.master.kv_punch_check_var.set(0)
        self.master.kv_sphere_check_var.set(0)
        self.master.maforce_cb_var.set(0)
        self.master.mass_model.set('point_mass')
        self.master.mdef_cb_var.set(0)
        self.master.mdef_x1_var.set(False)
        self.master.mdef_x2_var.set(False)
        self.master.mdef_x3_var.set(False)
        self.master.mdef_y1_var.set(False)
        self.master.mdef_y2_var.set(False)
        self.master.mdef_y3_var.set(False)
        self.master.mdis_cb_var.set(0)
        self.master.multifrequency_mode.set(True)
        self.master.mutip_sp.set('0.30')
        self.master.nper_sp.set('400')
        self.master.nperfin_sp.set('128')
        self.master.npp_sp.set('128')
        self.master.phase_x1_var.set(False)
        self.master.phase_x2_var.set(False)
        self.master.phase_x3_var.set(False)
        self.master.phase_y1_var.set(False)
        self.master.phase_y2_var.set(False)
        self.master.phase_y3_var.set(False)
        self.master.phi1a1_cb_var.set(False)
        self.master.phi1a2_cb_var.set(False)
        self.master.phi2a1_cb_var.set(False)
        self.master.phi2a2_cb_var.set(False)
        self.master.ps_cb_var.set(0)
        self.master.pts2_cb_var.set(0)
        self.master.pts_x1_var.set(False)
        self.master.pts_x2_var.set(False)
        self.master.pts_x3_var.set(False)
        self.master.pts_y1_var.set(False)
        self.master.pts_y2_var.set(False)
        self.master.pts_y3_var.set(False)
        self.master.ptst_cb_var.set(0)
        self.master.punch_check_var.set(0)
        self.master.q1_sp.set('200.00')
        self.master.q2_sp.set('1200.00')
        self.master.rescur_var.set(False)
        self.master.simulate_retrace_var.set(False)
        
        
        self.master.sneddon_check_var.set(0)
        self.master.tatara_check_var.set(0)
        self.master.tcd_x1_var.set(False)
        self.master.tcd_x2_var.set(False)
        self.master.tcd_x3_var.set(False)
        self.master.tcd_y1_var.set(False)
        self.master.tcd_y2_var.set(False)
        self.master.tcd_y3_var.set(False)
        self.master.tolerance_sp.set('9')
        self.master.tr_sp.set('5.00')
        self.master.vdw_check_var.set(0)
        self.master.virial2_cb_var.set(0)
        self.master.virialt_cb_var.set(0)
        self.master.viscosity_check_var.set(0)
        self.master.viscosity_cone_check_var.set(0)
        self.master.viscosity_punch_check_var.set(0)
        self.master.vt_cb_var.set(0)
        self.master.vts_x1_var.set(False)
        self.master.vts_x2_var.set(False)
        self.master.vts_x3_var.set(False)
        self.master.vts_y1_var.set(False)
        self.master.vts_y2_var.set(False)
        self.master.vts_y3_var.set(False)
        self.master.ymt_sp.set('160.00')
        
        self.master.z1_cb_var.set(0)
        self.master.z2_cb_var.set(0)
        self.master.z3_cb_var.set(0)
        self.master.z4_cb_var.set(0)
                
        
        self.master.zc_x1_var.set(False)
        self.master.zc_x2_var.set(False)
        self.master.zc_x3_var.set(False)
        self.master.zc_y1_var.set(False)
        self.master.zc_y2_var.set(False)
        self.master.zc_y3_var.set(False)
        self.master.zcmax_sp.set('22.00')
        self.master.zcmin_sp.set('2.5')
        self.master.zt_cb_var.set(0)
        
        self.master.imax_x1_var.set(False)
        self.master.imax_y1_var.set(False)
        self.master.imax_x2_var.set(False)
        self.master.imax_y2_var.set(False)
        self.master.imax_x3_var.set(False)
        self.master.imax_y3_var.set(False)
        
        self.master.kv_nanowire_check_var.set(0)
        self.master.viscosity_nanowire_check_var.set(0)
         
        #Variables
        self.master.toplevel1.sr_sp_var = 0.00
        self.master.toplevel1.saym_sp_var = 100.00
        self.master.toplevel1.saym2_sp_var = 100.0
        self.master.toplevel1.musam_sp_var = 0.3
        # Electrical
        self.master.toplevel1.epsilon_sp_var = 1.00
        self.master.toplevel1.sigmas_sp_var = 0.00
        self.master.toplevel1.sigmat_sp_var = 0.00
        self.master.toplevel1.landadeb_sp_var = 1.00
        # Other
        self.master.toplevel1.ham_sp_var = 1.00
        self.master.toplevel1.ham2_sp_var = 0.00
        self.master.toplevel1.a00_sp_var = 0.169
        self.master.toplevel1.angle_sp_var = 10.00
        self.master.toplevel1.sahe_sp_var = 10.00
        self.master.toplevel1.bec_bonded_var = True
        # Viscosity
        self.master.toplevel1.mvis_sp_var = 0.00
         
        # User customized
        self.master.toplevel1.custom1_sp_var = 0.00
        self.master.toplevel1.custom2_sp_var = 0.00
        self.master.toplevel1.custom3_sp_var = 0.00
        self.master.toplevel1.custom4_sp_var = 0.00
         
        # Magnetic
        self.master.toplevel1.magnetic_m_tip = 2.5 # SSS-MFMR
        self.master.toplevel1.magnetic_m_sample = 0.01 # 
        self.master.toplevel1.magnetic_B0_sample = 1.00 # 1 mT
        self.master.toplevel1.magnetic_k_sample = 0.1 # nota kself.appa = 2*np.pi/lambda
        
        self.master.on_activate_all_checks()
        self.master.point_mass_activated()
        self.master.single_mode_activated()

    def examples_menu_exit(self):
        self.examples_menu.destroy()