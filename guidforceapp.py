#!/usr/bin/python3
import tkinter as tk
import tkinter.ttk as ttk


class GuiDforceApp:
    def __init__(self, master=None):
        # build ui
        self.toplevel1 = tk.Tk() if master is None else tk.Toplevel(master)
        self.frame1 = ttk.Frame(self.toplevel1)
        self.label17 = ttk.Label(self.frame1)
        self.img_dForceLogo41 = tk.PhotoImage(file="icons/dForceLogo4-1.png")
        self.label17.configure(anchor="n", image=self.img_dForceLogo41)
        self.label17.pack(anchor="n")
        self.Toolbar_Frame = ttk.Labelframe(self.frame1)
        self.button_run = ttk.Button(self.Toolbar_Frame)
        self.img_run = tk.PhotoImage(file="icons/run.png")
        self.button_run.configure(
            compound="left", image=self.img_run, style="Toolbutton", text="Run"
        )
        self.button_run.configure(width=25)
        self.button_run.pack(side="top")
        self.button_run.configure(command=self.on_run_ex_clicked)
        self.button_abort = ttk.Button(self.Toolbar_Frame)
        self.img_abort = tk.PhotoImage(file="icons/abort.png")
        self.button_abort.configure(
            compound="left", image=self.img_abort, state="disabled", style="Toolbutton"
        )
        self.button_abort.configure(text="Abort!!!", width=25)
        self.button_abort.pack(side="top")
        self.button_abort.configure(command=self.on_abort_clicked)
        self.button_plot = ttk.Button(self.Toolbar_Frame)
        self.img_Plot = tk.PhotoImage(file="icons/Plot.png")
        self.button_plot.configure(
            compound="left", image=self.img_Plot, style="Toolbutton", text="Plot"
        )
        self.button_plot.configure(width=25)
        self.button_plot.pack(side="top")
        self.button_plot.configure(command=self.on_plot_ex_clicked)
        self.button_examples = ttk.Button(self.Toolbar_Frame)
        self.img_examples = tk.PhotoImage(file="icons/examples.png")
        self.button_examples.configure(
            compound="left",
            image=self.img_examples,
            style="Toolbutton",
            text="Examples",
        )
        self.button_examples.configure(width=25)
        self.button_examples.pack(side="top")
        self.button_examples.configure(command=self.on_examples_ex_clicked)
        self.button_about = ttk.Button(self.Toolbar_Frame)
        self.img_About = tk.PhotoImage(file="icons/About.png")
        self.button_about.configure(
            compound="left", image=self.img_About, style="Toolbutton", text="About"
        )
        self.button_about.configure(width=25)
        self.button_about.pack(side="top")
        self.button_about.configure(command=self.on_about_clicked)
        self.button_quit = ttk.Button(self.Toolbar_Frame)
        self.img_exit = tk.PhotoImage(file="icons/exit.png")
        self.button_quit.configure(
            compound="left", image=self.img_exit, style="Toolbutton", text="Quit"
        )
        self.button_quit.configure(width=25)
        self.button_quit.pack(side="top")
        self.button_quit.configure(command=self.gtk_main_quit)
        self.Toolbar_Frame.configure(text="Toolbar")
        self.Toolbar_Frame.pack(anchor="center", side="top")
        #======================PROGRESSBAR=====================================
        self.style = ttk.Style(master)
        # add label in the layout
        self.style.layout('text.Horizontal.TProgressbar', 
                     [('Horizontal.Progressbar.trough',
                       {'children': [('Horizontal.Progressbar.pbar',
                                      {'side': 'left', 'sticky': 'ns'})],
                        'sticky': 'nswe'}), 
                      ('Horizontal.Progressbar.label', {'sticky': 'nswe'})])
        # set initial text
        self.style.configure('text.Horizontal.TProgressbar', text='0 %', anchor='center')        
        
        
        self.progressbar1 = ttk.Progressbar(self.frame1,
                                            style='text.Horizontal.TProgressbar')
        self.progressbar1.configure(orient="horizontal")
        self.progressbar1.pack(anchor="s", expand="true", fill="x", side="bottom")
        self.label39 = ttk.Label(self.frame1)
        self.label39.configure(text="Advanced Force Microscopy Lab\ncsic.es")
        self.label39.pack(expand="false", fill="both", side="bottom")
        self.frame1.configure(height=200, width=200)
        self.frame1.pack(fill="y", side="right")
        self.tip_and_excitation = ttk.Notebook(self.toplevel1)
        self.labelframe1 = ttk.Labelframe(self.tip_and_excitation)
        self.label46 = ttk.Label(self.labelframe1)
        self.label46.configure(text="Driving frequencies")
        self.label46.grid(column=0, row=0, rowspan=2)
        self.separator16 = ttk.Separator(self.labelframe1)
        self.separator16.configure(orient="horizontal")
        self.separator16.grid(column=0, columnspan=2, row=2,sticky = "ew")
        self.label54 = ttk.Label(self.labelframe1)
        self.label54.configure(text="Dynamic System\t")
        self.label54.grid(column=0, row=3, rowspan=2)
        self.single_cb = ttk.Radiobutton(self.labelframe1)
        self.multifrequency_mode = tk.BooleanVar(value="1")
        self.single_cb.configure(
            text="Simple", value=1, variable=self.multifrequency_mode
        )
        self.single_cb.grid(column=1, row=0, sticky="ew")
        self.single_cb.configure(command=self.single_mode_activated)
        self.bmafm_cb = ttk.Radiobutton(self.labelframe1)
        self.bmafm_cb.configure(
            text="Bimodal AFM", value=0, variable=self.multifrequency_mode
        )
        self.bmafm_cb.grid(column=1, row=1, sticky="ew")
        self.bmafm_cb.configure(command=self.bimodal_mode_activated)
        self.pm_cb = ttk.Radiobutton(self.labelframe1)
        self.mass_model = tk.StringVar(value="point_mass")
        self.pm_cb.configure(
            text="Point-mass", value="point_mass", variable=self.mass_model
        )
        self.pm_cb.grid(column=1, row=3, sticky="ew")
        self.pm_cb.configure(command=self.point_mass_activated)
        self.eb_cb = ttk.Radiobutton(self.labelframe1)
        self.eb_cb.configure(
            text="Euler-Bernoulli", value="euler", variable=self.mass_model
        )
        self.eb_cb.grid(column=1, row=4, sticky="ew")
        self.eb_cb.configure(command=self.euler_activated)
        self.labelframe1.configure(height=200, text="Inputs", width=200)
        self.labelframe1.grid(column=0, row=0, sticky="ew")
        self.labelframe1.grid_anchor("center")
        self.labelframe1.rowconfigure(0, weight=1)
        self.labelframe1.rowconfigure("all", weight=1)
        self.labelframe1.columnconfigure(0, weight=1)
        self.labelframe1.columnconfigure("all", weight=1)
        self.tip_and_excitation.add(
            self.labelframe1, padding=10, text="Cantilever dynamics"
        )
        self.frame8 = ttk.Frame(self.tip_and_excitation)
        self.label24 = ttk.Label(self.frame8)
        self.label24.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label24.configure(text="Tip radius (nm)", width=30)
        self.label24.grid(column=1, row=0)
        self.tr_sp = ttk.Spinbox(self.frame8)
        self.tr_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "5.00"
        self.tr_sp.delete("0", "end")
        self.tr_sp.insert("0", _text_)
        self.tr_sp.grid(column=2, row=0)
        self.label26 = ttk.Label(self.frame8)
        self.label26.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label26.configure(text="Tip Young's modulus (GPa)", width=30)
        self.label26.grid(column=1, row=1)
        self.ymt_sp = ttk.Spinbox(self.frame8)
        self.ymt_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "160.00"
        self.ymt_sp.delete("0", "end")
        self.ymt_sp.insert("0", _text_)
        self.ymt_sp.grid(column=2, row=1)
        self.label28 = ttk.Label(self.frame8)
        self.label28.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label28.configure(text="Tip Poisson coeficient", width=30)
        self.label28.grid(column=1, row=2)
        self.mutip_sp = ttk.Spinbox(self.frame8)
        self.mutip_sp.configure(from_=0, increment=0.01, justify="right", to=1)
        _text_ = "0.30"
        self.mutip_sp.delete("0", "end")
        self.mutip_sp.insert("0", _text_)
        self.mutip_sp.grid(column=2, row=2)
        self.separator23 = ttk.Separator(self.frame8)
        self.separator23.configure(orient="horizontal")
        self.separator23.grid(column=0, columnspan=6, row=3,sticky = "ew")
        self.label6 = ttk.Label(self.frame8)
        self.label6.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label6.configure(text="Spring constant k\u2081 (N/m)", width=30)
        self.label6.grid(column=1, row=4)
        self.kc1_sp = ttk.Spinbox(self.frame8)
        self.kc1_sp.configure(
            font="TkTextFont", from_=0, increment=0.01, justify="right"
        )
        self.kc1_sp.configure(to=100000)
        _text_ = "0.90"
        self.kc1_sp.delete("0", "end")
        self.kc1_sp.insert("0", _text_)
        self.kc1_sp.grid(column=2, row=4)
        self.label32 = ttk.Label(self.frame8)
        self.label32.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label32.configure(text="Natural resonant frequency (kHz)", width=30)
        self.label32.grid(column=1, row=5)
        self.f01_sp = ttk.Spinbox(self.frame8)
        self.f01_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "100.00"
        self.f01_sp.delete("0", "end")
        self.f01_sp.insert("0", _text_)
        self.f01_sp.grid(column=2, row=5)
        self.fcorrection_q1 = ttk.Checkbutton(self.frame8)
        self.fcorrection_q1_var = tk.IntVar(value=0)
        self.fcorrection_q1.configure(
            offvalue=0,
            onvalue=1,
            text="Measured resonant frequency (low Q)",
            variable=self.fcorrection_q1_var,
        )
        self.fcorrection_q1.grid(column=1, columnspan=2, row=6, sticky="ew")
        self.label34 = ttk.Label(self.frame8)
        self.label34.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label34.configure(text="Quality factor Q\u2081", width=30)
        self.label34.grid(column=1, row=7)
        self.q1_sp = ttk.Spinbox(self.frame8)
        self.q1_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "200.00"
        self.q1_sp.delete("0", "end")
        self.q1_sp.insert("0", _text_)
        self.q1_sp.grid(column=2, row=7)
        self.label31 = ttk.Label(self.frame8)
        self.label31.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label31.configure(text="Mode 1", width=30)
        self.label31.grid(row=4, rowspan=4)
        self.separator24 = ttk.Separator(self.frame8)
        self.separator24.configure(orient="horizontal")
        self.separator24.grid(column=0, columnspan=6, row=8,sticky = "ew")
        self.label37 = ttk.Label(self.frame8)
        self.label37.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label37.configure(text="Mode 2", width=30)
        self.label37.grid(row=9, rowspan=4)
        self.label36 = ttk.Label(self.frame8)
        self.label36.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label36.configure(text="Spring constant k\u2082 (N/m)", width=30)
        self.label36.grid(column=1, row=9)
        self.kc2_sp = ttk.Spinbox(self.frame8)
        self.kc2_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "35.20"
        self.kc2_sp.delete("0", "end")
        self.kc2_sp.insert("0", _text_)
        self.kc2_sp.grid(column=2, row=9)
        self.label38 = ttk.Label(self.frame8)
        self.label38.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label38.configure(text="Natural resonant frequency (kHz)", width=30)
        self.label38.grid(column=1, row=10)
        self.f02_sp = ttk.Spinbox(self.frame8)
        self.f02_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "600"
        self.f02_sp.delete("0", "end")
        self.f02_sp.insert("0", _text_)
        self.f02_sp.grid(column=2, row=10)
        self.fcorrection_q2 = ttk.Checkbutton(self.frame8)
        self.fcorrection_q2_var = tk.IntVar(value=0)
        self.fcorrection_q2.configure(
            offvalue=0,
            onvalue=1,
            text="Measured resonant frequency (low Q)",
            variable=self.fcorrection_q2_var,
        )
        self.fcorrection_q2.grid(column=1, columnspan=2, row=11, sticky="ew")
        self.q2_sp = ttk.Spinbox(self.frame8)
        self.q2_sp.configure(from_=0, increment=0.01, justify="right", to=10000000)
        _text_ = "1200.00"
        self.q2_sp.delete("0", "end")
        self.q2_sp.insert("0", _text_)
        self.q2_sp.grid(column=2, row=12)
        self.label7 = ttk.Label(self.frame8)
        self.label7.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label7.configure(text="Quality factor Q\u2082", width=30)
        self.label7.grid(column=1, row=12)
        self.a01_sp = ttk.Spinbox(self.frame8)
        self.a01_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "0.00"
        self.a01_sp.delete("0", "end")
        self.a01_sp.insert("0", _text_)
        self.a01_sp.grid(column=5, row=4)
        self.label59 = ttk.Label(self.frame8)
        self.label59.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label59.configure(text="Free amplitude A\u2080\u2081 (nm)", width=30)
        self.label59.grid(column=4, row=4)
        self.label61 = ttk.Label(self.frame8)
        self.label61.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label61.configure(text="Driving frequency f\u2080\u2081 (KHz)", width=30)
        self.label61.grid(column=4, row=5)
        self.fd1_sp = ttk.Spinbox(self.frame8)
        self.fd1_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "0.00"
        self.fd1_sp.delete("0", "end")
        self.fd1_sp.insert("0", _text_)
        self.fd1_sp.grid(column=5, row=5)
        self.copy_fre1 = ttk.Checkbutton(self.frame8)
        self.copy_fre1_var = tk.IntVar(value=0)
        self.copy_fre1.configure(
            offvalue=0, onvalue=1, state="normal", text="Copy resonant frequency"
        )
        self.copy_fre1.configure(variable=self.copy_fre1_var)
        self.copy_fre1.grid(column=4, row=6)
        self.label66 = ttk.Label(self.frame8)
        self.label66.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label66.configure(text="Free amplitude A\u2080\u2082 (nm)", width=30)
        self.label66.grid(column=4, row=9)
        self.a02_sp = ttk.Spinbox(self.frame8)
        self.a02_sp.configure(from_=0, increment=0.00, justify="right", to=100000)
        _text_ = "0.00"
        self.a02_sp.delete("0", "end")
        self.a02_sp.insert("0", _text_)
        self.a02_sp.grid(column=5, row=9)
        self.label68 = ttk.Label(self.frame8)
        self.label68.configure(
            anchor="n", compound="bottom", justify="center", takefocus=True
        )
        self.label68.configure(text="Driving frequency f\u2080\u2082 (KHz)", width=30)
        self.label68.grid(column=4, row=10)
        self.fd2_sp = ttk.Spinbox(self.frame8)
        self.fd2_sp.configure(from_=0, increment=0.01, justify="right", to=100000)
        _text_ = "0.00"
        self.fd2_sp.delete("0", "end")
        self.fd2_sp.insert("0", _text_)
        self.fd2_sp.grid(column=5, row=10)
        self.copy_fre2 = ttk.Checkbutton(self.frame8)
        self.copy_fre2_var = tk.IntVar(value=0)
        self.copy_fre2.configure(
            offvalue=0,
            onvalue=1,
            text="Copy resonant frequency",
            variable=self.copy_fre2_var,
        )
        self.copy_fre2.grid(column=4, row=11)
        self.frame8.configure(height=200, width=1)
        self.frame8.pack(fill="both", side="top")
        self.frame8.grid_propagate(0)
        self.frame8.grid_anchor("sw")
        self.frame8.rowconfigure(0, weight=1)
        self.frame8.rowconfigure("all", weight=1)
        self.frame8.columnconfigure(0, weight=1)
        self.frame8.columnconfigure("all", weight=1)
        self.tip_and_excitation.add(self.frame8, text="Tip and excitation")
        self.frame9 = ttk.Frame(self.tip_and_excitation)
        self.vdw_check = ttk.Checkbutton(self.frame9)
        self.vdw_check_var = tk.IntVar(value=0)
        self.vdw_check.configure(offvalue=0, onvalue=1, state="normal", takefocus=False)
        self.vdw_check.configure(text="vdW (Sphere)", variable=self.vdw_check_var)
        self.vdw_check.grid(column=0, row=1, sticky="ew")
        self.vdw_check.configure(command=self.on_vdw_activated)
        
        # vdw cone check
        self.vdw_cone_check = ttk.Checkbutton(self.frame9)
        self.vdw_cone_check_var = tk.IntVar(value=0)
        self.vdw_cone_check.configure(offvalue=0, onvalue=1, state="normal", takefocus=False)
        self.vdw_cone_check.configure(text="vdW (Cone)", variable=self.vdw_cone_check_var)
        self.vdw_cone_check.grid(column=0, row=2, sticky="ew")
        self.vdw_cone_check.configure(command=self.on_vdw_cone_activated)
        
        # vdw punch check
        self.vdw_punch_check = ttk.Checkbutton(self.frame9)
        self.vdw_punch_check_var = tk.IntVar(value=0)
        self.vdw_punch_check.configure(offvalue=0, onvalue=1, state="normal", takefocus=False)
        self.vdw_punch_check.configure(text="vdW (Punch)", variable=self.vdw_punch_check_var)
        self.vdw_punch_check.grid(column=0, row=3, sticky="ew")
        self.vdw_punch_check.configure(command=self.on_vdw_punch_activated)
        
        
        
        self.dlvo_check = ttk.Checkbutton(self.frame9)
        self.dlvo_check_var = tk.IntVar(value=0)
        self.dlvo_check.configure(
            offvalue=0, onvalue=1, text="DLVO (Sphere)", variable=self.dlvo_check_var
        )
        self.dlvo_check.grid(column=0, row=4, sticky="ew")
        self.dlvo_check.configure(command=self.on_dlvo_activated)
        self.label73 = ttk.Label(self.frame9)
        self.label73.configure(text="Short range (elastic)")
        self.label73.grid(column=3, row=0)
        self.hertz_check = ttk.Checkbutton(self.frame9)
        self.hertz_check_var = tk.IntVar(value=0)
        self.hertz_check.configure(
            offvalue=0, onvalue=1, text="Sphere", variable=self.hertz_check_var
        )
        self.hertz_check.grid(column=3, row=1, sticky="ew")
        self.hertz_check.configure(command=self.on_hertz_activated)
        self.dmt_check = ttk.Checkbutton(self.frame9)
        self.dmt_check_var = tk.IntVar(value=0)
        self.dmt_check.configure(
            offvalue=0, onvalue=1, text="DMT (Sphere)", variable=self.dmt_check_var
        )
        self.dmt_check.grid(column=0, row=9, sticky="ew")
        self.dmt_check.configure(command=self.on_dmt_activated)
        self.jkr_check = ttk.Checkbutton(self.frame9)
        self.jkr_check_var = tk.IntVar(value=0)
        self.jkr_check.configure(
            offvalue=0, onvalue=1, text="JKR", variable=self.jkr_check_var
        )
        self.jkr_check.grid(column=3, row=5, sticky="ew")
        self.jkr_check.configure(command=self.on_jkr_activated)
        self.becc_check = ttk.Checkbutton(self.frame9)
        self.becc_check_var = tk.IntVar(value=0)
        self.becc_check.configure(
            offvalue=0, onvalue=1, text="Cone", variable=self.becc_check_var
        )
        self.becc_check.grid(column=3, row=10, sticky="ew")
        self.becc_check.configure(command=self.on_becc_activated)
        self.bepc_check = ttk.Checkbutton(self.frame9)
        self.bepc_check_var = tk.IntVar(value=0)
        self.bepc_check.configure(
            offvalue=0, onvalue=1, text="Flat cylinder", variable=self.bepc_check_var
        )
        self.bepc_check.grid(column=3, row=11, sticky="ew")
        self.bepc_check.configure(command=self.on_bepc_activated)
        self.tatara_check = ttk.Checkbutton(self.frame9)
        self.tatara_check_var = tk.IntVar(value=0)
        self.tatara_check.configure(
            offvalue=0, onvalue=1, text="Tatara", variable=self.tatara_check_var
        )
        self.tatara_check.grid(column=3, row=6, sticky="ew")
        self.tatara_check.configure(command=self.on_tatara_activated)
        self.label74 = ttk.Label(self.frame9)
        self.label74.configure(text="Short range (viscoelastic)")
        self.label74.grid(column=5, row=0)
        self.viscosity_check = ttk.Checkbutton(self.frame9)
        self.viscosity_check_var = tk.IntVar(value=0)
        self.viscosity_check.configure(
            offvalue=0, onvalue=1, text="Sphere", variable=self.viscosity_check_var
        )
        self.viscosity_check.grid(column=5, row=2, sticky="ew")
        self.viscosity_check.configure(command=self.on_viscosity_activated)
        self.adhy_check = ttk.Checkbutton(self.frame9)
        self.adhy_check_var = tk.IntVar(value=0)
        self.adhy_check.configure(
            offvalue=0,
            onvalue=1,
            text="Adhesion hysteresis",
            variable=self.adhy_check_var,
        )
        self.adhy_check.grid(column=0, row=7, sticky="ew")
        self.adhy_check.configure(command=self.on_adhesion_activated)
        self.label23 = ttk.Label(self.frame9)
        self.label23.configure(text="User Customized")
        self.label23.grid(column=7, row=0)
        self.custom1_check = ttk.Checkbutton(self.frame9)
        self.custom1_check_var = tk.IntVar(value=0)
        self.custom1_check.configure(
            offvalue=0,
            onvalue=1,
            text="User customized 1",
            variable=self.custom1_check_var,
        )
        self.custom1_check.grid(column=7, row=1, sticky="ew")
        self.custom1_check.configure(command=self.on_custom1_activated)
        self.custom2_check = ttk.Checkbutton(self.frame9)
        self.custom2_check_var = tk.IntVar(value=0)
        self.custom2_check.configure(
            offvalue=0,
            onvalue=1,
            text="User customized 2",
            variable=self.custom2_check_var,
        )
        self.custom2_check.grid(column=7, row=2, sticky="ew")
        self.custom2_check.configure(command=self.on_custom2_activated)
        self.button1 = ttk.Button(self.frame9)
        self.button1.configure(text="Show value table")
        self.button1.grid(column=7, row=10, sticky="ew")
        self.button1.configure(command=self.on_show_value_clicked)
        self.separator9 = ttk.Separator(self.frame9)
        self.separator9.configure(orient="vertical")
        self.separator9.grid(column=1, row=0, rowspan=11,sticky = "ns")
        self.separator10 = ttk.Separator(self.frame9)
        self.separator10.configure(orient="vertical")
        self.separator10.grid(column=4, row=0, rowspan=11,sticky = "ns")
        self.separator11 = ttk.Separator(self.frame9)
        self.separator11.configure(orient="vertical")
        self.separator11.grid(column=6, row=0, rowspan=11,sticky = "ns")
        self.magnetic_dipole_check = ttk.Checkbutton(self.frame9)
        self.magnetic_dipole_check_var = tk.IntVar(value=0)
        self.magnetic_dipole_check.configure(
            offvalue=0,
            onvalue=1,
            text="Dipole-dipole magnetic interaction",
            variable=self.magnetic_dipole_check_var,
        )
        self.magnetic_dipole_check.grid(column=0, row=5, sticky="ew")
        self.magnetic_dipole_check.configure(command=self.on_dipole_activated)
        self.magnetic_exponential_check = ttk.Checkbutton(self.frame9)
        self.magnetic_exponential_check_var = tk.IntVar(value=0)
        self.magnetic_exponential_check.configure(
            offvalue=0,
            onvalue=1,
            text="Transfer function magnetic interaction",
            variable=self.magnetic_exponential_check_var,
        )
        self.magnetic_exponential_check.grid(column=0, row=6, sticky="ew")
        self.magnetic_exponential_check.configure(command=self.on_exp_activated)
        self.sneddon_check = ttk.Checkbutton(self.frame9)
        self.sneddon_check_var = tk.IntVar(value=0)
        self.sneddon_check.configure(
            offvalue=0, onvalue=1, text="Cone", variable=self.sneddon_check_var
        )
        self.sneddon_check.grid(column=3, row=2, sticky="ew")
        self.sneddon_check.configure(command=self.on_sneddon_activated)
        self.kv_sphere_check = ttk.Checkbutton(self.frame9)
        self.kv_sphere_check_var = tk.IntVar(value=0)
        self.kv_sphere_check.configure(
            offvalue=0, onvalue=1, text="Sphere", variable=self.kv_sphere_check_var
        )
        self.kv_sphere_check.grid(column=5, row=7, sticky="ew")
        self.kv_sphere_check.configure(command=self.on_kv_sphere_activated)
        self.kv_cone_check = ttk.Checkbutton(self.frame9)
        self.kv_cone_check_var = tk.IntVar(value=0)
        self.kv_cone_check.configure(
            offvalue=0, onvalue=1, text="Cone", variable=self.kv_cone_check_var
        )
        self.kv_cone_check.grid(column=5, row=8, sticky="ew")
        self.kv_cone_check.configure(command=self.on_kv_cone_activated)
        self.punch_check = ttk.Checkbutton(self.frame9)
        self.punch_check_var = tk.IntVar(value=0)
        self.punch_check.configure(
            offvalue=0, onvalue=1, text="Flat cylinder", variable=self.punch_check_var
        )
        self.punch_check.grid(column=3, row=3, sticky="ew")
        self.punch_check.configure(command=self.on_punch_activated)
        self.viscosity_cone_check = ttk.Checkbutton(self.frame9)
        self.viscosity_cone_check_var = tk.IntVar(value=0)
        self.viscosity_cone_check.configure(
            offvalue=0, onvalue=1, text="Cone", variable=self.viscosity_cone_check_var
        )
        self.viscosity_cone_check.grid(column=5, row=3, sticky="ew")
        self.viscosity_cone_check.configure(command=self.on_viscosity_cone_activated)
        self.viscosity_punch_check = ttk.Checkbutton(self.frame9)
        self.viscosity_punch_check_var = tk.IntVar(value=0)
        self.viscosity_punch_check.configure(
            offvalue=0,
            onvalue=1,
            text="Flat cylinder",
            variable=self.viscosity_punch_check_var,
        )
        self.viscosity_punch_check.grid(column=5, row=4, sticky="ew")
        self.viscosity_punch_check.configure(command=self.on_viscosity_punch_activated)
        self.kv_punch_check = ttk.Checkbutton(self.frame9)
        self.kv_punch_check_var = tk.IntVar(value=0)
        self.kv_punch_check.configure(
            offvalue=0,
            onvalue=1,
            text="Flat cylinder",
            variable=self.kv_punch_check_var,
        )
        self.kv_punch_check.grid(column=5, row=9, sticky="ew")
        self.kv_punch_check.configure(command=self.on_kv_punch_activated)
        self.label35 = ttk.Label(self.frame9)
        self.label35.configure(text="Bottom Effect Correction")
        self.label35.grid(column=5, row=6, sticky="ew")
        self.behc_check = ttk.Checkbutton(self.frame9)
        self.behc_check_var = tk.IntVar(value=0)
        self.behc_check.configure(
            offvalue=0, onvalue=1, text="Sphere", variable=self.behc_check_var
        )
        self.behc_check.grid(column=3, row=9, sticky="ew")
        self.behc_check.configure(command=self.on_behc_activated)
        self.label40 = ttk.Label(self.frame9)
        self.label40.configure(text="Bottom Effect Correction (bec)")
        self.label40.grid(column=3, row=8, sticky="ew")
        self.label33 = ttk.Label(self.frame9)
        self.label33.configure(text="Long range")
        self.label33.grid(column=0, row=0)
        self.lennardJones_check = ttk.Checkbutton(self.frame9)
        self.lennardJones_check_var = tk.IntVar(value=0)
        self.lennardJones_check.configure(
            offvalue=0,
            onvalue=1,
            text="Lennard–Jones Potential",
            variable=self.lennardJones_check_var,
        )
        self.lennardJones_check.grid(column=0, row=10, sticky="ew")
        self.lennardJones_check.configure(command=self.on_lennardJones_activated)
        self.label48 = ttk.Label(self.frame9)
        self.label48.configure(text="Short Range - Long range ")
        self.label48.grid(column=0, row=8)
        self.benc_check = ttk.Checkbutton(self.frame9)
        self.benc_check_var = tk.IntVar(value=0)
        self.benc_check.configure(
            offvalue=0, onvalue=1, text="Nanowire", variable=self.benc_check_var
        )
        self.benc_check.grid(column=3, row=12, sticky="ew")
        self.benc_check.configure(command=self.on_benc_activated)
        self.nanowire_check = ttk.Checkbutton(self.frame9)
        self.nanowire_check_var = tk.IntVar(value=0)
        self.nanowire_check.configure(
            offvalue=0, onvalue=1, text="Nanowire", variable=self.nanowire_check_var
        )
        self.nanowire_check.grid(column=3, row=4, sticky="ew")
        self.nanowire_check.configure(command=self.on_nanowire_activated)
        self.label49 = ttk.Label(self.frame9)
        self.label49.configure(text="Kelvin-Voigt")
        self.label49.grid(column=5, row=1, sticky="ew")
        self.DELR_check = ttk.Checkbutton(self.frame9)
        self.DELR_check_var = tk.IntVar(value=0)
        self.DELR_check.configure(
            offvalue=0,
            onvalue=1,
            text="Two layered material",
            variable=self.DELR_check_var,
        )
        self.DELR_check.grid(column=3, row=7, sticky="ew")
        self.DELR_check.configure(command=self.on_DELR_activated)
        self.viscosity_nanowire_check = ttk.Checkbutton(self.frame9)
        self.viscosity_nanowire_check_var = tk.IntVar(value=0)
        self.viscosity_nanowire_check.configure(
            offvalue=0,
            onvalue=1,
            text="Nanowire",
            variable=self.viscosity_nanowire_check_var,
        )
        self.viscosity_nanowire_check.grid(column=5, row=5, sticky="ew")
        self.viscosity_nanowire_check.configure(
            command=self.on_viscosity_nanowire_activated
        )
        self.kv_nanowire_check = ttk.Checkbutton(self.frame9)
        self.kv_nanowire_check_var = tk.IntVar(value=0)
        self.kv_nanowire_check.configure(
            offvalue=0, onvalue=1, text="Nanowire", variable=self.kv_nanowire_check_var
        )
        self.kv_nanowire_check.grid(column=5, row=10, sticky="ew")
        self.kv_nanowire_check.configure(command=self.on_kv_nanowire_activated)
        self.frame9.configure(height=200, width=200)
        self.frame9.pack(expand="true", fill="both", side="top")
        self.frame9.rowconfigure("all", weight=1)
        self.frame9.columnconfigure("all", weight=1)
        self.tip_and_excitation.add(self.frame9, text="Forces")
        self.frame10 = ttk.Frame(self.tip_and_excitation)
        self.label147 = ttk.Label(self.frame10)
        self.label147.configure(anchor="w", text="Number of periods")
        self.label147.grid(column=0, row=1, sticky="ew")
        self.nper_sp = ttk.Spinbox(self.frame10)
        self.nper_sp.configure(
            font="TkDefaultFont", from_=1, increment=1, justify="right"
        )
        self.nper_sp.configure(to=100000, validate="none")
        _text_ = "400"
        self.nper_sp.delete("0", "end")
        self.nper_sp.insert("0", _text_)
        self.nper_sp.grid(column=1, row=1, sticky="ew")
        self.label150 = ttk.Label(self.frame10)
        self.label150.configure(anchor="nw", text="Number of points per period")
        self.label150.grid(column=0, row=2, sticky="ew")
        self.npp_sp = ttk.Spinbox(self.frame10)
        self.npp_sp.configure(from_=1, increment=1, justify="right", to=100000)
        self.npp_sp.configure(validate="none")
        _text_ = "128"
        self.npp_sp.delete("0", "end")
        self.npp_sp.insert("0", _text_)
        self.npp_sp.grid(column=1, row=2, sticky="ew")
        self.label151 = ttk.Label(self.frame10)
        self.label151.configure(
            anchor="nw", text="Interval of periods to calculate the steady state"
        )
        self.label151.grid(column=0, row=3, sticky="ew")
        self.nperfin_sp = ttk.Spinbox(self.frame10)
        self.nperfin_sp.configure(
            from_=0, increment=1, justify="right", to=10000000000000000
        )
        self.nperfin_sp.configure(validate="none")
        _text_ = "128"
        self.nperfin_sp.delete("0", "end")
        self.nperfin_sp.insert("0", _text_)
        self.nperfin_sp.grid(column=1, row=3, sticky="ew")
        self.label152 = ttk.Label(self.frame10)
        self.label152.configure(
            anchor="nw", text="Max average tip-surface distance (nm)"
        )
        self.label152.grid(column=0, row=4, sticky="ew")
        self.zcmax_sp = ttk.Spinbox(self.frame10)
        self.zcmax_sp.configure(
            from_=-10000000000000000000,
            increment=0.01,
            justify="right",
            to=1000000000000000000000,
        )
        self.zcmax_sp.configure(validate="none")
        _text_ = "22.00"
        self.zcmax_sp.delete("0", "end")
        self.zcmax_sp.insert("0", _text_)
        self.zcmax_sp.grid(column=1, row=4, sticky="ew")
        self.label153 = ttk.Label(self.frame10)
        self.label153.configure(
            anchor="nw", text="Min average tip-surface distance (nm)"
        )
        self.label153.grid(column=0, row=5, sticky="ew")
        self.zcmin_sp = ttk.Spinbox(self.frame10)
        self.zcmin_sp.configure(
            from_=-10000000000000000000000,
            increment=0.01,
            justify="right",
            to=10000000000,
        )
        self.zcmin_sp.configure(validate="none")
        _text_ = "2.50"
        self.zcmin_sp.delete("0", "end")
        self.zcmin_sp.insert("0", _text_)
        self.zcmin_sp.grid(column=1, row=5, sticky="ew")
        self.label154 = ttk.Label(self.frame10)
        self.label154.configure(
            anchor="nw", text="Average tip-surface distance step (nm)"
        )
        self.label154.grid(column=0, row=6, sticky="ew")
        self.deltazc_sp = ttk.Spinbox(self.frame10)
        self.deltazc_sp.configure(
            from_=0, increment=0.01, justify="right", to=10000000000000000000
        )
        self.deltazc_sp.configure(validate="none")
        _text_ = "2.50"
        self.deltazc_sp.delete("0", "end")
        self.deltazc_sp.insert("0", _text_)
        self.deltazc_sp.grid(column=1, row=6, sticky="ew")
        self.label155 = ttk.Label(self.frame10)
        self.label155.configure(
            anchor="nw", text="Distance for time domain analysis (nm)"
        )
        self.label155.grid(column=0, row=7, sticky="ew")
        self.fixedzc_sp = ttk.Spinbox(self.frame10)
        self.fixedzc_sp.configure(
            from_=-1000000000000000000000,
            increment=0.01,
            justify="right",
            to=100000000000000000000,
        )
        self.fixedzc_sp.configure(validate="none")
        _text_ = "10.00"
        self.fixedzc_sp.delete("0", "end")
        self.fixedzc_sp.insert("0", _text_)
        self.fixedzc_sp.grid(column=1, row=7, sticky="ew")
        self.label20 = ttk.Label(self.frame10)
        self.label20.configure(
            anchor="nw",
            foreground="#0000ff",
            text="Tolerance (10^-x, x recommended 9 or higher)",
        )
        self.label20.grid(column=0, row=8, sticky="ew")
        self.tolerance_sp = ttk.Spinbox(self.frame10)
        self.tolerance_sp.configure(
            font="TkDefaultFont", from_=1, increment=1, justify="right"
        )
        self.tolerance_sp.configure(
            state="normal", to=1000000000000000000000000, validate="focusin"
        )
        _text_ = "9"
        self.tolerance_sp.delete("0", "end")
        self.tolerance_sp.insert("0", _text_)
        self.tolerance_sp.grid(column=1, row=8, sticky="ew")
        self.simulate_retrace_cb = ttk.Checkbutton(self.frame10)
        self.simulate_retrace_var = tk.IntVar(value=0)
        self.simulate_retrace_cb.configure(
            offvalue=0, onvalue=1, text="Retraction", variable=self.simulate_retrace_var
        )
        self.simulate_retrace_cb.grid(column=1, row=9)
        self.button10 = ttk.Button(self.frame10)
        self.button10.configure(text="Optimized values")
        self.button10.grid(column=0, row=9)
        self.button10.configure(command=self.on_calculate_simulation_parameters)
        self.frame10.configure(height=1, width=1)
        self.frame10.pack(expand="true", fill="both", side="top")
        self.frame10.grid_anchor("center")
        self.frame10.rowconfigure("all", weight=1)
        self.frame10.columnconfigure("all", weight=1)
        self.tip_and_excitation.add(self.frame10, text="Numerical integration values")
        self.tip_and_excitation.pack(expand="true", fill="both", side="top")
        self.notebook4 = ttk.Notebook(self.toplevel1)
        self.frame5 = ttk.Frame(self.notebook4)
        self.ampst_cb = ttk.Checkbutton(self.frame5)
        self.ampst_cb_var = tk.IntVar(value=0)
        self.ampst_cb.configure(
            offvalue=0, onvalue=1, text="Amplitude (z_c)", variable=self.ampst_cb_var
        )
        self.ampst_cb.grid(column=0, row=0, sticky="ew")
        self.fasest_cb = ttk.Checkbutton(self.frame5)
        self.fasest_cb_var = tk.IntVar(value=0)
        self.fasest_cb.configure(
            offvalue=0, onvalue=1, text="Phase shift (z_c)", variable=self.fasest_cb_var
        )
        self.fasest_cb.grid(column=0, row=1, sticky="ew")
        self.etst_cb = ttk.Checkbutton(self.frame5)
        self.etst_cb_var = tk.IntVar(value=0)
        self.etst_cb.configure(
            offvalue=0,
            onvalue=1,
            text="Dissipated energy (z_c)",
            variable=self.etst_cb_var,
        )
        self.etst_cb.grid(column=0, row=2, sticky="ew")
        self.ptst_cb = ttk.Checkbutton(self.frame5)
        self.ptst_cb_var = tk.IntVar(value=0)
        self.ptst_cb.configure(
            offvalue=0,
            onvalue=1,
            text="Dissipated power (z_c)",
            variable=self.ptst_cb_var,
        )
        self.ptst_cb.grid(column=0, row=3, sticky="ew")
        self.amps2_cb = ttk.Checkbutton(self.frame5)
        self.amps2_cb_var = tk.IntVar(value=0)
        self.amps2_cb.configure(
            offvalue=0,
            onvalue=1,
            text="2nd Amplitude (z_c)",
            variable=self.amps2_cb_var,
        )
        self.amps2_cb.grid(column=1, row=0, sticky="ew")
        self.fases2_cb = ttk.Checkbutton(self.frame5)
        self.fases2_cb_var = tk.IntVar(value=0)
        self.fases2_cb.configure(
            offvalue=0,
            onvalue=1,
            text="2nd Phase shift (z_c)",
            variable=self.fases2_cb_var,
        )
        self.fases2_cb.grid(column=1, row=1, sticky="ew")
        self.ets2_cb = ttk.Checkbutton(self.frame5)
        self.ets2_cb_var = tk.IntVar(value=0)
        self.ets2_cb.configure(
            offvalue=0,
            onvalue=1,
            text="2nd Dissipated energy (z_c)",
            variable=self.ets2_cb_var,
        )
        self.ets2_cb.grid(column=1, row=2, sticky="ew")
        self.pts2_cb = ttk.Checkbutton(self.frame5)
        self.pts2_cb_var = tk.IntVar(value=0)
        self.pts2_cb.configure(
            offvalue=0,
            onvalue=1,
            text="2nd Dissipated power (z_c)",
            variable=self.pts2_cb_var,
        )
        self.pts2_cb.grid(column=1, row=3, sticky="ew")
        self.virial2_cb = ttk.Checkbutton(self.frame5)
        self.virial2_cb_var = tk.IntVar(value=0)
        self.virial2_cb.configure(
            offvalue=0, onvalue=1, text="2nd Virial (z_c)", variable=self.virial2_cb_var
        )
        self.virial2_cb.grid(column=1, row=4, sticky="ew")
        self.virialt_cb = ttk.Checkbutton(self.frame5)
        self.virialt_cb_var = tk.IntVar(value=0)
        self.virialt_cb.configure(
            offvalue=0, onvalue=1, text="Virial (z_c)", variable=self.virialt_cb_var
        )
        self.virialt_cb.grid(column=0, row=4, sticky="ew")
        self.mdef_cb = ttk.Checkbutton(self.frame5)
        self.mdef_cb_var = tk.IntVar(value=0)
        self.mdef_cb.configure(
            offvalue=0,
            onvalue=1,
            text="Mean deflection (z_c)",
            variable=self.mdef_cb_var,
        )
        self.mdef_cb.grid(column=0, row=5, sticky="ew")
        self.mdis_cb = ttk.Checkbutton(self.frame5)
        self.mdis_cb_var = tk.IntVar(value=0)
        self.mdis_cb.configure(
            offvalue=0,
            onvalue=1,
            text="Minimum distance (z_c)",
            variable=self.mdis_cb_var,
        )
        self.mdis_cb.grid(column=0, row=6, sticky="ew")
        self.maforce_cb = ttk.Checkbutton(self.frame5)
        self.maforce_cb_var = tk.IntVar(value=0)
        self.maforce_cb.configure(
            offvalue=0,
            onvalue=1,
            text="Maximum force (z_c)",
            variable=self.maforce_cb_var,
        )
        self.maforce_cb.grid(column=0, row=7, sticky="ew")
        self.ctime_cb = ttk.Checkbutton(self.frame5)
        self.ctime_cb_var = tk.IntVar(value=0)
        self.ctime_cb.configure(
            offvalue=0, onvalue=1, text="Contact time (z_c)", variable=self.ctime_cb_var
        )
        self.ctime_cb.grid(column=0, row=8, sticky="ew")
        self.indent_cb = ttk.Checkbutton(self.frame5)
        self.indent_cb_var = tk.IntVar(value=0)
        self.indent_cb.configure(
            offvalue=0, onvalue=1, state="normal", text="Indentation (z_c)"
        )
        self.indent_cb.configure(variable=self.indent_cb_var)
        self.indent_cb.grid(column=1, row=5, sticky="ew")
        self.lateralresolution_cb = ttk.Checkbutton(self.frame5)
        self.lateralresolution_cb_var = tk.IntVar(value=0)
        self.lateralresolution_cb.configure(
            offvalue=0, onvalue=1, state="normal", text="Lateral resolution (z_c)"
        )
        self.lateralresolution_cb.configure(variable=self.lateralresolution_cb_var)
        self.lateralresolution_cb.grid(column=1, row=6, sticky="ew")
        self.frame5.configure(borderwidth=10, height=1, width=1)
        self.frame5.grid(column=0, row=0, sticky="ns")
        self.frame5.rowconfigure("all", pad=10, weight=1)
        self.frame5.columnconfigure("all", pad=150, weight=1)
        self.notebook4.add(self.frame5, text="Distance domain")
        self.frame3 = ttk.Frame(self.notebook4)
        self.zt_cb = ttk.Checkbutton(self.frame3)
        self.zt_cb_var = tk.IntVar(value=0)
        self.zt_cb.configure(
            offvalue=0, onvalue=1, text="z(t)", variable=self.zt_cb_var
        )
        self.zt_cb.grid(column=0, row=0, sticky="ew")
        self.vt_cb = ttk.Checkbutton(self.frame3)
        self.vt_cb_var = tk.IntVar(value=0)
        self.vt_cb.configure(
            offvalue=0, onvalue=1, text="v(t)", variable=self.vt_cb_var
        )
        self.vt_cb.grid(column=0, row=1, sticky="ew")
        self.ft_cb = ttk.Checkbutton(self.frame3)
        self.ft_cb_var = tk.IntVar(value=0)
        self.ft_cb.configure(
            offvalue=0, onvalue=1, text="F(t)", variable=self.ft_cb_var
        )
        self.ft_cb.grid(column=0, row=2, sticky="ew")
        self.it_cb = ttk.Checkbutton(self.frame3)
        self.it_cb_var = tk.IntVar(value=0)
        self.it_cb.configure(
            offvalue=0, onvalue=1, text="Indentation (t)", variable=self.it_cb_var
        )
        self.it_cb.grid(column=0, row=3, sticky="ew")
        self.z1_cb = ttk.Checkbutton(self.frame3)
        self.z1_cb_var = tk.IntVar(value=0)
        self.z1_cb.configure(
            offvalue=0, onvalue=1, text="z\u2081(t)", variable=self.z1_cb_var
        )
        self.z1_cb.grid(column=1, row=0, sticky="ew")
        self.z2_cb = ttk.Checkbutton(self.frame3)
        self.z2_cb_var = tk.IntVar(value=0)
        self.z2_cb.configure(
            offvalue=0, onvalue=1, text="z\u2082(t)", variable=self.z2_cb_var
        )
        self.z2_cb.grid(column=1, row=1, sticky="ew")
        
        self.z3_cb = ttk.Checkbutton(self.frame3)
        self.z3_cb_var = tk.IntVar(value=0)
        #self.z3_cb.configure(state="disabled", text="z\u2083(t)", variable=self.z3_cb_var,onvalue = 1,offvalue = 0)
        self.z3_cb.configure(offvalue=0, onvalue=1, state="disabled", text="z\u2083(t)")
        self.z3_cb.configure(variable=self.z3_cb_var)
        self.z3_cb.grid(column=1, row=2, sticky="ew")
        
        self.z4_cb = ttk.Checkbutton(self.frame3)
        self.z4_cb_var = tk.IntVar(value=0)
        self.z4_cb.configure(offvalue=0, onvalue=1, state="disabled", text="z\u2084(t)")
        self.z4_cb.configure(variable=self.z4_cb_var)
        self.z4_cb.grid(column=1, row=3, sticky="ew")
        
        self.bias = ttk.Checkbutton(self.frame3)
        self.bias_var = tk.IntVar(value=0)
        self.bias.configure(
            offvalue=0, onvalue=1, text="Bias graph", variable=self.bias_var
        )
        self.bias.grid(column=1, row=7, sticky="ew")
        self.fd_cb = ttk.Checkbutton(self.frame3)
        self.fd_cb_var = tk.IntVar(value=0)
        self.fd_cb.configure(
            offvalue=0, onvalue=1, text="F(d)", variable=self.fd_cb_var
        )
        self.fd_cb.grid(column=0, row=4, sticky="ew")
        self.ps_cb = ttk.Checkbutton(self.frame3)
        self.ps_cb_var = tk.IntVar(value=0)
        self.ps_cb.configure(
            offvalue=0, onvalue=1, text="Phase space (z vs v)", variable=self.ps_cb_var
        )
        self.ps_cb.grid(column=0, row=5, sticky="ew")
        self.frame3.configure(borderwidth=10)
        self.frame3.grid(column=0, row=0)
        self.frame3.grid_anchor("nw")
        self.frame3.rowconfigure("all", pad=20, weight=1)
        self.frame3.columnconfigure("all", pad=300, weight=1)
        self.notebook4.add(self.frame3, text="Time domain")
        self.frame4 = ttk.Frame(self.notebook4)
        
        
        
        
        
        self.phi1a1_cb = ttk.Checkbutton(self.frame4)
        self.phi1a1_cb_var = tk.IntVar(value=0)
        self.phi1a1_cb.configure(text="ϕ\u2081 (A\u2081) ", variable=self.phi1a1_cb_var,offvalue=0, onvalue=1)
        self.phi1a1_cb.grid(column=0, row=0, sticky="ew")
        self.phi2a2_cb = ttk.Checkbutton(self.frame4)
        self.phi2a2_cb_var = tk.IntVar(value=0)
        self.phi2a2_cb.configure(text="ϕ\u2082 (A\u2082)", variable=self.phi2a2_cb_var,offvalue=0, onvalue=1)
        self.phi2a2_cb.grid(column=0, row=1, sticky="ew")
        self.phi2a1_cb = ttk.Checkbutton(self.frame4)
        self.phi2a1_cb_var = tk.IntVar(value=0)
        self.phi2a1_cb.configure(text="ϕ\u2082 (A\u2081)", variable=self.phi2a1_cb_var,offvalue=0, onvalue=1)
        self.phi2a1_cb.grid(column=1, row=0, sticky="ew")
        self.phi1a2_cb = ttk.Checkbutton(self.frame4)
        self.phi1a2_cb_var = tk.IntVar(value=0)
        self.phi1a2_cb.configure(text="ϕ\u2081 (A\u2082)", variable=self.phi1a2_cb_var,offvalue=0, onvalue=1)
        self.phi1a2_cb.grid(column=1, row=1, sticky="ew")
        self.frame4.configure(borderwidth=20)
        self.frame4.grid(column=0, row=0)
        self.frame4.grid_anchor("nw")
        self.frame4.rowconfigure("all", pad=80)
        self.frame4.columnconfigure("all", pad=300)
        
        self.notebook4.add(self.frame4, text="Cross-mode curves")
        self.frame2 = ttk.Frame(self.notebook4)
        self.fampt = ttk.Checkbutton(self.frame2)
        self.fampt_var = tk.IntVar(value=0)
        self.fampt.configure(text="Amp", variable=self.fampt_var,offvalue=0, onvalue=1)
        self.fampt.grid(column=0, row=0, sticky="ew")
        self.famp1 = ttk.Checkbutton(self.frame2)
        self.famp1_var = tk.IntVar(value=0)
        self.famp1.configure(text="Amp1", variable=self.famp1_var,offvalue=0, onvalue=1)
        self.famp1.grid(column=0, row=1, sticky="ew")
        self.famp2 = ttk.Checkbutton(self.frame2)
        self.famp2_var = tk.IntVar(value=0)
        self.famp2.configure(text="Amp2", variable=self.famp2_var,offvalue=0, onvalue=1)
        self.famp2.grid(column=0, row=2, sticky="ew")
        self.fphat = ttk.Checkbutton(self.frame2)
        self.fphat_var = tk.IntVar(value=0)
        self.fphat.configure(text="Phase", variable=self.fphat_var,offvalue=0, onvalue=1)
        self.fphat.grid(column=1, row=0, sticky="ew")
        self.fpha1 = ttk.Checkbutton(self.frame2)
        self.fpha1_var = tk.IntVar(value=0)
        self.fpha1.configure(text="Phase 1", variable=self.fpha1_var,offvalue=0, onvalue=1)
        self.fpha1.grid(column=1, row=1, sticky="ew")
        self.fpha2 = ttk.Checkbutton(self.frame2)
        self.fpha2_var = tk.IntVar(value=0)
        self.fpha2.configure(text="Phase 2", variable=self.fpha2_var,offvalue=0, onvalue=1)
        self.fpha2.grid(column=1, row=2, sticky="ew")
        self.rescur = ttk.Checkbutton(self.frame2)
        self.rescur_var = tk.IntVar(value=0)
        self.rescur.configure(text="Resonance curve", variable=self.rescur_var,offvalue = 0,onvalue = 1)
        self.rescur.grid(column=1, row=10, sticky="ew")
        self.frame2.configure(borderwidth=20)
        self.frame2.grid(column=0, row=0)
        self.frame2.grid_anchor("nw")
        self.frame2.rowconfigure("all", pad=50, weight=1)
        self.frame2.columnconfigure("all", pad=300, weight=1)
        self.notebook4.add(self.frame2, text="Frequency domain")
        self.frame6 = ttk.Frame(self.notebook4)
        self.label2 = ttk.Label(self.frame6)
        self.label2.grid(column=0, row=0)
        self.label3 = ttk.Label(self.frame6)
        self.label3.configure(text="Amplitude")
        self.label3.grid(column=0, row=2)
        self.label4 = ttk.Label(self.frame6)
        self.label4.configure(text="Phase")
        self.label4.grid(column=0, row=3)
        self.label5 = ttk.Label(self.frame6)
        self.label5.configure(text="z\u2080")
        self.label5.grid(column=0, row=4)
        self.label8 = ttk.Label(self.frame6)
        self.label8.configure(text="Dissipated power")
        self.label8.grid(column=0, row=5)
        self.label9 = ttk.Label(self.frame6)
        self.label9.configure(text="Dissipated energy")
        self.label9.grid(column=0, row=6)
        self.label10 = ttk.Label(self.frame6)
        self.label10.configure(text="Mean delfection")
        self.label10.grid(column=0, row=7)
        self.label11 = ttk.Label(self.frame6)
        self.label11.configure(text="Minimum distance")
        self.label11.grid(column=0, row=8)
        self.label12 = ttk.Label(self.frame6)
        self.label12.configure(text="Maximum force")
        self.label12.grid(column=0, row=9)
        self.label13 = ttk.Label(self.frame6)
        self.label13.configure(text="Contact time")
        self.label13.grid(column=0, row=10)
        self.label14 = ttk.Label(self.frame6)
        self.label14.configure(text="Virial")
        self.label14.grid(column=0, row=11)
        self.label1 = ttk.Label(self.frame6)
        self.label1.grid(column=0, row=1)
        self.label15 = ttk.Label(self.frame6)
        self.label15.configure(text="x")
        self.label15.grid(column=2, row=1)
        self.label16 = ttk.Label(self.frame6)
        self.label16.configure(text="y")
        self.label16.grid(column=3, row=1)
        self.amp_y1 = ttk.Checkbutton(self.frame6)
        self.amp_y1_var = tk.BooleanVar(value=0)
        self.amp_y1.configure(variable=self.amp_y1_var)
        self.amp_y1.grid(column=3, row=2)
        self.phase_y1 = ttk.Checkbutton(self.frame6)
        self.phase_y1_var = tk.BooleanVar(value=0)
        self.phase_y1.configure(variable=self.phase_y1_var)
        self.phase_y1.grid(column=3, row=3)
        self.zc_y1 = ttk.Checkbutton(self.frame6)
        self.zc_y1_var = tk.BooleanVar(value=0)
        self.zc_y1.configure(variable=self.zc_y1_var)
        self.zc_y1.grid(column=3, row=4)
        self.pts_y1 = ttk.Checkbutton(self.frame6)
        self.pts_y1_var = tk.BooleanVar(value=0)
        self.pts_y1.configure(variable=self.pts_y1_var)
        self.pts_y1.grid(column=3, row=5)
        self.ets_y1 = ttk.Checkbutton(self.frame6)
        self.ets_y1_var = tk.BooleanVar(value=0)
        self.ets_y1.configure(variable=self.ets_y1_var)
        self.ets_y1.grid(column=3, row=6)
        self.mdef_y1 = ttk.Checkbutton(self.frame6)
        self.mdef_y1_var = tk.BooleanVar(value=0)
        self.mdef_y1.configure(variable=self.mdef_y1_var)
        self.mdef_y1.grid(column=3, row=7)
        self.dmin_y1 = ttk.Checkbutton(self.frame6)
        self.dmin_y1_var = tk.BooleanVar(value=0)
        self.dmin_y1.configure(variable=self.dmin_y1_var)
        self.dmin_y1.grid(column=3, row=8)
        self.fmax_y1 = ttk.Checkbutton(self.frame6)
        self.fmax_y1_var = tk.BooleanVar(value=0)
        self.fmax_y1.configure(variable=self.fmax_y1_var)
        self.fmax_y1.grid(column=3, row=9)
        self.tcd_y1 = ttk.Checkbutton(self.frame6)
        self.tcd_y1_var = tk.BooleanVar(value=0)
        self.tcd_y1.configure(variable=self.tcd_y1_var)
        self.tcd_y1.grid(column=3, row=10)
        self.vts_y1 = ttk.Checkbutton(self.frame6)
        self.vts_y1_var = tk.BooleanVar(value=0)
        self.vts_y1.configure(variable=self.vts_y1_var)
        self.vts_y1.grid(column=3, row=11)
        self.amp_x1 = ttk.Checkbutton(self.frame6)
        self.amp_x1_var = tk.BooleanVar(value=0)
        self.amp_x1.configure(variable=self.amp_x1_var)
        self.amp_x1.grid(column=2, row=2)
        self.phase_x1 = ttk.Checkbutton(self.frame6)
        self.phase_x1_var = tk.BooleanVar(value=0)
        self.phase_x1.configure(variable=self.phase_x1_var)
        self.phase_x1.grid(column=2, row=3)
        self.zc_x1 = ttk.Checkbutton(self.frame6)
        self.zc_x1_var = tk.BooleanVar(value=0)
        self.zc_x1.configure(variable=self.zc_x1_var)
        self.zc_x1.grid(column=2, row=4)
        self.pts_x1 = ttk.Checkbutton(self.frame6)
        self.pts_x1_var = tk.BooleanVar(value=0)
        self.pts_x1.configure(variable=self.pts_x1_var)
        self.pts_x1.grid(column=2, row=5)
        self.ets_x1 = ttk.Checkbutton(self.frame6)
        self.ets_x1_var = tk.BooleanVar(value=0)
        self.ets_x1.configure(variable=self.ets_x1_var)
        self.ets_x1.grid(column=2, row=6)
        self.mdef_x1 = ttk.Checkbutton(self.frame6)
        self.mdef_x1_var = tk.BooleanVar(value=0)
        self.mdef_x1.configure(variable=self.mdef_x1_var)
        self.mdef_x1.grid(column=2, row=7)
        self.dmin_x1 = ttk.Checkbutton(self.frame6)
        self.dmin_x1_var = tk.BooleanVar(value=0)
        self.dmin_x1.configure(variable=self.dmin_x1_var)
        self.dmin_x1.grid(column=2, row=8)
        self.fmax_x1 = ttk.Checkbutton(self.frame6)
        self.fmax_x1_var = tk.BooleanVar(value=0)
        self.fmax_x1.configure(variable=self.fmax_x1_var)
        self.fmax_x1.grid(column=2, row=9)
        self.tcd_x1 = ttk.Checkbutton(self.frame6)
        self.tcd_x1_var = tk.BooleanVar(value=0)
        self.tcd_x1.configure(variable=self.tcd_x1_var)
        self.tcd_x1.grid(column=2, row=10)
        self.vts_x1 = ttk.Checkbutton(self.frame6)
        self.vts_x1_var = tk.BooleanVar(value=0)
        self.vts_x1.configure(variable=self.vts_x1_var)
        self.vts_x1.grid(column=2, row=11)
        self.label86 = ttk.Label(self.frame6)
        self.label86.configure(foreground="#0000ff", text="Plot 1")
        self.label86.grid(column=2, columnspan=2, row=0)
        self.amp_x2 = ttk.Checkbutton(self.frame6)
        self.amp_x2_var = tk.BooleanVar(value=0)
        self.amp_x2.configure(variable=self.amp_x2_var)
        self.amp_x2.grid(column=5, row=2)
        self.amp_y2 = ttk.Checkbutton(self.frame6)
        self.amp_y2_var = tk.BooleanVar(value=0)
        self.amp_y2.configure(variable=self.amp_y2_var)
        self.amp_y2.grid(column=6, row=2)
        self.phase_x2 = ttk.Checkbutton(self.frame6)
        self.phase_x2_var = tk.BooleanVar(value=0)
        self.phase_x2.configure(variable=self.phase_x2_var)
        self.phase_x2.grid(column=5, row=3)
        self.zc_x2 = ttk.Checkbutton(self.frame6)
        self.zc_x2_var = tk.BooleanVar(value=0)
        self.zc_x2.configure(variable=self.zc_x2_var)
        self.zc_x2.grid(column=5, row=4)
        self.pts_x2 = ttk.Checkbutton(self.frame6)
        self.pts_x2_var = tk.BooleanVar(value=0)
        self.pts_x2.configure(variable=self.pts_x2_var)
        self.pts_x2.grid(column=5, row=5)
        self.ets_x2 = ttk.Checkbutton(self.frame6)
        self.ets_x2_var = tk.BooleanVar(value=0)
        self.ets_x2.configure(variable=self.ets_x2_var)
        self.ets_x2.grid(column=5, row=6)
        self.mdef_x2 = ttk.Checkbutton(self.frame6)
        self.mdef_x2_var = tk.BooleanVar(value=0)
        self.mdef_x2.configure(variable=self.mdef_x2_var)
        self.mdef_x2.grid(column=5, row=7)
        self.dmin_x2 = ttk.Checkbutton(self.frame6)
        self.dmin_x2_var = tk.BooleanVar(value=0)
        self.dmin_x2.configure(variable=self.dmin_x2_var)
        self.dmin_x2.grid(column=5, row=8)
        self.fmax_x2 = ttk.Checkbutton(self.frame6)
        self.fmax_x2_var = tk.BooleanVar(value=0)
        self.fmax_x2.configure(variable=self.fmax_x2_var)
        self.fmax_x2.grid(column=5, row=9)
        self.tcd_x2 = ttk.Checkbutton(self.frame6)
        self.tcd_x2_var = tk.BooleanVar(value=0)
        self.tcd_x2.configure(variable=self.tcd_x2_var)
        self.tcd_x2.grid(column=5, row=10)
        self.vts_x2 = ttk.Checkbutton(self.frame6)
        self.vts_x2_var = tk.BooleanVar(value=0)
        self.vts_x2.configure(variable=self.vts_x2_var)
        self.vts_x2.grid(column=5, row=11)
        self.phase_y2 = ttk.Checkbutton(self.frame6)
        self.phase_y2_var = tk.BooleanVar(value=0)
        self.phase_y2.configure(variable=self.phase_y2_var)
        self.phase_y2.grid(column=6, row=3)
        self.zc_y2 = ttk.Checkbutton(self.frame6)
        self.zc_y2_var = tk.BooleanVar(value=0)
        self.zc_y2.configure(variable=self.zc_y2_var)
        self.zc_y2.grid(column=6, row=4)
        self.pts_y2 = ttk.Checkbutton(self.frame6)
        self.pts_y2_var = tk.BooleanVar(value=0)
        self.pts_y2.configure(variable=self.pts_y2_var)
        self.pts_y2.grid(column=6, row=5)
        self.ets_y2 = ttk.Checkbutton(self.frame6)
        self.ets_y2_var = tk.BooleanVar(value=0)
        self.ets_y2.configure(variable=self.ets_y2_var)
        self.ets_y2.grid(column=6, row=6)
        self.mdef_y2 = ttk.Checkbutton(self.frame6)
        self.mdef_y2_var = tk.BooleanVar(value=0)
        self.mdef_y2.configure(variable=self.mdef_y2_var)
        self.mdef_y2.grid(column=6, row=7)
        self.dmin_y2 = ttk.Checkbutton(self.frame6)
        self.dmin_y2_var = tk.BooleanVar(value=0)
        self.dmin_y2.configure(variable=self.dmin_y2_var)
        self.dmin_y2.grid(column=6, row=8)
        self.fmax_y2 = ttk.Checkbutton(self.frame6)
        self.fmax_y2_var = tk.BooleanVar(value=0)
        self.fmax_y2.configure(variable=self.fmax_y2_var)
        self.fmax_y2.grid(column=6, row=9)
        self.tcd_y2 = ttk.Checkbutton(self.frame6)
        self.tcd_y2_var = tk.BooleanVar(value=0)
        self.tcd_y2.configure(variable=self.tcd_y2_var)
        self.tcd_y2.grid(column=6, row=10)
        self.vts_y2 = ttk.Checkbutton(self.frame6)
        self.vts_y2_var = tk.BooleanVar(value=0)
        self.vts_y2.configure(variable=self.vts_y2_var)
        self.vts_y2.grid(column=6, row=11)
        self.amp_x3 = ttk.Checkbutton(self.frame6)
        self.amp_x3_var = tk.BooleanVar(value=0)
        self.amp_x3.configure(variable=self.amp_x3_var)
        self.amp_x3.grid(column=9, row=2)
        self.amp_y3 = ttk.Checkbutton(self.frame6)
        self.amp_y3_var = tk.BooleanVar(value=0)
        self.amp_y3.configure(variable=self.amp_y3_var)
        self.amp_y3.grid(column=10, row=2)
        self.phase_y3 = ttk.Checkbutton(self.frame6)
        self.phase_y3_var = tk.BooleanVar(value=0)
        self.phase_y3.configure(variable=self.phase_y3_var)
        self.phase_y3.grid(column=10, row=3)
        self.zc_y3 = ttk.Checkbutton(self.frame6)
        self.zc_y3_var = tk.BooleanVar(value=0)
        self.zc_y3.configure(variable=self.zc_y3_var)
        self.zc_y3.grid(column=10, row=4)
        self.pts_y3 = ttk.Checkbutton(self.frame6)
        self.pts_y3_var = tk.BooleanVar(value=0)
        self.pts_y3.configure(variable=self.pts_y3_var)
        self.pts_y3.grid(column=10, row=5)
        self.ets_y3 = ttk.Checkbutton(self.frame6)
        self.ets_y3_var = tk.BooleanVar(value=0)
        self.ets_y3.configure(variable=self.ets_y3_var)
        self.ets_y3.grid(column=10, row=6)
        self.mdef_y3 = ttk.Checkbutton(self.frame6)
        self.mdef_y3_var = tk.BooleanVar(value=0)
        self.mdef_y3.configure(variable=self.mdef_y3_var)
        self.mdef_y3.grid(column=10, row=7)
        self.dmin_y3 = ttk.Checkbutton(self.frame6)
        self.dmin_y3_var = tk.BooleanVar(value=0)
        self.dmin_y3.configure(variable=self.dmin_y3_var)
        self.dmin_y3.grid(column=10, row=8)
        self.fmax_y3 = ttk.Checkbutton(self.frame6)
        self.fmax_y3_var = tk.BooleanVar(value=0)
        self.fmax_y3.configure(variable=self.fmax_y3_var)
        self.fmax_y3.grid(column=10, row=9)
        self.tcd_y3 = ttk.Checkbutton(self.frame6)
        self.tcd_y3_var = tk.BooleanVar(value=0)
        self.tcd_y3.configure(variable=self.tcd_y3_var)
        self.tcd_y3.grid(column=10, row=10)
        self.vts_y3 = ttk.Checkbutton(self.frame6)
        self.vts_y3_var = tk.BooleanVar(value=0)
        self.vts_y3.configure(variable=self.vts_y3_var)
        self.vts_y3.grid(column=10, row=11)
        self.phase_x3 = ttk.Checkbutton(self.frame6)
        self.phase_x3_var = tk.BooleanVar(value=0)
        self.phase_x3.configure(variable=self.phase_x3_var)
        self.phase_x3.grid(column=9, row=3)
        self.zc_x3 = ttk.Checkbutton(self.frame6)
        self.zc_x3_var = tk.BooleanVar(value=0)
        self.zc_x3.configure(variable=self.zc_x3_var)
        self.zc_x3.grid(column=9, row=4)
        self.pts_x3 = ttk.Checkbutton(self.frame6)
        self.pts_x3_var = tk.BooleanVar(value=0)
        self.pts_x3.configure(variable=self.pts_x3_var)
        self.pts_x3.grid(column=9, row=5)
        self.ets_x3 = ttk.Checkbutton(self.frame6)
        self.ets_x3_var = tk.BooleanVar(value=0)
        self.ets_x3.configure(variable=self.ets_x3_var)
        self.ets_x3.grid(column=9, row=6)
        self.mdef_x3 = ttk.Checkbutton(self.frame6)
        self.mdef_x3_var = tk.BooleanVar(value=0)
        self.mdef_x3.configure(variable=self.mdef_x3_var)
        self.mdef_x3.grid(column=9, row=7)
        self.dmin_x3 = ttk.Checkbutton(self.frame6)
        self.dmin_x3_var = tk.BooleanVar(value=0)
        self.dmin_x3.configure(variable=self.dmin_x3_var)
        self.dmin_x3.grid(column=9, row=8)
        self.fmax_x3 = ttk.Checkbutton(self.frame6)
        self.fmax_x3_var = tk.BooleanVar(value=0)
        self.fmax_x3.configure(variable=self.fmax_x3_var)
        self.fmax_x3.grid(column=9, row=9)
        self.tcd_x3 = ttk.Checkbutton(self.frame6)
        self.tcd_x3_var = tk.BooleanVar(value=0)
        self.tcd_x3.configure(variable=self.tcd_x3_var)
        self.tcd_x3.grid(column=9, row=10)
        self.vts_x3 = ttk.Checkbutton(self.frame6)
        self.vts_x3_var = tk.BooleanVar(value=0)
        self.vts_x3.configure(variable=self.vts_x3_var)
        self.vts_x3.grid(column=9, row=11)
        self.label87 = ttk.Label(self.frame6)
        self.label87.configure(text="x")
        self.label87.grid(column=5, row=1)
        self.label88 = ttk.Label(self.frame6)
        self.label88.configure(text="y")
        self.label88.grid(column=6, row=1)
        self.label89 = ttk.Label(self.frame6)
        self.label89.configure(text="x")
        self.label89.grid(column=9, row=1)
        self.label90 = ttk.Label(self.frame6)
        self.label90.configure(text="y")
        self.label90.grid(column=10, row=1)
        self.label91 = ttk.Label(self.frame6)
        self.label91.configure(foreground="#0000ff", text="Plot 2")
        self.label91.grid(column=5, columnspan=2, row=0)
        self.label92 = ttk.Label(self.frame6)
        self.label92.configure(foreground="#0000ff", text="Plot 3")
        self.label92.grid(column=9, columnspan=2, row=0)
        self.separator40 = ttk.Separator(self.frame6)
        self.separator40.configure(orient="vertical")
        self.separator40.grid(column=1, ipady=160, row=0, rowspan=13, sticky="ns")
        self.separator42 = ttk.Separator(self.frame6)
        self.separator42.configure(orient="vertical")
        self.separator42.grid(column=4, ipady=160, row=0, rowspan=13, sticky="ns")
        self.separator43 = ttk.Separator(self.frame6)
        self.separator43.configure(orient="vertical")
        self.separator43.grid(column=7, ipady=160, row=0, rowspan=13, sticky="ns")
        self.separator3 = ttk.Separator(self.frame6)
        self.separator3.configure(orient="vertical")
        self.separator3.grid(column=11, ipady=160, row=0, rowspan=13, sticky="ns")
        self.labelcusmtaImax = ttk.Label(self.frame6)
        self.labelcusmtaImax.configure(text="Maximum Indentation")
        self.labelcusmtaImax.grid(column=0, row=12)
        self.imax_x1 = ttk.Checkbutton(self.frame6)
        self.imax_x1_var = tk.BooleanVar(value=0)
        self.imax_x1.configure(variable=self.imax_x1_var)
        self.imax_x1.grid(column=2, row=12)
        self.imax_y1 = ttk.Checkbutton(self.frame6)
        self.imax_y1_var = tk.BooleanVar(value=0)
        self.imax_y1.configure(variable=self.imax_y1_var)
        self.imax_y1.grid(column=3, row=12)
        self.imax_x2 = ttk.Checkbutton(self.frame6)
        self.imax_x2_var = tk.BooleanVar(value=0)
        self.imax_x2.configure(variable=self.imax_x2_var)
        self.imax_x2.grid(column=5, row=12)
        self.imax_y2 = ttk.Checkbutton(self.frame6)
        self.imax_y2_var = tk.BooleanVar(value=0)
        self.imax_y2.configure(variable=self.imax_y2_var)
        self.imax_y2.grid(column=6, row=12)
        self.imax_x3 = ttk.Checkbutton(self.frame6)
        self.imax_x3_var = tk.BooleanVar(value=0)
        self.imax_x3.configure(variable=self.imax_x3_var)
        self.imax_x3.grid(column=9, row=12)
        self.imax_y3 = ttk.Checkbutton(self.frame6)
        self.imax_y3_var = tk.BooleanVar(value=0)
        self.imax_y3.configure(variable=self.imax_y3_var)
        self.imax_y3.grid(column=10, row=12)
        self.frame6.configure(height=200, width=200)
        self.frame6.grid(column=0, row=0)
        self.frame6.rowconfigure("all", weight=1)
        self.frame6.columnconfigure("all", minsize=50, weight=1)
        self.notebook4.add(self.frame6, text="Custom graph")
        self.notebook4.pack(expand="true", fill="both", side="top")
        self.toplevel1.configure(height=200, width=200)

        # Main widget
        self.mainwindow = self.toplevel1

    def run(self):
        self.mainwindow.mainloop()

    def on_run_ex_clicked(self):
        pass

    def on_abort_clicked(self):
        pass

    def on_plot_ex_clicked(self):
        pass

    def on_examples_ex_clicked(self):
        pass

    def on_about_clicked(self):
        pass

    def gtk_main_quit(self):
        pass

    def single_mode_activated(self):
        pass

    def bimodal_mode_activated(self):
        pass

    def point_mass_activated(self):
        pass

    def euler_activated(self):
        pass

    def on_vdw_activated(self):
        pass
    
    def on_vdw_cone_activated(self):
        pass
    
    def on_vdw_punch_activated(self):
        pass

    def on_dlvo_activated(self):
        pass

    def on_hertz_activated(self):
        pass

    def on_dmt_activated(self):
        pass

    def on_jkr_activated(self):
        pass

    def on_becc_activated(self):
        pass

    def on_bepc_activated(self):
        pass

    def on_tatara_activated(self):
        pass

    def on_viscosity_activated(self):
        pass

    def on_adhesion_activated(self):
        pass

    def on_custom1_activated(self):
        pass

    def on_custom2_activated(self):
        pass

    def on_show_value_clicked(self):
        pass

    def on_dipole_activated(self):
        pass

    def on_exp_activated(self):
        pass

    def on_sneddon_activated(self):
        pass

    def on_kv_sphere_activated(self):
        pass

    def on_kv_cone_activated(self):
        pass

    def on_punch_activated(self):
        pass

    def on_viscosity_cone_activated(self):
        pass

    def on_viscosity_punch_activated(self):
        pass

    def on_kv_punch_activated(self):
        pass

    def on_behc_activated(self):
        pass

    def on_lennardJones_activated(self):
        pass

    def on_benc_activated(self):
        pass

    def on_nanowire_activated(self):
        pass

    def on_DELR_activated(self):
        pass

    def on_viscosity_nanowire_activated(self):
        pass

    def on_kv_nanowire_activated(self):
        pass

    def on_calculate_simulation_parameters(self):
        pass


if __name__ == "__main__":
    app = GuiDforceApp()
    app.run()
