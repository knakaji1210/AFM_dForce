
import tkinter as tk
import tkinter.ttk as ttk
import webbrowser

class credits_app():
    def __init__(self, master):         
        self.master = master           
        self.credits_menu = tk.Toplevel(master)
        self.credits_menu.wm_title("Credits")
        
        self.img = tk.PhotoImage(file="icons/icon.png")
        self.credits_menu.wm_iconphoto(False, self.img)
        
        
        self.frame = ttk.Frame(self.credits_menu)
        
        
        self.label1 = ttk.Label(self.frame)        
        self.img = tk.PhotoImage(file="icons/dForceLogo4-1.png")
        self.label1.configure(image=self.img, text="label2")
        self.label1.pack()
        self.label2 = ttk.Label(self.frame,
                                font='Helvetica 13 bold',                               
                                text="dForce 2.0")
        self.label2.pack()
        
        self.label3 = ttk.Label(self.frame,
                                text="Victor G. Gisbert & Ricardo Garcia",
                                padding = 10)
        self.label3.pack()
        
        
        self.label4 = ttk.Label(self.frame,
                                font='Helvetica 13 bold',
                                text="dForce 1.0")
        self.label4.pack()
        
        self.label5 = ttk.Label(self.frame,
                                text="Horacio Vargas, Pablo D. Garcia & Ricardo Garcia",
                                padding = 10)
        self.label5.pack()

        
        def callback(event):
            webbrowser.open_new(event.widget.cget("text"))
        
        
        lbl = ttk.Label(self.frame, 
                       text=r"http://www.icmm.csic.es/forcetool/", 
                       foreground="blue",
                       cursor="hand2",
                       padding = 10
                       )
        lbl.pack()
        lbl.bind("<Button-1>", callback)
        
        
        self.frame.configure(height=200, width=200)
        self.frame.rowconfigure("all", pad=20)
        self.frame.columnconfigure("all", pad=20)
        self.frame.pack(side="top")
        
        self.button = ttk.Button(self.credits_menu)
        self.button.configure(text="OK")
        self.button.pack(side="top")
        self.button.configure(command=self.credits_menu_exit)

        self.mainwindow = self.credits_menu

    def run(self):
        self.mainwindow.mainloop()

    def credits_menu_exit(self):
        self.credits_menu.destroy()