# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 12:27:41 2023

@author: user
"""



#!/usr/bin/python3
import tkinter as tk
import tkinter.ttk as ttk

from tkinter import font

# Widget which is a checkbutton where you can put any subscript
# It is a mixture of a checkbuttton widget and a text widget which 
# Can take more rich formatting
class CheckbuttonWithFormula(tk.Frame):
    def __init__(self, parent, label, var1, width = 30):
        
        tk.Frame.__init__(self, parent)
       
        self.var1 = var1
        #l1 = tk.Label(master,text="Hello")
        #fontname = font.nametofont(l1['font']).configure()["family"]
        #fontsize = font.nametofont(l1['font']).configure()["size"] 
        fontname,fontsize = "Segoe UI" , 9
        
        #test_font = font.Font(family=fontname, size=fontsize)
        
        #width = int(test_font.measure("Amplitude (zc)")/4 )
        
        
        text_entry_conf = {"height":1.5,
                           "width":width,
                           "borderwidth":0,
                           "font":(fontname,fontsize),
                           "background" : parent.winfo_toplevel().cget("background")} 
        
        
    
        fontsize_sub = int(fontsize*6.0/8.0)
        normal_conf = {"tagName":"normal", "offset":-7,"font":(fontname,fontsize) }
        subscript_conf = {"tagName":"subscript", "offset":-7-1,"font":(fontname,fontsize_sub) }
        superscript_conf = {"tagName":"superscript", "offset":-7+4,"font":(fontname,fontsize_sub) }

        
        self.checkbutton = ttk.Checkbutton(self,variable=self.var1, onvalue=1, offvalue=0)
        self.checkbutton.pack(side="left")       
        
             
               
        self.label = tk.Text(self,**text_entry_conf)
        self.label.bindtags((str(self.label), str(self), "all"))
        self.label.tag_configure(**normal_conf)
        self.label.tag_configure(**subscript_conf)
        self.label.tag_configure(**superscript_conf)
        self.label.insert("insert", label +" (z", "normal", "c", "subscript",")","normal")
        
        self.label.pack(side="left")    
        
        
    def disable(self):
        self.label["state"] = "disabled"
        self.checkbutton["state"] = "disabled"
    def enable(self):
        self.label["state"] = "normal"
        self.checkbutton["state"] = "normal"
        
            
# For testing
class CustomWidget(tk.Frame):
    def __init__(self, parent, label, default=""):
        tk.Frame.__init__(self, parent)

        self.label = tk.Label(self, text=label, anchor="w")
        self.entry = tk.Entry(self)
        self.entry.insert(0, default)

        self.label.pack(side="top", fill="x")
        self.entry.pack(side="bottom", fill="x", padx=4)
    


    def get(self):
        return self.entry.get()
    

        

# For testing
class NewProjectApp:
    def __init__(self, master=None):
        # build ui
        self.mainwindow = ttk.Frame(master)

        # Main widget
        self.ampst_cb_var = tk.IntVar(value=0)
        
        ch = CheckbuttonWithFormula(master,"Amplitude",self.ampst_cb_var)
        ch.pack(side = "top")
        
        
        
        ch = CustomWidget(master,"text")
        ch.pack(side = "top")
        
    def run(self):
        self.mainwindow.mainloop()



    def run(self):
        self.mainwindow.mainloop()

if __name__ == "__main__":
    root = tk.Tk()
    app = NewProjectApp(root)
    app.run()
    

