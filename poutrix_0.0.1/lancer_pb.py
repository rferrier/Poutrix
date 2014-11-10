# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 11:42:09 2014

@author: renaud
"""

#This code provides the graphical interface for Poutrix

#external importations :
import numpy as np
from scipy.linalg import *
import sp
import math as mh
from tkinter import *
import os as os

#internal importation :
from solveur_poutres_v4_4 import *
from parser_pout_v2 import parser
from assembly_v2 import *
    
class ScrolledCanvas(Frame):
    
    def __init__(self,parent,width,height,scrollregion):
        #I don't understand anything I am writing : thank you, Gerard Swinnen
        Frame.__init__(self,bd=2,relief=SUNKEN)
        self.can = Canvas(self,width=width-20,height=height-20,bg="white",scrollregion=scrollregion,bd=1)
        self.can.grid(row=0,column=0)
        vscrol = Scrollbar(self,orient=VERTICAL,command=self.can.yview,bd=1)
        hscrol = Scrollbar(self,orient=HORIZONTAL,command=self.can.xview,bd=1)
        self.can.configure(xscrollcommand=hscrol.set,yscrollcommand=vscrol.set)
        vscrol.grid(row=0,column=1,sticky=NS)
        hscrol.grid(row=1,column=0,sticky=EW)
        
        self.bind("<Configure>",self.redim)
        self.started=False
        
    def redim(self,event):
        if self.started == False:
            self.started = True
            return
        
        self.can.config(width=self.winfo_width()-30,height=self.winfo_height()-30)
        
class MainWindow(Tk):
    def __init__(self):
        Tk.__init__(self)
        
        #current directory :
        curDir = os.path.dirname(os.path.abspath(__file__))
        listDir = os.listdir(curDir)
        
        ButtonZone = ScrolledCanvas(self,500,300,(-200,0,200,20*len(curDir)))
        ButtonZone.pack(expand=YES,fill=BOTH,padx=6,pady=6)
        self.can = ButtonZone.can
        
        #buttons
        self.bu = []
        self.fb = []

        for i in range(len(listDir)):
            chcar = listDir[i]

            self.bu.append(Button(text=chcar,command = lambda chcar=chcar : self.start(chcar)))
            self.fb.append(self.can.create_window(0,30*i+20,window=self.bu[i]))

    def start(self,chcar):

        self.reader = Reader(chcar)

A = MainWindow()
A.title('  Poutrix')
#A.bind("<Button-1>", select)
A.mainloop()
#B = Reader("assembly2.txt")