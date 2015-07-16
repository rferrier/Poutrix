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
from solveur_poutres_v6 import *
from parser_pout_v2_2 import parser
from assembly_v3_1 import *
    
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
        
        Label(self,text="Choisir le fichier à ouvrir").pack(side=TOP)
        Button(self,text="actualiser avec le filtre txt",command = lambda txt=1 : self.actxt(txt,self.curDir)).pack()
        Button(self,text="actualiser sans le filtre txt",command = lambda txt=0 : self.actxt(txt,self.curDir)).pack()
        Button(self,text="retour au dossier de base",command = self.back).pack()

        #current directory :
        self.curDir = os.path.dirname(os.path.abspath(__file__))
        self.actxt(0,self.curDir)

    def start(self,chcar):
        try :
            #trying if it is an openable folder
            self.curDir1 = self.curDir + "\\" + chcar
            self.actxt(0,self.curDir1)
            self.curDir = self.curDir1
            
        except :
            chcar1 = self.curDir + "\\" + chcar
            self.reader = Reader(chcar1)
            
    def back(self):
        self.curDir = os.path.dirname(os.path.abspath(__file__))
        self.actxt(0,self.curDir)
    
    def actxt(self,txt,chcar):
        curDir = chcar
        listDir = os.listdir(curDir)
        if txt == 0:
            self.listDir = listDir
        else :
            self.listDir = []
            for i in listDir:
                if len(i) > 3:
                    if i[-4]+i[-3]+i[-2]+i[-1] == ".txt":
                        self.listDir.append(i)

        try :
            self.ButtonZone.destroy()
        except :
            pass
        #rebuilding of the buttons
        self.ButtonZone = ScrolledCanvas(self,500,300,(-200,0,200,40*len(self.listDir)))
        self.ButtonZone.pack(expand=YES,fill=BOTH,padx=6,pady=6)
        self.can = self.ButtonZone.can
        
        #buttons
        self.bu = []
        self.fb = []

        for i in range(len(self.listDir)):
            chcar = self.listDir[i]

            self.bu.append(Button(text=chcar,command = lambda chcar=chcar : self.start(chcar)))
            self.fb.append(self.can.create_window(0,30*i+20,window=self.bu[i]))

A = MainWindow()
A.title('  Poutrix')
#A.bind("<Button-1>", select)
A.mainloop()
#B = Reader("assembly2.txt")