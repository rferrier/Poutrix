#This code provides the finite elements solver for assemblies of beams

#Internal importations
import numpy as np
from scipy.linalg import *
from solveur_poutres_v4_4 import *

#External importations
import sp
import math as mh
from parser_pout_v2 import parser
import os as os

class Reader :
    #assembly data-process facility

    def __init__(self,data):
        
        #Recovering of the data
        self.title = data
        self.read_input(data)
        #print(len(self.input[0]))

        #distributing lists to readers
        self.link = []
        self.beam = []
        self.beamCut = [] #it will contain the points where the beam will be cut
        #print(self.beamCut)        
        self.beamName = {} #number : 'name' of the beam
        self.beamName1 = {} #the invert of beamName

        for i in range(len(self.input[0])) :
            #pre-processing of the beams
            self.beamCut.append([])
            self.beamName[self.input[0][i][0]] = i  
            self.beamName1[i] = self.input[0][i][0]
        #print(self.beamName)
        for i in range(len(self.input[1])) :
            #pre-processing of the links in order to mesh the beams
            #every beam becomes its links points (in order to place a node)
            self.beamCut[self.beamName[self.input[1][i][0]]].append(self.input[1][i][1])
            self.beamCut[self.beamName[self.input[1][i][2]]].append(self.input[1][i][3])
        
        for i in range(len(self.input[0])) :
            self.beam.append(BeamR(self.input[0][i],self))
            
        for i in range(len(self.input[1])) :
            self.link.append(LinkR(self.input[1][i],self))
            
        self.solr = SolvR(self.input[2],self)
        
        self.build_output()
        self.solve(self.assemlist)
            
    def read_input(self,data):
        #Reading of the imput file
        self.f = open(data,'r')
        self.text = ''

        for i in self.f:
            self.text = self.text + i

        self.f.close()

        #calling of the parser
        self.translate = parser()

        try :
            self.input = self.translate(self.text)
        except SyntaxError as error :
            print(" /!\ Error at[ligne,colonne]")
            print('    ',error)

        #print(self.input)

    def build_output(self):
      
        self.assemlist = []
        self.beamlist = []
        self.linklist = []
        
        for i in self.beam:
            i.build_output()
            self.beamlist.append(i.output1)
            
        for i in self.link :
            i.build_output()
            self.linklist.append(i.output)
        
        self.assemlist = [self.beamlist,self.linklist]
        #print(self.assemlist)
        #launching the calculus...

    def solve(self,data):

        #Launching of the solver
        if self.solr.type == 'static':
            #print(data)
            self.solv = Assembly_solv(data)
            self.solv.static_solv()

            #recovering of results
            self.u = self.solv.u
            #processing of results
            self.v = [0,0,0,0,0,0]
            #extracting theta and v from u
            self.ulen = int(len(self.u)/6)
            #print(self.ulen)
            for i in range(6): #empty lists
                self.v[i] = [0]*self.ulen
                for j in range(self.ulen):
                    self.v[i][j] = self.u[6*j+i][0]
                    
            #distributing v to beams :
            j = 0
            for i in self.beam:
                i.v = [0,0,0,0,0,0]
                for k in range(6):
                    i.v[k] = (self.v[k][j:j+len(i.elements)+1])
                j = j + len(i.elements) + 1
                #print(i.v)
            
            self.outputResults()

        if self.solr.type == 'modal':
            self.solv = Solver(data)
            self.solv.solve_modal_problem()
    
    def outputResults(self):
        #Building the output file
        self.outChaine = 'Resultats pour la simulation : ' + str(self.title) + '\n\nLes déplacements sont donnés dans la base globale\n\n'
        microlist = ['u_x','u_y','u_z','theta_x','theta_y','theta_z']
        self.outName = str(self.solr.outName) + '.txt'
        self.f = open(self.outName,'w')
        #len(self.beam[i].elements)
        alpha = 0
        for k in self.beam:

            self.outChaine = self.outChaine + '\n' + str(self.beamName1[alpha]) +'\n'
            for i in range(len(k.v)):
                self.outChaine = self.outChaine + '\n' + microlist[i] + '\n\n'
                for j in k.v[i]:
                    self.outChaine = self.outChaine + str(j) + '\n'  
            alpha = alpha + 1

        self.f.write(self.outChaine)
        self.f.close()
        os.startfile(self.outName)


class LinkR :
    #Link Reader

    def __init__(self,data,parent):
        
        self.parent = parent
        self.input = data
        #data processing
        self.leftB = parent.beamName[data[0]]
        self.rightB = parent.beamName[data[2]]
        
        simeq = False
        for i in parent.beam[self.leftB].ltelist:
            if data[1] < i+0.000001 and data[1] > i-0.000001:
                iretenu = i
                simeq = True
        if simeq == True:
            self.leftX = parent.beam[self.leftB].lte[iretenu]
            
        simeq = False
        for i in parent.beam[self.rightB].ltelist:
            if data[3] < i+0.000001 and data[3] > i-0.000001:
                iretenu = i
                simeq = True            
        if simeq == True:
            self.rightX = parent.beam[self.rightB].lte[iretenu]
            
        self.type = 'enca'
        self.repere = 'global'
        
    def build_output(self):
        self.output = [self.type,[self.leftB,self.leftX],[self.rightB,self.rightX],self.repere]

class BeamR :
    #Beam reader

    def __init__(self,data,parent):

        self.input = data
        #print(self.input)
        #data processing
        self.name = self.input[0]

        #characteristics
        self.leng = mh.sqrt((self.input[1][0][0]-self.input[1][1][0])**2 + (self.input[1][0][1]-self.input[1][1][1])**2 + (self.input[1][0][2]-self.input[1][1][2])**2)
        self.E = self.input[1][2]
        self.nu = self.input[1][3]
        #self.G = self.E/(2*(1+self.nu))
        self.S = self.input[1][4]
        self.Iy = self.input[1][5]
        self.Iz = self.input[1][6]
        #self.I = self.Iz + self.Iy
        self.rho = self.input[1][7]
        #self.mu = self.rho*self.S
        #self.j = self.rho*self.I
        
        #Links
        self.cut = parent.beamCut[parent.beamName[self.name]]
        self.points = self.cut[:]
        self.points.sort()
        #adding beginning and end
        if self.points == []:
            self.points.insert(0,0.0)
            self.points.append(self.leng)
        else :
            if self.points[0] > 0.001:
                self.points.insert(0,0.0)
            if self.points[-1] < self.leng - 0.001:
                self.points.append(self.leng)

        #bundaries
        self.uxg = self.input[2][0][0]
        self.uyg = self.input[2][0][1]
        self.uzg = self.input[2][0][2]
        self.txg = self.input[2][0][3]
        self.tyg = self.input[2][0][4]
        self.tzg = self.input[2][0][5]

        self.uxd = self.input[2][1][0]
        self.uyd = self.input[2][1][1]
        self.uzd = self.input[2][1][2]
        self.txd = self.input[2][1][3]
        self.tyd = self.input[2][1][4]
        self.tzd = self.input[2][1][5]

        self.fxg = self.input[3][0][0]
        self.fyg = self.input[3][0][1]
        self.fzg = self.input[3][0][2]
        self.mxg = self.input[3][0][3]
        self.myg = self.input[3][0][4]
        self.mzg = self.input[3][0][5]

        self.fxd = self.input[3][1][0]
        self.fyd = self.input[3][1][1]
        self.fzd = self.input[3][1][2]
        self.mxd = self.input[3][1][3]
        self.myd = self.input[3][1][4]
        self.mzd = self.input[3][1][5]

        #lineic efforts
        self.fx = self.input[3][2][0]
        self.fy = self.input[3][2][1]
        self.fz = self.input[3][2][2]
        self.mx = self.input[3][2][3]
        self.my = self.input[3][2][4]
        self.mz = self.input[3][2][5]

        self.list_char_mieux1 = [0]*6
        self.list_char = [self.fx,self.fy,self.fz,self.mx,self.my,self.mz]

        #transforming 1 or 0 elements liste into numbers

        for i in range(len(self.list_char)):
            if len(self.list_char[i]) == 0:
                self.list_char_mieux1[i] = 0.
                
            else :
                self.list_char_mieux1[i] = self.list_char[i][0]
        
        #Type of mesh
        self.elem_number = self.input[1][8]

        #Initialization of the elements
        self.elements = []

        #Processing of the variable characteristics
        self.list_char = [self.E,self.nu,self.S,self.Iy,self.Iz,self.rho,self.list_char_mieux1[0],self.list_char_mieux1[1],self.list_char_mieux1[2],self.list_char_mieux1[3],self.list_char_mieux1[4],self.list_char_mieux1[5]] 
        self.list_char_mieux2 = [0]*len(self.list_char)
        self.j = 0
        
        for i in range(len(self.list_char)) :
            #looking for a changing characteristic
            if type(self.list_char[i]) == str:#found !
                #building of the name of the file to open
                self.name1 = self.list_char[i] + '.dat'

                #building the mesh
                if self.j != 1:
                    self.copy_mesh(self.name1)
                    self.j = 1
                    #print(self.j)

                #writing the characteristic for every node
                self.list_char_mieux2[i] = []
                self.k = 1

                self.f = open(self.name1,'r')
                
                for n in self.f :
                #Taking of only one data for two
                    if self.k == 1:
                        self.k = 0
                        self.list_char_mieux2[i].append(float(n))
                    else:
                        self.k = 1

        if self.j == 0 :
            self.simple_mesh()

        for i in range(len(self.list_char)) :
            #Building of the list of constant char
            if type(self.list_char[i]) == float:
                self.list_char_mieux2[i] = [self.list_char[i]]*(len(self.elements) + 1)

        #print(self.list_char_mieux2) #from now, this list will be used
        #Building the lis of the no of nodes with links
        self.linkToElement()

    def simple_mesh(self):
        #A simple 1D meshing algorithm
        linit = 0.
        for j in self.points[1:] :   #cutting the beam between links
            self.elem_leng = (j-linit)/self.elem_number
            for i in range(self.elem_number):
                self.elements.append(self.elem_leng)
            linit = i

    def copy_mesh(self,name):
        #copying a mesh from a data file
        self.f = open(name,'r')

        self.k = 0
        for n in self.f:

            #Taking of only one data for two
            if self.k == 1:
                self.k = 0
                self.elements.append(float(n))
            else:
                self.k = 1

        self.f.close

        #Norming of the list
        self.len_befor_norm = sum(self.elements)

        for n in range(len(self.elements)) :
            self.elements[n] = self.elements[n]*self.leng/self.len_befor_norm

    def linkToElement(self):
        #building the dico containing the no of nodes and abscissa with links
        self.lte = {}
        self.ltelist = []
        for i in range(len(self.elements)+1) :
            alpha = sum(self.elements[:i])   #sum of len of elements before i = abscissa of link node of elem i
            #print(alpha)
            simeq = False
            #check if alpha ~= an elem of self.cut
            for j in self.cut : 
                if alpha < j+0.000001 and alpha > j-0.000001 :

                    simeq = True
            if simeq == True :
                self.lte[alpha] = i
                self.ltelist.append(alpha)
        #print(self.lte)
        #print(self.ltelist)

    def build_output(self):

        self.output = []

        #left boundaries
        self.output.append([0.]*6)
        self.list_char = [self.uxg,self.uyg,self.uzg,self.txg,self.tyg,self.tzg]
        for i in range(len(self.list_char)) :
            if len(self.list_char[i]) == 0:
                #free
                self.output[0][i] = 1
            else :
                self.output[0][i] = self.list_char[i][0]

        #nodes
        for i in range (len(self.list_char_mieux2[0])):
            self.output.append([[0.]*6,[0.]*7])

            self.output[i+1][0][0] = self.list_char_mieux2[6][i]#Fx
            self.output[i+1][0][1] = self.list_char_mieux2[7][i]#Fy
            self.output[i+1][0][2] = self.list_char_mieux2[8][i]#Fz
            self.output[i+1][0][3] = self.list_char_mieux2[9][i]#Mx
            self.output[i+1][0][4] = self.list_char_mieux2[10][i]#My
            self.output[i+1][0][5] = self.list_char_mieux2[11][i]#Mz
            
            if i == len(self.elements): #writing of the last and useless h
                self.output[i+1][1][0] = 1.#h
            else :
                self.output[i+1][1][0] = self.elements[i]#h

            self.output[i+1][1][1] = self.list_char_mieux2[5][i]*self.list_char_mieux2[2][i]#\mu = \rho*S
            self.output[i+1][1][2] = self.list_char_mieux2[0][i]*self.list_char_mieux2[2][i]#ES
            self.output[i+1][1][3] = self.list_char_mieux2[0][i]/(2*(1+self.list_char_mieux2[1][i]))*(self.list_char_mieux2[3][i] + self.list_char_mieux2[4][i])#GJ
            self.output[i+1][1][4] = self.list_char_mieux2[0][i]*self.list_char_mieux2[3][i]#EIy
            self.output[i+1][1][5] = self.list_char_mieux2[0][i]*self.list_char_mieux2[4][i]#EIz
            self.output[i+1][1][6] = self.list_char_mieux2[5][i]*(self.list_char_mieux2[3][i] + self.list_char_mieux2[4][i])#\rho*J
        #Local loading on the beam

        self.list_char_mieux3 = [0]*12
        self.list_char = [self.fxg,self.fyg,self.fzg,self.mxg,self.myg,self.mzg,self.fxd,self.fyd,self.fzd,self.mxd,self.myd,self.mzd]

        #transforming 1 or 0 elements liste into numbers

        for i in range(len(self.list_char)):
            if len(self.list_char[i]) == 0:
                self.list_char_mieux3[i] = 0.
                
            else :
                self.list_char_mieux3[i] = self.list_char[i][0]
                
        self.output[1][0][0] = self.output[1][0][0] + self.list_char_mieux3[0]#Fx
        self.output[1][0][1] = self.output[1][0][1] + self.list_char_mieux3[1]#Fy
        self.output[1][0][2] = self.output[1][0][2] + self.list_char_mieux3[2]#Fz
        self.output[1][0][3] = self.output[1][0][3] + self.list_char_mieux3[3]#Mx
        self.output[1][0][4] = self.output[1][0][4] + self.list_char_mieux3[4]#My
        self.output[1][0][5] = self.output[1][0][5] + self.list_char_mieux3[5]#Mz

        self.output[len(self.list_char_mieux2[0])][0][0] = self.output[len(self.list_char_mieux2[0])][0][0] + self.list_char_mieux3[6]#Fx
        self.output[len(self.list_char_mieux2[0])][0][1] = self.output[len(self.list_char_mieux2[0])][0][1] + self.list_char_mieux3[7]#Fy
        self.output[len(self.list_char_mieux2[0])][0][2] = self.output[len(self.list_char_mieux2[0])][0][2] + self.list_char_mieux3[8]#Fz
        self.output[len(self.list_char_mieux2[0])][0][3] = self.output[len(self.list_char_mieux2[0])][0][3] + self.list_char_mieux3[9]#Mx
        self.output[len(self.list_char_mieux2[0])][0][4] = self.output[len(self.list_char_mieux2[0])][0][4] + self.list_char_mieux3[10]#My
        self.output[len(self.list_char_mieux2[0])][0][5] = self.output[len(self.list_char_mieux2[0])][0][5] + self.list_char_mieux3[11]#Mz
        
        #right boundaries
        self.output.append([0.]*6)
        self.list_char = [self.uxd,self.uyd,self.uzd,self.txd,self.tyd,self.tzd]
        #print(self.list_char)
        for i in range(len(self.list_char)) :
            if len(self.list_char[i]) == 0:
                #free
                self.output[len(self.elements)+2][i] = 1
            else :
                self.output[len(self.elements)+2][i] = self.list_char[i][0]
        
        #dimensions of the beam :
        self.output1 = [[self.input[1][0][0],self.input[1][0][1],self.input[1][0][2],self.input[1][1][0],self.input[1][1][1],self.input[1][1][2]],self.output]

        #print(self.output1)  

class SolvR :
    #The reader for solver data

    def __init__(self,data,parent):
        self.type = data[0]
        self.outName = data[1]

class Assembly_solv :

    def __init__ (self,data):
        #recovery data
        self.beams = data[0]
        self.linksd = data[1]
        self.microsolv = []
        self.links = []
        
        self.create_microsolv()
        self.create_links()

        #computation of the theta

        self.theta1 = [0]*len(self.microsolv)
        self.theta2 = [0]*len(self.microsolv)
        for i in range(len(self.beams)):
            self.compute_angle(i)

    def create_microsolv (self):
        #creation of the local solvers for each beam
        for i in self.beams:
            self.microsolv.append(Solver(i))

    def create_links (self):
        #building of the links
        for i in self.linksd:
            #computation of the index of the two points linked :
            self.leftl = 0
            self.rightl = 0
            
            for j in range(i[1][0]):
                self.leftl = self.leftl + len(self.beams[j][1])-2

            self.leftl = self.leftl + i[1][1]

            for j in range(i[2][0]):
                self.rightl = self.rightl + len(self.beams[j][1])-2
                
            self.rightl = self.rightl + i[2][1]
            
            if i[0] == 'rot':
                self.links.append(Rotule([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self))

            if i[0] == 'enca':
                self.links.append(Encas([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self))

    def compute_mat(self):

        for i in self.microsolv:
            i.create_mesh()
            i.mesh.compute_leng()
            i.mesh.compute_A()
            i.mesh.compute_M()
            i.mesh.compute_b()

    def compute_angle(self,i):

        #computation of theta1 and theta2 for one beam
        self.x = self.beams[i][0][3] - self.beams[i][0][0]
        self.y = self.beams[i][0][4] - self.beams[i][0][1]
        self.z = self.beams[i][0][5] - self.beams[i][0][2]

        #avoiding of 1/0
        if self.x == 0.:
            if self.y > 0.:
                self.theta1[i] = pi/2
            if self.y < 0.:
                self.theta1[i] = -pi/2
            if self.y == 0.:
                if self.z > 0.:
                    self.theta2[i] = -pi/2
                else:
                    self.theta2[i] = pi/2

            else :
                self.theta2[i] = - mh.atan(self.z/mh.sqrt(self.x*self.x+self.y*self.y))

        else :
            self.theta1[i] = mh.atan(self.y/self.x)
            self.theta2[i] = - mh.atan(self.z/mh.sqrt(self.x*self.x+self.y*self.y))
            
        #print(self.x,self.y,self.z,self.theta1,self.theta2)

    def transform_mat (self):       
        #rotating of the matrix of a beam

        self.pas = [0.]*len(self.microsolv)
        self.pas1 = [0.]*len(self.microsolv)

        #initilaisation of the list of rigidity matrix befor assembly
        self.Mlist = [0.]*len(self.microsolv)
        self.Alist = [0.]*len(self.microsolv)
        self.blist = [0.]*len(self.microsolv)
        #computing of the elementary passage matrix
        for i in range(len(self.microsolv)) :
            for j in self.microsolv[i].mesh.elements:
                j.compute_P(self.theta1[i],self.theta2[i])

            #assembly of P in order to get the passage matrix for the beam

            self.pas[i] = np.mat([[0.]*self.microsolv[i].mesh.size]*self.microsolv[i].mesh.size)

            s = 0
            for ii in range(len(self.microsolv[i].mesh.elements)):
                for jj in range(len(self.microsolv[i].mesh.elements[ii].pas)):
                    for kk in range(len(self.microsolv[i].mesh.elements[ii].pas)):
                        self.pas[i][s+jj,s+kk] = self.pas[i][s+jj,s+kk] + self.microsolv[i].mesh.elements[ii].pas[jj,kk]

                #incrementation of the place for the pas
                s = s + len(self.microsolv[i].mesh.elements[ii].pas)-6

            #inverting the elementary matrix
            self.pas1[i] = np.transpose(self.pas[i])
            #print(self.pas[i])

            #writing of the rotation
            self.Mlist[i] = self.pas1[i] * self.microsolv[i].mesh.M * self.pas[i]
            self.Alist[i] = self.pas1[i] * self.microsolv[i].mesh.A * self.pas[i]
            self.blist[i] = self.microsolv[i].mesh.b

    def static_solv(self):
        self.compute_mat()
        self.compute_leng()

        self.transform_mat()
        self.assembly_A()
        self.compute_pen()
        for i in self.links:
            i.compute_matrix()
        self.assembly_b()
        self.u = solve(self.A,self.b)
        
    def modal_solv(self):
        self.compute_mat()
        self.compute_leng()
        
        self.transform_mat()
        self.assembly_A()
        self.assembly_M()
        self.compute_pen()
        for i in self.links:
            i.compute_matrix()
            
        self.ms = eig(self.A,self.M)
        self.omeg = self.ms[0]
        self.u = self.ms[1]

    def compute_leng(self):
        #computation of the size of the matrix
        self.leng = 0
        for i in self.microsolv:
            self.leng = self.leng + len(i.mesh.A)

    def compute_pen(self):
        #computation of the penalization matrix
        self.penlist = []
        for i in self.microsolv:
            self.penlist.append(i.mesh.elements[0].pen)
            self.penlist.append(i.mesh.elements[-1].pen)

        self.pen = max(self.penlist)

    def assembly_A(self):

        self.A = np.mat([[0.]*self.leng]*self.leng)
        s = 0
        for i in range(len(self.microsolv)):
            for ii in range(len(self.microsolv[i].mesh.A)):
                for jj in range(len(self.microsolv[i].mesh.A)):
                    self.A[s+ii,s+jj] = self.Alist[i][ii,jj]

            s = s + len(self.microsolv[i].mesh.A)
            
    def assembly_M(self):

        self.M = np.mat([[0.]*self.leng]*self.leng)
        s = 0
        for i in range(len(self.microsolv)):
            for ii in range(len(self.microsolv[i].mesh.M)):
                for jj in range(len(self.microsolv[i].mesh.M)):
                    self.M[s+ii,s+jj] = self.Mlist[i][ii,jj]

            s = s + len(self.microsolv[i].mesh.M)

    def assembly_b(self):

        self.b = np.mat([[0.]]*self.leng)
        s = 0
        for i in range(len(self.microsolv)):
            for ii in range(len(self.microsolv[i].mesh.b)):
                self.b[s+ii,0] = self.blist[i][ii,0]

            s = s + len(self.microsolv[i].mesh.b)

class Link:
    #Basic class for links

    def __init__(self,data,parent):
        #the two points linked
        self.left = data[0]
        self.right = data[1]
        self.rigidity = data[2]
        self.axis = data[3]
        self.parent = parent
        self.nbeaml = data[4]
        self.nbeamr = data[5]

    def compute_matrix(self):
        
        #updating parent.A
        for i in self.fixed :
            self.parent.A[i[0],i[0]] = self.parent.A[i[0],i[0]] + self.parent.pen
            self.parent.A[i[1],i[1]] = self.parent.A[i[1],i[1]] + self.parent.pen
            self.parent.A[i[1],i[0]] = self.parent.A[i[1],i[0]] - self.parent.pen
            self.parent.A[i[0],i[1]] = self.parent.A[i[0],i[1]] - self.parent.pen

class AppPlan (Link) :
    #appui plan
    def __init__(self,data,parents,norm):
        Link.__init__(self,data,parents)
        
class Rotule (Link) :
    #Rotule
    def __init__(self,data,parents):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2]]
        
class SphPlan (Link) :
    #sphère-plan
    def __init__(self,data,parents,norm):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6 + norm,self.right*6 + norm]]
        
class SphDoi (Link) :
    #sphérique à doigt
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)

class Pivot (Link) :
    #pivot
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)

class Gliss (Link) :
    #glissière
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)

class PivGliss (Link) :
    #Pivot-glissant
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)

class Helico (Link) :
    #hélicoïdale
    def __init__(self,data,parents,axis,step):
        Link.__init__(self,data,parents)

class SphCyl (Link) :
    #Sphère-cylindre
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)

class Encas (Link) :
    #encastrement
    def __init__(self,data,parents):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6+ 5]]


pi = 3.1415926535

#beam 1
"""
L1 = [[0.,0.,0.,0.,1.,0.],[[0.,0.,0.,0.,0.,0.],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,1.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[1,1,1,1,1,1]]] #[encastré,effort nodal nul *7,effort au bout : 10, libre]

L2 = [[0.,1.,0.,1.,1.,0.],[[1,1,1,1,1,1],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[0.,0.,0.,0.,0.,0.]]] #[encastré,effort nodal nul *7,effort au bout : 10, libre]

#A = Assembly_solv([[L1],[]])
A = Assembly_solv([[L1,L2],[['rot',[0,5],[1,0],'global']]])
#A.static_solv()
A.modal_solv()
print(A.u[:,0])"""

#[[[x0,y0,z0,xf,yf,zf],beam1],[beam2],[beam3]...,[[links]]]
#link : ['class',[n°beam,n°point],[n°beam,n°point],'axis system']
if __name__ == "__main__":
    B = Reader("assembly1.txt")
#print(B.v)