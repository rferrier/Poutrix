#This code provides the finite elements solver for assemblies of beams

#Internal importations
import numpy as np
from scipy.linalg import *
from solveur_poutres_v5 import *

#External importations
import sp
import math as mh
from parser_pout_v2_1 import parser
import os as os

class Reader :
    #assembly data-process facility

    def __init__(self,data):
        
        #Recovering of the data
        self.title = data
        self.read_input(data)
        #print(self.input[0])

        #distributing lists to readers
        self.link = []
        self.beam = []
        self.beamCut = [] #it will contain the points where the beam will be cut
        #print(self.beamCut)        
        self.beamName = {} #number : 'name' of the beam
        self.beamName1 = {} #the invert of beamName
        self.beamLeng = {} #name of the beam : leng of the beam

        for i in range(len(self.input[0])) :
            #pre-processing of the beams
            self.beamCut.append([])
            self.beamName[self.input[0][i][0]] = i  
            self.beamName1[i] = self.input[0][i][0]
            #pre-computation of the length :
            self.beamLeng[self.input[0][i][0]] = mh.sqrt((self.input[0][i][1][0][0]-self.input[0][i][1][1][0])**2 + (self.input[0][i][1][0][1]-self.input[0][i][1][1][1])**2 + (self.input[0][i][1][0][2]-self.input[0][i][1][1][2])**2)

        #print(self.beamName1)
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
            print('Assemblage de la matrice de rigidité')
            self.solv = Assembly_solv(data)
            print('Résolution du problème linéaire')
            self.solv.static_solv()

            #recovering of results
            print('Post-traitement et affichage des résultats')
            self.u = self.solv.u
            #processing of results
            self.cutVector([self.u])
            
            #printing in file :
            self.outChaine = 'Resultats pour la simulation statique : ' + str(self.title) + '\n\nLes déplacements sont donnés dans la base globale\net les efforts de cohésion dans les bases locales de la poutre\n'        

            self.outName = str(self.solr.outName) + '.txt'
            self.f = open(self.outName,'w')

            self.outputResults(self.outChaine,0)
            
            self.f.write(self.outChaine)
            self.f.close()
            os.startfile(self.outName)
            print('Fin de résolution du problème statique')

        if self.solr.type == 'modal':
            print('Assemblage des matrices de masse et de rigidité')
            self.solv = Assembly_solv(data)
            print('Résolution du problème aux valeurs propres généralisé')
            self.solv.modal_solv()
            
            #printing of the arg first modes :
            print('Post-traitement et affichage des résultats')
            
            self.outChaine = 'Resultats pour la simulation modale : ' + str(self.title) + '\n\nLes déplacements sont donnés dans la base globale\net les efforts de cohésion dans les bases locales de la poutre\n'        

            self.outName = str(self.solr.outName) + '.txt'
            self.f = open(self.outName,'w')

            listOfMod = list(self.solv.omeg)
            listOfRemoved = []
            for i in range(self.solr.arg) :
                a = max(listOfMod)
                for j in range(len(listOfMod)) :                #a = min(listOfMod)
                    if listOfMod[j] < a and listOfMod[j] not in listOfRemoved :
                        a = listOfMod[j]
                        indice = j
                        
                self.outChaine = self.outChaine + '\n' + 'MODE ' + str(i+1) + '\n' + '\omega = ' + str(abs(a)**(1/2)) + '\n'
                
                #extracting the good column in the matrix :
                vector = mat([[0]]*len(self.solv.u))
                vector[indice,0] = 1
                self.u = self.solv.u*vector
                
                self.cutVector([self.u])
                self.outputResults(self.outChaine,0)
                listOfRemoved.append(a)
            
            self.f.write(self.outChaine)
            self.f.close()
            os.startfile(self.outName)
            print('Fin de résolution du problème modal')
            
        if self.solr.type == 'dynamic':
            
            print('Assemblage des matrices de masse et de rigidité')
            self.solv = Assembly_solv(data)
            print('Inversion de la matrice et itérations de Newmark')
            self.solv.dyna_solv(self.solr.dt,self.solr.duree)
            
            self.u = self.solv.u
            
            #printing the results     
            print('Post-traitement et affichage des résultats')
            self.outChaine = 'Resultats pour la simulation dynamique : ' + str(self.title) + '\n\nLes déplacements sont donnés dans la base globale\net les efforts de cohésion dans les bases locales de la poutre\n'        

            self.outName = str(self.solr.outName) + '.txt'
            self.f = open(self.outName,'w')
            self.cutVector(self.u)
            if self.solr.postProI != {} :
                for i in range(len(self.solv.time)):
                    self.outChaine = self.outChaine + '\nPas de temps : ' + str(self.solv.time[i]) + '\n'
                    self.outputResults(self.outChaine,i)
                
            #Time depending results for one point (TIME_OUTPUT)
            for i in self.solr.postProIII :
                self.outChaine = self.outChaine + '\nRésultat dépendant du temps pour la poutre ' + str(i[0]) + ',\nà l\'abscisse curviligne ' + str(i[3]) + ', sur la composante ' + str(i[4]) + '\n'
                for j in range(len(self.solv.time)) :
                    self.outChaine = self.outChaine + '\n' + str(self.beam[self.beamName[i[0]]].v[j][i[2]][i[1]][0,0])
                
                self.outChaine = self.outChaine + '\n'
                
            self.f.write(self.outChaine)
            self.f.close()
            os.startfile(self.outName)
            print('Fin de résolution du problème dynamique')
                
    def cutVector(self,u):
        
        #Distributing u to beams and computing of F :
        j = 0
        
        for i1 in range(len(self.beam)):
            i = self.beam[i1]
            i.u = []
            i.uloc = []
            i.f = []
            i.v = []
            i.t = []
            
            for ti in range(len(u)) : #time
                i.u.append(u[ti][j:j+6*len(i.elements)+6])

                i.uloc.append(self.solv.pas[i1]*i.u[ti])           #base changing (local)               
                self.solv.microsolv[i1].mesh.compute_F1(i.uloc[ti])
                i.f.append(self.solv.microsolv[i1].mesh.F1)
                
                #Separating components
                
                i.v.append([[],[],[],[],[],[]])
                for k in range(6):
                    for k2 in range(int(len(i.u[ti])/6)):
                        i.v[ti][k].append(i.u[ti][6*k2+k][0])
                        
                i.t.append([[],[],[],[],[],[]])
                for k in range(6):
                    for k2 in range(int(len(i.f[ti])/6)):
                        i.t[ti][k].append(i.f[ti][6*k2+k][0])
                        
            j = j + 6*len(i.elements) + 6
        
    def outputResults(self,chaine,time):
        #Building the output file
        microlist = ['u_x','u_y','u_z','theta_x','theta_y','theta_z']
        microlist2 = ['F_x','F_y','F_z','M_x','M_y','M_z']

        alpha = 0
        for k in self.beam:
            
            if self.beamName1[alpha] in self.solr.postProI : #If this beam is asked

                self.outChaine = self.outChaine + '\n' + str(self.beamName1[alpha]) +'\n'

                for i1 in range(len(self.solr.postProI[self.beamName1[alpha]])) : #We take only the specified components
                    i = self.solr.postProI[self.beamName1[alpha]][i1] #the n° of the component
                    self.outChaine = self.outChaine + '\n' + microlist[i] + '\n\n'
                    for j in k.v[time][i]:
                        if type(j) == np.matrix :
                            self.outChaine = self.outChaine + str(abs(j[0,0])) + '\n'
                        else :
                            self.outChaine = self.outChaine + str(j) + '\n'
                                       
                for i1 in range(len(self.solr.postProII[self.beamName1[alpha]])) : #We take only the specified components
                    i = self.solr.postProII[self.beamName1[alpha]][i1] #the n° of the component
                    self.outChaine = self.outChaine + '\n' + microlist2[i] + '\n\n'
                    for j in k.t[time][i]:
                        self.outChaine = self.outChaine + str(j[0,0]) + '\n'
                    
            alpha = alpha + 1

class LinkR :
    #Link Reader

    def __init__(self,data,parent):
        
        self.parent = parent
        self.input = data
        #data processing
        self.leftB = parent.beamName[data[0]]
        self.leftL = parent.beamLeng[data[0]]
        self.rightB = parent.beamName[data[2]]
        self.rightL = parent.beamLeng[data[2]]
        epsilon = 0.000000001
        simeq = False
        for i in parent.beam[self.leftB].ltelist:
            if data[1] < i+epsilon/self.leftL and data[1] > i-epsilon/self.leftL:
                iretenu = i
                simeq = True
        if simeq == True:
            self.leftX = parent.beam[self.leftB].lte[iretenu]
            
        simeq = False
        for i in parent.beam[self.rightB].ltelist:
            if data[3] < i+epsilon/self.rightL and data[3] > i-epsilon/self.rightL:
                iretenu = i
                simeq = True            
        if simeq == True:
            self.rightX = parent.beam[self.rightB].lte[iretenu]
            
        self.type = self.input[4]  #class of link
        
        #axis/normal
        dico = {'x':0,'X':0,'y':1,'Y':1,'z':2,'Z':2}
        if len(self.input[5]) == 1 :
            self.axis = dico[self.input[5][0]]
        else :
            self.axis = 0
        
        self.repere = 'global'
        
    def build_output(self):
        self.output = [self.type,[self.leftB,self.leftX],[self.rightB,self.rightX],self.repere,self.axis]

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
        self.S = self.input[1][4]
        self.Iy = self.input[1][5]
        self.Iz = self.input[1][6]
        self.rho = self.input[1][7]
        
        #Links
        self.cut1 = parent.beamCut[parent.beamName[self.name]]
        self.cut = []
        #removing doubles elements :
        for i in self.cut1:
            if i not in self.cut :
                self.cut.append(i)
        #print(self.cut)
        self.points = self.cut[:]
        self.points.sort()
        #adding beginning and end
        if self.points == []:
            self.points.insert(0,0.0)
            self.points.append(self.leng)
        else :
            if self.points[0] > 0.000000001*self.leng:
                self.points.insert(0,0.0)
            if self.points[-1] < self.leng*(1-0.000000001):
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
            linit = linit + (i+1)*self.elem_leng    #increment

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
        epsilon = 0.000000001
        for i in range(len(self.elements)+1) :
            alpha = sum(self.elements[:i])   #sum of len of elements before i = abscissa of link node of elem i
            simeq = False
            #check if alpha ~= an elem of self.cut
            for j in self.cut : 
                if alpha < j+epsilon/self.leng and alpha > j-epsilon/self.leng :
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
            self.output.append([[0.]*6,[0.]*7,[0.]*6])

            self.output[i+1][2][0] = self.list_char_mieux2[6][i]#fx
            self.output[i+1][2][1] = self.list_char_mieux2[7][i]#fy
            self.output[i+1][2][2] = self.list_char_mieux2[8][i]#fz
            self.output[i+1][2][3] = self.list_char_mieux2[9][i]#mx
            self.output[i+1][2][4] = self.list_char_mieux2[10][i]#my
            self.output[i+1][2][5] = self.list_char_mieux2[11][i]#mz
            
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
        self.parent = parent
        self.type = data[0][0]
        
        #Processing of the optional arguments (ex : nb of modes)
        self.arg = 1
        if len(data[0][1]) == 1:
            self.arg = data[0][1][0]
        
        self.duree = 1
        if len(data[1]) == 1:
            self.duree = data[1][0]
            
        self.dt = 0.1
        if len(data[2]) == 1:
            self.dt = data[2][0]
            
        self.outName = data[3]
        
        #reading of postprocess instructions
        self.ppi = data[4]
        self.postProI = {} #'beamName':[components]
        self.postProII = {} #F
        dico = {'Ux':0,'Uy':1,'Uz':2,'Theta_x':3,'Theta_y':4,'Theta_z':5}
        for i in range(len(self.ppi)) :
            liste = []
            liste1 = []
            
            for j in range(len(self.ppi[i][1])) : #U
                liste.append(dico[self.ppi[i][1][j]])
            self.postProI[self.ppi[i][0]] = liste
            
            for j in range(len(self.ppi[i][2])) : #F
                liste1.append(dico[self.ppi[i][2][j]])
            self.postProII[self.ppi[i][0]] = liste1
        
        self.ppiii = data[5]
        self.postProIII = [] #('beamName',n° of node,n° of component,abscissa,'component')
        
        for i in self.ppiii :
            #print(self.parent.beamName)
            #Computing the n° of node associated with the given abscissa :
            abscissa = 0
            listOfEl = self.parent.beam[self.parent.beamName[i[0]]].elements

            for j in range(len(listOfEl)) : #the elements of the beam
                abscissa = abscissa + listOfEl[j] #increment of length
                
                if abscissa < i[1] + 0.000000001*self.parent.beam[self.parent.beamName[i[0]]].leng and abscissa > i[1] - 0.000000001*self.parent.beam[self.parent.beamName[i[0]]].leng :
                    theChosen = j #the good n° of node
                    
            self.postProIII.append((i[0],theChosen,dico[i[2]],i[1],i[2]))

class Assembly_solv :

    def __init__ (self,data):
        #data recovering
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
            
            #Link type management :
            
            if i[0][0] == 'rot':
                self.links.append(Rotule([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self))

            if i[0][0] == 'enca':
                self.links.append(Encas([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self))
            
            if i[0][0] == 'splan':
                self.links.append(SphPlan([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'aplan':
                self.links.append(AppPlan([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'pivot':
                self.links.append(Pivot([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'glis':
                self.links.append(Gliss([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'pglis':
                self.links.append(PivGliss([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'scyl':
                self.links.append(SphCyl([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))

            if i[0][0] == 'sdoi':
                self.links.append(SphDoi([self.leftl,self.rightl,0,i[3],i[1][0],i[2][0]],self,i[4]))


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
        #print(self.links)
        for i in self.links:
            i.compute_matrix()
        self.assembly_b()
        self.u = solve(self.A,self.b,sym_pos = True)
        
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
        
    def dyna_solv(self,dt,duree):
        self.compute_mat()
        self.compute_leng()

        self.transform_mat()
        self.assembly_A()
        self.assembly_M()
        self.compute_pen()
        for i in self.links:
            i.compute_matrix()
        self.assembly_b()

        #Newmark
        gamma = 0.5
        beta = 0.25
        nbt = int(duree/dt)

        self.A1=inv(self.M + beta*dt*dt*self.A);
        
        self.time = [0]                 #list of time
        self.acc = [self.b - self.b]    #list of acceleration
        self.u = [self.b - self.b]      #list of position
        self.v = [self.b - self.b]      #list of speed

        for i in range(1,nbt) :
            self.time.append(i*dt)
            self.acc.append(self.A1*(self.b - self.A * (self.u[i-1] + dt*self.v[i-1] + (dt*dt)/2*(1-2*beta)*self.acc[i-1])))
    
            self.u.append(self.u[i-1] + dt*self.v[i-1] + (dt*dt)/2*((1-2*beta)*self.acc[i-1] + 2*beta*self.acc[i]))
            self.v.append(self.v[i-1] + (1-gamma)*dt*self.acc[i-1] + gamma*dt*self.acc[i])


    def compute_leng(self):
        #computation of the size of the matrix
        self.leng = 0
        for i in self.microsolv:
            self.leng = self.leng + len(i.mesh.A)

    def compute_pen(self):
        #computation of the penalization matrix
        self.penlist = []
        self.pen = []
        
        for j in range(6):
            self.penlist.append([])
            for i in self.microsolv:
                #elem 0 is left linked while elem n is right linked
                self.penlist[j].append(i.mesh.elements[0].pen[j])
                self.penlist[j].append(i.mesh.elements[-1].pen[j+5])

        for j in range(6):
            self.pen.append(max(self.penlist[j]))

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
            #computing of the components by euclidian rest:
            fix0 = i[0] - 6*int(i[0]/6)
            fix1 = i[1] - 6*int(i[1]/6)
            fix = max(fix0,fix1)        #??
            #print(fix)
            #print(self.parent.pen)
            
            self.parent.A[i[0],i[0]] = self.parent.A[i[0],i[0]] + self.parent.pen[fix]
            self.parent.A[i[1],i[1]] = self.parent.A[i[1],i[1]] + self.parent.pen[fix]
            self.parent.A[i[1],i[0]] = self.parent.A[i[1],i[0]] - self.parent.pen[fix]
            self.parent.A[i[0],i[1]] = self.parent.A[i[0],i[1]] - self.parent.pen[fix]

class AppPlan (Link) :
    #appui plan
    def __init__(self,data,parents,norm):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6 + norm,self.right*6 + norm],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6 + 5]]
        self.fixed.remove([self.left*6 + 3 + norm,self.right*6 + 3 + norm])
        
        self.name = 'aplan' #name for the parser
        
class Rotule (Link) :
    #Rotule
    def __init__(self,data,parents):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2]]
        
        self.name = 'rot' #name for the parser        
        
class SphPlan (Link) :
    #sphère-plan
    def __init__(self,data,parents,norm):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6 + norm,self.right*6 + norm]]
        
        self.name = 'splan' #name for the parser        
        
class SphDoi (Link) :
    #sphérique à doigt
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3 + axis,self.right*6 + 3 + axis]]
        
        self.name = 'sdoi' #name for the parser        
        
class Pivot (Link) :
    #pivot
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6+ 5]]
        self.fixed.remove([self.left*6 + 3 + axis,self.right*6 + 3 + axis])

        self.name = 'pivot' #name for the parser

class Gliss (Link) :
    #glissière
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6+ 5]]
        self.fixed.remove([self.left*6 + axis,self.right*6 + axis])

        self.name = 'glis' #name for the parser

class PivGliss (Link) :
    #Pivot-glissant
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6+ 5]]
        self.fixed.remove([self.left*6 + axis,self.right*6 + axis])
        self.fixed.remove([self.left*6 + 3 + axis,self.right*6 + 3 + axis])

        self.name = 'pglis' #name for the parser

class Helico (Link) :
    #hélicoïdale
    def __init__(self,data,parents,axis,step):
        Link.__init__(self,data,parents)

        self.name = 'heli' #name for the parser

class SphCyl (Link) :
    #Sphère-cylindre
    def __init__(self,data,parents,axis):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2]]
        self.fixed.remove([self.left*6 + axis,self.right*6 + axis])

        self.name = 'scyl' #name for the parser

class Encas (Link) :
    #encastrement
    def __init__(self,data,parents):
        Link.__init__(self,data,parents)
        self.fixed = [[self.left*6,self.right*6],[self.left*6 + 1,self.right*6 + 1],[self.left*6 + 2,self.right*6 + 2],[self.left*6 + 3,self.right*6 + 3],[self.left*6 + 4,self.right*6 + 4],[self.left*6 + 5,self.right*6+ 5]]
        
        self.name = 'enca' #name for the parser

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
    B = Reader("dynatest.txt")
#print(B.v)