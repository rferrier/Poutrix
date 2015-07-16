#This code provides the finite elements solver for assemblies of beams

#Internal importations

from solveur_poutres_v6 import *
from parser_pout_v2_2 import parser

#External importations
import numpy as np
from scipy.linalg import *
import sp
import math as mh
import os as os

class Reader :
    #assembly data-process facility

    def __init__(self,data):
        
        #Recovering of the data
        self.title = data
        self.read_input(data)
        #print(self.input[0])

        #processing of the parameters            
        if len(self.input[5]) == 1 :
            self.param = ParamR(self.input[5][0],self) 

        #distributing lists to readers
        self.msys = []
        self.link = []
        self.beam = []
        self.beamCut = [] #it will contain the points where the beam will be cut
        #print(self.beamCut)   
        self.msysName = {} #name of the system : MobSystem object
        self.beamName = {} #number : 'name' of the beam
        self.beamName1 = {} #the invert of beamName
        self.beamLeng = {} #name of the beam : leng of the beam

        for i in range(len(self.input[1])) :
            #pre-processing of the beams
            self.beamCut.append([])
            self.beamName[self.input[1][i][0]] = i  
            self.beamName1[i] = self.input[1][i][0]
            #pre-computation of the length :
            self.beamLeng[self.input[1][i][0]] = mh.sqrt((eval(self.input[1][i][1][0][0])-eval(self.input[1][i][1][1][0]))**2 + (eval(self.input[1][i][1][0][1])-eval(self.input[1][i][1][1][1]))**2 + (eval(self.input[1][i][1][0][2])-eval(self.input[1][i][1][1][2]))**2)

        #print(self.beamName1)
        for i in range(len(self.input[3])) :
            #pre-processing of the links in order to mesh the beams
            #every beam becomes its links points (in order to place a node)
            self.beamCut[self.beamName[self.input[3][i][0]]].append(eval(self.input[3][i][1]))
            self.beamCut[self.beamName[self.input[3][i][2]]].append(eval(self.input[3][i][3]))
        
        for i in range(len(self.input[1])) :
            self.beam.append(BeamR(self.input[1][i],self))
            
        for i in range(len(self.input[2])):
            self.msys.append(MobSystem(self.input[2][i],self))
            
        for i in range(len(self.input[3])) :
            self.link.append(LinkR(self.input[3][i],self))
            
        self.solr = SolvR(self.input[4],self)
        
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
        if self.solr.type == 'static' or self.solr.type == 'statique':
            print('Assemblage de la matrice de rigidité')
            self.solv = Assembly_solv(data,self)
            print('Résolution du problème linéaire')
            self.solv.static_solv()

            #recovering of results
            print('Post-traitement et affichage des résultats')
            self.u = self.solv.u
            #processing of results
            self.cutVector([self.u])
            
            #printing in file :
            self.outChaine = 'Resultats pour la simulation statique : ' + str(self.title) + '\n\nLes déplacements et les efforts de cohésion sont donnés dans la base globale\n'        

            self.outName = str(self.solr.outName) + '.txt'
            self.f = open(self.outName,'w')

            self.outputResults(self.outChaine,0)
            
            self.f.write(self.outChaine)
            self.f.close()
            os.startfile(self.outName)
            print('Fin de résolution du problème statique')
            
        if self.solr.type == 'staticGD' or self.solr.type == 'statiqueGD':
            self.niter = self.solr.niter
            
            #first iteration
            self.solv = Assembly_solv(data,self)
            
            print('iteration : 1/'+str(self.niter))
            self.solv.static_solv()
            
            self.u = self.solv.u/self.niter
            
            for i in range(self.niter-1):
                
                #Updating matrix
                self.cutVector([self.u])
                self.updateKnF()
                
                print('iteration : '+str(i+2)+'/'+str(self.niter))
                
                self.solv.u = solve(self.solv.A,self.solv.b,sym_pos = True)
                
                self.u += self.solv.u/self.niter
            
            self.cutVector([self.u])
                
            #printing in file :
            self.outChaine = 'Resultats pour la simulation statique avec grands déplacements : ' + str(self.title) + '\n\nLes déplacements et les efforts de cohésion sont donnés dans la base globale\n'        

            self.outName = str(self.solr.outName) + '.txt'
            self.f = open(self.outName,'w')

            self.outputResults(self.outChaine,0)
            
            self.f.write(self.outChaine)
            self.f.close()
            os.startfile(self.outName)
            print('Fin de résolution du problème statique')
                

        if self.solr.type == 'modal':
            print('Assemblage des matrices de masse et de rigidité')
            self.solv = Assembly_solv(data,self)
            print('Résolution du problème aux valeurs propres généralisé')
            self.solv.modal_solv()
            
            #printing of the arg first modes :
            print('Post-traitement et affichage des résultats')
            
            self.outChaine = 'Resultats pour la simulation modale : ' + str(self.title) + '\n\nLes déplacements et les efforts de cohésion sont donnés dans la base globale\n'        

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
            
        if self.solr.type == 'dynamic' or self.solr.type == 'dynamique' :
            
            print('Assemblage des matrices de masse et de rigidité')
            self.solv = Assembly_solv(data,self)
            print('Inversion de la matrice et itérations de Newmark')
            self.solv.dyna_solv(self.solr.dt,self.solr.duree)
            
            self.u = self.solv.u
            
            #printing the results     
            print('Post-traitement et affichage des résultats')
            self.outChaine = 'Resultats pour la simulation dynamique : ' + str(self.title) + '\n\nLes déplacements et les efforts de cohésion sont donnés dans la base globale\n'        

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
                self.solv.microsolv[i1].mesh.compute_F1(i.u[ti])
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
            
    def updateKnF(self):
        #updating K and f for geometrical non-linear problems
        ti = 0
        for j in range(len(self.solv.microsolv)):
            j1 = self.solv.microsolv[j]
            #self.beam[i].straight = 'defo'
            #print(j1.mesh.elements)
            for i in range(len(j1.mesh.elements)) :
                i1 = j1.mesh.elements[i]
                #re-computing x,y and z
                xd = self.beam[j].v[ti][0][i]
                yd = self.beam[j].v[ti][1][i]
                zd = self.beam[j].v[ti][2][i]
                
                xf = i1.lenginit*mh.cos(i1.theta1init)*mh.cos(i1.theta2init) + self.beam[j].v[ti][0][i+1]
                yf = i1.lenginit*mh.sin(i1.theta1init)*mh.cos(i1.theta2init) + self.beam[j].v[ti][1][i+1]
                zf = -i1.lenginit*mh.sin(i1.theta2init) + self.beam[j].v[ti][2][i+1]
                #re-computation of h and theta
                h = mh.sqrt((xf-xd)**2 + (yf-yd)**2 + (zf-zd)**2)
                resuloc = compute_angle(xd,yd,zd,xf,yf,zf)
                
                #replacement in elemn caracteristics :
                i1.leng = h
                i1.theta1 = resuloc[0]
                i1.theta2 = resuloc[1]
                
                #re-computing of matrix :
                i1.compute_Ke()
                i1.compute_be()
                i1.compute_P()
            j1.mesh.compute_A()
            #print(j1.mesh.A)
            j1.mesh.compute_b()
            
        self.solv.assembly_A()
        for i in self.solv.links:
            i.update_matrix()

        self.solv.assembly_b()           
        
    def outputResults(self,chaine,time):
        #Building the output file
        microlist = ['u_x','u_y','u_z','theta_x','theta_y','theta_z']
        microlist2 = ['F_x','F_y','F_z','M_x','M_y','M_z']
        microlist3 = ['x','y','z']

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
       
       
            if self.beamName1[alpha] in self.solr.postProVI : #Plot Output

                self.outChaine = self.outChaine + '\n' + str(self.beamName1[alpha]) +'\n\n'

                #let's begin with the title
                for i in self.solr.postProIV[self.beamName1[alpha]] : #We take only the specified components
                    self.outChaine = self.outChaine + microlist3[i] + self.solr.separator
                for i in self.solr.postProV[self.beamName1[alpha]] :
                    self.outChaine = self.outChaine + microlist[i] + self.solr.separator
                for i in self.solr.postProVI[self.beamName1[alpha]] :
                    self.outChaine = self.outChaine + microlist2[i] + self.solr.separator
                    
                self.outChaine = self.outChaine + '\n'
                
                for j1 in range(len(self.solv.microsolv[alpha].mesh.nodes)) :
                    self.solv.microsolv[alpha].mesh.nodes[j1].compute_coords()
                    for i1 in range(len(self.solr.postProIV[self.beamName1[alpha]])) : #We take only the specified components
                        i = self.solr.postProIV[self.beamName1[alpha]][i1] #the n° of the component
                        j = self.solv.microsolv[alpha].mesh.nodes[j1].coords[self.solr.postProIV[self.beamName1[alpha]][i1]]
                        if type(j) == np.matrix :
                            self.outChaine = self.outChaine + str(abs(j[0,0])) + ';'
                        else :
                            self.outChaine = self.outChaine + str(j) + ';'

                    for i1 in range(len(self.solr.postProV[self.beamName1[alpha]])) : #We take only the specified components
                        i = self.solr.postProV[self.beamName1[alpha]][i1] #the n° of the component
                        
                        j = k.v[time][i][j1]
                        if type(j) == np.matrix :
                            self.outChaine = self.outChaine + str(abs(j[0,0])) + ';'
                        else :
                            self.outChaine = self.outChaine + str(j) + ';'
                            
                    for i1 in range(len(self.solr.postProVI[self.beamName1[alpha]])) : #We take only the specified components
                        i = self.solr.postProVI[self.beamName1[alpha]][i1] #the n° of the component
                        
                        #j = k.t[time][i][j1]
                        j = 'sorry:not_implemented'
                        if type(j) == np.matrix :
                            self.outChaine = self.outChaine + str(abs(j[0,0])) + ';'
                        else :
                            self.outChaine = self.outChaine + str(j) + ';'
                            
                    self.outChaine = self.outChaine + '\n'
                            
            alpha = alpha + 1

class ParamR :
    #Parameters reader
    def __init__(self,data,parent):
        
        self.parent = parent
        #processing of the parameters values
        for i in data:
            globals()[i[0]] = eval(i[1])
        
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
            if eval(data[1]) < i+epsilon/self.leftL and eval(data[1]) > i-epsilon/self.leftL:
                iretenu = i
                simeq = True
        if simeq == True:
            self.leftX = parent.beam[self.leftB].lte[iretenu]
            
        simeq = False
        for i in parent.beam[self.rightB].ltelist:
            if eval(data[3]) < i+epsilon/self.rightL and eval(data[3]) > i-epsilon/self.rightL:
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
        
        #axis system
        self.repere = 'global'
        if len(self.input[6]) == 1:
            self.repere = self.input[6][0]
        
        
    def build_output(self):
        self.output = [self.type,[self.leftB,self.leftX],[self.rightB,self.rightX],self.repere,self.axis]

class BeamR :
    #Beam reader

    def __init__(self,data,parent):

        self.straight = 'yes'    #by default, the beam is straight
        self.input = data
        #print(self.input)
        #data processing
        self.name = self.input[0]

        #characteristics
        self.leng = mh.sqrt((eval(self.input[1][0][0])-eval(self.input[1][1][0]))**2 + (eval(self.input[1][0][1])-eval(self.input[1][1][1]))**2 + (eval(self.input[1][0][2])-eval(self.input[1][1][2]))**2)
        self.origin = (eval(self.input[1][0][0]),eval(self.input[1][0][1]),eval(self.input[1][0][2])) #x,y,z at the left point of the beam
        self.E = self.input[1][2]
        self.nu = self.input[1][3]
        self.S = self.input[1][4]
        self.Iy = self.input[1][5]
        self.Iz = self.input[1][6]
        self.rho = self.input[1][7]
        
        #Links
        #print(parent.beamName)
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
                self.list_char_mieux1[i] = '0.'
                
            else :
                self.list_char_mieux1[i] = self.list_char[i][0]
        
        #Type of mesh
        self.elem_number = eval(self.input[1][8])
        
        #processing of initial deformation
        self.defo = self.input[1][9]
        if self.defo != [] :
            self.straight = 'defo'
            
            self.xdefo1 = self.defo[0][0]
            if self.xdefo1 == []:
                self.xdefo = '0.'
            else :
                self.xdefo  =self.xdefo1[0]
                
            self.ydefo1 = self.defo[0][1]
            if self.ydefo1 == []:
                self.ydefo = '0.'
            else :
                self.ydefo  =self.ydefo1[0]
                
            self.zdefo1 = self.defo[0][2]
            if self.zdefo1 == []:
                self.zdefo = '0.'
            else :
                self.zdefo  =self.zdefo1[0]


        #Initialization of the elements
        self.elements = []

        #Processing of the variable characteristics
        self.list_char = [self.E,self.nu,self.S,self.Iy,self.Iz,self.rho,self.list_char_mieux1[0],self.list_char_mieux1[1],self.list_char_mieux1[2],self.list_char_mieux1[3],self.list_char_mieux1[4],self.list_char_mieux1[5]] 
        self.list_char_mieux2 = [0]*len(self.list_char)
        self.j = 0

        for i in range(len(self.list_char)) :
            #looking for a changing characteristic
            if self.list_char[i][0] == '[':#found !
                #building of the name of the file to open
                self.name1 = self.list_char[i][1:-1] + '.dat'

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
                        
                self.f.close

        if self.j == 0 :
            self.simple_mesh()

        for i in range(len(self.list_char)) :
            #Building of the list of char if it is not in a file
            if self.list_char[i][0] != '[':
                X = 0  #init of curvilign abscissa
                self.list_char_mieux2[i] = []
                for j in range(len(self.elements)+1) :
                    self.list_char_mieux2[i].append(eval(self.list_char[i]))
                    if j<len(self.elements):
                        X = X + self.elements[j]

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
    
    def compute_angles(self):
        if self.straight == 'yes':
            #it's useless to compute every angle
            resu = compute_angle(eval(self.input[1][0][0]),eval(self.input[1][0][1]),eval(self.input[1][0][2]),eval(self.input[1][1][0]),eval(self.input[1][1][1]),eval(self.input[1][1][2]))
    
            self.theta1 = [resu[0]]*len(self.elements)
            self.theta2 = [resu[1]]*len(self.elements)
            
        if self.straight == 'defo':
            self.theta1 = []
            self.theta2 = []
            self.elements1 = []
            X = 0.
            
            #processing of the origin :
            self.origin = (self.origin[0]+eval(self.xdefo),self.origin[1]+eval(self.ydefo),self.origin[2]+eval(self.zdefo))
            #pre-computation of theta
            resu = compute_angle(eval(self.input[1][0][0]),eval(self.input[1][0][1]),eval(self.input[1][0][2]),eval(self.input[1][1][0]),eval(self.input[1][1][1]),eval(self.input[1][1][2]))
            for i in range(len(self.elements)) :
                #re-computing x,y and z
                xd = eval(self.xdefo)
                yd = eval(self.ydefo)
                zd = eval(self.zdefo)

                X = X + self.elements[i]
                
                xf = self.elements[i]*mh.cos(resu[0])*mh.cos(resu[1]) + eval(self.xdefo)
                yf = self.elements[i]*mh.sin(resu[0])*mh.cos(resu[1]) + eval(self.ydefo)
                zf = -self.elements[i]*mh.sin(resu[1]) + eval(self.zdefo)
                #re-computation of h and theta
                self.elements1.append(mh.sqrt((xf-xd)**2 + (yf-yd)**2 + (zf-zd)**2))
                resuloc = compute_angle(xd,yd,zd,xf,yf,zf)
                self.theta1.append(resuloc[0])
                self.theta2.append(resuloc[1])
            
            self.elements = self.elements1[:] #repalcement of h


    def build_output(self):

        self.compute_angles()

        self.output = []

        #left boundaries
        self.output.append([0.]*6)
        self.list_char = [self.uxg,self.uyg,self.uzg,self.txg,self.tyg,self.tzg]
        for i in range(len(self.list_char)) :
            if len(self.list_char[i]) == 0:
                #free
                self.output[0][i] = 1
            else :
                self.output[0][i] = eval(self.list_char[i][0])

        #nodes
        for i in range (len(self.list_char_mieux2[0])):
            self.output.append([[0.]*6,[0.]*9,[0.]*6])

            self.output[i+1][2][0] = self.list_char_mieux2[6][i]#fx
            self.output[i+1][2][1] = self.list_char_mieux2[7][i]#fy
            self.output[i+1][2][2] = self.list_char_mieux2[8][i]#fz
            self.output[i+1][2][3] = self.list_char_mieux2[9][i]#mx
            self.output[i+1][2][4] = self.list_char_mieux2[10][i]#my
            self.output[i+1][2][5] = self.list_char_mieux2[11][i]#mz
                
            if i < len(self.elements): #the last element of the list only represents a node
                self.output[i+1][1][0] = self.elements[i]#h
                self.output[i+1][1][7] = self.theta1[i]#theta1
                self.output[i+1][1][8] = self.theta2[i]#theta2

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
                self.list_char_mieux3[i] = eval(self.list_char[i][0])
                
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
                self.output[len(self.elements)+2][i] = eval(self.list_char[i][0])
        
        #dimensions of the beam :
        self.output1 = [[eval(self.input[1][0][0]),eval(self.input[1][0][1]),eval(self.input[1][0][2]),eval(self.input[1][1][0]),eval(self.input[1][1][1]),eval(self.input[1][1][2])],self.output,self.origin]

        #print(self.output1)

    def computeNoOfNode(self,givAbscissa) :
        """This function computes the n° of a node when it's abscissa is given"""
        abscissa = 0
        listOfEl = self.elements

        for j in range(len(listOfEl)) : #the elements of the beam
                abscissa = abscissa + listOfEl[j] #increment of length

                if abscissa < givAbscissa + 0.000000001*self.leng and abscissa > givAbscissa - 0.000000001*self.leng :
                    theChosen = j #the good n° of node
        
        return theChosen

class SolvR :
    #The reader for solver data

    def __init__(self,data,parent):
        self.parent = parent
        self.type = data[0][0]
        
        #Processing of the optional arguments (ex : nb of modes)
        self.arg = 1
        if len(data[0][1]) == 1:
            self.arg = eval(data[0][1][0])
        
        self.duree = 1
        if len(data[1]) == 1:
            self.duree = eval(data[1][0])
            
        self.dt = 0.1
        if len(data[2]) == 1:
            self.dt = eval(data[2][0])
        
        self.niter = 10
        if len(data[3]) == 1:
            self.niter = eval(data[3][0])
            
        self.outName = data[4]
        
        #reading of postprocess instructions
        self.ppi = data[5]
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
        
        self.ppiii = data[6]
        self.postProIII = [] #('beamName',n° of node,n° of component,abscissa,'component')
        
        for i in self.ppiii :
            #print(self.parent.beamName)
            #Computing the n° of node associated with the given abscissa :
                    
            theBeam = self.parent.beam[self.parent.beamName[i[0]]]
                    
            self.postProIII.append((i[0],theBeam.computeNoOfNode(eval(i[1])),dico[i[2]],i[1],i[2]))

        #reading of plot postprocess instructions
        self.ppiv = data[7]
        
        self.separator = ';'
        
        if self.ppiv!= []:
            if self.ppiv[0][4] != []:
                self.separator = self.ppiv[0][4][0]
                
        self.postProIV = {} #X,Y,Z'beamName':[components]
        self.postProV = {} #U
        self.postProVI = {} #F
        dicoxyz = {'X':0,'Y':1,'Z':2}
        for i in range(len(self.ppiv)) :
            liste = []
            liste1 = []
            liste2 = []
            
            for j in range(len(self.ppiv[i][1])) : #X
                liste.append(dicoxyz[self.ppiv[i][1][j]])
            self.postProIV[self.ppiv[i][0]] = liste
            
            for j in range(len(self.ppiv[i][2])) : #U
                liste1.append(dico[self.ppiv[i][2][j]])
            self.postProV[self.ppiv[i][0]] = liste1

            for j in range(len(self.ppiv[i][3])) : #F
                liste2.append(dico[self.ppiv[i][3][j]])
            self.postProVI[self.ppiv[i][0]] = liste2

class Assembly_solv :

    def __init__ (self,data,parent):
        self.parent = parent
        #data recovering
        self.beams = data[0]
        #print(self.beams)
        self.linksd = data[1]
        self.microsolv = []
        self.links = []
        
        self.create_microsolv()
        self.create_links()

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

    def static_solv(self):
        self.compute_mat()
        self.compute_leng()

        self.assembly_A()
        self.compute_pen()

        for i in self.links:
            i.compute_matrix()
        self.assembly_b()

        self.u = solve(self.A,self.b,sym_pos = True)
        
    def modal_solv(self):
        self.compute_mat()
        self.compute_leng()
        
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
                    self.A[s+ii,s+jj] = self.microsolv[i].mesh.A[ii,jj]

            s = s + len(self.microsolv[i].mesh.A)
            
    def assembly_M(self):

        self.M = np.mat([[0.]*self.leng]*self.leng)
        s = 0
        for i in range(len(self.microsolv)):
            for ii in range(len(self.microsolv[i].mesh.M)):
                for jj in range(len(self.microsolv[i].mesh.M)):
                    self.M[s+ii,s+jj] = self.microsolv[i].mesh.M[ii,jj]

            s = s + len(self.microsolv[i].mesh.M)

    def assembly_b(self):

        self.b = np.mat([[0.]]*self.leng)
        s = 0
        for i in range(len(self.microsolv)):
            for ii in range(len(self.microsolv[i].mesh.b)):
                self.b[s+ii,0] = self.microsolv[i].mesh.b[ii,0]

            s = s + len(self.microsolv[i].mesh.b)

class Link:
    #Basic class for links

    def __init__(self,data,parent):
        #the two points linked
        self.left = data[0]
        self.right = data[1]
        self.rigidity = data[2]
        self.systemID = data[3]
        self.parent = parent
        self.nbeaml = data[4]
        self.nbeamr = data[5]

    def compute_matrix(self):
        #creating the link's stiffness matrix
        if self.systemID == 'global' :
            for i in self.fixed :
                #computing of the components by euclidian rest:
                fix0 = i[0] - 6*int(i[0]/6)
                fix1 = i[1] - 6*int(i[1]/6)
                fix = max(fix0,fix1)        #??
                pen = self.parent.pen[fix]
                
                self.parent.A[i[0],i[0]] = self.parent.A[i[0],i[0]] + pen
                self.parent.A[i[1],i[1]] = self.parent.A[i[1],i[1]] + pen
                self.parent.A[i[1],i[0]] = self.parent.A[i[1],i[0]] - pen
                self.parent.A[i[0],i[1]] = self.parent.A[i[0],i[1]] - pen
        
        else :
        #we need to change base a matrix
            self.compute_Aloc()
            
            #base changing
            self.system = self.parent.parent.msysName[self.systemID]

            self.Aglob = np.transpose(self.system.pas) * self.Aloc * self.system.pas

            #updating parent.A (the global stiffness matrix of the linear problem)
            self.update_parent()

    def update_matrix(self) :
        """This function actualizes the matrix of a link. To be used in iterative processes"""
        #recovering the angles
        ti = 0
        beamm = self.parent.parent.beam[self.nbeaml] #the beam  
        
        thx = beamm.v[ti][3][self.left]
        thy = beamm.v[ti][4][self.left]
        thz = beamm.v[ti][5][self.left]

        pamat = compute_mat_zyx(thz,thy,thx)
        
        self.compute_Aloc()
        
        if self.systemID == 'global' :
            self.Aglob = np.transpose(pamat) * self.Aloc * pamat
            
        else :
            self.Aglob = np.transpose(self.system.pas) * np.transpose(pamat) * self.Aloc * pamat * self.system.pas
        
        self.update_parent()        
        
    def compute_Aloc(self) :
        """Computes the local stiffness matrix for the link"""
        self.Aloc = np.mat([[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0]])

        for i in self.fixed :
            #computing of the components by euclidian rest:
            fix0 = i[0] - 6*int(i[0]/6)
            fix1 = i[1] - 6*int(i[1]/6)
            fix = max(fix0,fix1)        #??
            pen = self.parent.pen[fix]
            
            #process of the link's stiffness
            self.Aloc[fix0,fix0] = pen/10000000   #to avoid error with length
            self.Aloc[6+fix1,6+fix1] = pen/10000000
            self.Aloc[fix0,6+fix1] = pen/10000000
            self.Aloc[6+fix1,fix0] = pen/10000000

    def update_parent(self) :
        """updates the global stiffness matrix"""
        for i in range(6) :

            self.parent.A[self.left*6+i,self.left*6+i] = self.parent.A[self.left*6+i,self.left*6+i] + self.Aglob[i,i]*10000000
            self.parent.A[self.right*6+i,self.right*6+i] = self.parent.A[self.right*6+i,self.right*6+i] + self.Aglob[i+6,i+6]*10000000
            self.parent.A[self.right*6+i,self.left*6+i] = self.parent.A[self.right*6+i,self.left*6+i] - self.Aglob[i+6,i]*10000000
            self.parent.A[self.left*6+i,self.right*6+i] = self.parent.A[self.left*6+i,self.right*6+i] - self.Aglob[i,i+6]*10000000

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

class MobSystem():
    #Mobile axis system (for the links)
    def __init__(self,data,parent):
        #print(data)
        self.parent = parent
        self.data = data
        
        self.name = self.data[0]
        
        self.parent.msysName[self.name] = self #put it in the dico
        
        #find the element used as reference
        self.beam = self.parent.beam[self.parent.beamName[self.data[1][0]]]
        self.ref = self.beam.computeNoOfNode(eval(self.data[1][1]))
        
        if self.data[2][0] == 'ZYX' :
            self.pas = compute_mat_zyx(eval(self.data[3]),eval(self.data[4]),eval(self.data[5]))


#Functions

def compute_angle(x1,y1,z1,x2,y2,z2):

    #computation of theta1 and theta2 for one beam
    x = x2-x1
    y = y2-y1
    z = z2-z1
    h = mh.sqrt(x**2+y**2)

    if x > 0.:
        if y > 0.:
            theta1 = mh.atan(y/x)
            theta2 = -mh.atan(z/h)
        if y < 0.:
            theta1 = mh.atan(y/x)
            theta2 = -mh.atan(z/h)
        if y == 0.:
            theta2 = -mh.atan(z/h)
            theta1 = 0.
            
    if x < 0.:
        if y > 0.:
            theta1 = pi-mh.atan(-y/x)
            theta2 = -mh.atan(z/h)
        if y < 0.:
            theta1 = -pi+mh.atan(y/x)
            theta2 = -mh.atan(z/h)
        if y == 0.:
            theta2 = -mh.atan(z/h)
            theta1 = 0.
        
    if x == 0.:
        if y > 0.:
            theta1 = pi/2
            theta2 = -mh.atan(z/h)
        if y < 0.:
            theta1 = -pi/2
            theta2 = -mh.atan(z/h)
        if y == 0.:
            if z > 0.:
                theta2 = -pi/2
            else:
                theta2 = pi/2
            theta1 = 0.
            
    #print(x,y,z,self.theta1,self.theta2)
    return theta1 , theta2

def compute_mat_zyx(th3,th2,th1):
    """Computation of the matrix using the zyx method"""
    c1 = mh.cos(th1)
    c2 = mh.cos(th2)
    c3 = mh.cos(th3)
    s1 = mh.sin(th1)
    s2 = mh.sin(th2)
    s3 = mh.sin(th3)
    
    return np.mat([[c1*c3,-s1*s2*c3+c1*s3,c1*s2*c3+s1*s3,0,0,0,0,0,0,0,0,0],[-c2*s3,s1*s2*s3+c1*c3,-c1*s2*s3+c3*s1,0,0,0,0,0,0,0,0,0],[-s2,-s1*c2,c1*c2,0,0,0,0,0,0,0,0,0],[0,0,0,c1*c3,-s1*s2*c3+c1*s3,c1*s2*c3+s1*s3,0,0,0,0,0,0],[0,0,0,-c2*s3,s1*s2*s3+c1*c3,-c1*s2*s3+c3*s1,0,0,0,0,0,0],[0,0,0,-s2,-s1*c2,c1*c2,0,0,0,0,0,0],[0,0,0,0,0,0,c1*c3,-s1*s2*c3+c1*s3,c1*s2*c3+s1*s3,0,0,0],[0,0,0,0,0,0,-c2*s3,s1*s2*s3+c1*c3,-c1*s2*s3+c3*s1,0,0,0],[0,0,0,0,0,0,-s2,-s1*c2,c1*c2,0,0,0],[0,0,0,0,0,0,0,0,0,c1*c3,-s1*s2*c3+c1*s3,c1*s2*c3+s1*s3],[0,0,0,0,0,0,0,0,0,-c2*s3,s1*s2*s3+c1*c3,-c1*s2*s3+c3*s1],[0,0,0,0,0,0,0,0,0,-s2,-s1*c2,c1*c2]])


if __name__ == "__main__":
    B = Reader("ressorts//ressort_berx1.txt")
    #B = Reader("dynatest.txt")
#print(B.v)