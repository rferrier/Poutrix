
#This code provides the finite elements solver for 1 dimension problems

#importations
import math as mh
import numpy as np
from scipy.linalg import *

class Solver :

    def __init__(self,data):
        #recovery of data from the object poutre
        self.data = data
        #print(self.data)
        #data processing

        self.left = self.data[1][0]                   #Boundaries conditions
        self.right = self.data[1][len(self.data[1])-1]#

    def __str__(self):
        #Function print

        return ('I am the solver for the following problem :\n'+str(self.data))
        
    def create_mesh(self):

        self.mesh = Mesh(self.data[1],self.left,self.right,self.data[2])

    def solve_static_problem(self):

        self.create_mesh()
        self.mesh.compute_leng()
        self.mesh.compute_A()
        self.mesh.compute_b()
        #print(len(self.mesh.A))
        #print(len(self.mesh.b))
        #print(self.mesh.elements)
        self.u = solve(self.mesh.A,self.mesh.b)

    def solve_modal_problem(self):

        #Computation of the modal solution for the beam

        self.create_mesh()
        self.mesh.compute_leng()
        self.mesh.compute_A()
        self.mesh.compute_M()
        #self.M1 = linalg.inv(self.mesh.M)
        self.ms = eig(self.mesh.A,self.mesh.M)
        self.omeg = self.ms[0]
        self.u = self.ms[1]

class Mesh :

    def __init__(self,data,left,right,origin):
        #recovery of data from object solver and processing
        self.data = data[1:-1]
        self.origin = origin
        self.left = left
        self.right = right
        self.numberOfLastElem = len(self.data)-2

        #creation of nodes and elements :
        self.nodes = [0]*len(self.data)
        self.elements = [0]*(len(self.data)-1)
        self.create_nodes()
        self.create_elements()

    def __str__(self):
        #Function print

        return ('I am the following mesh :\n'+str(self.data))

    def create_nodes(self):
        
        j=0
        for i in self.data:
            
            self.nodes[j] = Node(i[0],i[2],j,self)
            j = j+1
            
    def create_elements(self):

        #Re-writing of the bundaries
        self.left.extend(['free','free','free','free','free','free'])
        self.right2 = ['free','free','free','free','free','free']
        self.right2.extend(self.right)
        
        #the bundaries :
        self.elements[0] = Sing_Element(self.data[0][1],0,[0,1],self,self.left,'left')
        self.elements[self.numberOfLastElem] = Sing_Element(self.data[self.numberOfLastElem][1],self.numberOfLastElem,[self.numberOfLastElem,self.numberOfLastElem+1],self,self.right2,'right')
        
        #creation of the center
        j=1
        #print(self.data)
        for i in self.data[1:-2] :

            self.elements[j] = Reg_Element(i[1],j,[j,j+1],self)
            j = j+1
            #print(self.elements[j])

    def compute_leng(self):
        #Computation of the size of A and b

        self.size = 6
        
        for i in self.elements:
            self.size = self.size + len(i.Ke)-6
            #print(len(i.Ke))

        self.size = self.size
        #print(self.size)

    def compute_A(self):

        #assembly of the matrix A
        self.A = np.mat([[0.]*self.size]*self.size)

        s = 0
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].Ke)):
                for k in range(len(self.elements[i].Ke)):
                    self.A[s+j,s+k] = self.A[s+j,s+k] + self.elements[i].Ke[j,k]

            #incrementation of the place for the Ke
            s = s + len(self.elements[i].Ke)-6
            
    def compute_F1(self,U):
        #building of F1, vector of nodal forces
        self.F1 = np.mat([[0.]]*12*(self.numberOfLastElem+1))
        #print(self.F1)
        
        #rebuilding of boundaries elem
        self.foelem = self.elements[:]
        self.foelem[0] = Reg_Element(self.data[0][1],0,[0,1],self)
        self.foelem[self.numberOfLastElem] = Reg_Element(self.data[self.numberOfLastElem][1],self.numberOfLastElem,[self.numberOfLastElem,self.numberOfLastElem+1],self)
        
        for i in range(len(self.foelem)):
            Ui = np.mat([[U[6*i,0]],[U[6*i+1,0]],[U[6*i+2,0]],[U[6*i+3,0]],[U[6*i+4,0]],[U[6*i+5,0]],[U[6*i+6,0]],[U[6*i+7,0]],[U[6*i+8,0]],[U[6*i+9,0]],[U[6*i+10,0]],[U[6*i+11,0]]])          
            self.foelem[i].compute_Fe(Ui)
            for j in range(12):
                self.F1[12*i+j,0] = self.foelem[i].Fe[j,0]
        #print(self.F1)

    def compute_T(self,U):
        """This function builds the vector of nodal forces"""
        self.T = np.mat([[0.]]*6*len(self.nodes))
        #rebuilding of boundaries elem
        self.foelem = self.elements[:]
        self.foelem[0] = Reg_Element(self.data[0][1],0,[0,1],self)
        self.foelem[self.numberOfLastElem] = Reg_Element(self.data[self.numberOfLastElem][1],self.numberOfLastElem,[self.numberOfLastElem,self.numberOfLastElem+1],self)
        
        for i in range(len(self.foelem)):
            Ui = np.mat([[U[6*i,0]],[U[6*i+1,0]],[U[6*i+2,0]],[U[6*i+3,0]],[U[6*i+4,0]],[U[6*i+5,0]],[U[6*i+6,0]],[U[6*i+7,0]],[U[6*i+8,0]],[U[6*i+9,0]],[U[6*i+10,0]],[U[6*i+11,0]]])          
            self.foelem[i].compute_Fe(Ui)
        
        for k in range(len(self.nodes)):
            self.nodes[k].compute_T()
            
            self.T[k*6:(k+1)*6,0] = self.nodes[k].T

    def compute_M(self):

        #computation of the elementary mass

        for i in self.elements:
            i.compute_Me()
        
        #assembly of the matrix M
        self.M = np.mat([[0.]*self.size]*self.size)

        s = 0
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].Me)):
                for k in range(len(self.elements[i].Me)):
                    self.M[s+j,s+k] = self.M[s+j,s+k] + self.elements[i].Me[j,k]

            #incrementation of the place for the Me
            s = s + len(self.elements[i].Me)-6

    def compute_b(self):

        #computation of the elementary loads

        for i in self.elements:
            i.compute_be()

        #assembly of the vector b
        self.b = np.mat([[0.]]*self.size)

        s = 0
        for i in range(len(self.elements)):
            for j in range(len(self.elements[i].be)):
                self.b[s+j,0] = self.b[s+j,0] + self.elements[i].be[j,0]

            #incrementation of the size of the Ke(n-1)
            s = s + len(self.elements[i].Ke)-6

class Node :

    def __init__(self,data1,data2,number,parent):
        #recovery of data from object mesh
        self.effort = data1
        self.lineff = data2
        self.number = number #number is the global numerotation
        self.parent = parent
        
    def compute_coords(self): #computation of the coords of every node
        if self.number != 0:
            self.x = self.parent.nodes[self.number-1].x + self.parent.elements[self.number-1].leng*mh.cos(self.parent.elements[self.number-1].theta1init)*mh.cos(self.parent.elements[self.number-1].theta2init)
            self.y = self.parent.nodes[self.number-1].y + self.parent.elements[self.number-1].leng*mh.sin(self.parent.elements[self.number-1].theta1init)*mh.cos(self.parent.elements[self.number-1].theta2init)
            self.z = self.parent.nodes[self.number-1].z - self.parent.elements[self.number-1].leng*mh.sin(self.parent.elements[self.number-1].theta2init)
            
        else :
            self.x = self.parent.origin[0]
            self.y = self.parent.origin[1]
            self.z = self.parent.origin[2]
            
        self.coords = (self.x,self.y,self.z)
        
    def compute_T(self):
        """This function computes the strength seen at the node"""
        
        long = len(self.parent.nodes)

        if self.number != 0 and self.number != long-1 : #reg element
            forcelem = self.parent.foelem[0].Fe
            self.T = (forcelem[-6:len(forcelem),0] - self.parent.foelem[self.number].Fe[0:6,0])/2
        elif self.number == 0 :
            self.T = -self.parent.foelem[0].Fe[0:6,0]
        else : #self.number = long-1
            forcelem = self.parent.foelem[long-2].Fe
            self.T = forcelem[-6:len(forcelem),0]

class Element :
    #base class for elements
    def __init__(self,data,number,linked_nodes,parent):
        #recovery of data from object mesh
        self.leng = data[0]
        self.lenginit = data[0] #usefull for GD algorithm : will never move
        self.linmass = data[1]
        self.p = data[1:7]#Those are the mechanical properties
        self.number = number
        self.nodes = linked_nodes
        self.parent = parent
        self.theta1 = data[7]
        self.theta2 = data[8]
        self.theta1init = data[7]
        self.theta2init = data[8]
        self.compute_P()
        
    def compute_P(self):
        #Computation of the base changing matrix
        #suppression of the numerical error for trigonometric :
        if self.theta1 == 3.1415926535/2 or self.theta1 == 3*3.1415926535/2:
            ct1 = 0.
        else :
            ct1 = mh.cos(self.theta1)

        if self.theta2 == 3.1415926535/2 or self.theta1 == 3*3.1415926535/2:
            ct2 = 0.
        else :
            ct2 = mh.cos(self.theta2)

        if self.theta1 == 3.1415926535 :
            st1 = 0.
        else :
            st1 = mh.sin(self.theta1)

        if self.theta2 == 3.1415926535 :
            st2 = 0.
        else :
            st2 = mh.sin(self.theta2)
        
        self.pas = np.mat([[ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0,0,0,0],[-st1,ct1,0,0,0,0,0,0,0,0,0,0],[ct1*st2,st1*st2,ct2,0,0,0,0,0,0,0,0,0],[0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0],[0,0,0,-st1,ct1,0,0,0,0,0,0,0],[0,0,0,ct1*st2,st1*st2,ct2,0,0,0,0,0,0],[0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0],[0,0,0,0,0,0,-st1,ct1,0,0,0,0],[0,0,0,0,0,0,ct1*st2,st1*st2,ct2,0,0,0],[0,0,0,0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2],[0,0,0,0,0,0,0,0,0,-st1,ct1,0],[0,0,0,0,0,0,0,0,0,ct1*st2,st1*st2,ct2]])
        self.pas1 = np.transpose(self.pas)

class Reg_Element(Element) :
    #An element where u is free
    def __init__(self,data,number,linked_nodes,parent):
        Element.__init__(self,data,number,linked_nodes,parent)
        self.compute_Ke()

    def compute_Ke(self):
        #we compute the elementary rigidity matrix
        ES = self.p[1]
        GJ = self.p[2]
        EIz = self.p[4]
        EIy = self.p[3]
        h = self.leng
        self.Keloc = np.mat([[ES/h,0,0,0,0,0,-ES/h,0,0,0,0,0],[0,12*EIy/h/h/h,0,0,0,6*EIy/h/h,0,- 12*EIy/h/h/h,0,0,0,6*EIy/h/h],[0,0,12*EIz/h/h/h,0,-6*EIz/h/h,0,0,0,- 12*EIz/h/h/h,0,-6*EIz/h/h,0],[0,0,0,GJ/h,0,0,0,0,0,-GJ/h,0,0],[0,0,-6*EIz/h/h,0,4*EIz/h,0,0,0,6*EIz/h/h,0,2*EIz/h,0],[0,6*EIy/h/h,0,0,0,4*EIy/h,0,- 6*EIy/h/h,0,0,0,2*EIy/h],[-ES/h,0,0,0,0,0,ES/h,0,0,0,0,0],[0,- 12*EIy/h/h/h,0,0,0,- 6*EIy/h/h,0,12*EIy/h/h/h,0,0,0,- 6*EIy/h/h],[0,0,- 12*EIz/h/h/h,0,6*EIz/h/h,0,0,0,12*EIz/h/h/h,0,6*EIz/h/h,0],[0,0,0,-GJ/h,0,0,0,0,0,GJ/h,0,0],[0,0,-6*EIz/h/h,0,2*EIz/h,0,0,0,6*EIz/h/h,0,4*EIz/h,0],[0,6*EIy/h/h,0,0,0,2*EIy/h,0,- 6*EIy/h/h,0,0,0,4*EIy/h]])
        self.Ke = self.pas1 * self.Keloc * self.pas

    def compute_Me(self):
        #we compute the elementary rigidity matrix
        rS = self.p[0]
        rJ = self.p[5]
        h = self.leng
        self.Meloc = np.mat([[rS*h/3,0,0,0,0,0,rS*h/6,0,0,0,0,0],[0,13/35*h*rS,0,0,0,11/210*h*h*rS,0,9/70*h*rS,0,0,0,-13/420*h*h*rS],[0,0,13/35*h*rS,0,-11/210*h*h*rS,0,0,0,9/70*h*rS,0,13/420*h*h*rS,0],[0,0,0,rJ*h/3,0,0,0,0,0,rJ*h/6,0,0],[0,0,-11/210*h*h*rS,0,h*h*h*rS/105,0,0,0,-13/420*h*h*rS,0,-h*h*h/140*rS,0],[0,11/210*h*h*rS,0,0,0,h*h*h/105*rS,0,13/420*h*h*rS,0,0,0,-h*h*h/140*rS],[rS*h/6,0,0,0,0,0,rS*h/3,0,0,0,0,0],[0,9/70*h*rS,0,0,0,13*h*h/420*rS,0,13/35*h*rS,0,0,0,-11/210*rS*h*h],[0,0,9/70*h*rS,0,-13/420*h*h*rS,0,0,0,13/35*h*rS,0,11/210*h*h*rS,0],[0,0,0,rJ*h/6,0,0,0,0,0,rJ*h/3,0,0],[0,0,13/420*h*h*rS,0,-h*h*h/140*rS,0,0,0,11*h*h/210*rS,0,h*h*h/105*rS,0],[0,-13/420*h*h*rS,0,0,0,-h*h*h/140*rS,0,-11/210*h*h*rS,0,0,0,h*h*h/105*rS]])
        self.Me = self.pas1 * self.Meloc * self.pas

    def compute_be(self):
        h = self.leng
        F1 = self.parent.nodes[self.nodes[0]].effort
        fl1 = self.parent.nodes[self.nodes[0]].lineff
        F2 = self.parent.nodes[self.nodes[1]].effort
        fl2 = self.parent.nodes[self.nodes[1]].lineff
        self.be = np.mat([[F1[0]/2 + fl1[0]*h/2],[F1[1]/2 + fl1[1]*h/2],[F1[2]/2 + fl1[2]*h/2],[F1[3]/2 + fl1[3]*h/2],[F1[4]/2 + fl1[4]*h/2],[F1[5]/2 + fl1[5]*h/2],[F2[0]/2 + fl2[0]*h/2],[F2[1]/2 + fl2[1]*h/2],[F2[2]/2 + fl2[2]*h/2],[F2[3]/2 + fl2[3]*h/2],[F2[4]/2 + fl2[4]*h/2],[F2[5]/2 + fl2[5]*h/2]])

    def compute_Fe(self,Ue):
        #Computation of elementary nodal force
        self.Fe = self.Ke * Ue

class Sing_Element(Element) :
    #An element where u=0 on the left
    def __init__(self,data,number,linked_nodes,parent,compo,side):
        Element.__init__(self,data,number,linked_nodes,parent)
        self.compo = compo
        self.side = side
        self.compute_Ke()

    def compute_Ke(self):
        #we compute the elementary rigidity matrix
        ES = self.p[1]
        GJ = self.p[2]
        EIz = self.p[4]
        EIy = self.p[3]
        self.properties = [ES,GJ,EIy,EIz]
        #computation of the penalization factor
        
        h = self.leng
         #list of penalization factors :
        self.pen = [1000000.*ES/h,1000000.*12*EIy/h/h/h,1000000.*12*EIz/h/h/h,1000000.*GJ/h,1000000.*4*EIz/h,1000000.*4*EIy/h,1000000.*ES/h,1000000.*12*EIy/h/h/h,1000000.*12*EIz/h/h/h,1000000.*GJ/h,1000000.*4*EIz/h,1000000.*4*EIy/h]
        self.Keloc = np.mat([[ES/h,0,0,0,0,0,-ES/h,0,0,0,0,0],[0,12*EIy/h/h/h,0,0,0,6*EIy/h/h,0,- 12*EIy/h/h/h,0,0,0,6*EIy/h/h],[0,0,12*EIz/h/h/h,0,-6*EIz/h/h,0,0,0,- 12*EIz/h/h/h,0,-6*EIz/h/h,0],[0,0,0,GJ/h,0,0,0,0,0,-GJ/h,0,0],[0,0,-6*EIz/h/h,0,4*EIz/h,0,0,0,6*EIz/h/h,0,2*EIz/h,0],[0,6*EIy/h/h,0,0,0,4*EIy/h,0,- 6*EIy/h/h,0,0,0,2*EIy/h],[-ES/h,0,0,0,0,0,ES/h,0,0,0,0,0],[0,- 12*EIy/h/h/h,0,0,0,- 6*EIy/h/h,0,12*EIy/h/h/h,0,0,0,- 6*EIy/h/h],[0,0,- 12*EIz/h/h/h,0,6*EIz/h/h,0,0,0,12*EIz/h/h/h,0,6*EIz/h/h,0],[0,0,0,-GJ/h,0,0,0,0,0,GJ/h,0,0],[0,0,-6*EIz/h/h,0,2*EIz/h,0,0,0,6*EIz/h/h,0,4*EIz/h,0],[0,6*EIy/h/h,0,0,0,2*EIy/h,0,- 6*EIy/h/h,0,0,0,4*EIy/h]])
        #removing of the forbitten movements with the penalization method
        for i in range(12):
            if type(self.compo[i]) != str:
                self.Keloc[i,i] = self.Keloc[i,i] + self.pen[i]
        
        self.Ke = self.pas1 * self.Keloc * self.pas

    def compute_Me(self):
        #we compute the elementary rigidity matrix
        rS = self.p[0]
        rJ = self.p[5]
        h = self.leng
        
        self.Meloc = np.mat([[rS*h/3,0,0,0,0,0,rS*h/6,0,0,0,0,0],[0,13/35*h*rS,0,0,0,11/210*h*h*rS,0,9/70*h*rS,0,0,0,-13/420*h*h*rS],[0,0,13/35*h*rS,0,-11/210*h*h*rS,0,0,0,9/70*h*rS,0,13/420*h*h*rS,0],[0,0,0,rJ*h/3,0,0,0,0,0,rJ*h/6,0,0],[0,0,-11/210*h*h*rS,0,h*h*h*rS/105,0,0,0,-13/420*h*h*rS,0,-h*h*h/140*rS,0],[0,11/210*h*h*rS,0,0,0,h*h*h/105*rS,0,13/420*h*h*rS,0,0,0,-h*h*h/140*rS],[rS*h/6,0,0,0,0,0,rS*h/3,0,0,0,0,0],[0,9/70*h*rS,0,0,0,13*h*h/420*rS,0,13/35*h*rS,0,0,0,-11/210*rS*h*h],[0,0,9/70*h*rS,0,-13/420*h*h*rS,0,0,0,13/35*h*rS,0,11/210*h*h*rS,0],[0,0,0,rJ*h/6,0,0,0,0,0,rJ*h/3,0,0],[0,0,13/420*h*h*rS,0,-h*h*h/140*rS,0,0,0,11*h*h/210*rS,0,h*h*h/105*rS,0],[0,-13/420*h*h*rS,0,0,0,-h*h*h/140*rS,0,-11/210*h*h*rS,0,0,0,h*h*h/105*rS]])
        self.Me = self.pas1 * self.Meloc * self.pas
    def compute_be(self):

        h = self.leng
        F1 = self.parent.nodes[self.nodes[0]].effort
        fl1 = self.parent.nodes[self.nodes[0]].lineff
        F2 = self.parent.nodes[self.nodes[1]].effort
        fl2 = self.parent.nodes[self.nodes[1]].lineff
        
        if self.side == 'left' :
            self.be = np.mat([[F1[0] + fl1[0]*h],[F1[1] + fl1[1]*h],[F1[2] + fl1[2]*h],[F1[3] + fl1[3]*h],[F1[4] + fl1[4]*h],[F1[5] + fl1[5]*h],[F2[0]/2 + fl2[0]*h/2],[F2[1]/2 + fl2[1]*h/2],[F2[2]/2 + fl2[2]*h/2],[F2[3]/2 + fl2[3]*h/2],[F2[4]/2 + fl2[4]*h/2],[F2[5]/2 + fl2[5]*h/2]])
        else :
            self.be = np.mat([[F1[0]/2 + fl1[0]*h/2],[F1[1]/2 + fl1[1]*h/2],[F1[2]/2 + fl1[2]*h/2],[F1[3]/2 + fl1[3]*h/2],[F1[4]/2 + fl1[4]*h/2],[F1[5]/2 + fl1[5]*h/2],[F2[0] + fl2[0]*h],[F2[1] + fl2[1]*h],[F2[2] + fl2[2]*h],[F2[3] + fl2[3]*h],[F2[4] + fl2[4]*h],[F2[5] + fl2[5]*h]])

        #removing of the forbitten movements
        for i in range(12):
            if type(self.compo[i]) != str:
                self.be[i] = self.be[i] + self.pen[i]*self.compo[i]
                