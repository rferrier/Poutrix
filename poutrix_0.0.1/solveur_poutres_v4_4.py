
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

        self.left = self.data[1][0]                #Boundaries conditions
        self.right = self.data[1][len(self.data[1])-1]#

    def __str__(self):
        #Function print

        return ('I am the solver for the following problem :\n'+str(self.data))
        
    def create_mesh(self):

        self.mesh = Mesh(self.data[1],self.left,self.right)

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

    def __init__(self,data,left,right):
        #recovery of data from object solver and processing
        self.data = data[1:-1]
        #print(self.data)
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
            
            self.nodes[j] = Node(i[0],j)
            j = j+1
            
    def create_elements(self):

        #Re-writing of the bundaries
        self.left.extend([1,1,1,1,1,1])
        self.right2 = [1,1,1,1,1,1]
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

    def __init__(self,data,number):
        #recovery of data from object mesh
        self.effort = data
        self.number = number #number is the global numerotation

class Element :
    #base class for elements
    def __init__(self,data,number,linked_nodes,parent):
        #recovery of data from object mesh
        self.leng = data[0]
        self.linmass = data[1]
        self.p = data[1:]#Those are the mechanical properties
        self.number = number
        self.nodes = linked_nodes
        self.parent = parent     

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
        self.Ke = np.mat([[ES/h,0,0,0,0,0,-ES/h,0,0,0,0,0],[0,12*EIy/h/h/h,0,0,0,6*EIy/h/h,0,- 12*EIy/h/h/h,0,0,0,6*EIy/h/h],[0,0,12*EIz/h/h/h,0,-6*EIz/h/h,0,0,0,- 12*EIz/h/h/h,0,-6*EIz/h/h,0],[0,0,0,GJ/h,0,0,0,0,0,-GJ/h,0,0],[0,0,-6*EIz/h/h,0,4*EIz/h,0,0,0,6*EIz/h/h,0,2*EIz/h,0],[0,6*EIy/h/h,0,0,0,4*EIy/h,0,- 6*EIy/h/h,0,0,0,2*EIy/h],[-ES/h,0,0,0,0,0,ES/h,0,0,0,0,0],[0,- 12*EIy/h/h/h,0,0,0,- 6*EIy/h/h,0,12*EIy/h/h/h,0,0,0,- 6*EIy/h/h],[0,0,- 12*EIz/h/h/h,0,6*EIz/h/h,0,0,0,12*EIz/h/h/h,0,6*EIz/h/h,0],[0,0,0,-GJ/h,0,0,0,0,0,GJ/h,0,0],[0,0,-6*EIz/h/h,0,2*EIz/h,0,0,0,6*EIz/h/h,0,4*EIz/h,0],[0,6*EIy/h/h,0,0,0,2*EIy/h,0,- 6*EIy/h/h,0,0,0,4*EIy/h]])

    def compute_Me(self):
        #we compute the elementary rigidity matrix
        rS = self.p[0]
        rJ = self.p[5]
        h = self.leng
        self.Me = np.mat([[rS*h/3,0,0,0,0,0,rS*h/6,0,0,0,0,0],[0,13/35*h*rS,0,0,0,11/210*h*h*rS,0,9/70*h*rS,0,0,0,-13/420*h*h*rS],[0,0,13/35*h*rS,0,-11/210*h*h*rS,0,0,0,9/70*h*rS,0,13/420*h*h*rS,0],[0,0,0,rJ*h/3,0,0,0,0,0,rJ*h/6,0,0],[0,0,-11/210*h*h*rS,0,h*h*h*rS/105,0,0,0,-13/420*h*h*rS,0,-h*h*h/140*rS,0],[0,11/210*h*h*rS,0,0,0,h*h*h/105*rS,0,13/420*h*h*rS,0,0,0,-h*h*h/140*rS],[rS*h/6,0,0,0,0,0,rS*h/3,0,0,0,0,0],[0,9/70*h*rS,0,0,0,13*h*h/420*rS,0,13/35*h*rS,0,0,0,-11/210*rS*h*h],[0,0,9/70*h*rS,0,-13/420*h*h*rS,0,0,0,13/35*h*rS,0,11/210*h*h*rS,0],[0,0,0,rJ*h/6,0,0,0,0,0,rJ*h/3,0,0],[0,0,13/420*h*h*rS,0,-h*h*h/140*rS,0,0,0,11*h*h/210*rS,0,h*h*h/105*rS,0],[0,-13/420*h*h*rS,0,0,0,-h*h*h/140*rS,0,-11/210*h*h*rS,0,0,0,h*h*h/105*rS]])

    def compute_be(self):
        self.be = np.mat([[self.parent.nodes[self.nodes[0]].effort[0]/2],[self.parent.nodes[self.nodes[0]].effort[1]/2],[self.parent.nodes[self.nodes[0]].effort[2]/2],[self.parent.nodes[self.nodes[0]].effort[3]],[- self.parent.nodes[self.nodes[0]].effort[4]/2],[self.parent.nodes[self.nodes[0]].effort[5]/2],[self.parent.nodes[self.nodes[1]].effort[0]/2],[self.parent.nodes[self.nodes[1]].effort[1]/2],[self.parent.nodes[self.nodes[1]].effort[2]/2],[self.parent.nodes[self.nodes[1]].effort[3]/2],[- self.parent.nodes[self.nodes[1]].effort[4]/2],[self.parent.nodes[self.nodes[1]].effort[5]/2]])

    def compute_P(self,theta1,theta2):
        #Computation of the base changing matrix
        #suppression of the numerical error for trigonometric :
        if theta1 == 3.1415926535/2 or theta1 == 3*3.1415926535/2:
            ct1 = 0.
        else :
            ct1 = mh.cos(theta1)

        if theta2 == 3.1415926535/2 or theta1 == 3*3.1415926535/2:
            ct2 = 0.
        else :
            ct2 = mh.cos(theta2)

        if theta1 == 3.1415926535 :
            st1 = 0.
        else :
            st1 = mh.sin(theta1)

        if theta2 == 3.1415926535 :
            st2 = 0.
        else :
            st2 = mh.sin(theta2)
        
        self.pas = 1/2 * np.mat([[ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0,0,0,0],[-st1,ct1,0,0,0,0,0,0,0,0,0,0],[ct1*st2,st1*st2,ct2,0,0,0,0,0,0,0,0,0],[0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0],[0,0,0,-st1,ct1,0,0,0,0,0,0,0],[0,0,0,ct1*st2,st1*st2,ct2,0,0,0,0,0,0],[0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0],[0,0,0,0,0,0,-st1,ct1,0,0,0,0],[0,0,0,0,0,0,ct1*st2,st1*st2,ct2,0,0,0],[0,0,0,0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2],[0,0,0,0,0,0,0,0,0,-st1,ct1,0],[0,0,0,0,0,0,0,0,0,ct1*st2,st1*st2,ct2]])
  
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
        self.pen = 1000. * max(ES/h,12*EIy/h/h/h,12*EIz/h/h/h,GJ/h,4*EIz/h,4*EIy/h)
        self.Ke = np.mat([[ES/h,0,0,0,0,0,-ES/h,0,0,0,0,0],[0,12*EIy/h/h/h,0,0,0,6*EIy/h/h,0,- 12*EIy/h/h/h,0,0,0,6*EIy/h/h],[0,0,12*EIz/h/h/h,0,-6*EIz/h/h,0,0,0,- 12*EIz/h/h/h,0,-6*EIz/h/h,0],[0,0,0,GJ/h,0,0,0,0,0,-GJ/h,0,0],[0,0,-6*EIz/h/h,0,4*EIz/h,0,0,0,6*EIz/h/h,0,2*EIz/h,0],[0,6*EIy/h/h,0,0,0,4*EIy/h,0,- 6*EIy/h/h,0,0,0,2*EIy/h],[-ES/h,0,0,0,0,0,ES/h,0,0,0,0,0],[0,- 12*EIy/h/h/h,0,0,0,- 6*EIy/h/h,0,12*EIy/h/h/h,0,0,0,- 6*EIy/h/h],[0,0,- 12*EIz/h/h/h,0,6*EIz/h/h,0,0,0,12*EIz/h/h/h,0,6*EIz/h/h,0],[0,0,0,-GJ/h,0,0,0,0,0,GJ/h,0,0],[0,0,-6*EIz/h/h,0,2*EIz/h,0,0,0,6*EIz/h/h,0,4*EIz/h,0],[0,6*EIy/h/h,0,0,0,2*EIy/h,0,- 6*EIy/h/h,0,0,0,4*EIy/h]])
        #removing of the forbitten movements with the penalization method
        for i in range(12):
            if type(self.compo[i]) == float:
                self.Ke[i,i] = self.Ke[i,i] + self.pen
                #print(self.Ke[i,i])

    def compute_Me(self):
        #we compute the elementary rigidity matrix
        rS = self.p[0]
        rJ = self.p[5]
        h = self.leng
        
        self.Me = np.mat([[rS*h/3,0,0,0,0,0,rS*h/6,0,0,0,0,0],[0,13/35*h*rS,0,0,0,11/210*h*h*rS,0,9/70*h*rS,0,0,0,-13/420*h*h*rS],[0,0,13/35*h*rS,0,-11/210*h*h*rS,0,0,0,9/70*h*rS,0,13/420*h*h*rS,0],[0,0,0,rJ*h/3,0,0,0,0,0,rJ*h/6,0,0],[0,0,-11/210*h*h*rS,0,h*h*h*rS/105,0,0,0,-13/420*h*h*rS,0,-h*h*h/140*rS,0],[0,11/210*h*h*rS,0,0,0,h*h*h/105*rS,0,13/420*h*h*rS,0,0,0,-h*h*h/140*rS],[rS*h/6,0,0,0,0,0,rS*h/3,0,0,0,0,0],[0,9/70*h*rS,0,0,0,13*h*h/420*rS,0,13/35*h*rS,0,0,0,-11/210*rS*h*h],[0,0,9/70*h*rS,0,-13/420*h*h*rS,0,0,0,13/35*h*rS,0,11/210*h*h*rS,0],[0,0,0,rJ*h/6,0,0,0,0,0,rJ*h/3,0,0],[0,0,13/420*h*h*rS,0,-h*h*h/140*rS,0,0,0,11*h*h/210*rS,0,h*h*h/105*rS,0],[0,-13/420*h*h*rS,0,0,0,-h*h*h/140*rS,0,-11/210*h*h*rS,0,0,0,h*h*h/105*rS]])

    def compute_be(self):

        if self.side == 'left' :
            self.be = np.mat([[self.parent.nodes[self.nodes[0]].effort[0]],[self.parent.nodes[self.nodes[0]].effort[1]],[self.parent.nodes[self.nodes[0]].effort[2]],[self.parent.nodes[self.nodes[0]].effort[3]],[- self.parent.nodes[self.nodes[0]].effort[4]],[self.parent.nodes[self.nodes[0]].effort[5]],[self.parent.nodes[self.nodes[1]].effort[0]/2],[self.parent.nodes[self.nodes[1]].effort[1]/2],[self.parent.nodes[self.nodes[1]].effort[2]/2],[self.parent.nodes[self.nodes[1]].effort[3]/2],[- self.parent.nodes[self.nodes[1]].effort[4]/2],[self.parent.nodes[self.nodes[1]].effort[5]/2]])
        else :
            self.be = np.mat([[self.parent.nodes[self.nodes[0]].effort[0]/2],[self.parent.nodes[self.nodes[0]].effort[1]/2],[self.parent.nodes[self.nodes[0]].effort[2]/2],[self.parent.nodes[self.nodes[0]].effort[3]],[- self.parent.nodes[self.nodes[0]].effort[4]/2],[self.parent.nodes[self.nodes[0]].effort[5]/2],[self.parent.nodes[self.nodes[1]].effort[0]],[self.parent.nodes[self.nodes[1]].effort[1]],[self.parent.nodes[self.nodes[1]].effort[2]],[self.parent.nodes[self.nodes[1]].effort[3]],[- self.parent.nodes[self.nodes[1]].effort[4]],[self.parent.nodes[self.nodes[1]].effort[5]]])

            
        #removing of the forbitten movements
        for i in range(12):
            if type(self.compo[i]) == float:
                self.be[i] = self.be[i] + self.pen*self.compo[i]
                
    def compute_P(self,theta1,theta2):
        #Computation of the base changing matrix
        #suppression of the numerical error for trigonometric :
        if theta1 == 3.1415926535/2 or theta1 == -3.1415926535/2:
            ct1 = 0.
        else :
            ct1 = mh.cos(theta1)

        if theta2 == 3.1415926535/2 or theta1 == -3.1415926535/2:
            ct2 = 0.
        else :
            ct2 = mh.cos(theta2)

        if theta1 == 3.1415926535 :
            st1 = 0.
        else :
            st1 = mh.sin(theta1)

        if theta2 == 3.1415926535 :
            st2 = 0.
        else :
            st2 = mh.sin(theta2)

        if self.side == 'left' :
            self.pas = np.mat([[ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0,0,0,0],[-st1,ct1,0,0,0,0,0,0,0,0,0,0],[ct1*st2,st1*st2,ct2,0,0,0,0,0,0,0,0,0],[0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0,0,0,0],[0,0,0,-st1,ct1,0,0,0,0,0,0,0],[0,0,0,ct1*st2,st1*st2,ct2,0,0,0,0,0,0],[0,0,0,0,0,0,ct1*ct2/2,st1*ct2/2,-st2/2,0,0,0],[0,0,0,0,0,0,-st1/2,ct1/2,0,0,0,0],[0,0,0,0,0,0,ct1*st2/2,st1*st2/2,ct2/2,0,0,0],[0,0,0,0,0,0,0,0,0,ct1*ct2/2,st1*ct2/2,-st2/2],[0,0,0,0,0,0,0,0,0,-st1/2,ct1/2,0],[0,0,0,0,0,0,0,0,0,ct1*st2/2,st1*st2/2,ct2/2]])

        else :
            self.pas = np.mat([[ct1*ct2/2,st1*ct2/2,-st2/2,0,0,0,0,0,0,0,0,0],[-st1/2,ct1/2,0,0,0,0,0,0,0,0,0,0],[ct1/2*st2/2,st1*st2/2,ct2/2,0,0,0,0,0,0,0,0,0],[0,0,0,ct1*ct2/2,st1*ct2/2,-st2/2,0,0,0,0,0,0],[0,0,0,-st1/2,ct1/2,0,0,0,0,0,0,0],[0,0,0,ct1*st2/2,st1*st2/2,ct2/2,0,0,0,0,0,0],[0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2,0,0,0],[0,0,0,0,0,0,-st1,ct1,0,0,0,0],[0,0,0,0,0,0,ct1*st2,st1*st2,ct2,0,0,0],[0,0,0,0,0,0,0,0,0,ct1*ct2,st1*ct2,-st2],[0,0,0,0,0,0,0,0,0,-st1,ct1,0],[0,0,0,0,0,0,0,0,0,ct1*st2,st1*st2,ct2]])

#A = Solver([[0.,0.,0.,0.,0.,0.],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,0.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[[0.,0.,0.,1.,0.,0.],[1/2,1.,1.,1.,1.,1.,1.,1.,1.]],[1,1,1,1,1,1]]) #[encastré,effort nodal nul *7,effort au bout : 10, libre]

#[ux,uy,uz,mx,my,mz],[[Fx,Fy,Fz,Mx,My,Mz],[h,ES,GJ,EIy,EIz,rho S,rho J, rho Iy,rho Iz]]

"""A.solve_static_problem()
print(A.u)"""
#A.solve_modal_problem()
#print(A.u[:,0])
#for i in range(20,35):
#    print(A.u[:,i])
#    print(sqrt(A.omeg[i]))
#print(sqrt(A.omeg[0]))
#print(len(A.u))
#print(len(A.mesh.elements))

"""ES = 1
GJ = 1
EIz = 1
EIy = 1
h = 1
theta1 = 1
theta2 = 1
rS = 1
rJ = 1
rIz = 1
rIy = 1
Ke = np.mat([[ES/h,0,0,0,0,0,-ES/h,0,0,0,0,0],[0,12*EIy/h/h/h,0,0,0,6*EIy/h/h,0,- 12*EIy/h/h/h,0,0,0,6*EIy/h/h],[0,0,12*EIz/h/h/h,0,6*EIz/h/h,0,0,0,- 12*EIz/h/h/h,0,6*EIz/h/h,0],[0,0,0,GJ/h,0,0,0,0,0,-GJ/h,0,0],[0,0,6*EIz/h/h,0,4*EIz/h,0,0,0,- 6*EIz/h/h,0,2*EIz/h,0],[0,6*EIy/h/h,0,0,0,4*EIy/h,0,- 6*EIy/h/h,0,0,0,2*EIy/h],[-ES/h,0,0,0,0,0,ES/h,0,0,0,0,0],[0,- 12*EIy/h/h/h,0,0,0,- 6*EIy/h/h,0,12*EIy/h/h/h,0,0,0,- 6*EIy/h/h],[0,0,- 12*EIz/h/h/h,0,- 6*EIz/h/h,0,0,0,12*EIz/h/h/h,0,- 6*EIz/h/h,0],[0,0,0,-GJ/h,0,0,0,0,0,GJ/h,0,0],[0,0,6*EIz/h/h,0,2*EIz/h,0,0,0,- 6*EIz/h/h,0,4*EIz/h,0],[0,6*EIy/h/h,0,0,0,2*EIy/h,0,- 6*EIy/h/h,0,0,0,4*EIy/h]])
#pas = 1/2*np.mat([[cos(theta1)*cos(theta2),sin(theta1)*cos(theta2),-sin(theta2),0,0,0,0,0,0,0,0,0],[-sin(theta1),cos(theta1),0,0,0,0,0,0,0,0,0,0],[cos(theta1)*sin(theta2),sin(theta1)*sin(theta2),cos(theta2),0,0,0,0,0,0,0,0,0],[0,0,0,cos(theta1)*cos(theta2),sin(theta1)*cos(theta2),-sin(theta2),0,0,0,0,0,0],[0,0,0,-sin(theta1),cos(theta1),0,0,0,0,0,0,0],[0,0,0,cos(theta1)*sin(theta2),sin(theta1)*sin(theta2),cos(theta2),0,0,0,0,0,0],[0,0,0,0,0,0,cos(theta1)*cos(theta2),sin(theta1)*cos(theta2),-sin(theta2),0,0,0],[0,0,0,0,0,0,-sin(theta1),cos(theta1),0,0,0,0],[0,0,0,0,0,0,cos(theta1)*sin(theta2),sin(theta1)*sin(theta2),cos(theta2),0,0,0],[0,0,0,0,0,0,0,0,0,cos(theta1)*cos(theta2),sin(theta1)*cos(theta2),-sin(theta2)],[0,0,0,0,0,0,0,0,0,-sin(theta1),cos(theta1),0],[0,0,0,0,0,0,0,0,0,cos(theta1)*sin(theta2),sin(theta1)*sin(theta2),cos(theta2)]])
Me = np.mat([[rS/3,0,0,0,0,0,rS/6,0,0,0,0,0],[0,13/35*h*rIy,0,0,0,11/210*h*h*rIy,0,-9/70*h*rIy,0,0,0,-13/420*h*h*rIy],[0,0,13/35*h*rIz,0,11/210*h*h*rIz,0,0,0,-9/70*h*rIz,0,-13/420*h*h*rIz,0],[0,0,0,rJ/3,0,0,0,0,0,rJ/6,0,0],[0,0,11/210*h*h*rIz,0,h*h*h*rIz/105,0,0,0,13/420*h*h*rIz,0,-h*h*h/140*rIz,0],[0,11/210*h*h*rIy,0,0,0,h*h*h/105*rIy,0,13/420*h*h*rIy,0,0,0,-h*h*h/140*rIy],[rS/6,0,0,0,0,0,rS/3,0,0,0,0,0],[0,-9/70*h*rIy,0,0,0,13*h*h/420*rIy,0,13/35*h*rIy,0,0,0,-11/210*rIy*h*h],[0,0,-9/70*h*rIz,0,13/420*h*h*rIz,0,0,0,13/35*h*rIz,0,-11/210*h*h*rIz,0],[0,0,0,rJ/6,0,0,0,0,0,rJ/3,0,0],[0,0,-13/420*h*h*rIz,0,-h*h*h/140*rIz,0,0,0,-11*h*h/210*rIz,0,h*h*h/105*rIz,0],[0,-13/420*h*h*rIy,0,0,0,-h*h*h/140*rIy,0,-11/210*h*h*rIy,0,0,0,h*h*h/105*rIy]])
Mat = linalg.inv(Me)*Ke
#print(Mat)
print(eig(Ke,Me)[1][:,4])#[1][:,0]"""
