#This code provides the parser for simple beams

import sp

def parser():
    commentaire=sp.R(r'#.*')
    blancs=sp.R(r'\s+')
    nom=sp.R(r'[a-zA-Z_]')
    nomPout=sp.R(r'\w+')
    component=sp.R(r'Ux')|sp.R(r'Uy')|sp.R(r'Uz')|sp.R(r'Theta_x')|sp.R(r'Theta_y')|sp.R(r'Theta_z')
    nombre=sp.R(r'[0-9]+') / int  #convertit le nombre en entier
    nbVirgule=(sp.R(r'-?'+'[0-9]*'+'.'+'[0-9]+')|sp.R(r'-?'+'[0-9]+')|sp.R(r'-?'+'[0-9]+'+'.')) /float
    
    baliseDebutPoutre=sp.K(r'<')
    baliseDebutBoundaries=sp.K(r'<BOUNDARIES>')  
    baliseDebutDimensions=sp.K(r'<DIMENSIONS>')
    baliseDebutLeft=sp.K(r'<Left>')
    baliseDebutRight=sp.K(r'<Right>')
    baliseDebutLoading=sp.K(r'<LOADING>')
    baliseDebutLin=sp.K(r'<Lineic>')
    baliseDebutSolveur=sp.K(r'<SOLVER>')
    baliseDebutLink=sp.K(r'<LINK>')
    baliseDebutOutput=sp.K(r'<OUTPUT>')
    baliseDebutOutDyn=sp.K(r'<TIME_OUTPUT>')
    
    baliseFinPoutre=sp.K(r'</BEAM>')
    baliseFinBoundaries=sp.K(r'</BOUNDARIES>')
    baliseFinDimensions=sp.K(r'</DIMENSIONS>')
    baliseFinLeft=sp.K(r'</Left>')
    baliseFinRight=sp.K(r'</Right>')
    baliseFinLoading=sp.K(r'</LOADING>')
    baliseFinLin=sp.K(r'</Lineic>')
    baliseFinSolveur=sp.K(r'</SOLVER>')
    baliseFinLink=sp.K(r'</LINK>')
    baliseFinOutput=sp.K(r'</OUTPUT>')
    baliseFinOutDyn=sp.K(r'</TIME_OUTPUT>')
    
    egal=sp.K(r'=')
    depoin=sp.K(r':')
    
#########################################################################################################
    with sp.Separator(blancs|commentaire):  #on fait abstraction des blancs et des commentaires
        assembly=sp.Rule()
        
        beam=sp.Rule()
        link=sp.Rule()
        blocBoundaries=sp.Rule()
        blocDimensions=sp.Rule()
        blocLoading=sp.Rule()
        blocSolveur=sp.Rule()
        blocOutput=sp.Rule()
        blocOutDyn=sp.Rule()

        blocLeft=sp.Rule()
        blocRight=sp.Rule()
        blocLineic=sp.Rule()

        blocCoL=sp.Rule()
        blocCoR=sp.Rule()
        blocL=sp.Rule()
        blocE=sp.Rule()
        blocNu=sp.Rule()
        blocS=sp.Rule()
        blocIy=sp.Rule()
        blocIz=sp.Rule()
        blocRho=sp.Rule()
        blocNb=sp.Rule()
        blocU=sp.Rule()
        blocF=sp.Rule()
        
        blocUx=sp.Rule()
        blocUy=sp.Rule()
        blocUz=sp.Rule()
        blocTx=sp.Rule()
        blocTy=sp.Rule()
        blocTz=sp.Rule()

        blocFx=sp.Rule()
        blocFy=sp.Rule()
        blocFz=sp.Rule()
        blocMx=sp.Rule()
        blocMy=sp.Rule()
        blocMz=sp.Rule()

        blocSolve=sp.Rule()
        blocName1=sp.Rule()
        blocName2=sp.Rule()
        blocName=sp.Rule()
        blocx=sp.Rule()
        blocBNam=sp.Rule()
        
        blocTech=sp.Rule()
        blocDur=sp.Rule()
        blocComp=sp.Rule()
        
        blocAxis=sp.Rule()
        
        L=sp.Rule()
        E=sp.Rule()
        Type=sp.Rule()
        N=sp.Rule()

##########################################################################################################
        assembly|=beam[0:] & link[0:] & blocSolveur 

        link|=baliseDebutLink & blocName1 & blocx & blocName2 & blocx & blocSolve & blocAxis[0:1] & baliseFinLink
        beam|=baliseDebutPoutre & 'Beam Name' & egal & '"' & nomPout &'" >' & blocDimensions & blocBoundaries & blocLoading & baliseFinPoutre
        
        blocBoundaries|=baliseDebutBoundaries & blocLeft & blocRight & baliseFinBoundaries
        blocDimensions|=baliseDebutDimensions & blocCoL & blocCoR & blocE & blocNu & blocS & blocIy & blocIz & blocRho & blocNb & baliseFinDimensions
        blocLoading|=baliseDebutLoading & blocLeft & blocRight & blocLineic & baliseFinLoading
        blocSolveur|=baliseDebutSolveur & blocSolve & blocDur[0:1] & blocTech[0:1] & blocName & blocOutput[0:] & blocOutDyn[0:] & baliseFinSolveur
        blocOutput|=baliseDebutOutput & blocBNam & blocU & blocF & baliseFinOutput
        blocOutDyn|=baliseDebutOutDyn & blocBNam & blocx& blocComp & baliseFinOutDyn

        blocLeft|=baliseDebutLeft & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinLeft
        blocRight|=baliseDebutRight & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinRight
        blocLineic|=baliseDebutLin & blocFx[0:1] & blocFy[0:1] & blocFz[0:1] & blocMx[0:1] & blocMy[0:1] & blocMz[0:1] & baliseFinLin

        blocCoL|='Left' & depoin & L & ',' & L & ',' & L & ';'
        blocCoR|='Right' & depoin & L & ',' & L & ',' & L & ';'
        blocL|='L' & egal & L &';'
        blocE|='E' & egal & E &';'
        blocNu|='Nu' & egal & E &';'
        blocS|='S' & egal & E &';'
        blocIy|='Iy' & egal & E &';'
        blocIz|='Iz' & egal & E &';'
        blocRho|='Rho' & egal & E &';'
        blocNb|='n' & egal & N &';'
        blocU|='U' & egal & component[0:6] & ';'
        blocF|='F' & egal & component[0:6] & ';'
        
        blocUx|='Ux' & egal & L & ';'
        blocUy|='Uy' & egal & L & ';'
        blocUz|='Uz' & egal & L & ';'
        blocTx|='Theta_x' & egal & L & ';'
        blocTy|='Theta_y' & egal & L & ';'
        blocTz|='Theta_z' & egal & L & ';'

        blocFx|='fx' & egal & E & ';'
        blocFy|='fy' & egal & E & ';'
        blocFz|='fz' & egal & E & ';'
        blocMx|='mx' & egal & E & ';'
        blocMy|='my' & egal & E & ';'
        blocMz|='mz' & egal & E & ';'

        blocSolve|='Type' & egal & Type & N[0:1] & ';'
        blocName1|='Beam1' & egal & Type & ';'
        blocName2|='Beam2' & egal & Type & ';'
        blocName|='Output_name' & egal & Type & ';'
        blocx|='X' & egal & L & ';'
        blocBNam|='Beam' & egal & nomPout & ';'
        
        blocTech|='dt' & egal & L & ';'
        blocDur|='T' & egal & L & ';'
        blocComp|='Component' & egal & component & ';'
        
        blocAxis|='Axis' & egal & nom & ';'
        
        L|=nbVirgule
        E|=nbVirgule|nomPout
        Type|=nomPout
        N|=nombre

    return assembly

if __name__ == "__main__" :

    traduire=parser()

    texte1="""
<Beam Name="Poutre1" >

<DIMENSIONS>
Left : 0.0 , 0.0 , 0.0 ;#Coords of the beginning of the beam (m)
Right: 0.3 , 0.0 , 0.0 ;#Coords of the end of the beam (m)
E= 210000000000.0 ;#Young modulus (Pa)
Nu = 0.3 ;#Fish coefficient
S = 0.001 ;#Section(m^2)
Iy = 0.00000000001 ;#Quadratic moment (m^4)
Iz =dataiz ;#Quadratic moment (m^4)
Rho = 7800 ; #Volumic mass (kg.m-3)
n = 10 ; #Number of elements between 2 sides of a beam
</DIMENSIONS>

<BOUNDARIES> #boundaries conditions for displacement
<Left>
Ux = 0.0 ;
Uy = 0.0 ;
Uz = 0.0 ;
Theta_x = 0.0 ;
Theta_y=0.0 ;
Theta_z = 0.0 ;
</Left>
<Right>
Ux = 0.0 ;
</Right>
</BOUNDARIES>

<LOADING> #efforts on the beam
<Left> #ponctual effort on the left
Uy = 0.0 ;
Theta_x = 0.0 ;
</Left>
<Right> #ponctual effort on the right
Theta_z = 0.0 ;
</Right>
<Lineic> #Lineic effort
fy = -1.0;#Nm-1
mx = datamx;
</Lineic>
</LOADING>

</BEAM>


<Beam Name="Poutre2" >

<DIMENSIONS>
Left : 0.3 , 0.0 , 0.0 ;#Coords of the beginning of the beam (m)
Right : 0.4 , 0.0 , 0.0 ;#Coords of the end of the beam (m)
E = 210000000000.0 ;#Young modulus (Pa)
Nu = 0.3 ;#Fish coefficient
S = 0.001 ;#Section(m^2)
Iy = 0.00000000001 ;#Quadratic moment (m^4)
Iz = 0.0000000002 ;#Quadratic moment (m^4)
Rho = 7800 ; #Volumic mass (kg.m-3)
n = 10 ; #Number of elements between 2 sides of a beam
</DIMENSIONS>

<BOUNDARIES> #boundaries conditions for displacement
<Left>
Ux = 0.0 ;
Theta_y = 0.0 ;
</Left>
<Right>
Ux = 0 ;
Uy = 0.0 ;
Uz = 0.0 ;
Theta_x = 0.0 ;
Theta_y = 0.0 ;
Theta_z = 0.0 ;
</Right>
</BOUNDARIES>

<LOADING> #efforts on the beam
<Left> #ponctual effort on the left
Uy = 0.0 ;
Theta_x = 0.0 ;
</Left>
<Right> #ponctual effort on the right
Theta_z = 0.0 ;
</Right>
<Lineic> #Lineic effort
#nope
</Lineic>
</LOADING>

</BEAM>


<LINK>

Beam1 = Poutre1 ; #1st beam
X = 3.0 ; #curvilign absciss on beam 1
Beam2 = Poutre2 ; #2nd beam
X = 0.0 ;
Type = pivot ;
Axis = x ;

</LINK>

<SOLVER>

Type = static;
Output_name = toto;

<OUTPUT>
Beam = Poutre1 ;
U = Ux Uy Uz Theta_x Theta_y Theta_z ;
F = Ux Uy Uz Theta_x Theta_y Theta_z ;
</OUTPUT>

<OUTPUT>
Beam = Poutre2 ;
U = Ux Uy Uz Theta_x Theta_y Theta_z ;
F = Ux Uy Uz Theta_x Theta_y Theta_z ;
</OUTPUT>

<TIME_OUTPUT>
Beam = Poutre1 ;
X = 3.0 ;
Component = Ux ;
</TIME_OUTPUT>

</SOLVER>"""


    try :
        donnees=traduire(texte1)
    except SyntaxError as erreur :
        print("-Erreur à[ligne,colonne]")
        print('    ',erreur)

    print(donnees)