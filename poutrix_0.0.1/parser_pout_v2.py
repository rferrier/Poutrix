#This code provides the parser for simple beams

import sp

def parser():
    commentaire=sp.R(r'#.*')
    blancs=sp.R(r'\s+')
    nom=sp.R(r'[a-zA-Z]')
    nomPout=sp.R(r'\w+')
    nombre=sp.R(r'[0-9]+') / int  #convertit le nombre en entier
    nbVirgule=sp.R(r'[0-9]+'+'.'+'[0-9]+') /float
    baliseDebutPoutre=sp.K(r'<')
    baliseDebutBoundaries=sp.K(r'<BOUNDARIES>')  
    baliseDebutDimensions=sp.K(r'<DIMENSIONS>')
    baliseDebutLeft=sp.K(r'<Left>')
    baliseDebutRight=sp.K(r'<Right>')
    baliseDebutLoading=sp.K(r'<LOADING>')
    baliseDebutLin=sp.K(r'<Lineic>')
    baliseDebutSolveur=sp.K(r'<SOLVER>')
    baliseDebutLink=sp.K(r'<LINK>')
    
    baliseFinPoutre=sp.K(r'</BEAM>')
    baliseFinBoundaries=sp.K(r'</BOUNDARIES>')
    baliseFinDimensions=sp.K(r'</DIMENSIONS>')
    baliseFinLeft=sp.K(r'</Left>')
    baliseFinRight=sp.K(r'</Right>')
    baliseFinLoading=sp.K(r'</LOADING>')
    baliseFinLin=sp.K(r'</Lineic>')
    baliseFinSolveur=sp.K(r'</SOLVER>')
    baliseFinLink=sp.K(r'</LINK>')
    
#########################################################################################################
    with sp.Separator(blancs|commentaire):  #on fait abstraction des blancs et des commentaires
        assembly=sp.Rule()
        
        beam=sp.Rule()
        link=sp.Rule()
        blocBoundaries=sp.Rule()
        blocDimensions=sp.Rule()
        blocLoading=sp.Rule()
        blocSolveur=sp.Rule()

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
        
        L=sp.Rule()
        E=sp.Rule()
        Type=sp.Rule()
        N=sp.Rule()

##########################################################################################################
        assembly|=beam[0:] & link[0:] & blocSolveur 

        link|=baliseDebutLink & blocName1 & blocx & blocName2 & blocx & baliseFinLink
        beam|=baliseDebutPoutre & 'Beam Name="' & nomPout &'" >' & blocDimensions & blocBoundaries & blocLoading & baliseFinPoutre
        
        blocBoundaries|=baliseDebutBoundaries & blocLeft & blocRight & baliseFinBoundaries
        blocDimensions|=baliseDebutDimensions & blocCoL & blocCoR & blocE & blocNu & blocS & blocIy & blocIz & blocRho & blocNb & baliseFinDimensions
        blocLoading|=baliseDebutLoading & blocLeft & blocRight & blocLineic & baliseFinLoading
        blocSolveur|=baliseDebutSolveur & blocSolve & blocName & baliseFinSolveur

        blocLeft|=baliseDebutLeft & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinLeft
        blocRight|=baliseDebutRight & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinRight
        blocLineic|=baliseDebutLin & blocFx[0:1] & blocFy[0:1] & blocFz[0:1] & blocMx[0:1] & blocMy[0:1] & blocMz[0:1] & baliseFinLin

        blocCoL|='Left :' & L & ',' & L & ',' & L & ';'
        blocCoR|='Right :' & L & ',' & L & ',' & L & ';'
        blocL|='L =' & L &';'
        blocE|='E =' & E &';'
        blocNu|='Nu =' & E &';'
        blocS|='S =' & E &';'
        blocIy|='Iy =' & E &';'
        blocIz|='Iz =' & E &';'
        blocRho|='Rho =' & E &';'
        blocNb|='n ='& N &';'
        
        blocUx|='Ux =' & L & ';'
        blocUy|='Uy =' & L & ';'
        blocUz|='Uz =' & L & ';'
        blocTx|='Theta_x =' & L & ';'
        blocTy|='Theta_y =' & L & ';'
        blocTz|='Theta_z =' & L & ';'

        blocFx|='fx =' & E & ';'
        blocFy|='fy =' & E & ';'
        blocFz|='fz =' & E & ';'
        blocMx|='mx =' & E & ';'
        blocMy|='my =' & E & ';'
        blocMz|='mz =' & E & ';'

        blocSolve|='Type =' & Type & ';'
        blocName1|='Beam1 =' & Type & ';'
        blocName2|='Beam2 =' & Type & ';'
        blocName|='Output_name =' & Type & ';'
        blocx|='X =' & L & ';'
        L|=nbVirgule
        E|=nbVirgule|nomPout
        Type|=nomPout
        N|=nombre

    return assembly
"""
traduire=parser()

texte1=
<Beam Name="Poutre1" >

<DIMENSIONS>
Left : 0.0 , 0.0 , 0.0 ;#Coords of the beginning of the beam (m)
Right : 0.3 , 0.0 , 0.0 ;#Coords of the end of the beam (m)
E = 210000000000.0 ;#Young modulus (Pa)
Nu = 0.3 ;#Fish coefficient
S = 0.001 ;#Section(m^2)
Iy = 0.00000000001 ;#Quadratic moment (m^4)
Iz = dataiz ;#Quadratic moment (m^4)
Rho = 7800 ; #Volumic mass (kg.m-3)
n = 10 ; #Number of elements between 2 sides of a beam
</DIMENSIONS>

<BOUNDARIES> #boundaries conditions for displacement
<Left>
Ux = 0.0 ;
Uy = 0.0 ;
Uz = 0.0 ;
Theta_x = 0.0 ;
Theta_y = 0.0 ;
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
fy = 1.0;#Nm-1
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
Ux = 0.0 ;
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

</LINK>

<SOLVER>

Type = static;

</SOLVER>


try :
    donnees=traduire(texte1)
except SyntaxError as erreur :
    print("-Erreur à[ligne,colonne]")
    print('    ',erreur)

print(donnees)"""