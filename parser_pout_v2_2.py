#This code provides the parser for simple beams

import sp

def parser():
    commentaire=sp.R(r'#.*')
    blancs=sp.R(r'\s+')
    nom=sp.R(r'[a-zA-Z_]')
    nomPout=sp.R(r'\w+')
    nomFich=sp.R(r'\[\w+\]')
    component=sp.R(r'Ux')|sp.R(r'Uy')|sp.R(r'Uz')|sp.R(r'Theta_x')|sp.R(r'Theta_y')|sp.R(r'Theta_z')
    coord=sp.R(r'X')|sp.R(r'Y')|sp.R(r'Z')
    nombre=sp.R(r'[0-9]+') / int  #convertit le nombre en entier
    nbVirgule=(sp.R(r'-?'+'[0-9]*'+'.'+'[0-9]+')|sp.R(r'-?'+'[0-9]+')|sp.R(r'-?'+'[0-9]+'+'.')) /float
    express=sp.R(r'[\w|\+|\-|\*|\^|(|)|\.|/]+')
    separator=sp.R(r';|,|_|-')
    
    baliseDebutPoutre=sp.K(r'<')
    baliseDebutBoundaries=sp.K(r'<BOUNDARIES>')|sp.K(r'<LIMITES>')
    baliseDebutDimensions=sp.K(r'<DIMENSIONS>')
    baliseDebutLeft=sp.K(r'<Left>')|sp.K(r'<Gauche>')
    baliseDebutRight=sp.K(r'<Right>')|sp.K(r'<Droite>')
    baliseDebutLoading=sp.K(r'<LOADING>')|sp.K(r'<CHARGEMENT>')
    baliseDebutLin=sp.K(r'<Lineic>')|sp.K(r'<Lineique>')
    baliseDebutSolveur=sp.K(r'<SOLVER>')|sp.K(r'<SOLVEUR>')
    baliseDebutLink=sp.K(r'<LINK>')|sp.K(r'<LIAISON>')
    baliseDebutOutput=sp.K(r'<OUTPUT>')|sp.K(r'<SORTIE>')
    baliseDebutOutDyn=sp.K(r'<TIME_OUTPUT>')|sp.K(r'<SORTIE_TEMPS>')
    baliseDebutOutPlot=sp.K(r'<PLOT_OUTPUT>')|sp.K(r'<SORTIE_AFFI>')
    baliseDebutParam=sp.K(r'<PARAMETERS>')|sp.K(r'<PARAMETRES>')
    baliseDebutID=sp.K(r'<Init_defo>')|sp.K(r'<Defo_init>')
    baliseDebutRep=sp.K(r'<SYSTEM>')|sp.K(r'<REPERE>')
    baliseDebutRepM=sp.K(r'<MOB_SYSTEM>')|sp.K(r'<REPERE_MOB>')
    
    baliseFinPoutre=sp.K(r'</BEAM>')|sp.K(r'</POUTRE>')
    baliseFinBoundaries=sp.K(r'</BOUNDARIES>')|sp.K(r'</LIMITES>')
    baliseFinDimensions=sp.K(r'</DIMENSIONS>')
    baliseFinLeft=sp.K(r'</Left>')|sp.K(r'</Gauche>')
    baliseFinRight=sp.K(r'</Right>')|sp.K(r'</Droite>')
    baliseFinLoading=sp.K(r'</LOADING>')|sp.K(r'</CHARGEMENT>')
    baliseFinLin=sp.K(r'</Lineic>')|sp.K(r'</Lineique>')
    baliseFinSolveur=sp.K(r'</SOLVER>')|sp.K(r'</SOLVEUR>')
    baliseFinLink=sp.K(r'</LINK>')|sp.K(r'</LIAISON>')
    baliseFinOutput=sp.K(r'</OUTPUT>')|sp.K(r'</SORTIE>')
    baliseFinOutDyn=sp.K(r'</TIME_OUTPUT>')|sp.K(r'</SORTIE_TEMPS>')
    baliseFinOutPlot=sp.K(r'</PLOT_OUTPUT>')|sp.K(r'</SORTIE_AFFI>')
    baliseFinParam=sp.K(r'</PARAMETERS>')|sp.K(r'</PARAMETRES>')
    baliseFinID=sp.K(r'</Init_defo>')|sp.K(r'</Defo_init>')
    baliseFinRep=sp.K(r'</SYSTEM>')|sp.K(r'</REPERE>')
    baliseFinRepM=sp.K(r'</MOB_SYSTEM>')|sp.K(r'</REPERE_MOB>')
    
    gNomPoutre=sp.K(r'Beam Name')|sp.K(r'Nom Poutre')
    gGauche=sp.K(r'Left')|sp.K(r'Gauche')
    gDroite=sp.K(r'Right')|sp.K(r'Droite')
    gPout1=sp.K(r'Beam1')|sp.K(r'Poutre1')
    gPout2=sp.K(r'Beam2')|sp.K(r'Poutre2')
    gPout=sp.K(r'Beam')|sp.K(r'Poutre')
    gNomSort=sp.K(r'Output_name')|sp.K(r'Nom_sortie')
    gCompo=sp.K(r'Component')|sp.K(r'Composante')
    gAx=sp.K(r'Axis')|sp.K(r'Axe')
    gRepere=sp.K(r'System')|sp.K(r'Repere')
    gNom=sp.K(r'Name')|sp.K(r'Nom')
    
    egal=sp.K(r'=')
    depoin=sp.K(r':')
    
#########################################################################################################
    with sp.Separator(blancs|commentaire):  #on fait abstraction des blancs et des commentaires
        assembly=sp.Rule()
        
        beam=sp.Rule()
        link=sp.Rule()
        repe=sp.Rule()
        repemob=sp.Rule()
        
        blocBoundaries=sp.Rule()
        blocDimensions=sp.Rule()
        blocLoading=sp.Rule()
        blocSolveur=sp.Rule()
        blocOutput=sp.Rule()
        blocOutDyn=sp.Rule()
        blocOutPlot=sp.Rule()
        blocParam=sp.Rule()

        blocLeft=sp.Rule()
        blocRight=sp.Rule()
        blocLineic=sp.Rule()
        blocInDef=sp.Rule()   

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
        blocCo=sp.Rule()
        
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
        blocSeparator=sp.Rule()
        
        blocTech=sp.Rule()
        blocDur=sp.Rule()
        blocComp=sp.Rule()
        blocNiter=sp.Rule()
        blocNom=sp.Rule()
        
        blocAxis=sp.Rule()
        blocRepere=sp.Rule()
        affect=sp.Rule()
        blocRef=sp.Rule()
        
        L=sp.Rule()
        E=sp.Rule()
        Type=sp.Rule()
        N=sp.Rule()
    

##########################################################################################################
        assembly|=repe[0:] & beam[0:] & repemob[0:] & link[0:] & blocSolveur & blocParam[0:1]

        link|=baliseDebutLink & blocName1 & blocx & blocName2 & blocx & blocSolve & blocAxis[0:1] & blocRepere[0:1] & baliseFinLink
        beam|=baliseDebutPoutre & gNomPoutre & egal & '"' & nomPout &'" >' & blocDimensions & blocBoundaries & blocLoading & baliseFinPoutre
        repe|=baliseDebutRep & blocNom  & blocSolve & blocMz & blocMy & blocMx & baliseFinRep        
        repemob|=baliseDebutRepM & blocNom & blocRef & blocSolve & blocMz & blocMy & blocMx & baliseFinRepM
        
        blocBoundaries|=baliseDebutBoundaries & blocLeft & blocRight & baliseFinBoundaries
        blocDimensions|=baliseDebutDimensions & blocCoL & blocCoR & blocE & blocNu & blocS & blocIy & blocIz & blocRho & blocNb & blocInDef[0:1] & baliseFinDimensions
        blocLoading|=baliseDebutLoading & blocLeft[0:1] & blocRight[0:1] & blocLineic[0:1] & baliseFinLoading
        blocSolveur|=baliseDebutSolveur & blocSolve & blocDur[0:1] & blocTech[0:1] & blocNiter[0:1] & blocName & blocOutput[0:] & blocOutDyn[0:] & blocOutPlot[0:] & baliseFinSolveur
        blocOutput|=baliseDebutOutput & blocBNam & blocU & blocF & baliseFinOutput
        blocOutDyn|=baliseDebutOutDyn & blocBNam & blocx & blocComp & baliseFinOutDyn
        blocOutPlot|=baliseDebutOutPlot & blocBNam & blocCo & blocU & blocF & blocSeparator[0:1] & baliseFinOutPlot
        blocParam|=baliseDebutParam & affect[0:] & baliseFinParam

        blocLeft|=baliseDebutLeft & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinLeft
        blocRight|=baliseDebutRight & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & blocTx[0:1] & blocTy[0:1] & blocTz[0:1] & baliseFinRight
        blocLineic|=baliseDebutLin & blocFx[0:1] & blocFy[0:1] & blocFz[0:1] & blocMx[0:1] & blocMy[0:1] & blocMz[0:1] & baliseFinLin
        blocInDef|=baliseDebutID & blocUx[0:1] & blocUy[0:1] & blocUz[0:1] & baliseFinID 

        blocCoL|=gGauche & depoin & L & ',' & L & ',' & L & ';'
        blocCoR|=gDroite & depoin & L & ',' & L & ',' & L & ';'
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
        blocCo|='Coord' & egal & coord[0:3] & ';'
        
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
        blocName1|=gPout1 & egal & Type & ';'
        blocName2|=gPout2 & egal & Type & ';'
        blocName|=gNomSort & egal & Type & ';'
        blocx|='X' & egal & L & ';'
        blocBNam|=gPout & egal & nomPout & ';'
        blocSeparator|='Sep' & egal & separator & ';'
        
        blocTech|='dt' & egal & L & ';'
        blocDur|='T' & egal & L & ';'
        blocComp|=gCompo & egal & component & ';'
        blocNiter|='n_iter' & egal & N & ';'
        blocNom|=gNom & egal & nomPout & ';'
        
        blocAxis|=gAx & egal & nom & ';'
        blocRepere|=gRepere & egal & nomPout & ';'
        affect|=nomPout & egal & express & ';'
        blocRef|='ref' & egal & nomPout & ',' & express & ';'
        
        L|=express
        E|=express|nomFich#|express#|nomPout|nbVirgule
        Type|=nomPout
        N|=express
        
        
    return assembly

if __name__ == "__main__" :

    traduire=parser()

    texte1="""
<SYSTEM>
Name = Repere0 ;
Type = ZYX ;
mz = Pi1/16;
my = Pi1/8 ;
mx = -Pi1/16 ;
</SYSTEM>

<Beam Name="Poutre1" >

<DIMENSIONS>
Left : 0.0 , 0.0 , 0.0 ;#Coords of the beginning of the beam (m)
Right: 0.3 , 0.0 , 0.0 ;#Coords of the end of the beam (m)
E= 210000000000.0 ;#Young modulus (Pa)
Nu = 0.3 ;#Fish coefficient
S = 0.001 ;#Section(m^2)
Iy = 0.00000000001 ;#Quadratic moment (m^4)
Iz =[dataiz] ;#Quadratic moment (m^4)
Rho = 7800 ; #Volumic mass (kg.m-3)
n = 10 ; #Number of elements between 2 sides of a beam

<Init_defo>
Ux = 0.01*sin(X/0.3*pi) ;
</Init_defo>

</DIMENSIONS>

<LIMITES> #boundaries conditions for displacement
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
<Lineic> #Lineic effort
fy=-1.0;#Nm-1
mx = [datamx];
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

<MOB_SYSTEM>
Name = Repere0 ;
ref = Poutre2 , 0 ;
Type = ZYX ;
mz = Pi1/16;
my = Pi1/8 ;
mx = -Pi1/16 ;
</MOB_SYSTEM> 

<LINK>

Beam1 = Poutre1 ; #1st beam
X = 3.0 ; #curvilign absciss on beam 1
Beam2 = Poutre2 ; #2nd beam
X = 0.0 ;
Type = pivot ;
Axis = x ;
System = global ;

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

<PLOT_OUTPUT>
Beam = Poutre1 ;
Coord = X Y Z ;
U = Ux Uy Uz Theta_x Theta_y Theta_z ;
F = Ux Uy Uz Theta_x Theta_y Theta_z ;
Sep = ; ;
</PLOT_OUTPUT>

</SOLVER>

<PARAMETERS>

Pi1 = 3.1415926535 ;
L1 = cos(5.03) ;  #def of a param

</PARAMETERS>"""


    try :
        donnees=traduire(texte1)
    except SyntaxError as erreur :
        print("-Erreur à[ligne:colonne]")
        print('    ',erreur)

    print(donnees)