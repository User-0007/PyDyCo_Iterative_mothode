# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:55:48 2019

@author: physique
"""
#==================================================Si vous rencontrez un problème alrs RTFM ===============================================
import os
import time
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
import matplotlib.pyplot as plt
if os.name == 'tn' :
  from foncdyco import *
elif os.name =='posix':  
  from foncdyco.foncdyco import *
import numpy as np
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *
from decimal import Decimal
import matplotlib.gridspec as gridspec
from foncdyco.foncdyco import *
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++information su les potentiels+++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mu1=0.01                                                     
T1=300	
mu2=0 
T2=300
#pou une augmentation graduelle 
m=0 #le debut de l'incrémentation du potentiel électrique
#Preciser la valeur du pas :
#la valeur du pas pour le circuit électrique:
pasv=0.1
#la valeur du pas pour le circuit thermique :
past=10
    
#Si vous ne souhaitez pas augmenter le potentiel graduelement, décommentez cette instruction:
m=mu1 #<= cette instruction 
	

if T1>T2:
    T=T2
    #Si vous  souhaitez augmenter la temeprature graduelement, et T1>T2 commentez ces deux instructions:
    T=T1  #<=cette instruction 
    Tmax =T1 #<=cette instruction
elif T1<T2:
    T= T1
    #Si vous souhaitez  augmenter la temeprature graduelement, et T1<T2 commentez ces deux instructions:
    T=T2  #<=cette instruction
    Tmax=T2 #<=cette instruction
elif T1==T2:
    T=T1
    Tmax=T
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++Précision sur l'écart quadratique de la température et du potentiel électrique +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
precisionV=10e-15
precisionT=10e-5


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++information sur les dipoles+++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#RESISTANCE THERMIQUE DE {P,S,M,N,K,V}
RTP=(1/0.2)
RTS=(1/1e-2)
RTM=(1/0.2)
RTN=(1/0.2)
RTK=(1/10)
RTV=(1/10)
#RESISTANCE ELECTRIQUE DE {P,S,M,N,K,V}
REP=1/10
RES=1/(1.18e-15)
REM=1/0.001
REN=1/10
REK=0.8e-3
REV=1e-3
#COEFICIENT SEEBECK  DE {P,S,M,N,K,V}
SP=200e-6
SS=0  
SM=0
SN=-200e-6
SK=200e-6  
SV=-200e-6


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Ce programme est divisé en 3 parties 
#############################################################Debut de programme ####################################################################################################################################
#######################################################----ACTE 1--##########################################################-#########################################################################################################################################################################################################

mat = open("matériaux.txt","r")


f= mat.read()
c=f.index('\n')
ll= f.count('\n') 
if os.name == 'nt': 
  ll=ll+1
u1=[]
v1=[]
seeb=[]
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Recolte des informations 
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
for tip in f :
    
   if tip == 'p':
    u1.append(REP)
    v1.append(RTP)
    seeb.append(SP)
   elif tip == 's':
    u1.append(RES)
    seeb.append(SS)
    v1.append(RTS)
   elif tip == 'm':
    u1.append(REM)
    seeb.append(SM)
    v1.append(RTM)
   elif tip == 'n':
    u1.append(REN)
    seeb.append(SN)
    v1.append(RTN)
   elif tip == 'v':
    u1.append(REV)
    seeb.append(SV)
    v1.append(RTV)
   elif tip == 'k':
    u1.append(REK)
    seeb.append(SK)
    v1.append(RTK)

    
seeb=np.array(seeb)
seeb=seeb.reshape((ll,c))

u1= np.array(u1)
u1=u1.reshape((ll,c))

v1= np.array(v1)
v1=v1.reshape((ll,c))         

mat.close()


s=ll+1
#-----------------------------------------------------------------------
circuit = Circuit("circuit électrique")
res = tran_mat(ll,c,u1)
fseeb =[]
fseeb = tran_mat(ll,c,seeb)
circuit ,tp = Circuit_creation_R(circuit,res,ll,c)
Circuit_creation_V(circuit,ll,c,m,mu2)
#-----------------------------------------------------------------------  
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
sn=[]
fn=[]

#-----------------------------------------------------------------------
for node in analysis.nodes.values():
#     print('\033[34mNode {}: {:4.4f} V\033[0m'.format(str(node), float(node)))
     sn.append(str(node))
     fn.append(float(node))
#-----------------------------------------------------------------------
Vmoy=m_moy(ll,c,sn,fn)
VR  =delt_m(ll,c,sn,fn)    
                               
#-----------------------------------------------------------------------
circuiT = Circuit("circuit thermique")
ret=tran_mat(ll,c,v1)
circuiT,tpT=Circuit_creation_R(circuiT,ret,ll,c)
if T1>T2 :
               Circuit_creation_V(circuiT,ll,c,T,T2)
        
elif T1<T2 :
               Circuit_creation_V(circuiT,ll,c,T1,T)
        
else :
       
              Circuit_creation_V(circuiT,ll,c,T1,T2)   


#-----------------------------------------------------------------------
simulator = circuiT.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
snt=[]
fnt=[]



for node in analysis.nodes.values():
#     print('\033[31mNode {}: {:4.4f} K\033[0m'.format(str(node), float(node)))
     snt.append(str(node))
     fnt.append(float(node)) 
     
Tmoy=m_moy(ll,c,snt,fnt)
TR =delt_m(ll,c,snt,fnt)

print('-------------iteration---------------')

XVR=[]
XTR=[]
VR1=np.zeros(len(VR))
TR1=np.zeros(len(TR))
alpha=np.array(fseeb)
fois = 1
res=np.array(res)
VR=np.array(VR)
TR=np.array(TR)
#TR1=TR
#VR1=VR


#+++++++++++++++++++++++++++++++++++++++++++++++++Iterations+++++++++++++++++++++++++++Iteration
while m<=mu1 or  T<=Tmax :
  #       time.sleep(1)
         print ('--------iteration',fois) 
         if (m==mu1 and T ==Tmax and beta (VR1,VR)<precisionV and  beta (TR,TR1)<precisionT ):
                break 
         print('Ecart quadratique du potentiel électrique dans le réseau avec l''itération précedente',beta (VR1,VR),'V^2')
         print('Ecart quadratique de la température dans le réseau avec l''itération précedente',beta (TR1,TR),'K^2')               
   #      time.sleep(3)
	
         if (beta (VR1,VR)<10e-2 and  beta (TR,TR1)<10e-6 ):
    
             if (m<=mu1 and mu1!=0) :
                m=m+pasv
                if m>mu1:
                    m=mu1
                
             if (T1>T2 or T2<T1) and T<=Tmax:   
                T=T+past
                if T>Tmax:
                    T=Tmax
    #     time.sleep(3)
         #print('m,T',m,T,circuit)

         if (m<=mu1 or T <=Tmax ): 
          #print(fois ,'mu ',m ,'T',T)  
          fois=fois+1
          XTR.append(beta(TR1,TR))
          XVR.append(beta(VR1,VR))
          VR1=np.array(VR)
          TR1=np.array(TR)
          Tmoy1=np.array(Tmoy)
          Vmoy1=np.array(Vmoy)
          ina=In(alpha,1/res,TR)
          iea=iq(alpha,1/res,Tmoy,Vmoy,VR,TR)
          #print(iea,'Ie')
          print(ina,'IN')
          #time.sleep(10)
          #----------------------------------------------------
          #----------------------------------------------------
          del circuit
          circuit = Circuit("circuit électrique")
          circuit,tp = Circuit_creation_R(circuit,res,ll,c)
          Circuit_creation_I(circuit,ina,ll,c)
          Circuit_creation_V(circuit,ll,c,m,mu2)
#          print(circuit)
          #----------------------------------------------------  
          simulator = circuit.simulator(temperature=25, nominal_temperature=25)
          analysis = simulator.operating_point()
          sn=[]
          fn=[]
          #-----------------------------------------------------------------------
          for node in analysis.nodes.values():
              #     print('\033[34mNode {}: {:4.4f} V\033[0m'.format(str(node), float(node)))
              sn.append(str(node))
              fn.append(float(node))
          #-------------------------------------------
          #for resistance in tp :#Je suis la syntaxe A..
           #    resistance.minus.add_current_probe(circuit) 

          simulator = circuit.simulator(temperature=25, nominal_temperature=25)
          analysis = simulator.operating_point()                                                     

          #I0=[]          
          #for branche in analysis.branches.values():
#              print('branche {}: {:5.19f} A'.format(str(branche), float(branche))) 
              #I0.append(float(branche))
              
          
          Vmoy=m_moy(ll,c,sn,fn)
          VR  =delt_m(ll,c,sn,fn)  
          #----------------------------------------------------
          #----------------------------------------------------
          #----------------------------------------------------
          #----------------------------------------------------
          
          del circuiT
          circuiT = Circuit("circuit Thermique ")
          circuiT, tpT= Circuit_creation_R(circuiT,ret,ll,c)
          Circuit_creation_I(circuiT,iea,ll,c)
          if T1>T2 :
               Circuit_creation_V(circuiT,ll,c,T,T2)
         
          elif T1<T2 :
               Circuit_creation_V(circuiT,ll,c,T1,T)
         
          else :
         
              Circuit_creation_V(circuiT,ll,c,T1,T2)
#        
          #----------------------------------------------------  
          simulator = circuiT.simulator(temperature=25, nominal_temperature=25)
          analysis = simulator.operating_point()
          snt=[]
          fnt=[]
          #-----------------------------------------------------------------------
          for node in analysis.nodes.values():
         #     print('\033[34mNode {}: {:4.4f} V\033[0m'.format(str(node), float(node)))
               snt.append(str(node))
               fnt.append(float(node))
          #------------------------------------------
          IE=[] 
          for branche in analysis.branches.values():
 #              print('branche {}: {:5.19f} A'.format(str(branche), float(branche))) 
               IE.append(float(branche))
              
          bitch=IE[1]
          Tmoy=m_moy(ll,c,snt,fnt)
          TR=delt_m(ll,c,snt,fnt)
         # print(Tmoy,TR,'Tmoy,TR')
	  

print('Fin des iteration')
print('Création des fichiers VTK ...')

#+++++++++++++++++++++++++++++++++++++++Transformation en chaque matrice ++++++++++++++++++++++++++++++++++++++
global MatV
MatV= np.zeros(((ll+1)//2,c))
MatT= np.zeros(((ll+1)//2,c))
MatV=Netwk(ll,c,sn ,fn )
MatT=Netwk(ll,c,snt,fnt)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++Vtk pour les courants électrqiues+++++++++++++++++++++++++++++++++++++++++++++++++++++
zz=len(res)
I0=[]
rr=1
for resistance in tp :#Je suis la syntaxe A..
             resistance.minus.add_current_probe(circuit) 
simulator = circuit.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
for branche in analysis.branches.values():
     #xprint('branche {}: {:5.19f} A'.format(str(branche), float(branche)))
     if os.name=='nt':
       if str(branche)=='vr'+str(rr)+'_minus': 
           I0.append(float(branche))
           rr=rr+1
     elif os.name=='posix':
       if str(branche)=='VR'+str(zz)+'_minus': 
           I0.append(float(branche))
           zz=zz-1
if os.name=='posix':
  I0.reverse()
 
I0=I0-ina	

IN=carr(ll,c,I0)
vtk_IN(IN)#circuit simple
vtk_IN_vc(IN,MatV)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++Vtk pour les courants d'énergie et le courant de chaleur +++++++++++++++++++++++++++++++++++++++++++++++++++++
for resistance in tpT :#Je suis la syntaxe A..
         resistance.minus.add_current_probe(circuiT)
zz=len(res)
IE=[]
rr=1
simulator = circuiT.simulator(temperature=25, nominal_temperature=25)
analysis = simulator.operating_point()
for branche in analysis.branches.values():
    #print('branche {}: {:5.19f} A'.format(str(branche), float(branche))) 
    if os.name=='nt':
     if str(branche)=='vr'+str(rr)+'_minus': 
           IE.append(float(branche))
           rr=rr+1
    elif os.name=='posix':
     if str(branche)=='VR'+str(zz)+'_minus': 
           IE.append(float(branche))
           zz=zz-1
if os.name=='posix':     
 IE.reverse()	
IE=IE-iea 
IE0=carr(ll,c,IE)
vtk_IE(IE0)
vtk_IE_vc(IE0,MatV)
IE=np.array(IE)#---------
I0=np.array(I0)#---------
IQ=np.zeros(len(res))
IQ=IE-Vmoy*I0
IQ0=carr(ll,c,IQ)
vtk_IQ(IQ0)
vtk_IQ_vc(IQ0,MatV)
#+++++++++++++++++++++++++++++++++++++++++++++++Courant d'entropie et sistribution d'entropie+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
Is=IQ/Tmoy
Is0=carr(ll,c,Is)
vtk_Is(Is0) #courant d'entropie
vtk_Is_vc(Is0,MatV)
ITR =Invdelt_m(ll,c,snt,fnt)
ITR=np.array(ITR)#---------
R=np.array(res)#---------
S=(R*I0**2)/Tmoy+IQ*ITR
S=carr(ll,c,S)
vtk_S(S)  #Entropie dans le système 



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++equation électrique +++++++++++++++++ 
#cette équation est une tentative de retrouver une equation electrique qui suit la meme forme que l'equation de la chaleur sous les conditions de dirichlet  
def mufonc(x):
    sig=1/res[2]
    k=1/ret[2]
    a=alpha[2]
    l=len(snt)-1
    j=(T1-T2)/(9/1.48) 
    c=mu1
    f=(((-mu1+mu2)*x)/(l))+((a**2)*(j**2)*(l*x))/(sig*k)-((a**2)*(j**2)*x**2)/(sig*k)+mu1
    #f=(((-mu1+mu2)*x)/(l))-(2*a*(T2-T1)/(l))*x+(2*a*(T2-T1)/(l**2))*x**2+mu2
    return f
#                                               Ne pas en tenir compte 


#++++++++++++++++++++++++++++++++++++++++++++++++++++Equation de la chaleur++++++++++++++++++++++++++++++++++++++++++++++++++


def  chalfonc(x): 																
    sig=1/res[2]
    k=1/ret[2]
    l=len(snt)-1
    j=-I0[1] 
    c=T1
    f=(((T2-T1)*x)/(l))+(((j**2)*(l*x-x**2))/(2*sig*k))-(((j**2)*(x**2))/(2*sig*k))+c  
    return f
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

x= range (fois-1)

#ce graphe permet de voir l'evolution de l'écart quadratique du potentiel éléctrique et  de la température au cours du temps
#Si ça vous interesse déchochez toutes ces instructions 
#figure= pyplot.figure ()
#pyplot.plot(x,XVR,'blue')
#pyplot.plot(x,XTR,'red')
#pyplot.grid()
#pyplot.show()







#++++++++++++++++++++++++++++++++++++++++++++++++Tracer l'equation de la chaleur +++++++++++++++++++++++++++++++++++++++++++++++
#ATENTION :ces equations sont utilisées pour un réseau en une dimension 
yan=[]
Tlist=[]
mulist=[]
xlist =range(1,MatT.shape[0]+1)
for i in xlist :
    yan.append(chalfonc(i-1))  #les valeur de  l'equation de la chaleur 
#    mulist.append(mufonc(i-1)) # tentative de retrouver l'equation electrique  ne pas en tenir compte

ylist=[]
#fig= pyplot.figure ()
for x in xlist :
    i=str(x)
    ylist.append(fn [sn.index(i)])
    Tlist.append(fnt [snt.index(i)])
#pyplot.plot( xlist, yan,label='pink')#les valeurs sous l'equation  de la lachaleur
#pyplot.plot( xlist, ylist,label='Vij')#potentiel  electrique en chaque noeuds pour valider l'effet seebeck 
#plt.plot(xlist,mulist,label='equation')
#plt.plot(xlist,Tlist) # les valeurs retrouvé sur pydyco
#plt.grid()
#plt.legend()
#plt.xlabel("Noeuds")
#plt.ylabel("Potentiel électrique")
#pyplot.show() #pour valider l'effet seebeck 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##################
vtk_volt(MatV) #Fichier Vtk pour le potentiel électrique  
vtk_temp(MatT) #Fichier Vtk pour la temperature 


