# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:40:52 2019

@author: phys ue
"""#==================================================Si vous rencontrez un problème alors RTFM ===============================================
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from matplotlib import pyplot 
import numpy as np
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *
import matplotlib.gridspec as gridspec


#Fonction circuit ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    

def Circuit_creation_R(cir,dip,lin,col):
      tp=tuple() 
      h=0
      q=0
      k=1
      s=lin+1
      for l in range (0,s//2 ) :
              
        if  l%2 !=1 :  
          if q < (lin) :  
            
            for j   in  (range (1+l*col,col+l*col)) :   
                
                    re=(cir.R(k,j,j+1,dip[h]),)
                    k=k+1
                    
                    tp=tp+re
                    h=h+1  
            q=q+1
            
          if q < (lin) : 
              
            for v, j in zip ( range (col+l*col,l*col,-1) , range (col+1+col*l, 2*col+col*l+1 )):   
              
                 re=(cir.R(k,v,j,dip[h]),)  #je nomme mes resistance verticales 
                 k=k+1
                 tp=tp+re
                 h=h+1
             
            q=q+1    
              
           
            
        else:
              
            
          if q < (lin) :
              
            for j  in  (range (col+l*col,1+l*col,-1)) :   #si c'est paire je fait decrementation une  et je nomme mes resistances horizontamles  
              
                re=(cir.R(k,j,j-1,dip[h]),)
                k=k+1
                tp=tp+re
                h=h+1
                
            q=q+1
          if q < (lin) :
              
            for v, j in zip ( range (1+l*col,col+l*col+1) , range ( 2*col+col*l,col+col*l, -1 )):   
            
                re=(cir.R (k,v,j,dip[h]),) 
                k=k+1
                tp=tp+re
                h=h+1
               
            q=q+1
    
      return cir, tp
  
#circuit courent++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

def Circuit_creation_V(cir,lin,col,va,vb):
          s=lin+1
          for kk in range (1,col+1):
             cir.V(kk,kk,0 , va)
          for ss in range ((col)*(s//2)-col+1,col*(s//2)+1):  
             cir.V(ss,ss, 0,vb )
          return cir   
        
       
def Circuit_creation_I(cir,dip,lin,col):
      h=0
      q=0
      k=1
      s=lin+1
      for l in range (0,s//2 ) :
              
        if  l%2 !=1 :  
          if q < (lin) :  
            
            for j   in  (range (1+l*col,col+l*col)) :   
                
                    re=(cir.I(k,j,j+1,dip[h]),)
                    k=k+1
                    h=h+1  
            q=q+1
            
          if q < (lin) : 
              
            for v, j in zip ( range (col+l*col,l*col,-1) , range (col+1+col*l, 2*col+col*l+1 )):   
              
                 re=(cir.I(k,v,j,dip[h]),)  #je nomme mes resistance verticales 
                 k=k+1
                 h=h+1
             
            q=q+1    
              
           
            
        else:
              
            
          if q < (lin) :
              
            for j  in  (range (col+l*col,1+l*col,-1)) :   #si c'est paire je fait decrementation une  et je nomme mes resistances horizontamles  
              
                re=(cir.I(k,j,j-1,dip[h]),)
                k=k+1
                h=h+1
                
            q=q+1
          if q < (lin) :
              
            for v, j in zip ( range (1+l*col,col+l*col+1) , range ( 2*col+col*l,col+col*l, -1 )):   
            
                re=(cir.I (k,v,j,dip[h]),) 
                k=k+1
                h=h+1
               
            q=q+1
    
      return cir
  
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
#Fonction resi collone ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def tran_mat(lin,col,r ):
 mat=[] 
 s=lin+1
 q=0
 for l in range (0,s//2 ) :         
    if  l%2 !=1 :
      if q < (lin) :      
        for  i  in range (0,col-1) :        
                mat.append(r[q,i]) 
        q=q+1        
      if q < (lin) :      
        for i in range (col-1,-1,-1):   
             mat.append(r[q,i])
        q=q+1    
    else:            
      if q < (lin) :    
        for i  in range (0,col-1) :   
            mat.append(r[q,i])        
        q=q+1
      if q < (lin) :         
        for i in range (col-1,-1,-1):         
            mat.append(r[q,i])        
        q=q+1
 return mat 



#les moyennes ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def m_moy(lin,col,vn1,vn2):
    q=0 
    mmoy=[]  
    s=lin+1
    for l in range (0,s//2 ) : #♠je commence a faire ma boucle sur tout la toutes les lignes de ma matrice en prenanat une ligne sur deux  
        
      if  l%2 !=1 : # la je verfie si c'est paire ou imapaire 
          if q < (lin) : #on rajoute une condition pour ne pas depasser les valeurs de U  
            
            for j  in  (range (1+l*col,col+l*col)) :   #si c'est impaire je fait une incrementation et je nomme mes resistances horizontales     
            
                m1=str(j)
                m2=str(j+1)
                #print(m1,m2,'m1,m2 ')
                mmoy.append((vn2[vn1.index(m2)]+vn2[vn1.index(m1)])/2)
    #            Vmoy.append((fn[sn.index(m2)]*fn[sn.index(m1)])**0.5)
            q=q+1 
            
          if q < (lin) : 
              
            for v, j in zip ( range (col+l*col,l*col,-1) , range (col+1+col*l, 2*col+col*l+1 )):
                m1=str(v)
                m2=str(j)
                #print(m1,m2,'m1,m2')
                mmoy.append((vn2[vn1.index(m2)]+vn2[vn1.index(m1)])/2)
    #            Vmoy.append((fn[sn.index(m2)]*fn[sn.index(m1)])**0.5)
            q=q+1  
          
      else:
              if q < (lin): 
                  for j  in  (range (col+l*col,1+l*col,-1)) :
                    m1=str(j)
                    m2=str(j-1)
                    #print(m1,m2,q,'m1,m2')
                    mmoy.append((vn2[vn1.index(m2)]+vn2[vn1.index(m1)])/2)
    #                Vmoy.append((fn[sn.index(m2)]*fn[sn.index(m1)])**0.5)
                  q=q+1 
              
              if q < (lin) :
                   for v, j in zip(range (1+l*col,col+l*col+1) , range ( 2*col+col*l,col+col*l, -1 ) ):
                       m1=str(v)
                       m2=str(j)
                       #print(m1,m2,q,'m1,m2')
                       mmoy.append((vn2[vn1.index(m2)]+vn2[vn1.index(m1)])/2)
    #                   Vmoy.append((fn[sn.index(m2)]*fn[sn.index(m1)])**0.5)
                   q=q+1        
    return mmoy     

#Les differances ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def delt_m(lin,col,vn1,vn2):
    q=0 
    delm=[]  
    s=lin+1
    for l in range (0,s//2 ) : #♠je commence a faire ma boucle sur tout la toutes les lignes de ma matrice en prenanat une ligne sur deux  
        
      if  l%2 !=1 : # la je verfie si c'est paire ou imapaire 
          if q < (lin) : #on rajoute une condition pour ne pas depasser les valeurs de U  
            
            for j  in  (range (1+l*col,col+l*col)) :   #si c'est impaire je fait une incrementation et je nomme mes resistances horizontales     
            
                m1=str(j)
                m2=str(j+1)
                delm.append(vn2[vn1.index(m1)]-vn2[vn1.index(m2)])
            q=q+1 
            
          if q < (lin) : 
              
            for v, j in zip ( range (col+l*col,l*col,-1) , range (col+1+col*l, 2*col+col*l+1 )):
                m1=str(v)
                m2=str(j)
                delm.append(vn2[vn1.index(m1)]-vn2[vn1.index(m2)])
            q=q+1  
          
      else:
              if q < (lin): 
                  for j  in  (range (col+l*col,1+l*col,-1)) :
                    m1=str(j)
                    m2=str(j-1)
                    delm.append(vn2[vn1.index(m1)]-vn2[vn1.index(m2)])
                  q=q+1 
              
              if q < (lin) :
                   for v, j in zip(range (1+l*col,col+l*col+1) , range ( 2*col+col*l,col+col*l, -1 ) ):
                       m1=str(v)
                       m2=str(j)
                       delm.append(vn2[vn1.index(m1)]-vn2[vn1.index(m2)])
                   q=q+1        
    return delm     
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Invdelt_m(lin,col,vn1,vn2):
    q=0 
    delm=[]  
    s=lin+1
    for l in range (0,s//2 ) : #♠je commence a faire ma boucle sur tout la toutes les lignes de ma matrice en prenanat une ligne sur deux  
        
      if  l%2 !=1 : # la je verfie si c'est paire ou imapaire 
          if q < (lin) : #on rajoute une condition pour ne pas depasser les valeurs de U  
            
            for j  in  (range (1+l*col,col+l*col)) :   #si c'est impaire je fait une incrementation et je nomme mes resistances horizontales     
            
                m1=str(j)
                m2=str(j+1)
                delm.append(1/vn2[vn1.index(m1)]-1/vn2[vn1.index(m2)])
            q=q+1 
            
          if q < (lin) : 
              
            for v, j in zip ( range (col+l*col,l*col,-1) , range (col+1+col*l, 2*col+col*l+1 )):
                m1=str(v)
                m2=str(j)
                delm.append(1/vn2[vn1.index(m1)]-1/vn2[vn1.index(m2)])
            q=q+1  
          
      else:
              if q < (lin): 
                  for j  in  (range (col+l*col,1+l*col,-1)) :
                    m1=str(j)
                    m2=str(j-1)
                    delm.append(1/vn2[vn1.index(m1)]-1/vn2[vn1.index(m2)])
                  q=q+1 
              
              if q < (lin) :
                   for v, j in zip(range (1+l*col,col+l*col+1) , range ( 2*col+col*l,col+col*l, -1 ) ):
                       m1=str(v)
                       m2=str(j)
                       delm.append(1/vn2[vn1.index(m1)]-1/vn2[vn1.index(m2)])
                   q=q+1        
    return delm     






























def iq (alpha,sigma,tm,vm,delv,delt):
    i=(alpha*sigma*tm+vm*sigma)*delv+(alpha*vm*sigma+(alpha**2)*tm*sigma)*delt
    return i  

def In (alpha, sigma, delt) :
    i= alpha*sigma*delt 
    return i   
def beta (a,b): 
    f=sum ((b-a)**2)
    return f
#==================================================heat equation++++++++++++++++++++++++++++++++++     
def chalfonc (x):
    l=len(snt)-1
    j=-I0[1] 
    sig=1/ro[2]
    k=kp1[2]
#    c=T1+((T1-T2)/(l))+(-l*(j**2)+(j**2))/(2*sig*k)
    c=T1
    f=(((T2-T1)*x)/(l))+(((j**2)*l*x)/(2*sig*k))-(((j**2)*(x**2))/(2*sig*k))+c
    return f

#==================================================Fichier VTK scalars++++++++++++++++++++++++++++++++++ 
def Netwk(lin,col,si,fi):
    s=lin+1
    Mat=np.zeros(((s//2),col))
    it =0
    pom=1
    s= lin+1
    while (it < (s//2)) :
        for i in range (col): 
            j=str(pom)
            Mat[it,i]=fi[si.index(j)]  
            pom=pom+1
        it=it+1
 
        if  (it == (s//2)):
          break 
 
        for i in range (col-1,-1,-1):
            j=str(pom)
            Mat[it,i]=fi[si.index(j)]  
            pom=pom+1
        it=it+1
    return Mat 

def vtk_temp(mat):
  x=mat.shape[1]
  y=mat.shape[0]
  z=1
  c = open("temp.vtk", "w")	
  c.write('# vtk DataFile Version 3.0\nCube example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+'\nPOINTS  '+  str(x*y*z)+'  float')
  c.close()
  c = open("temp.vtk", "a")
  for k in range (z) :
     for j in range (y) :
        for i in range (x): 
             c.write('\n\t'+str(float(i))+'\t'+str(float(-j*x/y))+'\t'+str(float(k)))	
  c.write('\nPOINT_DATA \t' +str(x*y*z))
  c.write(' \n \nSCALARS Temp.(K) float\nLOOKUP_TABLE default')
  for i in range (y): 
      for j in range (x):
          c.write('\n\t'+str(mat[i,j]))
  c.close()
  return 

def vtk_volt(mat):
  x=mat.shape[1]
  y=mat.shape[0]
  z=1
  f = open("volt.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nCube example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("volt.vtk", "a")
  for k in range (z) :
     for j in range (y):
        for i in range (x): 
             f.write('\n\t'+str(float(i))+'\t'+str(float(-j*x/y))+'\t'+str(float(k)))	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nSCALARS Pot.(V) float\nLOOKUP_TABLE default')
  for i in range (y): 
        for j in range (x):
          f.write('\n\t'+str(mat[i,j]))
  f.close()
  return
 
#==================================================Fichier VTK Cuurent++++++++++++++++++++++++++++++++++ 
def carr(lin,col,vec):
    k=0
    j=0
    vm=np.zeros((lin,col))
    while j < lin :
     for i in range (col-1):
       vm[j,i]=vec[k]
       k=k+1
     j=j+1
     if j==lin:
       break 
     for i in range (col-1,-1,-1):
       vm[j,i]=vec[k]
       k=k+1
     j=j+1
    return vm

def vtk_IN(I):
  x=I.shape[1]
  y=I.shape[0]
  z=1
  f = open("IN.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IN.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i+0.5))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IN float\n')
  i=0
  while i< (y): 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
        i=i+1
        if i==y:
          break 
        for j in range (x):
          f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        i=i+1
  f.close()
  return 

def vtk_IE(I):
  x=I.shape[1]
  y=I.shape[0]
  z=1
  f = open("IE.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IE.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	

  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IE float\n')
  i=0
  while i< (y): 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
        i=i+1
        if i==y:
          break 
        for j in range (x):
          f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        i=i+1
  f.close()
  return 

def vtk_IQ(I):
  x=I.shape[1]
  y=I.shape[0]
  z=1
  f = open("IQ.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IQ.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i+0.5))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k)))
        j=j+1                  	

  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IQ float\n')
  i=0
  while i< (y): 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
        i=i+1
        if i==y:
          break 
        for j in range (x):
          f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        i=i+1
  f.close()
  return 



def vtk_Is(I):
  x=I.shape[1]
  y=I.shape[0]
  z=1
  f = open("Is.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("Is.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i+0.5))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS Is float\n')
  i=0
  while i< (y): 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
        i=i+1
        if i==y:
          break 
        for j in range (x):
          f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        i=i+1
  f.close()
  return 

def vtk_S(I):
  x=I.shape[1]
  y=I.shape[0]
  z=1
  f = open("S.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("S.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i+0.5))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(x-1)/y))+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nSCALARS S float\nLOOKUP_TABLE default')
  i=0
  while i< (y): 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t'+'\t\t')
        i=i+1
        if i==y:
          break 
        for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t'+'\t\t')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        i=i+1
  f.close()
  return 
#==================================================les nouveau courants changé==========================================
def vtk_IN_vc(I,Mat):
  x=Mat.shape[1]
  y=Mat.shape[0]
  lol=I.shape[0]
  z=1
  f = open("IN_vc.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IN_vc.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IN float\n')
  i=0 
  for j in range (x):
          f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
  i=i+1
  while i< (lol):
      for j in range (x):
          f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
      i=i+2
      if i==lol:
         break 
        #for j in range (x):
         # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        #i=i+1
  f.close()
  return 

#############################################################################################################################################"
###########################################pour le courant de chaleur #################################################################
def vtk_IQ_vc(I,Mat):
  x=Mat.shape[1]
  y=Mat.shape[0]
  lol=I.shape[0]
  z=1
  f = open("IQ_vc.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IQ_vc.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IQ float\n')
  i=0 
  for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
  i=i+1
  while i< (lol):
      for j in range (x):
          f.write('\n\t\t'+str(I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
      i=i+2
      if i==lol:
         break 
        #for j in range (x):
         # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        #i=i+1
  f.close()
  return 

###########################################################################################################################3######"
############"###################################Pour le  ################################################################################
def vtk_IE_vc(I,Mat):
  x=Mat.shape[1]
  y=Mat.shape[0]
  lol=I.shape[0]
  z=1
  f = open("IE_vc.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("IE_vc.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS IE float\n')
  i=0 
  for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
  i=i+1
  while i< (lol):
      for j in range (x):
          f.write('\n\t\t'+str(I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
      i=i+2
      if i==lol:
         break 
        #for j in range (x):
         # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        #i=i+1
  f.close()
  return 


#############################################################################################################################################
#################################################POur IS #################################################################################"
def vtk_Is_vc(I,Mat):
  x=Mat.shape[1]
  y=Mat.shape[0]
  lol=I.shape[0]
  z=1
  f = open("Is_vc.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("Is_vc.vtk", "a")
  j=0.0
  for k in range (z) :
     while j<(y):
        for i in range (x): 
          f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
        j=j+1
        if j==y:
           break
        for i in range (x): 
             f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
        j=j+1                  	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nVECTORS Is float\n')
  i=0 
  for j in range (x):
          f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0')
  i=i+1
  while i< (lol):
      for j in range (x):
          f.write('\n\t\t'+str(I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
      i=i+2
      if i==lol:
         break 
        #for j in range (x):
         # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
        #i=i+1
  f.close()
  return 





#==================================================faire apparaitre steamfield (mauvaise idée) ++++++++++++++++++++++++++++++++++ 
def steamflow(mat):
        w1=mat.shape[0]  
        w2=mat.shape[1]
        U=np.zeros((mat.shape[0],mat.shape[1]))
        V=np.zeros((mat.shape[0],mat.shape[1]))
        Y, X = np.mgrid[0:w1, 0:w2]
        j=0
        while j <w1:
           for i in range (w2):
             U[j,i]=mat[j,i]
             V[j,i]=0
           j=j+1
           if j==w1:
              break
           for i in range (w2):
             U[j,i]=0
             V[j,i]=mat[j,i] 
           j=j+1  
        fig = pyplot.figure(figsize=(20,10))
        gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 2])
        
        ax4 = fig.add_subplot()
        strm = ax4.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
        fig.colorbar(strm.lines)
        ax4.set_title('Varying Color')
        
        return pyplot.show
