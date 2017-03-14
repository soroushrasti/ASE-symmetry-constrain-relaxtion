# this only works for water molecule
from __future__ import print_function
import numpy as np
from numpy import *
from decimal import *

numberline = 1
f3 =  open('supercell.data')
for line in f3:
  numberline  = numberline +1
  if line.startswith("Atoms"):
     break

f = open("p", "w")       
with open("supercell.data") as f3:
    for line in f3:
        line  = line.rstrip('\n')
        if "atoms" in line:
             atoms=line.rstrip('\n').split(" ")
             print ("\n", file=f)
             print (line, file=f)
             print ("2 atom types", file=f)
             print (int(atoms[0])*2/3," bonds", file=f)
             print ("1 bond types", file=f)
             print (int(atoms[0])/3," angles", file=f)
             print ("1 angle types", file=f)
             print ("\n", file=f)
        if "xlo" in line:
             xpard = line.split(" ")
             print (line, file=f)
        if "ylo" in line:
             ypard = line.split(" ")
             print (line, file=f)
        if "zlo" in line:
             zpard = line.split(" ")
             print (line, file=f)
        if "xy" in line:
             trilin = line.split(" ")
             print (line, file=f)
             print ("\n", file=f)

atomnumber= int(atoms[0])
xpard= -float(xpard[0]) + float(xpard[1])
ypard= -float(ypard[0])+ float(ypard[1])
zpard= -float(zpard[0])+ float(zpard[1])
xypard = float(trilin[0])
xzpard = float(trilin[1])
yzpard = float(trilin[2])

bondlenght=1
oxygen=2
hydrogen=1
mol= 999*np.ones(atomnumber)
out = np.zeros(atomnumber)
charge= np.zeros(atomnumber)
data = np.genfromtxt("supercell.data", skip_header=numberline, max_rows=atomnumber)

typ,x,y,z = data[:,1], data[:,2], data[:,3], data[:,4]  
#typ,x,y,z = data[:,2], data[:,4], data[:,5], data[:,6] 


dis = 100*np.ones(shape=(atomnumber,atomnumber))
for i in range(0,atomnumber):
   for j in range(0,atomnumber):
        if i!=j:
             dis[i][j] = math.sqrt((x[j] - x[i])**2 + (y[j] - y[i])**2 + (z[j] - z[i])**2) 

kk=0
for i in range(0,atomnumber):
   if typ[i]==oxygen:
       kk= kk+1
       mol[i] = kk
       if typ[int(dis[:,i].argsort()[:2][0])]==hydrogen:
            if dis[int(dis[:,i].argsort()[:2][0])][i]<bondlenght:
                mol[int(dis[:,i].argsort()[:2][0])] = kk
            else:
                out[i]=1
       if typ[int(dis[:,i].argsort()[:2][1])]==hydrogen:
            if dis[int(dis[:,i].argsort()[:2][1])][i]<bondlenght:
                mol[int(dis[:,i].argsort()[:2][1])] = kk
            else:
                out[i]=1
       

for i in range(0,atomnumber):
   if mol[i]==999:
      for j in range(0,atomnumber):
           if out[j]==1:
                if y[i]-y[j]> ypard/2:
                      dis[i][j] = math.sqrt((y[j] - (y[i]-ypard))**2 + (x[j] - (x[i]-xypard))**2 + (z[j] - z[i])**2)
                      if dis[i][j] < bondlenght:
                           break
                if -y[i]+y[j]> ypard/2:
                      dis[i][j] = math.sqrt((y[j] - (y[i]+ypard))**2 + (x[j] - (x[i]+xypard))**2 + (z[j] - z[i])**2)
                      if dis[i][j] < bondlenght:
                           break
                if x[i]-x[j]> xpard/2:
                      dis[i][j] = math.sqrt((y[j] - y[i])**2 + (x[j] - (x[i]-xpard))**2 + (z[j] - z[i])**2)
                      if dis[i][j] < bondlenght:
                           break
                if -x[i]+x[j]> xpard/2:
                      dis[i][j] = math.sqrt((y[j] - y[i])**2 + (x[j] - (x[i]+xpard))**2 + (z[j] - z[i])**2)
                      if dis[i][j] < bondlenght:
                           break
                if -z[i]+z[j]> zpard/2:
                      dis[i][j] = math.sqrt((y[j] - y[i])**2 + (x[j] - x[i])**2 + (z[j] - (z[i]+zpard))**2)
                      if dis[i][j] < bondlenght:
                           break
                if z[i]-z[j]> zpard/2:
                      dis[i][j] = math.sqrt((y[j] - y[i])**2 + (x[j] - x[i])**2 + (z[j] - (z[i]-zpard))**2)
                      if dis[i][j] < bondlenght:
                           break

kk=0
for i in range(0,atomnumber):
   if typ[i]==oxygen:
       kk= kk+1
       mol[i] = kk
       if typ[int(dis[:,i].argsort()[:2][0])]==hydrogen:
            if dis[int(dis[:,i].argsort()[:2][0])][i]<bondlenght:
                mol[int(dis[:,i].argsort()[:2][0])] = kk
       if typ[int(dis[:,i].argsort()[:2][1])]==hydrogen:
            if dis[int(dis[:,i].argsort()[:2][1])][i]<bondlenght:
                mol[int(dis[:,i].argsort()[:2][1])] = kk



 # charge for each atom type
for i in range(0,atomnumber):
    if typ[i]==oxygen:
        charge[i]=-0.8476
    if typ[i]==hydrogen: 
        charge[i]=0.4238





    
# print coordinate  
print( '\nAtoms\n', file=f)               
for i in range(0,atomnumber):
   print( i+1, int(mol[i]), int(typ[i]),charge[i],x[i],y[i],z[i], file=f)
# print bonds
bb=1
print ('\nBonds\n', file=f)
for i in range(0,atomnumber):
    if typ[i]==oxygen:
        for j in range(0,atomnumber):
            if i!=j:
                if mol[i] == mol[j]:
                   print (bb,1,i+1,j+1, file=f) 
                   bb= bb+1
    
# print angles 
bb=1
ll=1
print ('\nAngles\n', file=f)
for i in range(0,atomnumber):
    if typ[i]==oxygen:
        for j in range(0,atomnumber):
            if i!=j:
                if mol[i] == mol[j]:
                    for k in range(0,atomnumber):
                        if k!=j:
                            if k!=i:
                               if mol[i] == mol[k]:
                                 if ll%2!=0: 
                                    print( bb,1,j+1,i+1, k+1, file=f) 
                                    bb= bb+1
                                 ll= ll+1
f.close()        
