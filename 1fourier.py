#!/usr/bin/env python3.8

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import csv
from scipy import stats
import scipy.fftpack

seuil_detection_fitness=0.0005


def analyse_ligne(row): 
    res=[]
    for i in range(int(len(row)/3-1)):
        if(row[3*i+2]!="NA"):
            nouv=float(row[3*i+2].replace(",","."))
            if nouv<seuil_detection_fitness: nouv=0
            res+=[nouv]
        else:
            res+=[-1] 
    return res

def analyse_fichier(): # ressort [W_t pour tout t] pour le fichier lu
    wt_tot=[]
    wt_cur=[]
    with open("doi_10/All_data/data_microMA/dataset_microMA_MutH.csv", newline='') as csvfile:
        reader=csv.reader(csvfile,delimiter=';')
        nrow=0
        for row in reader:
            if(nrow>=2):
                wt_cur=analyse_ligne(row)
                wt_tot+=[wt_cur]
            nrow+=1
    return wt_tot



# Lancement :

wt_real=analyse_fichier()

# PARAMÈTRE SUR LEQUEL IL FAUT JOUER : A ↓
A=5 # la fréquence maximale (le A de nos calculs)
nbreaks_xi=100 # précision sur xi

xmin=0
xmax=1
nbreaks_x=100 # précision sur x

x=np.linspace(xmin,xmax,nbreaks_x)
xi=np.linspace(0,A,nbreaks_xi)

m=[1,1.9,68,39,250,180,130,970,740,560,430] # liste des 10 premiers moments donnés dans l'article


def part_sum(m,xi):
    N=len(m)
    tm_sum=[(1j*xi)**k/np.math.factorial(k)*m[k] for k in range(N)]
    res=sum(tm_sum)
    return res

# transformée de fourier (pas inverse, du coup!) de la fonction caractéristique:
# x=valeurs en lesquelles on calcule Ff
# xi tel que: (fonction mathématiques Ff->) Ff(xi[i])=Ff[i] (<-tableau ff dans python)
def fourier_inverse(x,xi,Ff): 
    res=[]
    n=len(Ff)
    for xp in x:
        fx=0
        for i in range(n):
            fx+=Ff[i]*np.exp(-1j*xp*xi[i])
        res+=[fx/n]
    return res

phi_xi=[part_sum(m,xip) for xip in xi] # calcul de phi_chapeau
f=fourier_inverse(x,xi,phi_xi) # on en déduit f chapeau

# affichage: 
plt.plot(x,f)
plt.show()


