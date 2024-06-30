import numpy as np
import matplotlib.pyplot as plt


##################################################
def config_elec(Z):
    s1 = 0
    s2p2 = 0
    
    for i in range(1,Z+1):
        if i in range(1,3):
            s1+=1
        if i in range(3,10):
            s2p2+=1
    return s1,s2p2

def Z_2S2P(config):
    A=config[1]
    Z=sum(config)-(0.35*(A-1))-(0.85*config[0])
    return Z        

def Z_1S(config):
    A=config[0]
    Z=sum(config)-(0.31*(A-1))
    return Z
##################################################
def E_tot(config):
    E=0
    E=(config[0]*(-13.6)*(Z_1S(config)**2)/1)
    E+=(config[1]*(-13.6)*(Z_2S2P(config)**2)/(2**2))
    return E

def energies():
    E = []
    for k in range(3,10):
        config = config_elec(k)
        E.append(E_tot(config))
    return E
##################################################
def ion(z):
    res=config_elec(z)
    Z=0
    ev=res[1]
    e1=1
    e2=0.85
    e3=0.35
    e4=0.31
    E1=2*(-13.6*((z-0.31)**2))
    E2=8*(-13.6*((z-3.8)**2)/4)
    if res[0]==1:
        Z=z-((res[1]-2)*0.31)
        E=(z-1)*(-13.6*(Z**2))
        return E
    if res[0]==2:
        Z=z-1.7-((res[1]-2)*0.35)
        E=E1+(res[1]-1)*(-13.6*((Z**2)/4))
        return E
    if res[0]==3:
        Z=z-8.8-((res[1]-2)*0.35)
        E=E1+E2+(res[1]-1)*(-13.6*((Z**2)/9))
        return E

def listeion():
    L = []
    for k in range(3,10):
        L.append(ion(k))
    return L

def listeatome():
    L = []
    for k in range(3,10):
        L.append(E_tot(config_elec(k)))
    return L

def ionisation():
    L = []
    for k in range(7):
        L.append(listeion()[k] - listeatome()[k])
    return L


##################################################

def rayoncov(Z):
    return 0.215*(4/Z_2S2P(config_elec(Z))) + 0.148*2 +0.225

def listerayon():
    L = []
    for k in range(3,10):
        L.append(rayoncov(k))
    return L
R = [1.33,1.02,0.85,0.75,0.71,0.63,0.64]

def enega(Z):
    return 0.34*(Z_2S2P(config_elec(Z))/(rayoncov(Z)**2))+0.7

def listeenega():
    return [enega(k) for k in range(3,10)]
N = [0.98,1.57,2.04,2.55,3.04,3.44,3.98]

def graphnega():
    nega_list = listeenega()
    x = list(range(3,10))
    fig, ax = plt.subplots()

    ax.plot(x, nega_list, label='Méthode de Allred-Rochow', marker='o', color='red')

    ax.plot(x,N, label='Valeurs de Pauling', marker='o', color='black')

    ax.legend()
    ax.set_xlabel('Numero atomique')
    ax.set_ylabel('Electronégativité')
    ax.set_title('')
    plt.grid(True)
    plt.show()
    return

def graphrayon():
    rayon_list = listerayon()
    x = list(range(3, 10))
    fig, ax = plt.subplots()
    
    ax.plot(x, rayon_list, label='Rayons de covalence théoriques', marker='o', color='red')
    
    ax.plot(x, R, label='Valeurs expérimentales', marker='o', color='black')
    
    ax.legend()
    
    ax.set_xlabel('Numéro atomique')
    ax.set_ylabel('Rayon (A°)')
    ax.set_title("")
    plt.grid(True)  

    plt.show()
    return

C = [5.3917, 9.3227, 8.298, 11.2603, 14.5341, 13.6181, 17.4228]
def graphion():
    ionisation_vals = ionisation()
    x = list(range(3, 10))
    fig, ax = plt.subplots()
    
    ax.plot(x, ionisation_vals, label='Approx. de Slater', marker='o',color='red')
    
    ax.plot(x, C, label='Valeurs expérimentales', marker='o',color='black')
    
    ax.legend()
    
    ax.set_xlabel('Numéro atomique')
    ax.set_ylabel('Énergies (eV)')
    ax.set_title("Énergies d'ionisation pour les élements de la 2e ligne du tableau périodique.")
    plt.grid(True)  # Activer la grille

    plt.show()
    return

def graph():
    energies_list = energies()
    plt.figure()  
    plt.plot(energies_list, marker='o', linestyle='-', color='red') 
    plt.xlabel('Numéro atomique')  
    plt.ylabel('Énergie (eV)')  
    plt.title("Approximation de l'énergie des éléments selon la méthode de Slater.")
    plt.grid(True)  

    valeurs_x = range(len(energies_list))
    etiquettes_x = ['3', '4', '5', '6', '7', '8', '9']
    plt.xticks(valeurs_x, etiquettes_x)

    plt.show()
    return


graph()
graphrayon()
graphion()
graphnega()
