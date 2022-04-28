import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# Date de intrare

Ma = 22.8 # kg/s
D1m = 0.409 # m
A1 = 0.303 # m^2
Zc = 9 # nr. de trepte
nn = 8400 # rpm
alfa = 90 * np.pi/180 # rad
etac = 0.85
pic = 6.8

p1 = 101325 # Pa
T1 = 288 # K
k = 1.4 
R = 287 # J/kgK
cp = R*k/(k-1) # J/kgK

# Calculul parametrilor relevanti pentru regimul nominal

Kk = np.sqrt(k*((2/(k+1))**((k+1)/(k-1))))

lc = k/(k-1)*R*T1*(pic**((k-1)/k)-1)/etac

lctr1 = 0.9*lc/Zc

ql1n = Ma*np.sqrt(R*T1)/(p1*A1*np.sin(alfa)*Kk)

l1n = 2/np.pi*np.arcsin(1-(1-ql1n)**(np.pi/2*(2/(k+1))**(1/(k-1))))

C1an = l1n*np.sin(alfa)*np.sqrt(2*k/(k+1)*R*T1)

U1m = D1m*np.pi*nn/30/2

C1anb = C1an/U1m

lunb = lctr1/U1m**2

# Parametrii curgerii în miscarea relativa in sectiunea de intrare

W1 = np.sqrt(C1an**2+U1m**2+2*U1m*C1an*np.cos(alfa))

Tw1 = T1+(k-1)*(W1**2-C1an**2)/(2*k*R)

l1w = W1/np.sqrt(2*k*R*Tw1/(k+1))

# Trasarea liniei de pompaj

# Constantele specifice

K1 = 1-lunb

Keta = (1/lunb-1)**2

Xcp = (Keta*(1+2*K1)-np.sqrt(Keta**2*(1+2*K1)**2-3*K1*Keta*(2*Keta-K1+K1*Keta)))/(3*K1*Keta)

if Xcp < 0:
    Xcp = (2*K1+1-np.sqrt(4*K1-2*K1+1))/3/K1

C1apb = Xcp*C1anb

# Generarea liniei de pompaj

n = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1]) * nn

etacp = etac*(1-Keta*(1-Xcp)**2)*(2-n/nn)*n/nn

picp = (1+etacp*(n/np.sqrt(T1))**2*(1-K1*Xcp)*(np.pi*D1m/60)**2*Zc/0.9/cp)**(k/(k-1))

def q(l):
    return l*((k+1)/2*(1-(k-1)/(k+1)*l**2))**(1/(k-1))

l1p = Xcp*np.pi*D1m/60*C1anb/np.sin(alfa)*np.sqrt((k+1)/(2*k*R))*n/np.sqrt(T1)

Map = A1*np.sin(alfa)*q(l1p)*Kk/np.sqrt(R)


# Trasarea curbelor de n=ct

def plot_constant_n():
    
    def solve_equations(z, pict, nct):
        C1ab = z[0]
        etact = z[1]
        F = np.empty(2)
        F[0] = C1anb/(1-lunb)*(1-0.9*cp/(Zc*etact)*(60/(np.pi*D1m))**2*
                ((pict**((k-1)/k)-1)/((nct/np.sqrt(T1))**2))) - C1ab
        F[1] =  etac*(1-Keta*(C1ab/C1anb-1)**2)*(2-nct/nn)*nct/nn - etact
        return F

    picn = np.array([np.linspace(picp[0], picp[0], 10), np.linspace(picp[0], picp[1], 10), 
                    np.linspace(picp[1], picp[2], 10), np.linspace(picp[2], picp[3], 10), 
                    np.linspace(picp[3], picp[4], 10), np.linspace(picp[4], picp[5], 10), 
                    np.linspace(picp[5], picp[6], 10)])

    i, j = 0, 0
    MaT1p1 = []
    
    while(i < 7):
        if(j < 10):
            zGuess = np.array([1,1])
            z = fsolve(solve_equations, zGuess, (picn[i][j], n[i]))
            l1 = z[0]*np.sqrt((k+1)/2/k/R)*np.pi*D1m*n[i]/60/np.sqrt(T1)
            to_append = A1*np.sin(alfa)*q(l1)*Kk/np.sqrt(R)
            MaT1p1.append(to_append)
            j += 1
        else:
            plt.plot(MaT1p1, picn[i], label="n = " + str(int(n[i])) + " rpm")
            j = 0
            i += 1
            MaT1p1 = []
            

# Linia regimurilor optime

plt.figure(figsize=(8, 5))
plt.plot(Map, picp, label="Linia de pompaj")
plot_constant_n()
plt.xlabel("Map")
plt.ylabel("picp")
plt.title("Caracteristica compresorului axial")
plt.grid()
plt.legend(bbox_to_anchor = (1.01, 1))
plt.savefig("curbe_turatie_constanta_1.png", bbox_inches="tight", dpi=300)
plt.show()

