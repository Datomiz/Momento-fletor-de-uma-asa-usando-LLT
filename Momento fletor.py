# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:35:31 2023

@author: Davi
"""

'Código que faz o gráfico de esforço cortante e Momento Fletor da asa'

import numpy as np
from math import pi
import pylab as plt   
import pandas as pd
from sympy import symbols


'Número de sementos que a asa será separada em'

N = 1150 #Nº de segmentos
        
'Ângulo de torção da asa'

twist = 0

'Dados do perfil aerodinâmico'

Cla = 0.10606     # valor em 1/grau
alpha_0 = -10.5   # [graus]

'Ângulo de incidência/ataque'  

i_w = 15           #[graus]

'Constantes para o calculo da força'

phi = 1.20        #[kg/m³] Densidade do ar na temperatura analisada 
V = 12            #[m/s] Velocidade do vento



# lista_Cr = [0.54,0.54]    #asa de 2022
# lista_Ct = [0.54,0.24]    #asa de 2022
# lista_b  = [0.5,1.00]     #asa de 2022

# lista_Cr = [Cr,C1,C2]
# lista_Ct = [C1,C2,Ct]
# lista_b  = [b1,b2,b3]

lista_Cr = [  0.6 , 0.6 ]
lista_Ct = [  0.6 , 0.4 ]
lista_b = [  0.6 , 1.15 ]















'___________________________________________________________________________'

'Código em si, NÃO ALTERAR'

'___________________________________________________________________________'

'Calculo das linhas que representam a asa'

x=symbols('x')

segmentos=[]

for i in range(len(lista_b)):
    if i == 0:
        linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-0))*(x-0))
    if i > 0:     
        linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-lista_b[i-1]))*(x-lista_b[i-1]))
    
    segmentos.append(linha)

listax = np.linspace(0,lista_b[len(lista_b)-1],N) #divido N vezes
listay = []
listacont = np.arange(1,len(lista_b),1)

for x in listax:
    
    if x < lista_b[0]:
        y = eval(str(segmentos[0]))
    
    if x > lista_b[0]:
        for i in listacont:
            if x > lista_b[i-1] and x < lista_b[i]:
                y = eval(str(segmentos[i]))
                
    listay.append(y)

'___________________________________________________________________________'

'calculo da corda média'

lista_S = []
lista_C = []

for i in range(len(lista_b)):
    
    TR_i = lista_Ct[i]/lista_Cr[i]
    cmac=(2/3)*lista_Cr[i]*((1+TR_i+(TR_i**2))/(1+TR_i))
    
    if i == 0:
        S_i = lista_b[i]*cmac
    if i > 0:
        S_i = (lista_b[i]-lista_b[i-1])*cmac
    
    S_i = 2*S_i
    
    lista_C.append(cmac)
    lista_S.append(S_i)

S = sum(lista_S)
soma_cima = 0
for i in range(len(lista_C)):
    
    soma_cima = soma_cima + (lista_C[i]*lista_S[i])
    

MAC = soma_cima/S


'___________________________________________________________________________'


'Calculo do TR efetivo'

R_TR = []
R_S  = []

lista_contar= range(len(lista_Cr))

for i in lista_contar:
    
    Cri = lista_Cr[i]
    Cti = lista_Ct[i]
    
    if i > 0:
        bi = lista_b[i] - lista_b[i-1]
    else:
        bi  = lista_b[i]
    
    TRi=Cti/Cri
    cmaci=(2/3)*Cri*((1+TRi+(TRi**2))/(1+TRi))
    Si=bi*cmaci
    
    Si = Si*2 #isso pq os valores de b são pra só um lado
    
    R_TR.append(TRi)
    R_S.append(Si)

soma = sum(R_S)
medias = []

for i in lista_contar:
    me = R_S[i]*R_TR[i]
    medias.append(me)

cima = sum(medias)

taper = cima/soma

'___________________________________________________________________________'


     #%% Cálculos ###


b = lista_b[len(lista_b)-1]*2

Cla = Cla * 180/pi

alpha_twist = twist         
a_2d = Cla

theta = np.linspace((pi / (2 * N)), (pi / 2), N, endpoint=True)

alpha = np.linspace(i_w + alpha_twist, i_w, N) 
z = (b / 2) * np.cos(theta)

listay = listay[::-1]

c = np.array(listay)


mu = c * a_2d / (4 * b)

LHS = mu * (np.array(alpha) - alpha_0) / (180/pi)

RHS = []
for i in range(1, 2 * N + 1, 2):
    RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))
    RHS.append(RHS_iter)
    


test = np.asarray(RHS)
x = np.transpose(test)
inv_RHS = np.linalg.inv(x)
    
A = np.matmul(inv_RHS, LHS) 

AR = b/MAC        #Aspect ratio

CL_wing = (pi * AR * A[0]) 

listay = listay[::-1]
    
"_____________________________________________________________________"
        
        # Plot CL distribution
        

k = np.divide((4 * b), c) #4b/c_i
    
soma = np.zeros((N,N))
          
for j in range(1,N+1,1):

             row_to_be_added = (np.sin((2*j-1)*theta)) * A[j-1] * k    
             soma[j-1] =  row_to_be_added         
        
CL = soma.sum(axis=0)
CLvec = np.append(0.001, CL)# Colocar o zero
y_sr=[]
y_sr.append(b/2)


lista_bs = y_sr
lista_CLs = CLvec[::-1]


for i in range(0,N,1):
    z = i*(b/2)/N
    y_sr.append(z)

y_sr.pop(0)
y_sr.append(b/2)


y_sr = y_sr[::-1]

print('\nCL = ',round(CL_wing,5))


'Esforço cortante'

area = 0

for i in range(len(lista_bs)-1):
    
    h_tri = abs(lista_CLs[i] - lista_CLs[i+1])
    
    base  =  lista_bs[i+1] - lista_bs[i]
    
    h_ret = min([abs(lista_CLs[i]),abs(lista_CLs[i+1])])

    a_tri =  (base * h_tri)/2
    
    a_ret = base * h_ret

    area = area + a_tri + a_ret

resultante_real = area*S*(V**2)*phi/2

Resultante_LLT = sum(lista_CLs)/N

fator = resultante_real/(Resultante_LLT*phi*(V**2)*S/2)

cortante = []
vx = []

for i in range((len(lista_CLs)*2)):
    
    if i == 0:
        
        L = fator*Resultante_LLT*phi*(V**2)*S/2
        
        cortante.append(L)
        
        vx.append(lista_bs[int(i/2)])
        
    else:
        
        if i % 2 != 0: #impar
        
            L = L - (fator*lista_CLs[int(i/2)]*phi*(V**2)*S/(2*N))
            
            cortante.append(L)
            vx.append(lista_bs[int(i/2)])
        
        else: #par
            
            cortante.append(L)
            vx.append(lista_bs[int(i/2)])
    
'Momento Fletor'

momento = []

for i in range(len(cortante)):
    
    cort = list(cortante[::-1])
    
    if i == 0:
        
        momento.append(cort[0]*b/2)
        momento.append(cort[0]*b/2)
    
    if i % 2 == 0 and i !=0: #par e diferente de 0

        M = momento[i-1] - (cort[i]*(vx[i]-vx[i-1]))  
        
        momento.append(M)
        momento.append(M)

momento = momento[::-1]

        
print('\nEsforço máximo =',round(max(cortante),3),'N')
print('Momento máximo =',round(abs(min(momento)),3),'N*m')


try:
    plt.style.use('extensys-gd')
except:
    pass



plt.figure(num=None, figsize=(7, 5), dpi=200, facecolor='w', edgecolor='k')
plt.plot(y_sr,CLvec)
plt.title("Distribuição elíptica de sustentação\nLifting Line Theory")
plt.xlabel("Envergadura [m]")
plt.ylabel("Coeficiente de sustentação [CL]")
plt.show()


plt.figure(num=None, figsize=(7, 5), dpi=200, facecolor='w', edgecolor='k')
plt.plot(vx,cortante,color='b')
plt.title("Esforço Cortante na asa V ="+str(V))
plt.xlabel("Envergadura [m]")
plt.ylabel("Força [N]")
# plt.axis([0,0.1,400,470]) #para V = 24 m/s
plt.show()

plt.figure(num=None, figsize=(7, 5), dpi=200, facecolor='w', edgecolor='k')
plt.plot(vx,momento,color='b')
plt.title("Momento Fletor na asa V ="+str(V))
plt.xlabel("Envergadura [m]")
plt.ylabel("Momento Fletor [N*m]")
plt.show()


arquivo = {'V [m/s]' : V,
           'Rho [Kg/m³]': phi,
           'b [m]': vx,
           'Cortante [N]': cortante,
           'Fletor [N*m]': momento}

df = pd.DataFrame(arquivo)
df.to_excel('Esforço e Momento da asa.xlsx',index=False)



print('\nTudo feito! Cheque a pasta onde está o programa, lá tem um arquivo de Excel com os resultados')


'Desenho da asa'

x=symbols('x')

segmentos=[]

for i in range(len(lista_b)):
    if i == 0:
        linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-0))*(x-0))
    if i > 0:     
        linha = lista_Cr[i] + (((lista_Ct[i]-lista_Cr[i])/(lista_b[i]-lista_b[i-1]))*(x-lista_b[i-1]))
    
    segmentos.append(linha)

listax = np.linspace(0,lista_b[len(lista_b)-1],100)
listay = []
listacont = np.arange(1,len(lista_b),1)

for x in listax:
    
    if x < lista_b[0]:
        y = eval(str(segmentos[0]))
    
    if x > lista_b[0]:
        for i in listacont:
            if x > lista_b[i-1] and x < lista_b[i]:
                y = eval(str(segmentos[i]))
                
    listay.append(y)

listax2 = []

for i in range(len(listay)):
    listax2_v = -listax[i]      
    listax2.append(listax2_v)

for i in range(len(listay)):
    listay[i] = -listay[i]

px = [-lista_b[len(lista_b)-1],-lista_b[len(lista_b)-1]]
py = [0,-lista_Ct[len(lista_Ct)-1]]

bx = [lista_b[len(lista_b)-1],lista_b[len(lista_b)-1]]
by = [0,-lista_Ct[len(lista_Ct)-1]]

cx = [-lista_b[len(lista_b)-1],lista_b[len(lista_b)-1]]
cy = [0,0]

plt.figure(num=None, figsize=(10.5, 7.5), dpi=200, facecolor='w', edgecolor='k')
plt.axis([-0.1-lista_b[len(lista_b)-1],lista_b[len(lista_Ct)-1]+0.1,-max(lista_Cr)-0.05,0.05])
plt.title('Desenho da asa')
plt.plot(listax,listay,color='blue')
plt.plot(listax2,listay,color='blue')
plt.plot(cx,cy,color='b')
plt.plot(px,py,color='b')
plt.plot(bx,by,color='b')
plt.xlabel('b[m]')
plt.ylabel('C[m]')
ax = plt.gca()
ax.set_aspect('equal',adjustable = 'box')
plt.show()





