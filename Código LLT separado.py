# -*- coding: utf-8 -*-
"""
Created on Tue May 18 17:50:35 2021

@author: davim
"""

import numpy as np
from math import pi
import matplotlib.pylab as plt   
import pandas as pd
from sympy import symbols


N = 1000 #Nº de segmentos
        
#Ângulos

twist = 0

#perfil

Cla = 0.104     # valor em 1/grau
alpha_0 = -10.66

#Configuação da Aeronave  

i_w = 5

Tabela = False #Se você quiser que o programa crie uma tabela de CL vs alpha

# lista_Cr = [0.54,0.54]    #asa de 2022
# lista_Ct = [0.54,0.24]    #asa de 2022
# lista_b  = [0.5,1.00]     #asa de 2022

# lista_Cr = [Cr,C1,C2]
# lista_Ct = [C1,C2,Ct]
# lista_b  = [b1,b2,b3]

lista_Cr = [0.6,0.6]
lista_Ct = [0.6,0.4]
lista_b  = [0.6,1.15]


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


# Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # root chord (m)  

b = lista_b[len(lista_b)-1]*2

Cla = Cla * 180/pi

alpha_twist = twist         
a_2d = Cla

theta = np.linspace((pi / (2 * N)), (pi / 2), N, endpoint=True)


alpha = np.linspace(i_w + alpha_twist, i_w, N) 
z = (b / 2) * np.cos(theta)

listay = listay[::-1]

c = np.array(listay)
# print(c)
# c = Croot * (1 - (1 - taper) * np.cos(theta))  # Mean Aerodynamics
# print(c)

mu = c * a_2d / (4 * b)

LHS = mu * (np.array(alpha) - alpha_0) / (180/pi)  # .reshape((self.N-1),1)# Left Hand Side

RHS = []
for i in range(1, 2 * N + 1, 2):
    RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))  # .reshape(1,self.N)
    # print(RHS_iter,"RHS_iter shape")
    RHS.append(RHS_iter)
    


test = np.asarray(RHS)
x = np.transpose(test)
inv_RHS = np.linalg.inv(x)
    
A = np.matmul(inv_RHS, LHS) #An

# A=inv_RHS*LHS
    
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
for i in range(0,N,1):
    z = i*(b/2)/N
    y_sr.append(z)

y_sr.pop(0)
y_sr.append(b/2)

y_sr = y_sr[::-1]

print('\nCL = ',round(CL_wing,5))

b=b/2

try:
    plt.style.use('extensys-gd')
except:
    pass

plt.figure(num=None, figsize=(10.5, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(y_sr,CLvec,color='r')
# plt.title("Distribuição elíptica de sustentação\nLifting Line Theory")
plt.xlabel("b [m]")
plt.ylabel("CL")
# plt.axis([0,b,0,max(CLvec)])
plt.show()

px=[0,0]
py=[0,lista_Cr[0]]

bx = [lista_b[len(lista_b)-1],lista_b[len(lista_b)-1]]
by = [0,lista_Ct[len(lista_Ct)-1]]

cx = [0,lista_b[len(lista_b)-1]]
cy = [0,0]
        
plt.figure(num=None, figsize=(10.5, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(listax,listay,color='blue')
plt.plot(cx,cy,color='b')
plt.plot(px,py,color='b')
plt.plot(bx,by,color='b')
plt.title('Desenho da Superfície')
plt.axis([-0.1,lista_b[len(lista_Ct)-1]+0.1,-0.1,max(lista_Cr)+0.1])
plt.xlabel('b[m]')
plt.ylabel('C[m]')
plt.show()


valor_max = max(listay)

valor_porc_corda = []

for i in listay:
    porcentagem = i*100/valor_max
    valor_porc_corda.append(porcentagem)
    

valor_porc_CL = []

valor_max_CL = max(CLvec)

for i in CLvec:
    porcentagem = i*100/valor_max_CL
    valor_porc_CL.append(porcentagem)
    
valor_porc_CL.pop(len(valor_porc_CL)-1)

valor_porc_CL = list(reversed(valor_porc_CL))

    
px=[0,0]
py=[0,lista_Cr[0]*100/valor_max]

bx = [lista_b[len(lista_b)-1],lista_b[len(lista_b)-1]]
by = [0,lista_Ct[len(lista_Ct)-1]*100/valor_max]

cx = [0,lista_b[len(lista_b)-1]]
cy = [0,0]

plt.figure(num=None, figsize=(10.5, 7), dpi=200, facecolor='w', edgecolor='k')
plt.plot(listax,valor_porc_corda,color='blue',label='Asa')
plt.plot(listax,valor_porc_CL,color='r',label='CL')
plt.plot(cx,cy,color='b')
plt.plot(px,py,color='b')
plt.plot(bx,by,color='b')
# plt.title('Desenho da Superfície')
plt.axis([-0.1,lista_b[len(lista_Ct)-1]+0.1,-10,110])
plt.xlabel('b[m]')
plt.ylabel('% C[m] ou % CL')
plt.legend(bbox_to_anchor=(0.9, 0.85, 0.25, 0.1), ncol=1, frameon=True, shadow=True)
plt.show()

dados_dist_CL={
    'b local':y_sr[::-1],
    'CL':CLvec[::-1],
}


df=pd.DataFrame(dados_dist_CL)
df.to_excel('Distribuição_de_sustentação.xlsx',index=False)


'_____________________________________________________________________________________________'
"Para achar curva de cl vs alpha"

if Tabela == True:

    vetorcl=[]
    vetoral=[]
    vetor1=[]
    vetor2=[]
    vetor3=[]
    vetorp=[]
    listaal=np.arange(0.1,12,0.1)
    
    print('Valores de Cl vs alpha')
    for i_w in listaal:
    
                
                
                
                #%% Cálculos ###
        
        
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

        listax = np.linspace(0,1,N) #divido N vezes
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
            
            lista_C.append(cmac)
            lista_S.append(S_i)

        soma = sum(lista_S)
        soma_cima = 0
        for i in range(len(lista_C)):
            
            soma_cima = soma_cima + (lista_C[i]*lista_S[i])
            

        MAC = soma_cima/soma


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


        # Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # root chord (m)     

        Cla = Cla * 180/pi

        alpha_twist = twist         
        a_2d = Cla

        theta = np.linspace((pi / (2 * N)), (pi / 2), N, endpoint=True)


        alpha = np.linspace(i_w + alpha_twist, i_w, N) 
        z = (b / 2) * np.cos(theta)

        c = np.array(listay)
        # print(c)
        # c = Croot * (1 - (1 - taper) * np.cos(theta))  # Mean Aerodynamics
        # print(c)

        mu = c * a_2d / (4 * b)

        LHS = mu * (np.array(alpha) - alpha_0) / 57.3  # .reshape((self.N-1),1)# Left Hand Side

        RHS = []
        for i in range(1, 2 * N + 1, 2):
            RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))  # .reshape(1,self.N)
            # print(RHS_iter,"RHS_iter shape")
            RHS.append(RHS_iter)
            


        test = np.asarray(RHS)
        x = np.transpose(test)
        inv_RHS = np.linalg.inv(x)
            
        A = np.matmul(inv_RHS, LHS) #An

        # A=inv_RHS*LHS
            
        AR = b/MAC         #Aspect ratio

        CL_wing = (pi * AR * A[0])
            
            
        "_____________________________________________________________________"
                

        
        vetorcl.append(round(CL_wing,5))
        # print('CL = ',round(CL_wing,5))
        # print(round(CL_wing,4))
    
        vetoral.append(round(i_w,1))
        vetor1.append('(')
        vetor2.append(')')
        vetor3.append(',')
        vetorp.append('pontuar')
    
    
    dados_Cl_VT={
        'p':vetorp,
        'a':vetor1,
        'Alpha':vetoral,
        'b':vetor3,
        'CL':vetorcl,
        'c':vetor2,
    }
    
    
    df=pd.DataFrame(dados_Cl_VT)
    df.to_csv('Curva de CL.txt',index=False,sep='\t')
    df.to_excel('Curva de CL.xlsx',index=False)
    
    print('Feito!')
    
    # print(vetoral)
    # print(vetorcl)

#EXEMPLO USADO NO LIVRO:
    

# N = 10 #Nº de segmentos
        
# #Geometria
# S = 25             
# b=14.1421
# MAC=1.76777
 

# #Ângulos
# taper=0.6
# twist = -1
# sweep = 0

# #perfil
# Cla = 6.3
# alpha_0 = -1.5

#         #Configuração da Aeronave  
# i_w = 2
    
    
