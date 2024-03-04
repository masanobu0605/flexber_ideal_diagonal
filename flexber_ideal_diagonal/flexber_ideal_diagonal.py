import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
import sympy as sp
from scipy.special import *
'''
##--変数定義．大文字のみ定義すればいい--

L = 2.54
E = 70*10**9
D = 4.2*10**(-3) #円柱断面形状を想定している．
alpha = 60 #degで入力していいよ
sigma_y = 210 * 10**6

H,B = 4.2*10**(-3),4.2*10**(-3) #Hは高さ，Bは幅
I = B*H**3/12
Zp = B * H**2/6

#I = np.pi*(D**4)/64
#Zp = D**3/6 #全塑性モーメントの計算に使う．円柱断面形状です
alpha = alpha * np.pi /180
ei = E * I


##--先端たわみ角からpの算出--
phi_0_data = np.loadtxt("phi_0.csv" ,delimiter=",", encoding="shift-jis") #先端たわみ角のimport
phi_0_data = phi_0_data * np.pi / 180 #先端たわみ角をラジアン表記に直す
p = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    p[i] = np.sin((alpha + phi_0_data[i]) * 0.5)
print("進捗度10%|pの導出完了")

##!--mの算出(mは行う積分の下限になる．0 -> m -> pi/2 の順番で大きい)---
m = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    m[i] = np.arcsin(np.sin(alpha* 0.5)  / p[i])
print("進捗度20%|m の導出完了")

##--関数の準備F(start_point,stop_point,m)------------------
def Finte(a1,a2,p):
    n = 500
    d_theta = (a2 - a1) / n
    theta_before = a1
    Fanswer = 0
    
    for i in range(n):
        theta_after = theta_before + d_theta
        Fanswer = Fanswer  + (((1 - (p**2)*(np.sin(theta_before))**2)**(-0.5)) + ((1 - (p**2)*(np.sin(theta_after))**2)**(-0.5)) )*d_theta*0.5
        theta_before = theta_after
    return Fanswer
##--関数の準備E(start_point,stop_point,m)------------------
def Einte(a1,a2,p):
    n = 500
    d_theta = (a2 - a1) / n
    theta_before = a1
    Eanswer = 0
    
    for i in range(n):
        theta_after = theta_before + d_theta
        Eanswer = Eanswer + (((1 - (p**2)*(np.sin(theta_before))**2)**(0.5)) + ((1 - (p**2)*(np.sin(theta_after))**2)**(0.5)))*d_theta*0.5
        theta_before = theta_after
    return Eanswer

##--kの算出---
k = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    k[i] = (Finte(m[i],np.pi * 0.5,p[i])/L)
print("進捗度25%|k の導出完了")

##--nの算出---
n = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    n[i] = np.arcsin(np.sin((phi_0_data[i] + alpha) * 0.5)/p[i])
print("進捗度25%|k の導出完了")

##--Pの導出
P_load = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    P_load[i] = (k[i]**2 * ei)
print("進捗度50%|P_loadの導出完了")

##--xc水平方向変位の導出の導出/
xc = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    xc[i] = (1/k[i]) * (np.cos(alpha) * (Finte(0,m[i],p[i]) - Finte(0,np.pi * 0.5,p[i]) + 2 * Einte(m[i],np.pi * 0.5,p[i])) 
                        + 2 * p[i] *np.sin(alpha) * np.cos(m[i]))
print("進捗度70%|水平方向変位の導出完了")

##--yc垂直方向変位
yc = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    yc[i] = (1/k[i]) * (2*p[i]* np.cos(alpha) * np.cos(m[i]) - 
                        np.sin(alpha) * (Finte(0,m[i],p[i]) - Finte(0,np.pi * 0.5 ,p[i]) + 2*Einte(0,np.pi *0.5, p[i]) - 2*Einte(0,m[i], p[i])))
print("進捗度90%|垂直方向変位の導出完了")


##--微小変形理論に基づく垂直方向変位--
delta_v2 = np.zeros(len(phi_0_data))
for i in range(len(phi_0_data)):
    delta_v2[i] = (P_load[i] * np.sin(alpha) * L**3)/(3 * ei)

##--弾塑性判定についての計算--
Mp = Zp*sigma_y
safety_ratio = np.zeros(len(phi_0_data))
Mx= np.zeros(len(phi_0_data))
P_safety_ratio = 0

for i in range(len(phi_0_data)):
    Mx[i] = P_load[i] * np.sin(alpha) *xc[i] + P_load[i] * np.cos(alpha) *yc[i]
    safety_ratio[i] = Mp / Mx[i]
    if safety_ratio[i] < 1 and safety_ratio[i-1] > 1:
        P_safety_ratio = P_load[i]


##--データ整理・出力
datasheet = np.ndarray((len(phi_0_data),8))
phi_0_data = np.loadtxt("phi_0.csv" ,delimiter=",", encoding="shift-jis")
for i in range(len(phi_0_data)):
    datasheet[i,0] = phi_0_data[i]
    datasheet[i,1] = p[i]
    datasheet[i,2] = m[i]
    datasheet[i,3] = P_load[i]
    datasheet[i,4] = xc[i]
    datasheet[i,5] = yc[i]
    datasheet[i,6] = delta_v2[i]
    datasheet[i,7] = safety_ratio[i]
#---データ出力のフォーマットは
'''|phi_0_data(deg)|p|m|P_load[N]|xc|yc|small_deflection|安全率|'''
def fig():
    ##---荷重による変位量を3種類プロット
    fig, ax1 = plt.subplots()
    x = datasheet[:,3]
    y1 = datasheet[:,4]
    y2 = datasheet[:,5]
    y3 = datasheet[:,6]

    ax1.plot(x, y1, 'o', color='r', markersize=4)
    ax1.plot(x, y2, '*', color='b', markersize=4)
    #ax1.plot(x, y3, '+', color='k', markersize=4)
    ax1.axvline(x = P_safety_ratio)

    #ax1.legend(["H_Displacement","V_deflection","V_small_deflection","弾塑性の境界"],prop = {"family" : "MS Gothic"})
    ax1.legend(["x_position","y_position","弾塑性の境界"],prop = {"family" : "MS Gothic"})


    ax1.set_title("荷重による梁先端の推移 | L = " + str(L) + "[m]",fontname = 'MS Gothic')
    ax1.set_xlabel("Load[N] | P")
    ax1.set_ylabel("Displacement[m]")
    plt.savefig("荷重による変位量の推移")


    ##---安全率をプロット
    fig, ax2 = plt.subplots()
    x = datasheet[:,3]
    y1 = datasheet[:,7]
    ax2.axhline(y=1)
    ax2.plot(x, y1, 'o', color='r', markersize=4)

    ax2.set_title("荷重による安全率の推移",fontname = 'MS Gothic')
    ax2.set_title("安全率の推移",fontname = 'MS Gothic')
    ax2.set_xlabel("Load[N] | P")
    ax2.set_ylabel("safety_ratio[-]")
    plt.savefig("荷重による安全率の推移")



    ##---自由端位置
    fig, ax3 = plt.subplots()
    x = datasheet[:,4]
    y1 = -datasheet[:,5]
    ax3.plot(x, y1, 'o', color='r', markersize=4)
    ax3.set_xlim(0,L)
    ax3.set_title("自由端位置",fontname = 'MS Gothic')
    ax3.set_xlabel("x [m]")
    ax3.set_ylabel("y [m]")
    plt.savefig("自由端位置")



    print("進捗度100%|計算完了")
    print("画像のプロットを行います")

datasheet_pd = pd.DataFrame(datasheet,columns=["phi0","p","m","P[N]","xc","yc","small_delta_v","safety_ratio"])
datasheet_pd.to_csv('flexber_ideal_diagonal_cluclate_data.csv')
fig()

plt.show()
