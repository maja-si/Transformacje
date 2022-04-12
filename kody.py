#metoda Hirvonena

import math
import numpy as np

#Hirvonen XYZ na fi, lam, h
def hirvonen (X, Y, Z, a, e2):
    r = math.sqrt(X**2 + Y**2)
    fi_n = math.atan(Z/(r*(1-e2)))
    eps = 0.000001/3600 *math.pi/180 # radiany
    fi = fi_n*2
    while np.abs(fi_n - fi) > eps:
        fi = fi_n
        N = a/np.sqrt(1-e2*np.sin(fi_n)**2)
        h = r/np.cos(fi_n) -N
        fi_n = math.atan(Z/(r*(1-e2*(N/(N + h)))))
    lam = math.atan(Y/X)
    h = r/np.cos(fi_n) -N
    return fi, lam, h


#radiany na stopnie minuty sekundy do 5 miejsc po przecinku
def s_m_s5(fi):
    fi = fi*180/math.pi # stopnie
    d =np.floor(fi)
    m = np.floor((fi - d)*60)
    s = round((fi -d -m/60)*3600, 5)
    print(d, 'st', m, 'min', s, 'sek')


#s, A, z(w radianach) do n, e, u
def s_A_z2neu(s, A, z):
   
    n = s*np.sin(z)*np.cos(A)
    e = s*np.sin(z)*np.sin(A)
    u = s*np.cos(z)
    return n, e, u

#n, e, u(w radianach) do X, Y, Z
def neu2dXYZ(n, e, u, F1, L1):
    
    R = np.array([[-np.sin(F1) * np.cos(L1), -np.sin(L1), np.cos(F1)*np.cos(L1)],
              [-np.sin(F1)*np.sin(L1), np.cos(L1), np.cos(F1)*np.sin(L1)],
              [np.cos(F1), 0, np.sin(F1)]])
    dx = np.linalg.inv(R.transpose()) @ np.array([n, e, u])
    return dx

#dx = neu2dXYZ(n, e, u, fi, lam)

#Xb = X + dx[0]
#Yb = Y + dx[1]
#Zb = Z + dx[2]


# Kivioji fi_A, lam_A, A_AB(w radianach), S_AB na fi_B, lam_B, A_BA
def Kivioj(fi, lam, A_AB, S_AB, e2, a):
    n = round(S_AB/1000)
    ds = S_AB/n #s/n
    for i in range(n):
        M = (a*(1-e2))/math.sqrt((1-e2*(math.sin(fi))**2)**3)
        N = a/math.sqrt(1-e2*(math.sin(fi))**2)
        
        dfi = (ds*math.cos(A_AB))/M
        dA = ((math.sin(A_AB) * math.tan(fi))/N)*ds
        
        fi_m = fi + dfi/2
        A_m = A_AB + dA/2
        M_m = (a*(1-e2))/(math.sqrt(1-e2*np.sin(fi_m)**2)**3)
        N_m = a/(math.sqrt(1-e2*np.sin(fi_m)**2))
        
        dfi_p = (ds*math.cos(A_m))/M_m
        dlam_p = (ds*math.sin(A_m))/(N_m*math.cos(fi_m))
        dA_p = ((math.sin(A_m) * math.tan(fi_m))/N_m)*ds
        fi = fi +dfi_p
        lam = lam + dlam_p
        A_AB = A_AB + dA_p
    fi_k = fi
    A_k = A_AB 
    lam_k = lam 
    
    if A_k < math.pi:
        A_k = A_k + math.pi
    elif A_k >= math.pi:
        A_k = A_k - math.pi
    return fi_k, lam_k, A_k

#Vincent fi_A, lam_A, fi_B, lam_B (w radianach) na S_AB, A_AB, A_BA
def Vincent(fi, lam, fi_k, lam_k, a, e2):
    b = a*np.sqrt(1-e2)
    f = 1-(b/a)
    d_lam = lam_k -lam
    U_A = np.arctan((1-f)*np.tan(fi))
    U_B = np.arctan((1-f)*np.tan(fi_k))
    eps = 0.000001/3600 *math.pi/180 # radiany
    L=d_lam
 
    while True:
        
        sin_sgm = np.sqrt((np.cos(U_B)*np.sin(L))**2 + (np.cos(U_A)*np.sin(U_B) - np.sin(U_A)*np.cos(U_B)*np.cos(L))**2)
        cos_sgm = np.sin(U_A)*np.sin(U_B) + np.cos(U_A)*np.cos(U_B)*np.cos(L)
        sgm = np.arctan(sin_sgm/cos_sgm)
        
       
        sin_alfa = (np.cos(U_A)*np.cos(U_B)*np.sin(L))/sin_sgm
        cos2_alfa = 1 - sin_alfa**2
       
        cos_2sgm_m = cos_sgm - (2*np.sin(U_A)*np.sin(U_B))/cos2_alfa
      
        C = (f/16)*cos2_alfa*(4 + f*(4 - 3*cos2_alfa))
        
        L1 = L
        L= d_lam + (1 - C)*f*sin_alfa*(sgm + C*sin_sgm*(cos_2sgm_m + C*cos_sgm*((-1)+2*(cos_2sgm_m)**2)))
        
        if np.abs(L1 - L)<eps:
            break
    
    u2 = ((a**2 - b**2)/b**2)*cos2_alfa
    A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320-175*u2)))
    B = (u2/1024)*(256 + u2*(-128 + u2*(74-(47*u2))))
    
    d_sgm = B*sin_sgm*(cos_2sgm_m + (1/4)*B*(cos_sgm*(-1+2*(cos_2sgm_m**2)) - (1/6)*B*cos_2sgm_m*(-3+4*(sin_sgm**2))*(-3+4*(cos_2sgm_m**2))))
   
    s_AB = b*A*(sgm - d_sgm)
    A_AB = np.arctan((np.cos(U_B)*np.sin(L))/(np.cos(U_A)*np.sin(U_B) - np.sin(U_A)*np.cos(U_B)*np.cos(L)))
    A_BA = np.arctan((np.cos(U_A)*np.sin(L))/(-np.sin(U_A)*np.cos(U_B) + np.cos(U_A)*np.sin(U_B)*np.cos(L)))+np.pi
    return s_AB, A_AB, A_BA

#stopnie minuty sekundy z radianów do 3 miejsc po przecinku
def s_m_s3(fi):
    fi = fi*180/math.pi # stopnie
    d =np.floor(fi)
    m = np.floor((fi - d)*60)
    s = round((fi -d -m/60)*3600, 3)
    print(d, 'st', m, 'min', s, 'sek')
  
    
  
    
#uklad 1992 fi, lam(w radianach) na x92, y92
def u1992(fi, lam, e2, m_0, a):
    N = a/(np.sqrt(1-e2 * np.sin(fi)**2))
    t = np.tan(fi)
    e_2 = e2/(1-e2)
    n2 = e_2 * (np.cos(fi))**2
    
    lam_0 = math.radians(19)
    l = lam - lam_0
    
    A0 = 1 - (e2/4) - ((3*(e2**2))/64) - ((5*(e2**3))/256)   
    A2 = (3/8) * (e2 + ((e2**2)/4) + ((15 * (e2**3))/128))
    A4 = (15/256) * (e2**2 + ((3*(e2**3))/4))
    A6 = (35 * (e2**3))/3072 
    

    sig = a * ((A0*fi) - (A2*np.sin(2*fi)) + (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
    
    x = sig + ((l**2)/2) * N *np.sin(fi) * np.cos(fi) * (1 + ((l**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
    y = l * (N*math.cos(fi)) * (1 + ((((l**2)/6) * (math.cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
    x92 = round(m_0*x - 5300000, 3)
    y92 = round(m_0*y + 500000, 3)
    return x92, y92 



#uklad 2000 fi, lam(w radianach) na x00, y00
def u2000(fi, lam, e2, m, a):
    N = a/math.sqrt(1-e2*math.sin(fi)**2)
    t = np.tan(fi)
    e_2 = e2/(1-e2)
    n2 = e_2 * (np.cos(fi))**2
    
    lam = math.degrees(lam)
    if lam>13.5 and lam <16.5:
        s = 5
        lam0 = 15
    elif lam>16.5 and lam <19.5:
        s = 6
        lam0 = 18
    elif lam>19.5 and lam <22.5:
        s = 7
        lam0 = 21
    elif lam>22.5 and lam <25.5:
        s = 8
        lam0 = 24
        
    lam = math.radians(lam)
    lam0 = math.radians(lam0)
    l = lam - lam0

    A0 = 1 - (e2/4) - ((3*(e2**2))/64) - ((5*(e2**3))/256)   
    A2 = (3/8) * (e2 + ((e2**2)/4) + ((15 * (e2**3))/128))
    A4 = (15/256) * (e2**2 + ((3*(e2**3))/4))
    A6 = (35 * (e2**3))/3072 
    

    sig = a * ((A0*fi) - (A2*np.sin(2*fi)) + (A4*np.sin(4*fi)) - (A6*np.sin(6*fi)))
    
    x = sig + ((l**2)/2) * N *np.sin(fi) * np.cos(fi) * (1 + ((l**2)/12) * ((math.cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((math.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
    y = l * (N*math.cos(fi)) * (1 + ((((l**2)/6) * (math.cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (math.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
    x00 = round(m * x, 3) 
    y00 = round(m * y + (s*1000000) + 500000, 3)
     
    return(x00, y00)

#z 2000/1992 na X_GK, Y_GK, lam_0, m
def odcechowanie(X, Y, m_92, m_00):
    if X < 1000000 and Y < 1000000:
        print('zamiana z 1992 na G-K')
        X_GK = (X + 5300000)/m_92
        Y_GK = (Y - 500000)/m_92
        lam_0 = math.radians(19)
        m = m_92
    elif X > 1000000 and Y > 1000000:
        if Y > 5000000 and Y < 6000000:
            s = 5
            lam_0 = math.radians(15)
        elif Y > 6000000 and Y < 7000000:
            s = 6
            lam_0 =  math.radians(18)
        elif Y > 7000000 and Y < 8000000:
            s = 7
            lam_0 =  math.radians(21)
        elif Y > 8000000 and Y < 9000000:
            s = 8
            lam_0 =  math.radians(24)
        print('zamiana z 2000 na G-K')
        X_GK = X/m_00
        Y_GK = (Y - (s * 1000000) - 500000)/m_00
        m = m_00
    return(X_GK, Y_GK, lam_0, m)


#z X_GK, Y_GK, lam_0, m na fi, lam
def geodezyjne(X_GK, Y_GK, e2, a, lam_0, m):
    A0 = 1 - (e2/4) - (3*(e2**2))/64 - (5*(e2**3))/256
    A2 = 3/8 * (e2 + ((e2**2)/4) + ((15*(e2**3))/128))
    A4 = 15/256 * ((e2**2) + (3*(e2**3))/4)
    A6 = (35*(e2**3))/3072
    eps = 0.000001/3600 *math.pi/180
    b1 = X_GK/(a * A0)
    while True:
        b0 = b1
        b1= (X_GK/(a * A0 )) + ((A2/A0) * np.sin(2*b0)) - ((A4/A0) * np.sin(4*b0)) + ((A6/A0) * np.sin(6*b0))
        b = b1
        if np.abs(b1 - b0) <= eps:
            break
        
    
    e_2 = e2/(1-e2)
    N = a/math.sqrt(1 - e2 * np.sin(b)**2)
    t = np.tan(b)
    n2 = e_2 * np.cos(b)**2
    
    fi = b - (t/2) * (((Y_GK/(N))**2) * (1 + n2) - (1/12) * ((Y_GK/( N))**4) * (5 + (3 * t**2) + (6*n2) - (6 * n2 * t**2) - (3 * n2**2) - (9 * t**2 * n2**2)) + (1/360) * ((Y_GK/(N))**6) * (61 + (90 * t**2) + (45 * t**4) + (107 * n2) - (162 * t**2 * n2) - (45 * (t**4) * n2)))
    l = (1/np.cos(b)) * ((Y_GK/(N)) - ((1/6) * (Y_GK/(N))**3 * (1 + 2 * t**2 + n2)) + ((1/120) * (Y_GK/( N))**5 * (5 + (28 * t**2) + (24 * t**4) + (6 * n2) + (8 * n2 * t**2))))
    lam = lam_0 + l
    return(fi, lam)

#geodezyjne fi, lam(w radianach), h na XYZ
def geod2XYZ(fi, lam, h, a, e2):
    N = a/math.sqrt(1-e2*math.sin(fi)**2)
    X = (N + h) * math.cos(fi) * math.cos(lam)
    Y = (N + h) * math.cos(fi) * math.sin(lam)
    Z = (N*(1-e2) + h) * math.sin(fi)
    return(X, Y, Z)


#M fi w radianach
def M(fi, a, e2):
    M = (a*(1-e2))/((np.sqrt(1 - e2 * (np.sin(fi)**2)))**3)
    return(M)

#N fi w radianach
def N(fi, a, e2):
    N = a/(np.sqrt(1 - e2 * (np.sin(fi)**2)))
    return(N)
