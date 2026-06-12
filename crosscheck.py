from mpmath import mp, mpf
mp.dps = 60

N = 600
t = [0]*(N+1); t[0] = 1
P2 = [0]*(N+1)
for n in range(1, N+1):
    m = n-1
    P2[m] = sum(t[i]*t[m-i] for i in range(m+1))
    c3   = sum(t[d]*P2[m-d] for d in range(m+1))
    cTT2 = sum(t[j//2]*t[m-j] for j in range(0, m+1, 2))
    cT3  = t[m//3] if m % 3 == 0 else 0
    t[n] = (c3 + 3*cTT2 + 2*cT3)//6

rho = mpf('0.355181742314377392882244473647632636708747')

def extrap(fn, ns):
    us = [mpf(1)/n for n in ns]
    tab = [fn(n) for n in ns]
    for k in range(1, len(ns)):
        for i in range(len(ns)-k):
            tab[i] = (us[i]*tab[i+1] - us[i+k]*tab[i])/(us[i]-us[i+k])
    return tab[0]

CR = extrap(lambda n: mpf(t[n])*rho**n*mpf(n)**mpf('1.5'), list(range(240, 601, 30)))
print("my rooted constant  C_R =", mp.nstr(CR, 30))
print("Kotesovec (OEIS 2015)   = 0.517875906458893536993162356993")
