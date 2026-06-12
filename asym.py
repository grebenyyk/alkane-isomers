import sys

N = 620

# ---------- exact integer series ----------
# T(x): g.f. of A000598 (alkyl radicals; rooted trees, <=3 children per node), t[0]=1
# T = 1 + (x/6)(T(x)^3 + 3 T(x) T(x^2) + 2 T(x^3))
t = [0]*(N+1); t[0] = 1
P2 = [0]*(N+1)                       # running coefficients of T^2
for n in range(1, N+1):
    m = n-1
    P2[m] = sum(t[i]*t[m-i] for i in range(m+1))
    c3   = sum(t[dd]*P2[m-dd] for dd in range(m+1))       # [x^m] T^3
    cTT2 = sum(t[j//2]*t[m-j] for j in range(0, m+1, 2))  # [x^m] T(x)T(x^2)
    cT3  = t[m//3] if m % 3 == 0 else 0                   # [x^m] T(x^3)
    s = c3 + 3*cTT2 + 2*cT3
    assert s % 6 == 0, n
    t[n] = s // 6
P2[N] = sum(t[i]*t[N-i] for i in range(N+1))

assert t[:15] == [1,1,1,2,4,8,17,39,89,211,507,1238,3057,7639,19241], "A000598 mismatch"
print("A000598 OK")

def mul(A, B):
    C = [0]*(N+1)
    for i, a in enumerate(A):
        if a:
            for j in range(N-i+1):
                if B[j]:
                    C[i+j] += a*B[j]
    return C

T2 = [0]*(N+1); T3 = [0]*(N+1); T4 = [0]*(N+1)
for k in range(N+1):
    if 2*k <= N: T2[2*k] = t[k]
    if 3*k <= N: T3[3*k] = t[k]
    if 4*k <= N: T4[4*k] = t[k]

P4   = mul(P2, P2)     # T^4
P2T2 = mul(P2, T2)     # T^2 * T(x^2)
T2sq = mul(T2, T2)     # T(x^2)^2
TT3  = mul(t,  T3)     # T * T(x^3)

# Otter/dissymmetry:  I(x) = x*Z(S4;T) - ((T-1)^2 - (T(x^2)-1))/2
I = [0]*(N+1)
for n in range(1, N+1):
    znum = P4[n-1] + 6*P2T2[n-1] + 3*T2sq[n-1] + 8*TT3[n-1] + 6*T4[n-1]
    assert znum % 24 == 0, n
    fnum = P2[n] - 2*t[n] - (t[n//2] if n % 2 == 0 else 0)
    assert fnum % 2 == 0, n
    I[n] = znum//24 - fnum//2

known = [1,1,1,2,3,5,9,18,35,75,159,355,802,1858,4347,10359,24894,60523,148284,366319]
print("I(1..20) match A000602:", I[1:21] == known)
print("digits of I(600):", len(str(I[600])))

# ---------- high-precision singularity constants ----------
from mpmath import mp, mpf
mp.dps = 60
K = 350

def Tval(x):
    s = mpf(0); xk = mpf(1)
    for k in range(K+1):
        if t[k]: s += t[k]*xk
        xk *= x
    return s

def Tder(x):
    s = mpf(0); xk = mpf(1)   # xk = x^(k-1)
    for k in range(1, K+1):
        s += k*t[k]*xk
        xk *= x
    return s

# Solve: a = 1 + (rho/6)(a^3 + 3a*T(rho^2) + 2T(rho^3)),  (rho/2)(a^2 + T(rho^2)) = 1
rho = mpf('0.3551817'); a = mpf('2.117')
for it in range(60):
    x2 = rho*rho; x3 = x2*rho
    E2, E3 = Tval(x2), Tval(x3)
    dE2 = Tder(x2)*2*rho
    dE3 = Tder(x3)*3*x2
    F1 = 1 + rho/6*(a**3 + 3*a*E2 + 2*E3) - a
    F2 = rho/2*(a*a + E2) - 1
    J11 = (a**3 + 3*a*E2 + 2*E3)/6 + rho/6*(3*a*dE2 + 2*dE3)
    J12 = rho/2*(a*a + E2) - 1
    J21 = (a*a + E2)/2 + rho/2*dE2
    J22 = rho*a
    det = J11*J22 - J12*J21
    drho = (F1*J22 - F2*J12)/det
    da   = (J11*F2 - J21*F1)/det
    rho -= drho; a -= da
    if abs(drho) < mpf(10)**-55 and abs(da) < mpf(10)**-55:
        break

alpha = 1/rho
print("iters:", it)
print("rho   =", mp.nstr(rho, 42))
print("alpha =", mp.nstr(alpha, 42))
print("T(rho)=", mp.nstr(a, 42))

# ---------- Richardson/Neville extrapolation in 1/n ----------
def extrap(fn, ns):
    us = [mpf(1)/n for n in ns]
    tab = [fn(n) for n in ns]
    m = len(ns)
    for k in range(1, m):
        for i in range(m-k):
            tab[i] = (us[i]*tab[i+1] - us[i+k]*tab[i])/(us[i]-us[i+k])
    return tab[0]

def s_of(n):
    return mpf(I[n]) * rho**n * mpf(n)**mpf('2.5')

ns1 = list(range(240, 601, 30))
ns2 = list(range(200, 601, 40))
C   = extrap(s_of, ns1)
Cb  = extrap(s_of, ns2)
print("C  =", mp.nstr(C, 30))
print("C' =", mp.nstr(Cb, 30))
print("|C-C'| =", mp.nstr(abs(C-Cb), 5))

c1  = extrap(lambda n: (s_of(n)/C - 1)*n, ns1)
c1b = extrap(lambda n: (s_of(n)/Cb - 1)*n, ns2)
print("c1 =", mp.nstr(c1, 20), "/", mp.nstr(c1b, 20))

c2  = extrap(lambda n: ((s_of(n)/C - 1)*n - c1)*n, ns1)
c2b = extrap(lambda n: ((s_of(n)/Cb - 1)*n - c1b)*n, ns2)
print("c2 =", mp.nstr(c2, 15), "/", mp.nstr(c2b, 15))

c3  = extrap(lambda n: (((s_of(n)/C - 1)*n - c1)*n - c2)*n, ns1)
print("c3 =", mp.nstr(c3, 10))

# ---------- error table ----------
print()
print("  n    relerr 1-term   relerr 2-term   relerr 3-term")
for n in [25, 50, 100, 200, 400, 600]:
    base = C * alpha**n / mpf(n)**mpf('2.5')
    f1 = base
    f2 = base*(1 + c1/n)
    f3 = base*(1 + c1/n + c2/n**2)
    e = lambda f: float(f/mpf(I[n]) - 1)
    print(f"{n:4d}   {e(f1):+.3e}      {e(f2):+.3e}      {e(f3):+.3e}")

print()
print("I(25) =", I[25])
print("I(50) =", I[50])
print("I(100)=", I[100])
