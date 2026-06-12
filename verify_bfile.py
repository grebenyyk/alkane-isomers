N = 620
t = [0]*(N+1); t[0] = 1
P2 = [0]*(N+1)
for n in range(1, N+1):
    m = n-1
    P2[m] = sum(t[i]*t[m-i] for i in range(m+1))
    c3   = sum(t[d]*P2[m-d] for d in range(m+1))
    cTT2 = sum(t[j//2]*t[m-j] for j in range(0, m+1, 2))
    cT3  = t[m//3] if m % 3 == 0 else 0
    t[n] = (c3 + 3*cTT2 + 2*cT3)//6
P2[N] = sum(t[i]*t[N-i] for i in range(N+1))

def mul(Aa, Bb):
    C = [0]*(N+1)
    for i, a in enumerate(Aa):
        if a:
            for j in range(N-i+1):
                if Bb[j]:
                    C[i+j] += a*Bb[j]
    return C

T2 = [0]*(N+1); T3 = [0]*(N+1); T4 = [0]*(N+1)
for k in range(N+1):
    if 2*k <= N: T2[2*k] = t[k]
    if 3*k <= N: T3[3*k] = t[k]
    if 4*k <= N: T4[4*k] = t[k]
P4 = mul(P2,P2); P2T2 = mul(P2,T2); T2sq = mul(T2,T2); TT3 = mul(t,T3)
I = [0]*(N+1)
for n in range(1, N+1):
    znum = P4[n-1] + 6*P2T2[n-1] + 3*T2sq[n-1] + 8*TT3[n-1] + 6*T4[n-1]
    fnum = P2[n] - 2*t[n] - (t[n//2] if n % 2 == 0 else 0)
    I[n] = znum//24 - fnum//2
A = [1] + I[1:]

import os
_here = os.path.dirname(os.path.abspath(__file__))
bf = {}
for line in open(os.path.join(_here, 'b000602.txt')):
    line = line.strip()
    if not line or line.startswith('#'): continue
    n_, v_ = line.split()
    bf[int(n_)] = int(v_)

bad = [n for n in range(0, 621) if bf.get(n) != A[n]]
print("b-file entries:", len(bf), " | mismatches in 0..620:", bad if bad else "NONE — all 621 values agree")
