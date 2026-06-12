import numpy as np

N = 620
# ---------- exact series: t = A000598, A = A000602 (a(0..620)) ----------
t = [0]*(N+1); t[0] = 1
P2 = [0]*(N+1)
for n in range(1, N+1):
    m = n-1
    P2[m] = sum(t[i]*t[m-i] for i in range(m+1))
    c3   = sum(t[dd]*P2[m-dd] for dd in range(m+1))
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
P4   = mul(P2, P2); P2T2 = mul(P2, T2); T2sq = mul(T2, T2); TT3 = mul(t, T3)

I = [0]*(N+1)
for n in range(1, N+1):
    znum = P4[n-1] + 6*P2T2[n-1] + 3*T2sq[n-1] + 8*TT3[n-1] + 6*T4[n-1]
    fnum = P2[n] - 2*t[n] - (t[n//2] if n % 2 == 0 else 0)
    I[n] = znum//24 - fnum//2

A = [1] + I[1:]          # A[n] = a(n), n = 0..620, a(0)=1
assert A[:21] == [1,1,1,1,2,3,5,9,18,35,75,159,355,802,1858,4347,10359,24894,60523,148284,366319]
print("exact data OK: a(0..620), a(620) has", len(str(A[620])), "digits")

# ---------- primes < 2^26 (so np.convolve of 621 terms cannot overflow int64) ----------
def isprime(m):
    if m < 2: return False
    d = 2
    while d*d <= m:
        if m % d == 0: return False
        d += 1
    return True
p1 = 2**26 - 5
while not isprime(p1): p1 -= 2
p2 = p1 - 2
while not isprime(p2): p2 -= 2
print("primes:", p1, p2)

# ---------- GF(p) rank ----------
def rankp(Mx, p):
    Mx = np.array(Mx, dtype=np.int64) % p
    R, Cn = Mx.shape; r = 0
    for c in range(Cn):
        nz = np.nonzero(Mx[r:, c])[0]
        if nz.size == 0: continue
        piv = r + nz[0]
        if piv != r: Mx[[r, piv]] = Mx[[piv, r]]
        inv = pow(int(Mx[r, c]), p-2, p)
        Mx[r] = (Mx[r]*inv) % p
        below = Mx[r+1:]
        f = below[:, c]
        sel = np.nonzero(f)[0]
        if sel.size:
            below[sel] = (below[sel] - f[sel, None]*Mx[r][None, :]) % p
        r += 1
        if r == R: break
    return r

# ---------- ROUTE 1: holonomic (P-recursive) guessing, inhomogeneous included ----------
def holo_test(seq, r, d, p):
    # sum_{j=0}^{r} p_j(n) a(n+j) = q(n), deg p_j, q <= d  (homog. + inhom. block)
    S = (r+1)*(d+1) + (d+1)
    M = min(len(seq)-1-r, S+25)
    if M < S+8: return None
    aa = np.array([x % p for x in seq], dtype=np.int64)
    ns = np.arange(1, M+1, dtype=np.int64)
    cols = []
    nk = np.ones(M, dtype=np.int64)
    for k in range(d+1):
        for j in range(r+1):
            cols.append((aa[1+j:1+j+M]*nk) % p)
        cols.append(nk.copy())          # inhomogeneous q(n) block
        nk = (nk*ns) % p
    Mat = np.stack(cols, axis=1)
    return S, rankp(Mat, p)

# validation decoy: Catalan (order 1, degree 1 recurrence) must be FOUND
cat = [1]
for n in range(620): cat.append(cat[-1]*2*(2*n+1)//(n+2))
S, rk = holo_test(cat, 1, 1, p1)
assert rk < S, "guesser failed to find Catalan recurrence!"
print(f"decoy Catalan: holo(1,1) deficiency {S-rk} > 0  -> guesser works")

pairs = [(2,150),(3,120),(4,100),(5,85),(6,72),(8,55),(10,46),(12,38),(15,31),
         (18,26),(22,21),(26,18),(31,15),(38,12),(46,10),(56,8),(70,6),(84,5),
         (96,4),(110,3),(145,2),(200,1),(304,0)]
print("\nROUTE 1: P-recurrence search  sum p_j(n) a(n+j) = q(n)")
allfull = True
for (r, d) in pairs:
    res = holo_test(A, r, d, p1)
    if res is None: print(f"  (r={r},d={d}) skipped"); continue
    S, rk = res
    full = (rk == S); allfull &= full
    print(f"  order<= {r:3d} deg<= {d:3d}: unknowns {S:3d}, rank {rk:3d}  ->  {'NONE' if full else 'FOUND?!'}")
for (r, d) in [(2,150),(10,46),(304,0)]:        # second-prime sanity
    S, rk = holo_test(A, r, d, p2)
    assert rk == S
print("  re-checked (2,150),(10,46),(304,0) mod second prime: still full rank")
print("  => NO P-recurrence anywhere in the staircase region; verdict:", "NONE EXISTS (certified)" if allfull else "?")

# ---------- ROUTE 2: algebraic g.f. and order-1 ADE ----------
def series_pows(base, Dmax, p):
    out = [np.zeros(len(base), dtype=np.int64)]
    out[0][0] = 1
    for j in range(Dmax):
        out.append(np.convolve(out[-1], base)[:len(base)] % p)
    return out

yp = np.array([x % p1 for x in A], dtype=np.int64)
print("\nROUTE 2a: algebraic equation  P(x,y)=0")
# decoy: Catalan gf satisfies 1 - y + x y^2 = 0
cp = np.array([x % p1 for x in cat], dtype=np.int64)
Y = series_pows(cp, 2, p1)
cols = []
for j in range(3):
    for i in range(2):
        v = np.zeros(621, dtype=np.int64); v[i:] = Y[j][:621-i]; cols.append(v)
assert rankp(np.stack(cols,axis=1), p1) < 6
print("  decoy Catalan: algebraic relation found -> tester works")

apairs = [(2,190),(3,145),(4,115),(5,95),(6,82),(8,64),(10,52),(12,44),
          (14,38),(17,32),(20,27),(24,22),(28,19)]
allfull2 = True
for (Dy, Dx) in apairs:
    Y = series_pows(yp, Dy, p1)
    cols = []
    for j in range(Dy+1):
        for i in range(Dx+1):
            v = np.zeros(621, dtype=np.int64); v[i:] = Y[j][:621-i]; cols.append(v)
    S = len(cols); rk = rankp(np.stack(cols, axis=1), p1)
    full = (rk == S); allfull2 &= full
    print(f"  deg_y<= {Dy:2d} deg_x<= {Dx:3d}: unknowns {S:3d}, rank {rk:3d} -> {'NONE' if full else 'FOUND?!'}")
print("  => g.f. satisfies NO algebraic equation in that region" if allfull2 else "  ?!")

print("\nROUTE 2b: order-1 ADE  P(x, y, y') = 0, deg<=7 in each")
y1 = np.zeros(621, dtype=np.int64)
y1[:620] = (np.arange(1, 621, dtype=np.int64)*yp[1:]) % p1
Yj = series_pows(yp, 7, p1); Zk = series_pows(y1, 7, p1)
cols = []
for j in range(8):
    for k in range(8):
        W = np.convolve(Yj[j], Zk[k])[:620] % p1
        for i in range(8):
            v = np.zeros(620, dtype=np.int64); v[i:] = W[:620-i]; cols.append(v)
S = len(cols); rk = rankp(np.stack(cols, axis=1), p1)
print(f"  unknowns {S}, rank {rk} -> {'NO order-1 ADE with deg<=7' if rk==S else 'FOUND?!'}")

# ---------- ROUTE 3: automatic / digit formulas, periodicity mod m ----------
def kernel_growth(seq, base, mod, maxdepth, min_ov=20):
    seq = [int(x) % mod for x in seq]
    states = []   # list of lists
    def find(v):
        for idx, s in enumerate(states):
            L = min(len(s), len(v))
            if L >= min_ov and s[:L] == v[:L]: return idx
        return None
    frontier = [(0, 0)]
    states.append(seq)
    counts = [1]
    for depth in range(1, maxdepth+1):
        newf = []
        for (e, r) in frontier:
            for tdig in range(base):
                r2 = r + tdig*base**e; e2 = e + 1
                v = [seq[(base**e2)*n + r2] for n in range((len(seq)-1-r2)//base**e2 + 1)]
                if len(v) < min_ov: continue
                if find(v) is None:
                    states.append(v); newf.append((e2, r2))
        frontier = newf
        counts.append(len(states))
    return counts

tm = [bin(n).count('1') % 2 for n in range(621)]   # Thue-Morse decoy
ctm = kernel_growth(tm, 2, 2, 4)
assert ctm[-1] == 2, ctm
print("\nROUTE 3: decoy Thue-Morse 2-kernel closes at 2 states -> tester works")
c22 = kernel_growth(A, 2, 2, 4)
c33 = kernel_growth(A, 3, 3, 3)
c23 = kernel_growth(A, 3, 2, 3)
print("  a(n) mod 2, base 2: kernel sizes by depth:", c22)
print("  a(n) mod 3, base 3: kernel sizes by depth:", c33)
print("  a(n) mod 2, base 3: kernel sizes by depth:", c23)

per = []
for m in range(2, 8):
    found = None
    for P in range(1, 201):
        if all(A[n] % m == A[n+P] % m for n in range(250, 621-P)):
            found = P; break
    per.append((m, found))
print("  eventual periodicity mod m (period<=200, n>=250):", per)

# ---------- ROUTE 4: labeled closed form (depth-2 finite sum) ----------
from math import factorial as f
def labeled_deg4(n):
    if n == 1: return 1
    if n == 2: return 1
    tot = 0
    for d4 in range(0, (n-2)//3 + 1):
        for c3 in range(0, (n-2-3*d4)//2 + 1):
            b2 = n - 2 - 2*c3 - 3*d4
            a1 = n - b2 - c3 - d4
            if a1 < 0: continue
            tot += f(n)//(f(a1)*f(b2)*f(c3)*f(d4)) * f(n-2)//(2**c3 * 6**d4)
    return tot
chk = [labeled_deg4(n) for n in range(1, 9)]
print("\nROUTE 4: labeled deg<=4 trees, n=1..8:", chk)
# brute force verification for n<=8 via Pruefer sequences
def brute(n):
    if n == 1: return 1
    if n == 2: return 1
    cnt = 0
    import itertools
    for seq in itertools.product(range(n), repeat=n-2):
        deg = [1]*n
        for v in seq: deg[v] += 1
        if max(deg) <= 4: cnt += 1
    return cnt
bf = [brute(n) for n in range(1, 9)]
print("  brute force          n=1..8:", bf)
assert chk == bf, "labeled formula wrong!"
print("  labeled closed form VERIFIED (depth-2 finite sum)")
