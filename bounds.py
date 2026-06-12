import os
_here = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(_here, 'asym.py')).read())   # reuses I[], rho, alpha, C, c1, c2, c3

print("\n---- certified bound scan over 25 <= n <= 600 ----")
m1 = m2 = m3 = 0
for n in range(25, 601):
    base = C * alpha**n / mpf(n)**mpf('2.5')
    e1 = abs(base/mpf(I[n]) - 1)
    e2 = abs(base*(1 + c1/n)/mpf(I[n]) - 1)
    e3 = abs(base*(1 + c1/n + c2/n**2)/mpf(I[n]) - 1)
    m1 = max(m1, e1*n)
    m2 = max(m2, e2*n*n)
    m3 = max(m3, e3*n**3)
print("max n  * |err(1-term)| =", float(m1))
print("max n^2* |err(2-term)| =", float(m2))
print("max n^3* |err(3-term)| =", float(m3))

# rounding demo
print("\n---- where rounding the 3-term formula still reproduces I(n) exactly ----")
ok = 0
for n in range(1, 81):
    f3 = C * alpha**n / mpf(n)**mpf('2.5') * (1 + c1/n + c2/n**2)
    if int(mp.nint(f3)) == I[n]:
        ok = n if ok == n-1 or ok == 0 else ok
print("3-term formula, n, round(F), I(n) for n=10..16:")
for n in range(10, 17):
    f3 = C * alpha**n / mpf(n)**mpf('2.5') * (1 + c1/n + c2/n**2)
    print(n, int(mp.nint(f3)), I[n])
