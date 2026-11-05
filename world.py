"""One kind of cell. Every body has four heritable genes: light, absorb, engulf, kin range, each in [0,1]; pass 9 adds kin range and an energy store.
Metabolism is Troost, Kooi & Kooijman 2005 (eqs 1-11, Table 3), individual-based, with a third route added.
No plant, animal or fungus type anywhere. No fitness function: a body divides when it has the mass.
Mass is closed: every unit is in C (inorganic), D (dead matter) or a body. Light is open."""
import os
import numpy as np
from numba import njit

L = 40                # lattice side; bodies carry continuous positions
DT = 0.1              # days per step
KC, KD, KE = 500.0, 2500.0, 2500.0   # KE: not in the paper, set equal to KD
ZC, ZN, ZCH = 10.0, 10.0, 10.0
JL = 5.0 / 25.0
YDV, KM, H = 1.3, 0.1, 0.1
JM = 2.6              # all three routes' maximum flux
Y0 = 1.1
TOTAL = 4000.0        # mass per lattice cell, the paper's is 800: see plan.md
V0 = 300.0            # mass of a newborn; a body divides at 2*V0
VMIN = 0.5 * V0       # below this a body dies; not in the paper, which has no bodies
MU, SIGMA = 0.02, 0.02  # per gene per birth
DIFF = 0.2            # field diffusion per step
SPEED, TURN = 0.3, 0.6  # lattice units per step, radians sd per step
NMAX = 60000
NOHUNT = os.environ.get("NOHUNT", "0") == "1"   # paper arm N: QE is 0 in founders and every child, so nobody can engulf
NF = 17               # float columns per body below
X, Y_, HD, V, QA, QH, QE, COST, CL, QK, EN, QS, DM, HS, PT = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14   # pass 13: QS dormancy gene, DM dormant flag
CMAX = 2000000         # clade table size; at the cap no new clade is made
CDIST = 0.25          # a body whose genes sit farther than this from its clade's founding genes founds a new clade
MAINT = 5.5632 * float(os.environ.get("MAINT_MULT", "1"))   # pass 9 cost x multiple (plan.md ladder)
EMAX = 20 * MAINT     # energy store cap; a newborn gets half its parent's
NCA = 10              # accumulator columns per clade (9 = kin-spared)
LITE = 0.40            # pass 21: light at or above which a column wakes its sleepers and admits no sleep
DCOST, DABS = 0.2, 0.5   # pass 13: dormant maintenance and absorption, as shares of awake
WAKE = 0.2            # a dormant body wakes at (QS + WAKE) x EMAX, capped at EMAX
PMAX = 2000000        # engulfing pairs held between censuses

@njit(cache=True)
def seed_numba(s):
    np.random.seed(s)

@njit(cache=True)
def cost(qa, qh, qe, c, z, yint):
    """Eq 11, extended symmetrically to three genes; yint multiplies every pair's product."""
    return (Y0 + c * (1 - (1 - qa) ** z) + c * (1 - (1 - qh) ** z) + c * (1 - (1 - qe) ** z)
            + yint * (qa * qh + qa * qe + qh * qe))

@njit(cache=True)
def f_a(xc):
    """Eq 9 with X_N >> K_N (DIN not a state here), so f_N = 1."""
    fc = xc / (KC + xc)
    a = ZC * fc + 1e-12
    fch = (1 + 1 / ZC) / (1 + 1 / a + 1 / JL - 1 / (a + JL))
    a2 = ZN
    b2 = ZCH * fch + 1e-12
    return (1 + 1 / ZN + 1 / ZCH - 1 / (ZN + ZCH)) / (1 + 1 / a2 + 1 / b2 - 1 / (a2 + b2))

@njit(cache=True)
def mutate(q, mu):
    if np.random.random() < mu:
        q = q + np.random.normal() * SIGMA
        if q < 0.0: q = 0.0
        if q > 1.0: q = 1.0
    return q

@njit(cache=True)
def step(B, n, C, D, c, z, yint, led, flash, light, patch, mu, CT, cacc, cs, pairs, pc, stepno):
    """Advance one step. led: mass by route [auto, absorb, engulf], deaths [hazard, starve, engulfed], births.
    light: global multiplier on the light route's flux (the sky); patch[a,b]: static per-cell factor on it (pass 11). mu: mutation chance per gene per birth.
    led[7:10], led[10:13]: income by route of bodies with light gene under 0.5 / at or above 0.5.
    CT[k] = founding genes (3), parent, born step; cs[0] = clades made. cacc[k] = income by route (3), births, hunger deaths,
    hazard deaths, engulfed deaths, engulfings made. pairs/pc: (victim clade, killer clade) per engulfing since the last census.
    flash: positions of this step's engulfings (film only), flash[0,0] holds the count."""
    fA = np.empty((L, L)); fH = np.empty((L, L)); M = np.zeros((L, L))
    for a in range(L):
        for b in range(L):
            fA[a, b] = f_a(C[a, b]); fH[a, b] = D[a, b] / (KD + D[a, b])
    head = -np.ones((L, L), np.int64); nxt = -np.ones(NMAX, np.int64); cnt = np.zeros((L, L), np.int64)
    alive = np.ones(NMAX, np.bool_)
    for i in range(n):
        a = int(B[i, X]); b = int(B[i, Y_])
        nxt[i] = head[a, b]; head[a, b] = i; cnt[a, b] += 1; M[a, b] += B[i, V]
    nf = 0
    n0 = n
    for i in range(n0):
        if not alive[i]: continue
        a = int(B[i, X]); b = int(B[i, Y_])
        v = B[i, V]; y = B[i, COST]
        # light: eq 6, 2
        grp = 7 if B[i, QA] < 0.5 else 10
        dm = B[i, DM] > 0.5
        ab = DABS if dm else 1.0
        mc = (DCOST + 0.3 * B[i, HS]) if dm else 1.0   # pass 15: sleeping safely costs up to half an awake upkeep
        g = ab * v * B[i, QA] * JM * fA[a, b] * light[a] * patch[a, b] / y * DT
        if g > C[a, b]: g = C[a, b]
        cl = int(B[i, CL])
        C[a, b] -= g; led[0] += g; led[grp] += g; vnew = v + g; cacc[cl, 0] += g
        # absorb: eq 7, 3, 4
        du = ab * v * B[i, QH] * JM * fH[a, b] * YDV * DT
        if du > D[a, b]: du = D[a, b]
        g = du / (YDV * y)
        D[a, b] -= du; C[a, b] += du - g; led[1] += g; led[grp + 1] += g; vnew += g; cacc[cl, 1] += g
        # engulf: same form as absorb with living neighbours as the food, taken whole
        qe = B[i, QE]
        others = M[a, b] - v
        if (not dm) and qe > 0.0 and cnt[a, b] > 1 and others > 0.0:
            want = v * qe * JM * others / (KE + others) * YDV * DT
            k = np.random.randint(cnt[a, b]); t = head[a, b]
            for _ in range(k): t = nxt[t]
            if t != i and t >= 0 and t < n0 and alive[t] and np.random.random() < want / B[t, V]:
                dk = ((B[t, QA] - B[i, QA]) ** 2 + (B[t, QH] - B[i, QH]) ** 2 + (B[t, QE] - B[i, QE]) ** 2) ** 0.5
                if dk < B[i, QK]:
                    cacc[cl, 9] += 1   # refused: the target's genes sit within the eater's kin range; the meal is lost
                else:
                    m = B[t, V]
                    g = m / (YDV * y)
                    C[a, b] += m - g; led[2] += g; led[grp + 2] += g; vnew += g; cacc[cl, 2] += g
                    alive[t] = False; M[a, b] -= m; led[4] += 1
                    cv = int(B[t, CL]); cacc[cv, 6] += 1; cacc[cl, 8] += 1
                    if pc[0] < PMAX:
                        pairs[pc[0], 0] = cv; pairs[pc[0], 1] = cl; pc[0] += 1
                    else:
                        pc[1] += 1
                    if nf < flash.shape[0] - 1:
                        nf += 1; flash[nf, 0] = B[i, X]; flash[nf, 1] = B[i, Y_]
        # maintenance: eq 1
        en = B[i, EN] + (vnew - v) - MAINT * mc   # energy: this step's income (mass gained) less the per-step cost
        if en > EMAX: en = EMAX
        mm = mc * vnew * KM * DT; vnew -= mm; C[a, b] += mm
        if en <= 0.0:
            D[a, b] += vnew; alive[i] = False; M[a, b] -= v; led[5] += 1; cacc[cl, 4] += 1; continue
        B[i, EN] = en
        # death: eq 1c, as a hazard per body
        if np.random.random() < H * DT * ((1.0 - 0.9999 * B[i, HS]) if dm else 1.0):   # pass 15: shield cuts a sleeper's hazard
            D[a, b] += vnew; alive[i] = False; M[a, b] -= v; led[3] += 1; cacc[cl, 5] += 1; continue
        if vnew < VMIN:
            D[a, b] += vnew; alive[i] = False; M[a, b] -= v; led[5] += 1; cacc[cl, 4] += 1; continue
        # dormancy: under QS x EMAX a body goes dormant; it wakes at min(QS + WAKE, 1) x EMAX
        qs = B[i, QS]
        lit = light[a] >= LITE   # pass 21: sleep is a thing of the dark
        if dm:
            if lit or en >= min(qs + WAKE, 1.0) * EMAX: B[i, DM] = 0.0; dm = False; led[13 + a // 10] += 1.0   # pass 21: a lit column wakes its sleeper whatever its store
        elif (not lit) and en < qs * EMAX:
            B[i, DM] = 1.0; dm = True; led[17 + a // 10] += 1.0   # pass 21: no sleep in a lit column
        # division
        while (not dm) and vnew >= 2 * V0 and n < NMAX:
            vnew *= 0.5
            B[n, X] = B[i, X]; B[n, Y_] = B[i, Y_]; B[n, HD] = np.random.random() * 6.2831853
            B[n, V] = vnew
            B[n, EN] = 0.5 * en; en *= 0.5
            qa2 = mutate(B[i, QA], mu); qh2 = mutate(B[i, QH], mu); qe2 = 0.0 if NOHUNT else mutate(B[i, QE], mu)
            B[n, QK] = mutate(B[i, QK], mu); B[n, QS] = mutate(B[i, QS], mu); B[n, HS] = mutate(B[i, HS], mu); B[n, DM] = 0.0; B[n, PT] = mutate(B[i, PT], mu)   # pass 23: heritable pull toward light
            B[n, QA] = qa2; B[n, QH] = qh2; B[n, QE] = qe2; B[n, COST] = cost(qa2, qh2, qe2, c, z, yint)
            d2 = (qa2 - CT[cl, 0]) ** 2 + (qh2 - CT[cl, 1]) ** 2 + (qe2 - CT[cl, 2]) ** 2
            if d2 > CDIST * CDIST and cs[0] < CMAX:
                k = cs[0]; cs[0] += 1
                CT[k, 0] = qa2; CT[k, 1] = qh2; CT[k, 2] = qe2; CT[k, 3] = cl; CT[k, 4] = stepno
                B[n, CL] = k; cacc[k, 3] += 1
            else:
                B[n, CL] = cl; cacc[cl, 3] += 1
            alive[n] = True; led[6] += 1
            n += 1
        B[i, V] = vnew; B[i, EN] = en
    # moving and compaction
    w = 0
    for i in range(n):
        if not alive[i]: continue
        if B[i, DM] > 0.5:
            if w != i: B[w] = B[i]
            w += 1
            continue
        h = B[i, HD] + np.random.normal() * TURN
        pt = B[i, PT]
        if pt > 0.0 and np.random.random() < pt:   # pass 23: face the brighter neighbouring column
            a = int(B[i, X]); lr = light[(a + 1) % L]; ll = light[(a - 1) % L]
            if lr > ll: h = 0.0
            elif ll > lr: h = 3.1415927
        x = (B[i, X] + SPEED * np.cos(h)) % L; yy = (B[i, Y_] + SPEED * np.sin(h)) % L
        if x >= L: x = 0.0
        if yy >= L: yy = 0.0
        B[i, X] = x; B[i, Y_] = yy; B[i, HD] = h
        if w != i: B[w] = B[i]
        w += 1
    # diffusion on the torus
    for F in (C, D):
        G = F.copy()
        for a in range(L):
            an = (a + 1) % L; ap = (a - 1) % L
            for b in range(L):
                bn = (b + 1) % L; bp = (b - 1) % L
                F[a, b] = G[a, b] + DIFF * (G[an, b] + G[ap, b] + G[a, bn] + G[a, bp] - 4 * G[a, b])
    flash[0, 0] = nf
    return w

def start(seed, c, z, yint, founders=(0.5, 0.1, 0.05), n0=None):  # kin range 0, energy EMAX/2
    """Clones of one founder: the paper's starting mixotroph, with a little engulf so the gene is not at its wall."""
    rng = np.random.default_rng(seed); seed_numba(seed)
    n0 = n0 or int(0.25 * TOTAL * L * L / V0)
    B = np.zeros((NMAX, NF))
    B[:n0, X] = rng.random(n0) * L; B[:n0, Y_] = rng.random(n0) * L; B[:n0, HD] = rng.random(n0) * 6.2831853
    B[:n0, V] = V0 * (1 + rng.random(n0))
    B[:n0, QA], B[:n0, QH], B[:n0, QE] = founders
    if NOHUNT: B[:n0, QE] = 0.0
    B[:n0, COST] = float(cost(founders[0], founders[1], founders[2], c, z, yint))
    B[:n0, PT] = 0.5; B[:n0, CL] = 0.0; B[:n0, QK] = 0.0; B[:n0, EN] = 0.5 * EMAX
    rest = TOTAL - B[:n0, V].sum() / (L * L)
    C = np.full((L, L), rest / 2); D = np.full((L, L), rest / 2)
    return B, n0, C, D

def total_mass(B, n, C, D):
    return B[:n, V].sum() + C.sum() + D.sum()

def new_record():
    """Clade table with the founder as clade 0, accumulators, and the engulfing pair buffer."""
    CT = np.zeros((CMAX, 5)); cs = np.zeros(1, np.int64); cacc = np.zeros((CMAX, NCA))
    pairs = np.zeros((PMAX, 2), np.int32); pc = np.zeros(2, np.int64)
    return CT, cs, cacc, pairs, pc

def reseed(seed, step):
    """Exact replay: the generator is set from (seed, step) at every checkpoint."""
    seed_numba((seed * 1000003 + step) % 2147483647)
