"""usage: run.py SKY SEED STEPS OUTDIR [RESUME]. Sky A: light 0.30; B: 0.24; C: sine 0.23..0.60, period 3,000,000, starting at 0.60 (a typed input).
Writes OUTDIR/census.jsonl (per clade, every 3,000 steps), clades.jsonl (one row per new clade), ck/step_N.npz every 150,000 steps
(generator reseeded from (seed, step) there), status.json. Curve A at c=1 (z=1.3, yint=c), pass 7's constants."""
import sys, json, os, time, math, glob, numpy as np
import world as W
PATCH = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'patch.npy'))   # static light factor per cell, mean 1

CENSUS, CKPT = 3000, 150000
C_COST, Z, FOUNDER = 1.0, 1.3, (0.5, 0.1, 0.05)

def sky_light(sky, s):
    if sky.startswith("L"): return float(sky[1:])   # smoke only: constant light
    if sky == "A": return 0.30
    if sky == "B": return 0.24
    return 0.415 + 0.185 * math.cos(2 * math.pi * s / 3_000_000)   # sky C; sky D reports column 0, which is exactly C

LATTICE = 40

def sky_arr(sky, s):
    """Per-column light, length LATTICE, float64. Sky D: sky C at step s - 750,000 x a / (L - 1) for column a. Others: the scalar everywhere."""
    if sky in ("G", "U"):   # G: sky D frozen at step 0; U: the mean of G's columns, everywhere
        g = np.array([0.415 + 0.185 * math.cos(2 * math.pi * (0 - 750_000 * a / (LATTICE - 1)) / 3_000_000) for a in range(LATTICE)])
        return g if sky == "G" else np.full(LATTICE, g.mean())
    if sky == "D":
        return np.array([0.415 + 0.185 * math.cos(2 * math.pi * (s - 750_000 * a / (LATTICE - 1)) / 3_000_000) for a in range(LATTICE)])
    return np.full(LATTICE, sky_light(sky, s))

def save_ck(path, B, n, C, D, CT, cs, s):
    ids = np.unique(B[:n, W.CL].astype(np.int64))   # only clades with living bodies are ever read again
    np.savez_compressed(path, B=B[:n], C=C, D=D, ids=ids, CT=CT[ids], ncl=int(cs[0]), step=s)

def load_ck(path):
    z = np.load(path)
    B = np.zeros((W.NMAX, W.NF)); n = len(z["B"]); B[:n] = z["B"]
    CT, cs, cacc, pairs, pc = W.new_record(); CT[z["ids"]] = z["CT"]; cs[0] = int(z["ncl"])
    return B, n, z["C"].copy(), z["D"].copy(), CT, cs, cacc, pairs, pc, int(z["step"])

def census(s, B, n, CT, cs, cacc, pairs, pc, C, D, clades_f, led):
    cl = B[:n, W.CL].astype(np.int64); k = int(cs[0])
    cnt = np.bincount(cl, minlength=k)[:k]; mass = np.bincount(cl, weights=B[:n, W.V], minlength=k)[:k]
    gm = np.stack([np.bincount(cl, weights=B[:n, W.QA + j], minlength=k)[:k] for j in range(3)], 1)
    gk = np.bincount(cl, weights=B[:n, W.QK], minlength=k)[:k]
    gq = np.bincount(cl, weights=B[:n, W.QS], minlength=k)[:k]; dmf = (B[:n, W.DM] > 0.5).astype(np.float64); gd = np.bincount(cl, weights=dmf, minlength=k)[:k]
    gh = np.bincount(cl, weights=B[:n, W.HS], minlength=k)[:k]; gdh = np.bincount(cl, weights=B[:n, W.HS] * dmf, minlength=k)[:k]; gdn = gd
    pf = PATCH[B[:n, W.X].astype(np.int64), B[:n, W.Y_].astype(np.int64)]
    gp = np.bincount(cl, weights=pf, minlength=k)[:k]
    p = pairs[:pc[0]]
    byv = {}
    if len(p):
        u, c = np.unique(p.astype(np.int64), axis=0, return_counts=True)
        for (v, kk), cc in zip(u, c): byv.setdefault(int(v), {})[int(kk)] = int(cc)
    act = np.where(cacc[:k].any(1) | (cnt > 0))[0]
    keep = [int(a) for a in act if cnt[a] >= 0.002 * max(n, 1)]
    ks = set(keep); thin = np.array([a for a in act if int(a) not in ks], dtype=np.int64)
    def row(a):
        return dict(id=a, bodies=int(cnt[a]), mass=round(float(mass[a]), 1),
                    genes=[round(float(x), 3) for x in (gm[a] / cnt[a] if cnt[a] else [0, 0, 0])],
                    income=[round(float(x), 1) for x in cacc[a, :3]], births=int(cacc[a, 3]),
                    d_hunger=int(cacc[a, 4]), d_hazard=int(cacc[a, 5]), d_engulfed=int(cacc[a, 6]),
                    engulfed_by=byv.get(a, {}), engulfings_made=int(cacc[a, 8]),
                    k=round(float(gk[a] / cnt[a]), 3) if cnt[a] else 0.0, kin_spared=int(cacc[a, 9]),
                    patch=round(float(gp[a] / cnt[a]), 3) if cnt[a] else 0.0,
                    qs=round(float(gq[a] / cnt[a]), 4) if cnt[a] else 0.0, dormant=round(float(gd[a] / cnt[a]), 4) if cnt[a] else 0.0,
                    hs=round(float(gh[a] / cnt[a]), 4) if cnt[a] else 0.0, hs_dorm=round(float(gdh[a] / gd[a]), 4) if gd[a] > 0 else 0.0)
    rows = [row(a) for a in keep]
    if len(thin):
        t = thin; tc = int(cnt[t].sum())
        o = dict(id="other", clades=len(t), bodies=tc, mass=round(float(mass[t].sum()), 1),
                 genes=[round(float(x), 3) for x in (gm[t].sum(0) / tc if tc else [0, 0, 0])],
                 income=[round(float(x), 1) for x in cacc[t, :3].sum(0)], births=int(cacc[t, 3].sum()),
                 d_hunger=int(cacc[t, 4].sum()), d_hazard=int(cacc[t, 5].sum()), d_engulfed=int(cacc[t, 6].sum()),
                 engulfed_by={}, engulfings_made=int(cacc[t, 8].sum()),
                 k=round(float(gk[t].sum() / tc), 3) if tc else 0.0, kin_spared=int(cacc[t, 9].sum()),
                 patch=round(float(gp[t].sum() / tc), 3) if tc else 0.0,
                 qs=round(float(gq[t].sum() / tc), 4) if tc else 0.0, dormant=round(float(gd[t].sum() / tc), 4) if tc else 0.0,
                 hs=round(float(gh[t].sum() / tc), 4) if tc else 0.0, hs_dorm=round(float(gdh[t].sum() / gd[t].sum()), 4) if gd[t].sum() > 0 else 0.0)
        rows.append(o)
    bands = []
    xa = B[:n, W.X].astype(np.int64); hun = (gm[:, 2] / np.maximum(cnt, 1) > 0.15)[cl] if n else np.zeros(0, bool)
    for lo in (0, 10, 20, 30):
        m = (xa >= lo) & (xa < lo + 10); nb = int(m.sum()); nd = int((B[:n, W.DM][m] > 0.5).sum())
        bands.append(dict(cols=[lo, lo + 9], n=nb, dormant=nd, share=round(nd / nb, 4) if nb else 0.0, hunters=int((m & hun).sum()), pt=round(float(B[:n, W.PT][m].mean()), 4) if nb else 0.0))
    for j in range(4): bands[j]['woke'] = int(led[13 + j]); bands[j]['slept'] = int(led[17 + j])   # pass 20: transitions since the last census (reset there)
    return dict(step=s, n=int(n), top5=round(float((B[:n, W.X] < 5).mean()), 4) if n else 0.0, pt_mean=round(float(B[:n, W.PT].mean()), 4) if n else 0.0, pt_gt05=round(float((B[:n, W.PT] > 0.5).mean()), 4) if n else 0.0, pt_lt02=round(float((B[:n, W.PT] < 0.2).mean()), 4) if n else 0.0, light39=round(float(sky_arr(SKY, s)[-1]), 4), bands=bands, C=round(float(C.mean()), 1), D=round(float(D.mean()), 1), light=round(sky_light(SKY, s), 4),
                clades_made=k, pairs_dropped=int(pc[1]), clades=rows)

def write_new_clades(f, CT, cs, done):
    for k in range(done, int(cs[0])):
        f.write(json.dumps(dict(id=k, step=int(CT[k, 4]), parent=int(CT[k, 3]), genes=[round(float(x), 4) for x in CT[k, :3]])) + "\n")
    f.flush(); return int(cs[0])

if __name__ == "__main__":
    SKY, seed, steps, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    resume = len(sys.argv) > 5
    os.makedirs(f"{out}/ck", exist_ok=True)
    if resume:
        cks = sorted(glob.glob(f"{out}/ck/step_*.npz")); last = cks[-1]
        B, n, C, D, CT, cs, cacc, pairs, pc, s0 = load_ck(last)
        for fn, key in (("census", "step"), ("clades", "step")):
            rows = [l for l in open(f"{out}/{fn}.jsonl") if json.loads(l)[key] <= s0]
            open(f"{out}/{fn}.jsonl", "w").writelines(rows)
        W.reseed(seed, s0)
    else:
        B, n, C, D = W.start(seed, C_COST, Z, C_COST, founders=FOUNDER)
        CT, cs, cacc, pairs, pc = W.new_record(); CT[0, :3] = FOUNDER; CT[0, 3] = -1; cs[0] = 1; s0 = 0
        W.reseed(seed, 0); save_ck(f"{out}/ck/step_{0:09d}.npz", B, n, C, D, CT, cs, 0)
        open(f"{out}/census.jsonl", "w").close(); open(f"{out}/clades.jsonl", "w").close()
    led = np.zeros(21); flash = np.zeros((2, 2)); t0 = time.time(); end = "DONE"; s = s0
    cf = open(f"{out}/census.jsonl", "a"); kf = open(f"{out}/clades.jsonl", "a")
    if not resume:
        kf.write(json.dumps(dict(id=0, step=0, parent=-1, genes=list(FOUNDER))) + "\n"); done = 1
    else:
        done = int(cs[0])
    for s in range(s0 + 1, steps + 1):
        n = W.step(B, n, C, D, C_COST, Z, C_COST, led, flash, sky_arr(SKY, s), PATCH, W.MU, CT, cacc, cs, pairs, pc, s)
        if s % CENSUS == 0 or n == 0:
            cf.write(json.dumps(census(s, B, n, CT, cs, cacc, pairs, pc, C, D, kf, led)) + "\n"); cf.flush()
            done = write_new_clades(kf, CT, cs, done)
            cacc[:int(cs[0])] = 0; pc[:] = 0; led[13:21] = 0
            json.dump(dict(sky=SKY, seed=seed, step=s, n=int(n), wall=round(time.time() - t0)), open(f"{out}/progress.json", "w"))
        if n == 0:
            end = "EXTINCT"; break
        if s % CKPT == 0:
            save_ck(f"{out}/ck/step_{s:09d}.npz", B, n, C, D, CT, cs, s); W.reseed(seed, s)
    json.dump(dict(sky=SKY, seed=seed, end=end, last_step=s, steps=steps, wall=round(time.time() - t0), n=int(n), clades=int(cs[0])), open(f"{out}/status.json", "w"))
