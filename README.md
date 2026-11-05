# When a world breeds out a useful instinct

Every founder in this small world is born wanting to move toward the light, and the light is where food is. When does the world take that want away with nobody removing it? We test two candidate causes: crowding at the light, and predators that evolve to hunt where prey gathers.

The paper is [PAPER.md](PAPER.md).

## The world (five rules, see `world.py`)
1. Light is strongest at one edge of a 40-column glass.
2. Every founder carries a gene PT, the chance of stepping toward brighter light; founders start at PT 0.5.
3. Bodies in one cell share one food supply, so a crowd goes hungry; energy at zero kills.
4. Some bodies can engulf a neighbour (gene QE); hunters are bodies that evolved it.
5. Children copy their parents' genes with small changes (mutation rate and size 0.02).

Upkeep (`MAINT_MULT`) is the energy cost of staying alive each step.

## Arms
`run.py SKY SEED STEPS OUTDIR` takes the sky as its first argument. Environment: `MAINT_MULT` (upkeep multiple) and `NOHUNT=1` (QE is zero in founders and children, so nobody can engulf). Sky `G` is a fixed gradient, `U` an even sky at the gradient's mean, `D` the moving sky, `C` the whole glass brightening and dimming together. Install with `uv sync`; the run needs numba and numpy.

| arm | sky | hunters | upkeep | steps | seeds | command (from this folder) |
|---|---|---|---|---|---|---|
| H | D moving | on | x2 | 3M | 1-3 | `MAINT_MULT=2 python run.py D SEED 3000000 runs/H_SEED` |
| N | D moving | off | x2 | 3M | 1-3 | `NOHUNT=1 MAINT_MULT=2 python run.py D SEED 3000000 runs/N_SEED` |
| E | C whole-glass | on | x2 | 3M | 1-3 | `MAINT_MULT=2 python run.py C SEED 3000000 runs/E_SEED` |
| GH | G fixed gradient | on | x2 | 1.5M | 1-3 | `MAINT_MULT=2 python run.py G SEED 1500000 runs/GH_SEED` |
| GN | G fixed gradient | off | x2 | 1.5M | 1-3 | `NOHUNT=1 MAINT_MULT=2 python run.py G SEED 1500000 runs/GN_SEED` |
| UH | U even | on | x2 | 1.5M | 1-3 | `MAINT_MULT=2 python run.py U SEED 1500000 runs/UH_SEED` |
| UN | U even | off | x2 | 1.5M | 1-3 | `NOHUNT=1 MAINT_MULT=2 python run.py U SEED 1500000 runs/UN_SEED` |
| GD | G fixed gradient | off | x3 | 1.5M | 1-3 | `NOHUNT=1 MAINT_MULT=3 python run.py G SEED 1500000 runs/GD_SEED` |

Each run writes `census.jsonl` every 3,000 steps (the csv files in `data/` are per-checkpoint summaries of mean PT, bodies, and share in the five brightest columns). The moving-sky arms (H, N, E) were first runs, sampled every 150,000 steps; the 'moving sky' file `data/moving_sky.csv` holds them. The moving sky was discarded as the main design: its no-hunter runs moved in lockstep with the sky's wave, so the sky drove part of the curve. That is why the sky is fixed in the other arms.

One copy of `run.py` and `world.py` serves every arm. The three scratch versions differed only in small ways: the earliest `run.py` lacked the `G`/`U` skies, the second lacked the top-five-column share in the census, and the earliest `world.py` lacked `NOHUNT`. The versions here are the newest.

## Redraw the figures
```
python -m venv venv && venv/bin/pip install matplotlib numpy pandas
venv/bin/python figures.py
```
`figures.py` reads only `data/*.csv` and writes `fig1_moving_sky.png`, `fig2_fixed_sky.png`, `fig3_density.png`.

## What the runs show
On a fixed light gradient, the want toward light is held near 0.4 without hunters and is driven to about 0.02 within 30,000 steps once hunters evolve. A hunterless world denser and more crowded at the light than the hunter world kept the want and raised it. Crowding does not breed the want out; hunters do.

## Limits
1. Counts are not matched: GH is the less dense arm (peaks 9,000-10,800 bodies; GD 11,600-12,300), and no upkeep setting reached 7,000-11,000 without extinction. So "fewer bodies lowers PT" is not excluded by GD.
2. Upkeep itself moved PT (x2: 0.41; x3: 0.65-0.84), untested why, so GN and GD differ in two ways.
3. One world, one gradient, one length, three seeds per arm.
4. The predictions made before the runs were wrong twice (density arm PT 0.3-0.45; even sky holding at 0.5).

## Licence
MIT. Copyright (c) 2026 mycelium.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated files, to deal in the software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies, subject to the inclusion of this notice. The software is provided "as is", without warranty of any kind.
