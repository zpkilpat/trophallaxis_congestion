# Chemotactic recruitment offsets congestion in honeybee trophallaxis

MATLAB code that reproduces every figure of the paper, in the main text and the appendices.

Honeybees share food mouth to mouth (trophallaxis), and a hungry bee must first find a donor. In the model, searchers are velocity-jump walkers on a ring that bias their turning toward a scent cue, and fed bees that linger at a donor sharpen that cue. Because a donor serves the bees around it in random order, the same lingering bees that guide searchers also crowd the donor. The code computes exact first-passage statistics for cue-guided search, the recruitment gain a colony draws from its scenting bees, and a reduced model of colony feeding with one and two donors, and it checks the reduction against individual-based simulations. Feeding is fastest at an intermediate abandonment rate, and with two donors the cue reports which donor is serving fastest, because only served bees scent.

## Citation

William Magrogan, Orit Peleg, and Zachary P. Kilpatrick. Chemotactic recruitment offsets congestion in honeybee trophallaxis. Submitted to PRX Life (2026).

```bibtex
@article{magrogan2026chemotactic,
  author  = {Magrogan, William and Peleg, Orit and Kilpatrick, Zachary P.},
  title   = {Chemotactic recruitment offsets congestion in honeybee trophallaxis},
  journal = {PRX Life, submitted},
  year    = {2026}
}
```

## Requirements

MATLAB R2023b or later. No toolboxes are needed. The figures are written with `exportgraphics`.

## Quick start

Clone the repository, `trophallaxis_congestion`, from https://github.com/zpkilpat/trophallaxis_congestion. Then start MATLAB in the repository root, since outputs are written to paths relative to it.

```matlab
run_all                   % all nine article figures, then the regression checks
figure3                   % a single figure
figure3('fast', true)     % coarse grids and few realizations, a quick smoke test (figures 3 to 8; `figure9` reads the production two-donor table that `figure5` caches, so run `figure5` first)
```

From a shell:

```bash
matlab -batch "run_all"
```

`figureN` draws Fig. N of the article. It writes the panels to `output/figureN/` and the lettered figure to `output/figs/figN.pdf`. Figs. 3 and 4 share the single-donor gain tables, and Figs. 5 and 6 share the two-donor arrival table; each is cached in the figure's own folder. `figure9` reads the two-donor table that `figure5` caches, and `run_all` runs them in order.

## Figures

| file | article figure | contents |
| --- | --- | --- |
| `figure1.m` | 1, main text | model schematic, search geometry, population flow, and the two abandonment regimes |
| `figure2.m` | 2, main text | search for one donor: mean first-passage time, reliance, gain surface |
| `figure3.m` | 3, main text | recruitment gain, feeding-time landscape, lingering boundary at one donor |
| `figure4.m` | 4, main text | reduced model against simulation, crowd composition, optimal abandonment against fed fraction |
| `figure5.m` | 5, main text | two competing donors: total benefit, cue contrast, optimal abandonment |
| `figure6.m` | 6, main text | gain elasticity, reduction against simulation, stochastic symmetry at identical donors |
| `figure7.m` | 7, appendix | search for one donor across parameters |
| `figure8.m` | 8, appendix | the averaged recruitment gain and its limits |
| `figure9.m` | 9, appendix | abandonment limits, reduction error, non-exponential service |

Each figure file holds its panels as local functions (`panel_a`, `panel_b`, ...) below the driver.

## Caches and runtime

Monte Carlo ensembles, gain tables, and parameter sweeps are cached as `.mat` files under `output/figureN/`. Each cache stores a key built from the parameters, the grids, and a tag, and a cache whose key does not match is rebuilt, so a cache can never silently feed a figure computed under other settings. Delete a cache file to force a rebuild.

Caches are not tracked in this repository. The cached ensembles used for the published figures are available from Zenodo at DOI `10.5281/zenodo.XXXXXXX`. Unpack the archive in the repository root so that the files land in `output/figure1/` to `output/figure9/` (and `output/linger_tables_cache.mat`, used by `selftest`).

| | runtime |
| --- | --- |
| `run_all` with the Zenodo caches | about 10 minutes |
| `run_all` from scratch | roughly 2 to 3 hours |

Times are for MATLAB R2023b on a laptop with an 8-core Intel Core i9. Most of a rebuild goes to the individual-based ensembles of Figs. 4, 6, and 9. Random streams are seeded per realization, so a rebuild reproduces the cached ensembles.

## Repository layout

```
.
├── run_all.m               regenerate every figure, then run selftest
├── figure1.m ... figure9.m one file per figure, panels as local functions
├── src/
│   ├── model/              search, gain, reduced dynamics, individual-based models
│   │   ├── afpt.m          exact first passage on the ring with a displaced tent cue
│   │   ├── lingermodel.m   single donor, gain tables, reduced dynamics, individual-based model
│   │   ├── linger_tables.m cached ridge table and gain surface
│   │   ├── atwo.m          two-target splitting probabilities and arrival rates
│   │   ├── twogain.m       exact two-donor arrival table (key from twogain_key.m)
│   │   ├── twodonor.m      reduced two-donor system and its service-asymmetry sweep
│   │   ├── twoasym.m       large abandonment-rate expansion at two donors
│   │   ├── abm2.m          two-donor individual-based model
│   │   └── abm_svc.m       single donor with a general service-time law
│   ├── plotting/           figure style, export, assembly, schematic glyphs
│   └── checks/             selftest and the printed checks quoted in the text
├── output/
│   ├── figs/               the article's figures fig1 to fig9 (PDF)
│   └── figureN/            panel exports and caches (not tracked)
├── LICENSE
└── README.md
```

## Reference values

For the representative parameters `N = 120`, `mu = 1`, `v = 1`, `beta0 = 1`, `sigma = 0.8`, `s_m = 1`, `kbar = 0.0314`, `p0 = 0.75`. The drivers print these and warn if a value moves.

| quantity | value |
| --- | --- |
| critical offset `p*/sigma` | 0.878 |
| gain ceiling `G(0)` | 3.434 |
| first-scenter value `g1` at `p0 = 0.75, 1.13, 1.52` | 1.552, 1.297, 1.157 |
| simulated optimum | `delta = 0.200`, `T_1/2 = 90.2 ± 0.7` |
| reduced optimum | `delta = 0.248`, `T_1/2 = 85.7` |

## License

MIT, see `LICENSE`.

## Funding

This work was supported by the National Science Foundation through grants DMS-2527338 and IIS-2616531. ZPK acknowledges support from the National Institute for Theory and Mathematics in Biology through the National Science Foundation (grant number DMS-2235451) and the Simons Foundation (grant number MP-TMPS-00005320).
