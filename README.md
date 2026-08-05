# Heat Waves Dashboard

An interactive version of this repository's heat-metric figures, laid out like the
heat-waves page of the [US Extreme Weather and Climate Change
Dashboard](https://extremeweather.thehonestbroker.org/heat-waves/) and using that site's
published change-detection rule.

**Live: https://aedessler.github.io/extreme-heat-dashboard/**

Open `index.html` in a browser. That is the whole deployment: one static page plus a
generated data file. No build step, no server, no dependencies, no network requests at
view time. Anything that can serve static files can host it; GitHub Pages serves this
repository's `main` branch directly.

## What it shows

Four metrics across three dataset/region cells — twelve series, switchable from the top of
the page:

| | CONUS | NH land, 24–50°N |
| --- | --- | --- |
| **GHCN-Daily** | 1,266 stations, FLs.52j-adjusted, area-weighted | *not offered — see below* |
| **Berkeley Earth** | 281 land grid points | 1,140 land grid points |

- **Heat Wave Index** — events per site per year (EPA/Kunkel)
- **TXx** — annual hottest daily high, °F (ETCCDI)
- **WSDI** — warm spell duration, days per year (ETCCDI)
- **TN90p** — warm nights, % of days (ETCCDI)

"NH land" is all Northern Hemisphere land between 24°N and 50°N at every longitude,
CONUS included. That is the same latitude band as the CONUS box, so the comparison holds
latitude fixed and varies only longitude — and it is exactly the band the Berkeley Earth
`preprocessed_nh_*` files were built for.

### Why there is no hemispheric GHCN series

GHCN aggregation *is* area-weighted, by the repository's own
`temperature_data.area_weights` — `cos(lat) / stations-per-2°-cell`, which averages
stations within a cell and then combines cells in proportion to area. That is
algebraically the same two-stage scheme Berkeley uses; Berkeley skips the first stage only
because it is already a grid.

But weighting can only redistribute among cells that contain stations, and in this band
GHCN barely has any outside the US. Using Berkeley's land mask as the denominator on the
same 2° grid:

| Cell | Land cells reached |
| --- | --- |
| GHCN CONUS | 220 / 281 (78% of land area) |
| GHCN 24–50°N band | 274 / 1,140 (23%) |
| — outside CONUS | **60 / 873 (6.9%)** |
| Berkeley CONUS | 281 / 281 |
| Berkeley 24–50°N band | 1,140 / 1,140 |

A hemispheric GHCN series built this way keeps 78% of its area weight in the US and
correlates **0.962** with the CONUS box of the same raw stations — it is CONUS wearing a
hemisphere label. Coarsening the weighting grid moves the number without fixing the
problem (10° cells drop the US share to 50% and the correlation to 0.65, but then a single
station stands in for a 10°×10° region). So the cell is not built at all, and the page
disables the control and says why. Berkeley Earth is the hemispheric option.

One more thing worth stating plainly: the FLs.52j homogeneity offsets exist only for US
stations, so even a raw hemispheric GHCN series could never have been homogenized. That is
a second, independent reason Berkeley Earth is the right instrument for the band.

## What the data says

Over the full 1900–2024 record, with Mann–Kendall at p < 0.10:

| Metric | GHCN CONUS | Berkeley CONUS | Berkeley NH |
| --- | --- | --- | --- |
| Heat Wave Index | **detected** (p=0.041) | **detected** (p<0.001) | **detected** (p<0.001) |
| TXx | no (p=0.483) | **detected** (p=0.014) | **detected** (p<0.001) |
| WSDI | no (p=0.103) | **detected** (p<0.001) | **detected** (p<0.001) |
| TN90p | **detected** (p<0.001) | **detected** (p<0.001) | **detected** (p<0.001) |

Berkeley Earth detects an increase in all eight of its cells. GHCN CONUS detects two of
four: TXx shows essentially no monotonic trend across the full record, and WSDI lands just
the wrong side of the threshold at p=0.103. Those two therefore read **dataset dependent**
in the summary table, which is the reference site's own label for exactly this situation.

Set the window to 1970–2024 and **every cell detects an increase**, all at p ≤ 0.003. The
disagreement is not about whether recent decades have warmed; it is about how much a
full-record linear fit is diluted by the 1930s, which are prominent in the GHCN station
record and much weaker in the Berkeley grid.

## Change detection

A change is called **detected** when the Mann–Kendall two-sided p < 0.10 — the IPCC's
example threshold rather than the conventional 0.05. The reported rate is the Theil–Sen
slope. The test recomputes for whatever window you select, so narrowing the window changes
both the slope and the verdict.

The [reference dashboard](https://extremeweather.thehonestbroker.org/methodology) pairs
that likelihood test with a second, magnitude-versus-variability criterion — the fitted
change must reach 25% of the variable's 5th-to-95th percentile spread. **That second test
is deliberately not applied here.** Detection on this page is significance alone.

Significance is not the same as importance, so each chart still draws the 66% and 90%
historical range bands: a trend can clear p < 0.10 while remaining small against ordinary
year-to-year spread, and the bands let you see that directly. Bear in mind too that a
linear trend is a poor summary of a series whose largest feature is a single decade — for
CONUS the 1930s dominate.

## Metric definitions

All four run through one implementation, `scripts/metrics.py`, for every dataset and
region, so a difference between two panels is a difference in the data rather than in the
arithmetic. A *site* is a station for GHCN and a land grid point for Berkeley Earth.

Site screen (GHCN): within 1900–2024, at least 100 years with observations and valid
unflagged TMAX and TMIN on at least 80% of all possible days. Berkeley grid points qualify
wherever the land mask leaves paired data. Only site-years at least 90% complete are
scored.

- **Heat Wave Index** — rolling 4-day means of (TMAX+TMIN)/2 within each year; an event is
  a continuous run of overlapping windows above the site's 1-in-10-year threshold, itself
  the 90th percentile of that site's annual maximum 4-day means. The threshold is taken
  over every year's maximum while only complete years are scored; that is the parent
  project's convention, kept so the numbers match.
- **TXx** — the hottest daily maximum at each site each year, averaged across sites, °F.
- **WSDI** — days inside runs of at least 6 consecutive days whose TMAX exceeds the site's
  calendar-day 90th percentile, built from 1961–1990 with a 5-day centered circular window.
- **TN90p** — the percentage of each year's days whose TMIN exceeds the same style of
  calendar-day 90th percentile threshold.

Spatial mean: stations by `cos(latitude) / stations-in-2°-cell`, grid points by
`cos(latitude)`, renormalized over the sites valid in each year.

### Levels differ between datasets; compare the trends

A thermometer and a 2° grid cell do not measure the same thing. Averaging over a cell damps
daily extremes, so gridded TXx sits about 7°F below station TXx over CONUS, while gridded
WSDI runs about 4 days a year above station WSDI because a smoother series crosses and stays
above its own calendar-day threshold more persistently. This is a property of spatial
support, not an error in either dataset, and the affected panels say so. It is also why the
detection rule — which is scale-free, comparing a trend with that series' own variability —
is the right basis for comparing the two.

### Agreement with the repository's published results

Three cells can be checked against files this repository already publishes, and all three
match to the precision those files are stored at:

| Series | Reference | Agreement |
| --- | --- | --- |
| `ghcn_conus_hwi` | `data/heat_wave_index_raw_adj_wtd.csv` | exact |
| `ghcn_conus_txx` | `data/annual_hottest_daily_high_raw_adj_wtd.csv` | exact |
| `berkeley_conus_hwi` | `data/heat_wave_index_berkeley.csv` | exact |

WSDI additionally reproduces `plot_modules/wsdi.py` bit for bit. TN90p has no existing
reference implementation here, but shares the verified threshold machinery with WSDI.

### Not reproduced

The reference page also carries WSDI and TXx from NOAA nClimGrid-Daily. That homogenized
gridded product is not part of this repository, so those two panels are absent; Berkeley
Earth plays the equivalent role as the homogenized, gridded comparison.

## Rebuilding the data

```bash
python scripts/build_dashboard_data.py
```

Writes `dashboard-series.json` (canonical, diffable) and `dashboard-series.js` (what the
page loads). Both come from the same payload and a test asserts they agree. A full rebuild
takes roughly 20 minutes.

### Requirements

Viewing needs none of this — the generated data is committed. These are only for
regenerating it, which also requires checking this directory out inside the parent
project ([GHCN-daily-US-corrected](https://github.com/aedessler/GHCN-daily-US-corrected)),
since the builder reads its checkpoints and the tests compare against its published CSVs.

- **Python with NumPy, pandas and xarray.**
- **The parent project's checkpoints**, about 830 MB in `../data/`:
  `checkpoint_tmax.dat`, `checkpoint_tmin.dat`, `checkpoint_meta.npz`.
- **Berkeley Earth preprocessed daily files.** The builder looks in
  `$BERKELEY_DAILY_DIR`, defaulting to the DOE 6.3 `processed_data` directory in iCloud.
  It needs `preprocessed_us_TMAX_data.nc`, `preprocessed_us_TMIN_data.nc` (about 750 MB
  each) and `preprocessed_nh_TMAX_data.nc`, `preprocessed_nh_TMIN_data.nc` (about 4.5 GB
  each). Provenance and the generating code are documented in
  [`README_processed_data.md`](https://github.com/aedessler/DOE_report_section_6_3/blob/main/README_processed_data.md).
- **About 2 GB of disk** for `cache/ghcn/`, where NOAA's per-station files are cached. The
  cache is git-ignored and reused on later runs.

### Options

| Flag | Purpose |
| --- | --- |
| `--cells ghcn:conus,berkeley:nh` | Rebuild only these cells; the rest are kept from the existing output |
| `--berkeley-dir PATH` | Where the `preprocessed_*.nc` files live |
| `--cache-dir PATH` | Where NOAA station downloads are cached, default `cache/ghcn` |
| `--workers N` | Parallel station workers, default 10 |
| `--limit N` | Process only the first N GHCN candidates (development) |
| `--output PATH` | JSON output path; the `.js` companion is written alongside it |

Valid cells are `ghcn:conus`, `berkeley:conus` and `berkeley:nh`; asking for `ghcn:nh`
fails with the reason. Because `--cells` merges into the existing output,
`--cells berkeley:conus` is a fast way to re-derive one panel without touching the other
eight series.

## Verifying

```bash
python -m unittest discover tests
```

Checks that all twelve series exist and no others, that detection is Mann–Kendall alone
with no magnitude parameter left in the payload to imply otherwise, that the excluded GHCN
hemispheric cell is genuinely absent *and* carries an explanation the page can show, that
every region still has at least one dataset, that the published CONUS cells reproduce their CSVs, that every
field the page reads is present, that all four metrics in a cell share one site set, that
gridded cells reach all of their land, and that values stay physically plausible.

To view the page over HTTP rather than `file://`:

```bash
python -m http.server 8123
```

## Files

| Path | |
| --- | --- |
| `index.html` | The entire page: markup, styling, statistics, charts and exports |
| `dashboard-series.js` | Generated data the page loads |
| `dashboard-series.json` | The same payload, for reuse elsewhere |
| `scripts/metrics.py` | The four metric kernels, shared by every dataset and region |
| `scripts/build_dashboard_data.py` | Builder and orchestration |
| `tests/test_dashboard_data.py` | Checks on the generated data |
| `cache/` | Git-ignored NOAA download cache |
