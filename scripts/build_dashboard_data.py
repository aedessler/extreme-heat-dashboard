#!/usr/bin/env python3
"""Build every series the dashboard charts: 2 datasets x 2 regions x 4 metrics.

    dataset  ghcn      station observations, NOAA GHCN-Daily
             berkeley  gridded land data, Berkeley Earth daily

    region   conus     24.5-49.5N, 125-66W -- the parent project's station set
             nh        24-50N, all longitudes (includes CONUS)

    metric   hwi, txx, wsdi, tn90p -- see scripts/metrics.py

Every cell runs through the same kernels in ``metrics.py``, so a difference
between two cells is a difference in the data, never in the arithmetic. The
CONUS GHCN cell reproduces the parent project's published figures exactly;
tests assert this.

GHCN CONUS uses the FLs.52j homogeneity-adjusted checkpoints because those
offsets exist for US stations. GHCN NH uses raw observations because they do
not exist outside the US, so an adjusted hemispheric series is not available at
any price. That asymmetry is recorded on each series and shown on the page.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics as M  # noqa: E402


DASHBOARD_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = DASHBOARD_DIR.parent
DEFAULT_CACHE = DASHBOARD_DIR / "cache" / "ghcn"
DEFAULT_OUTPUT = DASHBOARD_DIR / "dashboard-series.json"
DEFAULT_BERKELEY_DIR = Path(
    "/Users/adessler/Library/Mobile Documents/iCloud~md~obsidian/Documents"
    "/BulletVault/plotting code/DOE 6.3 extreme heat figures/processed_data"
)

STATIONS_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
INVENTORY_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"
STATION_DATA_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_station/{station}.csv.gz"

MIN_YEARS = 100
MIN_DATA_FRAC = 0.80
GHCN_SOURCE_URL = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/"
BERKELEY_SOURCE_URL = "https://berkeleyearth.org/data/"

# (label, lat_min, lat_max, lon_min, lon_max)
REGIONS = {
    "conus": {
        "label": "CONUS",
        "blurb": "The contiguous United States.",
        "box": ("the CONUS box (24.5-49.5N, 125-66W)", 24.5, 49.5, -125.0, -66.0),
    },
    "nh": {
        "label": "NH land, 24-50°N",
        "blurb": "All Northern Hemisphere land between 24°N and 50°N, CONUS included.",
        "box": ("the 24-50N land band, all longitudes", 24.0, 50.0, -180.0, 180.0),
    },
}
DATASETS = {
    "ghcn": {"label": "GHCN-Daily", "siteLabel": "stations"},
    "berkeley": {"label": "Berkeley Earth", "siteLabel": "land grid points"},
}

# GHCN is deliberately not offered for the hemispheric band. Its long-record
# stations occupy only 274 of that band's 1,140 land cells, and just 60 of the
# 873 outside CONUS -- under 7%. Area weighting redistributes among cells that
# hold stations and cannot give weight to empty ones, so the "hemispheric"
# station series correlates 0.96 with the CONUS box of the same raw stations.
# It would be CONUS wearing a hemisphere label, so it is not built at all.
AVAILABLE_CELLS = [("ghcn", "conus"), ("berkeley", "conus"), ("berkeley", "nh")]
EXCLUDED_CELLS = {
    "ghcn_nh": (
        "its qualifying long-record stations reach under 7% of the land cells in this band "
        "outside CONUS, so a hemispheric station estimate is not supportable — it would track "
        "CONUS almost exactly (r = 0.96) whatever weighting is applied. Berkeley Earth covers "
        "every land cell in the band and is the hemispheric option here."
    ),
}


@dataclass(frozen=True)
class Candidate:
    station: str
    latitude: float
    longitude: float


@dataclass
class SiteSet:
    """Per-site metric values plus the coordinates needed to weight them."""

    values: dict[str, np.ndarray]  # metric -> (sites, years)
    latitudes: np.ndarray
    longitudes: np.ndarray
    weights: np.ndarray
    identifiers: list[str] | None = None


# ---------------------------------------------------------------
# GHCN raw archive
# ---------------------------------------------------------------
def _download(url: str, destination: Path) -> None:
    """Download atomically, leaving an existing non-empty cache untouched."""
    if destination.exists() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "ghcn-dashboard/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=90) as response, temporary.open("wb") as output:
            while block := response.read(1024 * 1024):
                output.write(block)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def _candidate_stations(inventory_path: Path, lat_min: float, lat_max: float) -> list[Candidate]:
    """Stations in a latitude band whose TMAX/TMIN inventories overlap >=100 years.

    This is only a cheap pre-filter against the inventory so the builder does not
    download the whole archive; the real completeness test needs the data itself.
    """
    elements: dict[str, dict[str, tuple[int, int, float, float]]] = {}
    with inventory_path.open(encoding="ascii") as source:
        for line in source:
            parts = line.split()
            if len(parts) < 6 or parts[3] not in {"TMAX", "TMIN"}:
                continue
            station, latitude, longitude, element, first, last = parts[:6]
            lat = float(latitude)
            if not lat_min <= lat <= lat_max:
                continue
            elements.setdefault(station, {})[element] = (
                max(int(first), M.STUDY_START),
                min(int(last), M.STUDY_END),
                lat,
                float(longitude),
            )

    candidates = []
    for station, records in elements.items():
        if not {"TMAX", "TMIN"}.issubset(records):
            continue
        overlap_start = max(records["TMAX"][0], records["TMIN"][0])
        overlap_end = min(records["TMAX"][1], records["TMIN"][1])
        if overlap_end - overlap_start + 1 < MIN_YEARS:
            continue
        candidates.append(Candidate(station, records["TMAX"][2], records["TMAX"][3]))
    return sorted(candidates, key=lambda item: item.station)


def _read_station(path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Parse one cached station file into (1, years, 366) TMAX/TMIN arrays."""
    shape = (1, M.N_YEARS, M.N_DOY)
    tmax = np.full(shape, np.nan, dtype=np.float32)
    tmin = np.full(shape, np.nan, dtype=np.float32)
    arrays = {"TMAX": tmax, "TMIN": tmin}
    try:
        with gzip.open(path, "rt", encoding="ascii", newline="") as source:
            for row in csv.reader(source):
                # station, YYYYMMDD, element, value (tenths degC), m-flag, q-flag, s-flag
                # A non-blank q-flag means the observation failed NOAA quality control.
                if len(row) < 6 or row[2] not in arrays or row[5]:
                    continue
                stamp = row[1]
                year = int(stamp[:4])
                if year < M.STUDY_START or year > M.STUDY_END:
                    continue
                day_index = date(year, int(stamp[4:6]), int(stamp[6:8])).timetuple().tm_yday - 1
                arrays[row[2]][0, year - M.STUDY_START, day_index] = float(row[3]) / 10.0
    except (OSError, EOFError, ValueError, csv.Error):
        return None

    possible = M.days_in_year().sum()
    for element in (tmax, tmin):
        if np.isfinite(element).sum() / possible < MIN_DATA_FRAC:
            return None
        if np.count_nonzero(np.any(np.isfinite(element), axis=2)) < MIN_YEARS:
            return None
    return tmax, tmin


def _station_metrics(job: tuple[str, Path]) -> tuple[str, dict[str, list[float]] | None]:
    """Worker: every metric for one station, or None if it fails the screen."""
    station, path = job
    pair = _read_station(path)
    if pair is None:
        return station, None
    try:
        values = M.compute_all(*pair)
    except ValueError:
        return station, None
    return station, {name: array[0].tolist() for name, array in values.items()}


def build_ghcn_stations(cache_dir: Path, region: str, workers: int,
                        limit: int | None = None) -> SiteSet:
    """Raw, area-weighted GHCN for every qualifying station in a region."""
    inventory = cache_dir / "ghcnd-inventory.txt"
    _download(STATIONS_URL, cache_dir / "ghcnd-stations.txt")
    _download(INVENTORY_URL, inventory)

    _, lat_min, lat_max, _, _ = REGIONS[region]["box"]
    candidates = _candidate_stations(inventory, lat_min, lat_max)
    if limit is not None:
        candidates = candidates[:limit]
    print(f"  candidates in {lat_min:g}-{lat_max:g}N: {len(candidates):,}", flush=True)

    missing = [c for c in candidates
               if not (cache_dir / "stations" / f"{c.station}.csv.gz").exists()]
    if missing:
        print(f"  downloading {len(missing):,} station files...", flush=True)
        for done, candidate in enumerate(missing, start=1):
            _download(STATION_DATA_URL.format(station=candidate.station),
                      cache_dir / "stations" / f"{candidate.station}.csv.gz")
            if done % 100 == 0 or done == len(missing):
                print(f"    {done:,}/{len(missing):,}", flush=True)

    jobs = [(c.station, cache_dir / "stations" / f"{c.station}.csv.gz") for c in candidates]
    by_station = {c.station: c for c in candidates}
    kept: list[Candidate] = []
    collected: dict[str, list[list[float]]] = {name: [] for name in M.METRIC_ORDER}

    with ProcessPoolExecutor(max_workers=workers) as executor:
        for done, (station, result) in enumerate(
                executor.map(_station_metrics, jobs, chunksize=8), start=1):
            if result is not None:
                kept.append(by_station[station])
                for name in M.METRIC_ORDER:
                    collected[name].append(result[name])
            if done % 250 == 0 or done == len(jobs):
                print(f"    processed {done:,}/{len(jobs):,}; retained {len(kept):,}",
                      flush=True)

    if not kept:
        raise RuntimeError(f"No stations in region {region!r} passed the completeness screen")

    order = np.argsort([c.station for c in kept])
    kept = [kept[int(i)] for i in order]
    latitudes = np.asarray([c.latitude for c in kept], dtype=float)
    longitudes = np.asarray([c.longitude for c in kept], dtype=float)
    values = {name: np.asarray(collected[name], dtype=float)[order] for name in M.METRIC_ORDER}
    return SiteSet(values, latitudes, longitudes,
                   M.station_area_weights(latitudes, longitudes),
                   [c.station for c in kept])


def build_ghcn_conus() -> SiteSet:
    """Adjusted, area-weighted CONUS from the parent project's checkpoints."""
    sys.path.insert(0, str(PROJECT_DIR))
    from temperature_data import (  # noqa: PLC0415
        area_weights, load_good_stations, station_coords,
    )

    identifiers, tmax, tmin = load_good_stations()
    print(f"  stations: {len(identifiers):,}", flush=True)
    latitudes, longitudes = station_coords(identifiers)
    values = M.compute_all(tmax, tmin)
    return SiteSet(values, np.asarray(latitudes), np.asarray(longitudes),
                   area_weights(identifiers, M.GRID_DEGREES), list(identifiers))


# ---------------------------------------------------------------
# Berkeley Earth gridded land data
# ---------------------------------------------------------------
def _read_berkeley(path: Path):
    import pandas as pd  # noqa: PLC0415
    import xarray as xr  # noqa: PLC0415

    if not path.exists():
        raise FileNotFoundError(f"Berkeley Earth input not found: {path}")
    with xr.open_dataset(path) as dataset:
        subset = dataset["temperature"].sel(
            time=slice(f"{M.STUDY_START}-01-01", f"{M.STUDY_END}-12-31")
        ).transpose("time", "latitude", "longitude")
        return (
            pd.DatetimeIndex(subset.time.values),
            np.asarray(subset.latitude.values, dtype=float),
            np.asarray(subset.longitude.values, dtype=float),
            np.asarray(subset.values, dtype=np.float32),
        )


def build_berkeley(data_dir: Path, region: str, block: int = 400) -> SiteSet:
    """Every metric on the Berkeley Earth land grid for one region.

    Sites are processed in blocks because the hemispheric grid is large enough
    that holding every intermediate for all points at once is wasteful.
    """
    prefix = "us" if region == "conus" else "nh"
    tmax_time, latitudes, longitudes, tmax = _read_berkeley(
        data_dir / f"preprocessed_{prefix}_TMAX_data.nc")
    tmin_time, tmin_lat, tmin_lon, tmin = _read_berkeley(
        data_dir / f"preprocessed_{prefix}_TMIN_data.nc")
    if not (np.array_equal(latitudes, tmin_lat) and np.array_equal(longitudes, tmin_lon)):
        raise ValueError("TMAX and TMIN grids differ")

    shared = tmax_time.intersection(tmin_time)
    if len(shared) == 0:
        raise ValueError("TMAX and TMIN share no dates")
    if not tmax_time.equals(tmin_time):
        print(f"  aligning to {len(shared):,} shared dates "
              f"({shared[0].date()} through {shared[-1].date()})", flush=True)
        tmax = tmax[tmax_time.get_indexer(shared)]
        tmin = tmin[tmin_time.get_indexer(shared)]

    grid_lat, grid_lon = np.meshgrid(latitudes, longitudes, indexing="ij")
    point_lat, point_lon = grid_lat.ravel(), grid_lon.ravel()
    flat_tmax = tmax.reshape(len(shared), -1)
    flat_tmin = tmin.reshape(len(shared), -1)
    del tmax, tmin

    # Keep grid points with at least one paired observation. That removes ocean
    # cells without relying on the very large, time-repeated land mask variable.
    land = np.any(np.isfinite(flat_tmax) & np.isfinite(flat_tmin), axis=0)
    point_lat, point_lon = point_lat[land], point_lon[land]
    flat_tmax, flat_tmin = flat_tmax[:, land], flat_tmin[:, land]
    n_sites = int(land.sum())
    print(f"  land grid points: {n_sites:,}", flush=True)

    year_index = shared.year.to_numpy() - M.STUDY_START
    doy_index = shared.dayofyear.to_numpy() - 1
    collected: dict[str, list[np.ndarray]] = {name: [] for name in M.METRIC_ORDER}
    for start in range(0, n_sites, block):
        stop = min(start + block, n_sites)
        shape = (stop - start, M.N_YEARS, M.N_DOY)
        site_tmax = np.full(shape, np.nan, dtype=np.float32)
        site_tmin = np.full(shape, np.nan, dtype=np.float32)
        site_tmax[:, year_index, doy_index] = flat_tmax[:, start:stop].T
        site_tmin[:, year_index, doy_index] = flat_tmin[:, start:stop].T
        for name, array in M.compute_all(site_tmax, site_tmin).items():
            collected[name].append(array)
        print(f"    metrics for grid points {stop:,}/{n_sites:,}", flush=True)

    values = {name: np.concatenate(blocks, axis=0) for name, blocks in collected.items()}
    return SiteSet(values, point_lat, point_lon, M.grid_area_weights(point_lat))


# ---------------------------------------------------------------
# Coverage accounting
# ---------------------------------------------------------------
def _domain_area(box) -> float:
    _, lat_min, lat_max, lon_min, lon_max = box
    rows = range(math.floor(lat_min / M.GRID_DEGREES), math.ceil(lat_max / M.GRID_DEGREES))
    columns = math.ceil(lon_max / M.GRID_DEGREES) - math.floor(lon_min / M.GRID_DEGREES)
    return sum(math.cos(math.radians((row + 0.5) * M.GRID_DEGREES)) * columns for row in rows)


def _cell_set(latitudes: np.ndarray, longitudes: np.ndarray) -> set[tuple[int, int]]:
    ilat = np.floor(np.asarray(latitudes) / M.GRID_DEGREES).astype(np.int64)
    ilon = np.floor(np.asarray(longitudes) / M.GRID_DEGREES).astype(np.int64)
    return set(zip(ilat.tolist(), ilon.tolist()))


def _cos_area(cells) -> float:
    return sum(math.cos(math.radians((row + 0.5) * M.GRID_DEGREES)) for row, _ in cells)


def _coverage(site_set: SiteSet, region: str, dataset: str,
              land_cells: list | None = None) -> dict:
    """How much of the nominal domain this site set actually samples.

    ``areaFraction`` counts every occupied 2-degree cell in full, land or not,
    so it is an upper bound on coverage rather than a flattering estimate.

    ``landCells*`` is the number that actually matters for a station network.
    Area weighting can only redistribute weight among cells that have data, so
    the honest measure of a station series is what fraction of the domain's
    *land* cells it occupies at all. Berkeley Earth's land mask supplies the
    denominator, on the same 2-degree grid.
    """
    box = REGIONS[region]["box"]
    occupied = _cell_set(site_set.latitudes, site_set.longitudes)

    coverage = {
        "sites": int(len(site_set.latitudes)),
        "siteLabel": DATASETS[dataset]["siteLabel"],
        "occupiedCells": len(occupied),
        "domain": box[0],
        "areaFraction": round(_cos_area(occupied) / _domain_area(box), 4),
    }

    if land_cells:
        land = {tuple(cell) for cell in land_cells}
        sampled = occupied & land
        coverage["landCellsSampled"] = len(sampled)
        coverage["landCellsTotal"] = len(land)
        coverage["landAreaFraction"] = round(_cos_area(sampled) / _cos_area(land), 4)

    if site_set.identifiers is not None:
        share = site_set.weights / site_set.weights.sum()
        by_country: Counter[str] = Counter()
        for identifier, value in zip(site_set.identifiers, share):
            by_country[identifier[:2]] += float(value)
        coverage["weightShareByCountry"] = [
            {"country": code, "share": round(value, 4)}
            for code, value in by_country.most_common(6)
        ]

    bands: Counter[str] = Counter()
    for value in site_set.latitudes:
        low = int(value // 10) * 10
        bands[f"{low}-{low + 10}N"] += 1
    coverage["sitesByLatitudeBand"] = [
        {"band": band, "sites": count}
        for band, count in sorted(bands.items(), key=lambda item: int(item[0].split("-")[0]))
    ]

    longitudes: Counter[str] = Counter()
    for value in site_set.longitudes:
        low = int(value // 60) * 60
        longitudes[f"{low}..{low + 60}"] += 1
    coverage["sitesByLongitudeBand"] = [
        {"band": band, "sites": count}
        for band, count in sorted(longitudes.items(), key=lambda item: int(item[0].split("..")[0]))
    ]
    return coverage


# ---------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------
def _points(values: np.ndarray, decimals: int) -> list[dict]:
    years = M.years_axis()
    return [
        {"x": int(year), "y": round(float(value), decimals + 3)}
        for year, value in zip(years, values)
        if np.isfinite(value)
    ]


def _series_for(dataset: str, region: str, site_set: SiteSet,
                land_cells: list | None = None) -> dict[str, dict]:
    coverage = _coverage(site_set, region, dataset, land_cells)
    adjustment = ("FLs.52j-adjusted" if dataset == "ghcn" and region == "conus"
                  else "Raw (unadjusted)" if dataset == "ghcn"
                  else "Berkeley Earth homogenized")
    if dataset == "ghcn":
        source = ("NOAA GHCN-Daily with FLs.52j homogenization offsets"
                  if region == "conus" else "NOAA GHCN-Daily raw daily observations")
        source_url = GHCN_SOURCE_URL
        weighting = "2-degree station-density and cosine-latitude weighted"
    else:
        source = "Berkeley Earth daily gridded land data"
        source_url = BERKELEY_SOURCE_URL
        weighting = "cosine-latitude weighted on the native 2-degree grid"

    out = {}
    for metric in M.METRIC_ORDER:
        annual = M.regional_mean(site_set.values[metric], site_set.weights)
        points = _points(annual, M.METRICS[metric]["decimals"])
        if not points:
            continue
        out[f"{dataset}_{region}_{metric}"] = {
            "id": f"{dataset}_{region}_{metric}",
            "dataset": dataset,
            "region": region,
            "metric": metric,
            "adjustment": adjustment,
            "weighting": weighting,
            "source": source,
            "sourceUrl": source_url,
            "coverage": coverage,
            "reliableWindow": [points[0]["x"], points[-1]["x"]],
            "points": points,
        }
    return out


def _payload_skeleton() -> dict:
    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "method": {
            "durationDays": M.DURATION,
            "returnYears": M.RETURN_YEARS,
            "minimumAnnualCompleteness": M.MIN_ANNUAL_FRAC,
            "minimumRecordYears": MIN_YEARS,
            "minimumOverallCompleteness": MIN_DATA_FRAC,
            "gridDegrees": M.GRID_DEGREES,
            "startYear": M.STUDY_START,
            "endYear": M.STUDY_END,
            "referencePeriod": [M.REF_START, M.REF_END],
            "percentile": M.PCTILE,
            "percentileWindowDays": M.WINDOW,
            "minimumSpellDays": M.MIN_RUN,
        },
        # Detection is Mann-Kendall significance alone: no magnitude-versus-
        # variability check, so no magnitude parameter belongs in here to imply
        # otherwise. A test asserts it stays absent.
        "detection": {
            "pThreshold": 0.10,
            "test": "Mann-Kendall",
            "slope": "Theil-Sen",
        },
        "datasets": {key: dict(value) for key, value in DATASETS.items()},
        "regions": {
            key: {"label": value["label"], "blurb": value["blurb"], "domain": value["box"][0]}
            for key, value in REGIONS.items()
        },
        "metrics": M.METRICS,
        "metricOrder": M.METRIC_ORDER,
        "datasetOrder": ["ghcn", "berkeley"],
        "regionOrder": ["conus", "nh"],
        "availableCells": [f"{dataset}_{region}" for dataset, region in AVAILABLE_CELLS],
        "excludedCells": EXCLUDED_CELLS,
        "series": {},
    }


def _write(payload: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    companion = output.with_suffix(".js")
    companion.write_text(
        "/* Generated by scripts/build_dashboard_data.py -- do not edit. */\n"
        f"window.HEAT_WAVE_DATA = {json.dumps(payload, separators=(',', ':'))};\n",
        encoding="utf-8",
    )
    print(f"Wrote {output}")
    print(f"Wrote {companion}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cells", default="all",
                        help="comma-separated dataset:region cells to rebuild, or 'all' "
                             "(e.g. ghcn:conus,berkeley:nh). Cells not rebuilt are kept "
                             "from the existing output.")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--berkeley-dir", type=Path,
                        default=Path(os.environ.get("BERKELEY_DAILY_DIR", DEFAULT_BERKELEY_DIR)))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=min(10, os.cpu_count() or 4))
    parser.add_argument("--limit", type=int, help="development only: first N GHCN candidates")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    if args.cells == "all":
        cells = list(AVAILABLE_CELLS)
    else:
        cells = []
        for token in args.cells.split(","):
            dataset, _, region = token.strip().partition(":")
            if (dataset, region) not in AVAILABLE_CELLS:
                excluded = EXCLUDED_CELLS.get(f"{dataset}_{region}")
                parser.error(
                    f"{token!r} is not built: {excluded}" if excluded else
                    f"unknown cell {token!r}; choose from "
                    f"{', '.join(f'{d}:{r}' for d, r in AVAILABLE_CELLS)}")
            cells.append((dataset, region))

    payload = _payload_skeleton()
    previous_regions: dict[str, dict] = {}
    if args.output.exists():
        previous = json.loads(args.output.read_text(encoding="utf-8"))
        previous_regions = previous.get("regions", {})
        rebuilding = {f"{d}_{r}_" for d, r in cells}
        payload["series"] = {
            key: value for key, value in previous.get("series", {}).items()
            if not any(key.startswith(prefix) for prefix in rebuilding)
        }
        if payload["series"]:
            print(f"Keeping {len(payload['series'])} series from the existing output")

    # Berkeley's land mask is the denominator for "how much land does this
    # station network actually reach", so build gridded cells before station
    # cells when both are requested, and carry the mask across partial rebuilds.
    land_cells = {region: previous_regions.get(region, {}).get("landCells")
                  for region in REGIONS}
    for dataset, region in sorted(cells, key=lambda cell: cell[0] != "berkeley"):
        print(f"\n=== {DATASETS[dataset]['label']} / {REGIONS[region]['label']} ===", flush=True)
        if dataset == "ghcn":
            site_set = (build_ghcn_conus() if region == "conus"
                        else build_ghcn_stations(args.cache_dir, region, args.workers, args.limit))
        else:
            site_set = build_berkeley(args.berkeley_dir, region)
            land_cells[region] = sorted(
                list(cell) for cell in _cell_set(site_set.latitudes, site_set.longitudes))
        payload["series"].update(
            _series_for(dataset, region, site_set, land_cells[region]))
        for metric in M.METRIC_ORDER:
            key = f"{dataset}_{region}_{metric}"
            if key in payload["series"]:
                values = [p["y"] for p in payload["series"][key]["points"]]
                print(f"  {metric:6s} n={len(values):3d}  mean={np.mean(values):.3f}  "
                      f"max={np.max(values):.3f}")
        coverage = payload["series"][f"{dataset}_{region}_hwi"]["coverage"]
        if "landCellsSampled" in coverage:
            print(f"  land cells sampled: {coverage['landCellsSampled']:,}"
                  f"/{coverage['landCellsTotal']:,} "
                  f"({coverage['landAreaFraction'] * 100:.1f}% of land area)")

    for region, cells_for_region in land_cells.items():
        if cells_for_region:
            payload["regions"][region]["landCells"] = cells_for_region
    _write(payload, args.output)


if __name__ == "__main__":
    main()
