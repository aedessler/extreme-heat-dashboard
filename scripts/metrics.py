#!/usr/bin/env python3
"""Heat-metric kernels shared by every dataset and region on the dashboard.

Each kernel takes site-major arrays shaped ``(sites, years, 366)`` in degrees
Celsius and returns a ``(sites, years)`` value per site-year, with incomplete
site-years set to NaN. "Site" is a station for GHCN and a land grid point for
Berkeley Earth, so the same code produces all four dataset/region combinations
and no metric can drift between them.

The definitions follow the parent project so the CONUS numbers reproduce the
published figures exactly:

  hwi    EPA/Kunkel Annual Heat Wave Index -- plot_heat_wave_index.py
  txx    annual maximum daily high         -- plot_modules/annual_hottest.py
  wsdi   ETCCDI warm spell duration index  -- plot_modules/wsdi.py
  tn90p  ETCCDI warm nights                -- same conventions as wsdi

Kernels are vectorised over the site axis, so passing a single site as
``(1, years, 366)`` gives the same answer as passing the whole set at once.
That is what lets the per-station streaming path and the whole-array gridded
path share one implementation.
"""

from __future__ import annotations

import warnings

import numpy as np


STUDY_START, STUDY_END = 1900, 2024
N_YEARS = STUDY_END - STUDY_START + 1
N_DOY = 366

# Heat Wave Index (EPA/Kunkel)
DURATION = 4
RETURN_YEARS = 10.0

# ETCCDI percentile indices
REF_START, REF_END = 1961, 1990
PCTILE = 90
WINDOW = 5
MIN_RUN = 6

MIN_ANNUAL_FRAC = 0.90
GRID_DEGREES = 2.0


def is_leap(year: int) -> bool:
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


def years_axis() -> np.ndarray:
    return np.arange(STUDY_START, STUDY_END + 1)


def days_in_year() -> np.ndarray:
    return np.asarray([366 if is_leap(int(y)) else 365 for y in years_axis()])


def _completeness(data: np.ndarray) -> np.ndarray:
    """Fraction of each site-year's calendar days that carry a value."""
    return np.isfinite(data).sum(axis=2) / days_in_year()[None, :]


def _mask_incomplete(values: np.ndarray, source: np.ndarray,
                     min_annual_frac: float) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    out[_completeness(source) < min_annual_frac] = np.nan
    return out


# ---------------------------------------------------------------
# Heat Wave Index
# ---------------------------------------------------------------
def rolling_mean_valid(data: np.ndarray, duration: int) -> np.ndarray:
    """Rolling mean within each year; a window needs all `duration` days."""
    if duration < 1 or duration > data.shape[2]:
        raise ValueError("duration must be between 1 and the number of day slots")
    n_windows = data.shape[2] - duration + 1
    total = np.zeros((*data.shape[:2], n_windows), dtype=np.float32)
    count = np.zeros((*data.shape[:2], n_windows), dtype=np.int16)
    for shift in range(duration):
        block = data[:, :, shift:shift + n_windows]
        total += np.nan_to_num(block, nan=0.0)
        count += np.isfinite(block)
    return np.where(count == duration, total / duration, np.nan)


def heat_wave_index(tmax: np.ndarray, tmin: np.ndarray,
                    duration: int = DURATION,
                    return_years: float = RETURN_YEARS,
                    min_annual_frac: float = MIN_ANNUAL_FRAC) -> np.ndarray:
    """Events per site-year: runs of overlapping above-threshold windows.

    The threshold is the ``100 * (1 - 1/return_years)`` percentile of the site's
    annual maximum rolling means -- the empirical return level. It is taken over
    every year's maximum while only complete years are scored, which is the
    parent project's convention.
    """
    if return_years <= 1:
        raise ValueError("return_years must be greater than 1")
    daily_mean = (tmax + tmin) / 2.0
    roll = rolling_mean_valid(daily_mean, duration)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        annual_max = np.nanmax(roll, axis=2)
        threshold = np.nanpercentile(
            annual_max, 100.0 * (1.0 - 1.0 / return_years), axis=1)

    exceeds = np.isfinite(roll) & (roll > threshold[:, None, None])
    event_start = exceeds.copy()
    event_start[:, :, 1:] &= ~exceeds[:, :, :-1]
    return _mask_incomplete(event_start.sum(axis=2), daily_mean, min_annual_frac)


# ---------------------------------------------------------------
# TXx -- annual maximum daily high
# ---------------------------------------------------------------
def annual_hottest(tmax: np.ndarray,
                   min_annual_frac: float = MIN_ANNUAL_FRAC) -> np.ndarray:
    """Hottest daily high per site-year, in degrees Celsius."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        hottest = np.nanmax(tmax, axis=2)
    return _mask_incomplete(hottest, tmax, min_annual_frac)


def c_to_f(values: np.ndarray) -> np.ndarray:
    return values * 9.0 / 5.0 + 32.0


# ---------------------------------------------------------------
# ETCCDI percentile indices: WSDI and TN90p
# ---------------------------------------------------------------
def calendar_thresholds(data: np.ndarray, pctile: int = PCTILE,
                        window: int = WINDOW, ref_start: int = REF_START,
                        ref_end: int = REF_END, chunk: int = 300) -> np.ndarray:
    """Per-site calendar-day percentile, shape (sites, 366).

    For each calendar day the percentile is taken over all reference-period
    values inside a centered `window`-day circular window, which enlarges the
    sample and smooths the seasonal cycle. Calendar days with no reference data
    become NaN and therefore never qualify as hot.
    """
    sites, _, doy = data.shape
    years = years_axis()
    in_reference = (years >= ref_start) & (years <= ref_end)
    if not in_reference.any():
        raise ValueError("reference period lies outside the study period")
    half = window // 2
    thresholds = np.empty((sites, doy), np.float32)
    for start in range(0, sites, chunk):
        stop = min(start + chunk, sites)
        reference = data[start:stop][:, in_reference, :]
        stacked = np.concatenate(
            [np.roll(reference, shift, axis=2) for shift in range(-half, half + 1)],
            axis=1,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            thresholds[start:stop] = np.nanpercentile(stacked, pctile, axis=1)
    return thresholds


def days_in_runs(hot: np.ndarray, min_run: int) -> np.ndarray:
    """Days lying in runs of at least `min_run` consecutive True, per row."""
    doy = hot.shape[1]
    position = np.broadcast_to(np.arange(doy, dtype=np.int32), hot.shape)
    last_false = np.maximum.accumulate(
        np.where(hot, np.int32(-1), position), axis=1)
    run_length = np.where(hot, position - last_false, 0)
    following = np.zeros_like(hot)
    following[:, :-1] = hot[:, 1:]
    at_run_end = hot & ~following
    return (run_length * (at_run_end & (run_length >= min_run))).sum(axis=1)


def _hot_day_mask(data: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    """Days strictly above the site's calendar-day threshold."""
    return np.isfinite(data) & (data > thresholds[:, None, :])


def wsdi(tmax: np.ndarray, min_annual_frac: float = MIN_ANNUAL_FRAC,
         min_run: int = MIN_RUN, **threshold_kwargs) -> np.ndarray:
    """Warm Spell Duration Index: days per site-year inside qualifying runs.

    Missing days both shorten and break warm spells, so incomplete site-years
    are dropped rather than scored low.
    """
    thresholds = calendar_thresholds(tmax, **threshold_kwargs)
    hot = _hot_day_mask(tmax, thresholds)
    counts = np.stack(
        [days_in_runs(hot[:, index, :], min_run) for index in range(hot.shape[1])],
        axis=1,
    )
    return _mask_incomplete(counts, tmax, min_annual_frac)


def tn90p(tmin: np.ndarray, min_annual_frac: float = MIN_ANNUAL_FRAC,
          **threshold_kwargs) -> np.ndarray:
    """Warm nights: percentage of each site-year's days above the TN90 threshold."""
    thresholds = calendar_thresholds(tmin, **threshold_kwargs)
    hot = _hot_day_mask(tmin, thresholds)
    valid = np.isfinite(tmin).sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        percent = np.where(valid > 0, 100.0 * hot.sum(axis=2) / valid, np.nan)
    return _mask_incomplete(percent, tmin, min_annual_frac)


# ---------------------------------------------------------------
# Spatial aggregation
# ---------------------------------------------------------------
def station_area_weights(latitudes: np.ndarray, longitudes: np.ndarray,
                         grid_deg: float = GRID_DEGREES) -> np.ndarray:
    """cos(lat) / stations-in-cell, normalised to mean one.

    Matches ``temperature_data.area_weights``: point networks are denser in some
    regions than others, so each occupied cell contributes in proportion to its
    area rather than its station count.
    """
    if grid_deg <= 0:
        raise ValueError("grid size must be positive")
    ilat = np.floor(np.asarray(latitudes) / grid_deg).astype(np.int64)
    ilon = np.floor(np.asarray(longitudes) / grid_deg).astype(np.int64)
    cell_area = np.cos(np.radians((ilat + 0.5) * grid_deg))
    _, inverse, counts = np.unique(
        ilat * np.int64(100000) + ilon, return_inverse=True, return_counts=True)
    weights = cell_area / counts[inverse]
    return weights * len(ilat) / weights.sum()


def grid_area_weights(latitudes: np.ndarray) -> np.ndarray:
    """cos(lat) for an already-regular grid, normalised to mean one."""
    weights = np.cos(np.radians(np.asarray(latitudes, dtype=float)))
    return weights * len(weights) / weights.sum()


def regional_mean(site_year: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """Annual mean over sites, renormalised over the sites valid in each year."""
    if weights is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmean(site_year, axis=0)
    valid = np.isfinite(site_year)
    column = np.asarray(weights, dtype=float)[:, None]
    numerator = np.where(valid, site_year * column, 0.0).sum(axis=0)
    denominator = np.where(valid, column, 0.0).sum(axis=0)
    return numerator / np.where(denominator > 0, denominator, np.nan)


# ---------------------------------------------------------------
# Metric registry -- the single place a metric is described
# ---------------------------------------------------------------
METRICS = {
    "hwi": {
        "label": "Heat Wave Index",
        "unit": "events per site per year",
        "axisLabel": "Events per site per year",
        "summary": (
            "Mean number of rare heat-wave events per site per year. An event is a "
            "continuous run of overlapping 4-day mean-temperature windows above that "
            "site's 1-in-10-year threshold."
        ),
        "reference": "Kunkel et al. (1999); CCSP (2008); EPA (2021)",
        "decimals": 3,
    },
    "txx": {
        "label": "TXx (annual hottest daily high)",
        "unit": "°F",
        "axisLabel": "Hottest daily high (°F)",
        "summary": (
            "The single hottest daily maximum temperature observed at each site each "
            "year, averaged across sites."
        ),
        "reference": "ETCCDI TXx",
        "decimals": 2,
    },
    "wsdi": {
        "label": "WSDI (warm spell duration)",
        "unit": "days per year",
        "axisLabel": "Days in warm spells",
        "summary": (
            "Days per year inside runs of at least 6 consecutive days whose daily high "
            "exceeds the site's calendar-day 90th percentile for 1961–1990."
        ),
        "reference": "ETCCDI WSDI",
        "decimals": 2,
    },
    "tn90p": {
        "label": "TN90p (warm nights)",
        "unit": "% of days",
        "axisLabel": "Warm nights (% of days)",
        "summary": (
            "Percentage of days each year whose daily low exceeds the site's "
            "calendar-day 90th percentile for 1961–1990."
        ),
        "reference": "ETCCDI TN90p",
        "decimals": 2,
    },
}
METRIC_ORDER = ["hwi", "txx", "wsdi", "tn90p"]


def compute_all(tmax: np.ndarray, tmin: np.ndarray) -> dict[str, np.ndarray]:
    """Every metric for one block of sites, as (sites, years) arrays.

    TXx is returned in Fahrenheit to match the parent project's figures; the
    other metrics are unit-free counts or percentages.
    """
    return {
        "hwi": heat_wave_index(tmax, tmin),
        "txx": c_to_f(annual_hottest(tmax)),
        "wsdi": wsdi(tmax),
        "tn90p": tn90p(tmin),
    }
