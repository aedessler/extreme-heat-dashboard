"""Checks on the generated dashboard series.

Two things matter here. First, the CONUS cells must be the same numbers the
repository publishes elsewhere, since the whole point of the page is that it
shows the project's own results. Second, every field the page reads must exist,
because a static page has no server to fall back on when a key is missing.

Run from the dashboard directory:  python -m unittest discover tests
"""

import csv
import itertools
import json
import unittest
from pathlib import Path

DASHBOARD_DIR = Path(__file__).resolve().parents[1]
PROJECT_DIR = DASHBOARD_DIR.parent
SERIES_JSON = DASHBOARD_DIR / "dashboard-series.json"
SERIES_JS = DASHBOARD_DIR / "dashboard-series.js"
DATA_DIR = PROJECT_DIR / "data"

# series key -> (published csv, column, decimals stored in that csv). These are
# the cells the parent project publishes, so they are the ones that can be
# checked against something. The tolerance is one unit in the file's last stored
# place: anything larger is a real disagreement, anything smaller is its rounding.
PUBLISHED = {
    "ghcn_conus_hwi": ("heat_wave_index_raw_adj_wtd.csv", "weighted_heat_wave_index", 6),
    "ghcn_conus_txx": ("annual_hottest_daily_high_raw_adj_wtd.csv",
                       "weighted_annual_hottest_daily_high_f", 4),
    "berkeley_conus_hwi": ("heat_wave_index_berkeley.csv", "berkeley_heat_wave_index", 6),
}


def read_published(name, column):
    with (DATA_DIR / name).open(encoding="utf-8") as handle:
        return {int(row["year"]): float(row[column])
                for row in csv.DictReader(handle) if row[column]}


@unittest.skipUnless(SERIES_JSON.exists(), f"{SERIES_JSON} has not been generated")
class DashboardDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads(SERIES_JSON.read_text(encoding="utf-8"))

    def cells(self):
        return [tuple(cell.split("_")) for cell in self.data["availableCells"]]

    # -- structure -----------------------------------------------------
    def test_every_available_cell_has_all_metrics(self):
        for dataset, region in self.cells():
            for metric in self.data["metricOrder"]:
                with self.subTest(cell=f"{dataset}/{region}/{metric}"):
                    self.assertIn(f"{dataset}_{region}_{metric}", self.data["series"])

    def test_no_unexpected_series(self):
        expected = {
            f"{dataset}_{region}_{metric}"
            for dataset, region in self.cells()
            for metric in self.data["metricOrder"]
        }
        self.assertEqual(set(self.data["series"]), expected)

    def test_excluded_cells_are_absent_and_explained(self):
        # GHCN cannot support a hemispheric estimate, so it must not be built
        # for that region and the page must be able to say why.
        every = set(itertools.product(self.data["datasetOrder"], self.data["regionOrder"]))
        for dataset, region in every - set(self.cells()):
            with self.subTest(cell=f"{dataset}/{region}"):
                self.assertIn(f"{dataset}_{region}", self.data["excludedCells"])
                self.assertTrue(self.data["excludedCells"][f"{dataset}_{region}"].strip())
                for metric in self.data["metricOrder"]:
                    self.assertNotIn(f"{dataset}_{region}_{metric}", self.data["series"])

    def test_every_region_has_at_least_one_dataset(self):
        for region in self.data["regionOrder"]:
            with self.subTest(region=region):
                self.assertTrue(any(r == region for _, r in self.cells()))

    def test_the_page_can_render_every_field_it_reads(self):
        for series_id, series in self.data["series"].items():
            with self.subTest(series=series_id):
                for field in ("dataset", "region", "metric", "adjustment", "weighting",
                              "source", "sourceUrl", "coverage", "reliableWindow", "points"):
                    self.assertIn(field, series)
                self.assertIn(series["dataset"], self.data["datasets"])
                self.assertIn(series["region"], self.data["regions"])
                self.assertIn(series["metric"], self.data["metrics"])
                coverage = series["coverage"]
                for field in ("sites", "siteLabel", "occupiedCells", "domain", "areaFraction",
                              "sitesByLatitudeBand", "sitesByLongitudeBand"):
                    self.assertIn(field, coverage)
                self.assertGreater(coverage["sites"], 0)
                self.assertTrue(0.0 < coverage["areaFraction"] <= 1.0)

    def test_metric_descriptions_are_complete(self):
        for metric in self.data["metricOrder"]:
            with self.subTest(metric=metric):
                meta = self.data["metrics"][metric]
                for field in ("label", "unit", "axisLabel", "summary", "reference", "decimals"):
                    self.assertIn(field, meta)

    # -- series shape --------------------------------------------------
    def test_years_are_sorted_unique_and_inside_the_study_period(self):
        method = self.data["method"]
        for series_id, series in self.data["series"].items():
            years = [point["x"] for point in series["points"]]
            with self.subTest(series=series_id):
                self.assertEqual(years, sorted(years))
                self.assertEqual(len(years), len(set(years)))
                self.assertGreaterEqual(min(years), method["startYear"])
                self.assertLessEqual(max(years), method["endYear"])
                self.assertEqual(series["reliableWindow"], [years[0], years[-1]])
                self.assertGreaterEqual(len(years), 4)

    def test_values_are_finite_and_physically_plausible(self):
        limits = {"hwi": (0.0, 400.0), "txx": (-40.0, 150.0),
                  "wsdi": (0.0, 366.0), "tn90p": (0.0, 100.0)}
        for series_id, series in self.data["series"].items():
            low, high = limits[series["metric"]]
            for point in series["points"]:
                with self.subTest(series=series_id, year=point["x"]):
                    self.assertIsInstance(point["y"], (int, float))
                    self.assertTrue(low <= point["y"] <= high,
                                    f"{point['y']} outside [{low}, {high}]")

    def test_all_metrics_in_a_cell_share_one_site_set(self):
        for dataset, region in self.cells():
            coverages = [self.data["series"][f"{dataset}_{region}_{metric}"]["coverage"]
                         for metric in self.data["metricOrder"]]
            with self.subTest(cell=f"{dataset}/{region}"):
                for coverage in coverages[1:]:
                    self.assertEqual(coverage, coverages[0])

    def test_conus_ghcn_is_the_adjusted_station_set(self):
        self.assertEqual(self.data["series"]["ghcn_conus_hwi"]["adjustment"], "FLs.52j-adjusted")

    def test_hemispheric_cells_cover_more_ground_than_conus(self):
        for dataset, region in self.cells():
            if region != "nh":
                continue
            conus = self.data["series"][f"{dataset}_conus_hwi"]["coverage"]
            hemisphere = self.data["series"][f"{dataset}_nh_hwi"]["coverage"]
            with self.subTest(dataset=dataset):
                self.assertGreater(hemisphere["occupiedCells"], conus["occupiedCells"])

    def test_gridded_cells_reach_all_their_land(self):
        # Berkeley Earth is a grid, so it should occupy every land cell in its
        # domain. A station network will not, which is the point of the field.
        for dataset, region in self.cells():
            coverage = self.data["series"][f"{dataset}_{region}_hwi"]["coverage"]
            if "landCellsSampled" not in coverage:
                continue
            with self.subTest(cell=f"{dataset}/{region}"):
                self.assertLessEqual(coverage["landCellsSampled"], coverage["landCellsTotal"])
                if dataset == "berkeley":
                    self.assertEqual(coverage["landCellsSampled"], coverage["landCellsTotal"])

    # -- agreement with the parent project's published results -----------
    def test_detection_is_mann_kendall_significance_alone(self):
        # A 10% likelihood threshold. Detection is significance alone, with no
        # magnitude-versus-variability check, so no magnitude parameter should
        # be present to imply otherwise.
        self.assertEqual(self.data["detection"]["pThreshold"], 0.10)
        self.assertEqual(self.data["detection"]["test"], "Mann-Kendall")
        self.assertNotIn("magnitudeFraction", self.data["detection"])

    def test_published_cells_reproduce_their_csvs(self):
        for series_id, (name, column, decimals) in PUBLISHED.items():
            path = DATA_DIR / name
            with self.subTest(series=series_id):
                if not path.exists():
                    self.skipTest(f"{path} is not present")
                published = read_published(name, column)
                points = self.data["series"][series_id]["points"]
                self.assertEqual(len(points), len(published))
                tolerance = 10.0 ** -decimals
                for point in points:
                    self.assertAlmostEqual(point["y"], published[point["x"]], delta=tolerance,
                                           msg=f"{series_id} year {point['x']}")

    @unittest.skipUnless(SERIES_JS.exists(), f"{SERIES_JS} has not been generated")
    def test_javascript_companion_carries_the_same_payload(self):
        text = SERIES_JS.read_text(encoding="utf-8")
        prefix = "window.HEAT_WAVE_DATA = "
        start = text.index(prefix) + len(prefix)
        self.assertEqual(json.loads(text[start:text.rindex(";")]), self.data)


if __name__ == "__main__":
    unittest.main()
