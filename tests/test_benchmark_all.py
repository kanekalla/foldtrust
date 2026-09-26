"""Test benchmark all command."""

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from foldtrust.benchmark.runner import run_all_benchmarks


def test_benchmark_all_dispatches_layers(tmp_path, monkeypatch):
    """Test that benchmark all dispatches all layers correctly."""
    # Mock all layer functions
    with (
        patch("foldtrust.benchmark.runner.run_layer1_tests") as mock_layer1,
        patch("foldtrust.benchmark.runner.run_layer2_benchmark") as mock_layer2,
        patch("foldtrust.benchmark.runner.run_layer3_calibration") as mock_layer3,
        patch("foldtrust.benchmark.runner.run_layer4_shape_analysis") as mock_layer4,
        patch("foldtrust.benchmark.runner.run_layer5_analysis") as mock_layer5,
        patch("foldtrust.benchmark.runner.subprocess.run") as mock_subprocess,
        patch(
            "foldtrust.benchmark.layer5_figures.generate_all_layer5_figures"
        ) as mock_layer5_figures,
    ):
        # Setup return values
        mock_layer1.return_value = {"tests_passed": 5, "tests_total": 5}
        mock_layer2.return_value = {"n_structures": 200}
        mock_layer3.return_value = {"ece": 0.05, "auroc": 0.88}
        mock_layer4.return_value = {
            "fse_start": 13462,
            "fse_end": 13542,
            "datasets": {"a": {}, "b": {}},
        }
        mock_layer5.return_value = [
            {"case": "case1"},
            {"case": "case2"},
            {"case": "case3"},
            {"case": "case4"},
            {"case": "case5"},
        ]

        # Mock cache and cases directories
        cache_dir = tmp_path / "data" / "_cache"
        cache_dir.mkdir(parents=True)
        cases_dir = tmp_path / "data" / "cases"
        cases_dir.mkdir(parents=True)

        monkeypatch.chdir(tmp_path)

        output_dir = tmp_path / "output"

        # Run benchmark all
        results = run_all_benchmarks(output_dir, verbose=False)

        # Verify all layers were called
        assert mock_layer1.called
        assert mock_layer2.called
        assert mock_layer3.called
        assert mock_layer4.called
        assert mock_layer5.called
        assert mock_layer5_figures.called

        # Verify results structure
        assert "layers" in results
        assert "layer1" in results["layers"]
        assert "layer2" in results["layers"]
        assert "layer3" in results["layers"]
        assert "layer4_shape" in results["layers"]
        assert "layer5" in results["layers"]

        # Verify summary file was written
        summary_file = output_dir / "benchmark_summary.json"
        assert summary_file.exists()

        with open(summary_file) as f:
            summary = json.load(f)
            assert summary["layers"]["layer1"]["tests_passed"] == 5
            assert summary["layers"]["layer2"]["n_structures"] == 200
            assert summary["layers"]["layer3"]["ece"] == 0.05
            assert summary["layers"]["layer4_shape"]["n_datasets"] == 2
            assert summary["layers"]["layer5"]["n_cases"] == 5
