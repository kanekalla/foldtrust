"""Test that benchmark all propagates errors correctly."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from foldtrust.benchmark.runner import run_all_benchmarks


def test_runner_collects_layer_errors():
    """Test that layer errors are collected in results."""
    with patch("foldtrust.benchmark.runner.run_layer1_tests") as mock_layer1:
        mock_layer1.side_effect = RuntimeError("Layer 1 failed")

        results = run_all_benchmarks(Path("/tmp/test_output"), verbose=False)

        assert "errors" in results
        assert len(results["errors"]) > 0
        assert any("Layer 1" in err for err in results["errors"])


def test_runner_collects_render_errors():
    """Test that render errors are collected."""
    # Mock all layers to succeed
    with (
        patch("foldtrust.benchmark.runner.run_layer1_tests", return_value={}),
        patch("foldtrust.benchmark.runner.run_layer2_benchmark", return_value={}),
        patch("foldtrust.benchmark.runner.run_layer3_calibration", return_value={}),
        patch("foldtrust.benchmark.runner.run_layer4_shape_analysis", return_value={}),
        patch("foldtrust.benchmark.runner.run_layer5_analysis", return_value={}),
        patch(
            "foldtrust.benchmark.layer5_figures.generate_all_layer5_figures",
            return_value=None,
        ),
        patch("foldtrust.benchmark.runner.subprocess.run") as mock_subprocess,
    ):
        # Make render script fail
        result = MagicMock()
        result.returncode = 1
        result.stderr = "Missing file"
        mock_subprocess.return_value = result

        results = run_all_benchmarks(Path("/tmp/test_output"), verbose=False)

        assert "errors" in results
        assert any("Table rendering" in err for err in results["errors"])
