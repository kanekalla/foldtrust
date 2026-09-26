"""Test that benchmark all propagates errors correctly."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import typer.testing

from foldtrust.benchmark.runner import run_all_benchmarks
from foldtrust.cli import app


def test_runner_collects_layer_errors(tmp_path):
    """Test that layer errors are collected in results."""
    with (
        patch("foldtrust.benchmark.runner.run_layer1_tests") as mock_layer1,
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
        mock_layer1.side_effect = RuntimeError("Layer 1 failed")
        
        # Mock subprocess to return success for render scripts
        result = MagicMock()
        result.returncode = 0
        mock_subprocess.return_value = result

        results = run_all_benchmarks(tmp_path, verbose=False)

        assert "errors" in results
        assert len(results["errors"]) > 0
        assert any("Layer 1" in err for err in results["errors"])


def test_runner_collects_render_errors(tmp_path):
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
        result.stdout = ""
        result.stderr = "Missing file"
        mock_subprocess.return_value = result

        results = run_all_benchmarks(tmp_path, verbose=False)

        assert "errors" in results
        assert any("Table rendering" in err for err in results["errors"])


def test_cli_all_exits_nonzero_on_error(tmp_path):
    """Test that CLI exits 1 when a layer fails."""
    runner = typer.testing.CliRunner()
    
    with (
        patch("foldtrust.benchmark.runner.run_layer1_tests") as mock_layer1,
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
        mock_layer1.side_effect = RuntimeError("Layer 1 failed")
        result_obj = MagicMock()
        result_obj.returncode = 0
        mock_subprocess.return_value = result_obj
        
        result = runner.invoke(app, ["benchmark", "all", "-o", str(tmp_path)])
        
        assert result.exit_code == 1
        assert "finished with errors" in result.stdout


def test_cli_all_exits_zero_on_success(tmp_path):
    """Test that CLI exits 0 when all layers succeed."""
    runner = typer.testing.CliRunner()
    
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
        result_obj = MagicMock()
        result_obj.returncode = 0
        mock_subprocess.return_value = result_obj
        
        result = runner.invoke(app, ["benchmark", "all", "-o", str(tmp_path)])
        
        assert result.exit_code == 0
        assert "Benchmark complete" in result.stdout
