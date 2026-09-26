"""FoldTrust CLI interface."""

from pathlib import Path

import typer
from rich.console import Console

from foldtrust.core import process_sequence
from foldtrust.utils import find_case_directories

app = typer.Typer(help="FoldTrust: RNA structure reliability reports")
console = Console()


@app.command()
def report(
    sequence: Path = typer.Argument(..., help="Path to FASTA file"),
    output: Path = typer.Option("out", "-o", "--output", help="Output directory"),
    html: bool = typer.Option(True, help="Generate HTML report"),
    markdown: bool = typer.Option(True, help="Generate Markdown report"),
):
    """Generate structure reliability report for a single sequence."""
    if not sequence.exists():
        console.print(f"[red]Error: {sequence} not found[/red]")
        raise typer.Exit(1)

    console.print(f"[bold]Processing {sequence.name}...[/bold]")

    output.mkdir(parents=True, exist_ok=True)

    try:
        result = process_sequence(sequence, output, generate_html=html, generate_markdown=markdown)
        console.print(f"[green]✓ Report generated in {output}[/green]")
        console.print(f"  Verdict: [bold]{result['verdict']}[/bold]")
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def batch(
    cases_dir: Path = typer.Argument(..., help="Directory containing case folders"),
    output: Path = typer.Option("out/batch", "-o", "--output", help="Output directory"),
):
    """Process all cases in a directory."""
    if not cases_dir.exists():
        console.print(f"[red]Error: {cases_dir} not found[/red]")
        raise typer.Exit(1)

    case_dirs = find_case_directories(cases_dir)

    if not case_dirs:
        console.print(f"[yellow]No case directories found in {cases_dir}[/yellow]")
        raise typer.Exit(1)

    console.print(f"[bold]Processing {len(case_dirs)} cases...[/bold]")

    output.mkdir(parents=True, exist_ok=True)

    success_count = 0
    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"

        if not sequence_file.exists():
            console.print(f"[yellow]⚠ Skipping {case_name}: no sequence.fa[/yellow]")
            continue

        case_output = output / case_name
        console.print(f"  Processing {case_name}...")

        try:
            process_sequence(sequence_file, case_output)
            console.print(f"  [green]✓ {case_name}[/green]")
            success_count += 1
        except Exception as e:
            console.print(f"  [red]✗ {case_name}: {e}[/red]")

    console.print(f"\n[bold]Completed: {success_count}/{len(case_dirs)} cases[/bold]")


@app.command()
def demo():
    """Run demo on all MVP cases in data/cases."""
    cases_dir = Path("data/cases")
    output_dir = Path("examples/out")

    if not cases_dir.exists():
        console.print(f"[red]Error: {cases_dir} not found[/red]")
        console.print("Run this command from the repository root.")
        raise typer.Exit(1)

    console.print("[bold cyan]FoldTrust Demo: MVP Disease Cases[/bold cyan]\n")

    case_dirs = find_case_directories(cases_dir)

    if not case_dirs:
        console.print("[yellow]No cases found in data/cases[/yellow]")
        raise typer.Exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    for case_dir in case_dirs:
        case_name = case_dir.name
        sequence_file = case_dir / "sequence.fa"

        if not sequence_file.exists():
            continue

        case_output = output_dir / case_name
        console.print(f"[bold]{case_name}[/bold]")

        try:
            result = process_sequence(sequence_file, case_output)
            console.print(f"  Verdict: [bold]{result['verdict']}[/bold]")
            console.print(f"  Output: {case_output}")
        except Exception as e:
            console.print(f"  [red]Error: {e}[/red]")

        console.print()

    console.print("[green]✓ Demo complete[/green]")
    console.print(f"Reports available in: {output_dir}")


@app.command()
def benchmark(
    layer: str = typer.Argument(
        ..., help="Layer to run: 'layer1', 'layer2', 'layer3', 'layer5', 'all', or 'layer4_shape'"
    ),
    output: Path = typer.Option("benchmarks/outputs", "-o", "--output", help="Output directory"),
    full: bool = typer.Option(False, "--full", help="Run full benchmark (no subsampling)"),
):
    """Run benchmark analyses to validate FoldTrust predictions."""
    output.mkdir(parents=True, exist_ok=True)

    # Layer 1: Scoring correctness
    if layer == "layer1":
        console.print("[bold cyan]Running Layer 1: Scoring Correctness Tests...[/bold cyan]")
        from foldtrust.benchmark.scoring import run_layer1_tests

        layer1_output = output / "layer1"
        layer1_output.mkdir(parents=True, exist_ok=True)

        results = run_layer1_tests(layer1_output)
        console.print(
            f"\n[green]✓ Layer 1 complete: {results['tests_passed']}/{results['tests_total']} tests passed[/green]"
        )
        console.print(f"Results in {layer1_output}")
        return

    # Layer 2: Structure accuracy
    if layer == "layer2":
        console.print("[bold cyan]Running Layer 2: Structure Accuracy Benchmark...[/bold cyan]")
        from foldtrust.benchmark.layer2_accuracy import run_layer2_benchmark

        cache_dir = Path("data/_cache")
        if not cache_dir.exists():
            console.print(
                "[red]Error: data/_cache not found. Extract benchmark data bundle first.[/red]"
            )
            raise typer.Exit(1)

        layer2_output = output / "layer2"
        layer2_output.mkdir(parents=True, exist_ok=True)

        results = run_layer2_benchmark(
            cache_dir=cache_dir,
            output_dir=layer2_output,
            use_full=full,
            sample_size=None if full else 200,
        )
        console.print(
            f"\n[green]✓ Layer 2 complete: {results['n_structures']} structures analyzed[/green]"
        )
        console.print(f"Results in {layer2_output}")
        return

    # Layer 3: Calibration analysis
    if layer == "layer3":
        console.print("[bold cyan]Running Layer 3: Calibration Analysis...[/bold cyan]")
        from foldtrust.benchmark.layer3_calibration import run_layer3_calibration

        cache_dir = Path("data/_cache")
        layer2_output = output / "layer2"

        if not cache_dir.exists():
            console.print(
                "[red]Error: data/_cache not found. Extract benchmark data bundle first.[/red]"
            )
            raise typer.Exit(1)

        if not (layer2_output / "sample_ids.csv").exists():
            console.print(
                "[yellow]Warning: sample_ids.csv not found. Run layer2 first for consistent sampling.[/yellow]"
            )

        layer3_output = output / "layer3"
        layer3_output.mkdir(parents=True, exist_ok=True)

        results = run_layer3_calibration(
            cache_dir=cache_dir,
            output_dir=layer3_output,
            use_full=full,
            sample_size=None if full else 200,
        )
        console.print(
            f"\n[green]✓ Layer 3 complete: ECE={results['ece']:.4f}, AUROC={results['auroc']:.4f}[/green]"
        )
        console.print(f"Results in {layer3_output}")
        return

    # Layer 5: Robustness analysis
    if layer == "layer5":
        console.print("[bold cyan]Running Layer 5: Robustness Analysis...[/bold cyan]")
        from foldtrust.benchmark.layer5_robustness import run_layer5_analysis

        cases_dir = Path("data/cases")
        if not cases_dir.exists():
            console.print("[red]Error: data/cases not found[/red]")
            raise typer.Exit(1)

        cache_dir = Path("data/_cache")
        layer5_output = output / "layer5"
        layer5_output.mkdir(parents=True, exist_ok=True)

        results = run_layer5_analysis(cases_dir, layer5_output, cache_dir)
        console.print(f"\n[green]✓ Layer 5 complete: {len(results)} cases analyzed[/green]")
        console.print(f"Results in {layer5_output}")

        # Generate figures
        console.print("\nGenerating figures...")
        from foldtrust.benchmark.layer5_figures import generate_all_layer5_figures

        figures_dir = output / "figures"
        generate_all_layer5_figures(layer5_output, figures_dir)

        return

    # Layer 4: SHAPE analysis (standalone)
    if layer == "layer4_shape":
        console.print("[bold cyan]Running Layer 4: SHAPE Agreement Analysis...[/bold cyan]")
        from foldtrust.benchmark.shape import run_layer4_shape_analysis

        cache_dir = Path("data/_cache")
        if not cache_dir.exists():
            console.print(
                "[red]Error: data/_cache not found. Extract benchmark data bundle first.[/red]"
            )
            raise typer.Exit(1)

        layer4_output = output / "layer4_shape"
        run_layer4_shape_analysis(cache_dir, layer4_output)
        console.print(f"\n[green]✓ Layer 4 complete. Results in {layer4_output}[/green]")
        return

    # Comprehensive benchmark runner (all layers)
    if layer == "all":
        try:
            from foldtrust.benchmark.runner import run_all_benchmarks
        except ImportError:
            console.print("[red]Error: benchmark runner not found[/red]")
            console.print("[yellow]For now, use 'layer4_shape' for SHAPE analysis[/yellow]")
            raise typer.Exit(1)

        console.print("[bold cyan]Running complete 6-layer benchmark analysis...[/bold cyan]")
        console.print("[dim]This will take several minutes...[/dim]\n")

        try:
            results = run_all_benchmarks(output, verbose=True, use_full=full)

            if results.get("errors"):
                console.print("\n[bold red]✗ Benchmark finished with errors:[/bold red]")
                for error in results["errors"]:
                    console.print(f"  - {error}")
                raise typer.Exit(1)

            console.print("\n[bold green]✓ Benchmark complete![/bold green]")
            console.print(f"Results saved to: [bold]{output}[/bold]")
            console.print(f"Summary: {output / 'benchmark_summary.json'}")
            console.print(f"Figures: {output / 'figures'}")

        except Exception as e:
            console.print(f"[red]Error running benchmark: {e}[/red]")
            raise typer.Exit(1)
    else:
        console.print(
            "[yellow]Supported layers: 'layer1', 'layer2', 'layer3', 'layer5', 'all', or 'layer4_shape'[/yellow]"
        )
        console.print("[dim]Examples:[/dim]")
        console.print("[dim]  foldtrust benchmark layer1[/dim]")
        console.print("[dim]  foldtrust benchmark layer2 --full[/dim]")
        console.print("[dim]  foldtrust benchmark layer3[/dim]")
        console.print("[dim]  foldtrust benchmark layer5[/dim]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
