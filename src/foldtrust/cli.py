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
    analysis: str = typer.Argument(
        ...,
        help="Analysis type: 'all', 'layer4_shape', or individual layers"
    ),
    output: Path = typer.Option("benchmarks/outputs", "-o", "--output", help="Output directory"),
    n_windows: int = typer.Option(50, help="Number of genome control windows (layer4_shape)"),
    seed: int = typer.Option(42, help="Random seed"),
):
    """Run benchmark analyses to validate FoldTrust predictions."""
    output.mkdir(parents=True, exist_ok=True)

    # Layer 4: SHAPE analysis (standalone)
    if analysis == "layer4_shape":
        console.print("[bold cyan]Running Layer 4: SHAPE Agreement Analysis...[/bold cyan]")
        from foldtrust.benchmark.shape import run_layer4_shape_analysis
        
        cache_dir = Path("data/_cache")
        if not cache_dir.exists():
            console.print("[red]Error: data/_cache not found. Extract benchmark data bundle first.[/red]")
            raise typer.Exit(1)
        
        layer4_output = output / "layer4_shape"
        run_layer4_shape_analysis(
            cache_dir,
            layer4_output,
            n_genome_windows=n_windows,
            seed=seed
        )
        console.print(f"\n[green]✓ Layer 4 complete. Results in {layer4_output}[/green]")
        return
    
    # Comprehensive benchmark runner (all layers)
    if analysis == "all":
        try:
            from foldtrust.benchmark.runner import run_all_benchmarks
        except ImportError:
            console.print("[red]Error: benchmark runner not found[/red]")
            console.print("[yellow]For now, use 'layer4_shape' for SHAPE analysis[/yellow]")
            raise typer.Exit(1)
        
        console.print("[bold cyan]Running complete 6-layer benchmark analysis...[/bold cyan]")
        console.print("[dim]This will take several minutes...[/dim]\n")
        
        try:
            results = run_all_benchmarks(output, verbose=True)
            
            console.print("\n[bold green]✓ Benchmark complete![/bold green]")
            console.print(f"Results saved to: [bold]{output}[/bold]")
            console.print(f"Summary: {output / 'benchmark_summary.json'}")
            console.print(f"Figures: {output / 'figures'}")
            
        except Exception as e:
            console.print(f"[red]Error running benchmark: {e}[/red]")
            raise typer.Exit(1)
    else:
        console.print(f"[yellow]Supported: 'all' (complete benchmark) or 'layer4_shape' (SHAPE analysis)[/yellow]")
        console.print("[dim]Use: foldtrust benchmark layer4_shape[/dim]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
