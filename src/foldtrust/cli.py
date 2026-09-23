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


if __name__ == "__main__":
    app()
