"""Report generation for HTML and Markdown formats."""

from pathlib import Path
from typing import Dict

from foldtrust.utils import format_structure_ascii


def generate_html_report(result: Dict, output_path: Path):
    """Generate HTML report from analysis results."""

    metadata = result.get("metadata", {})
    disease = metadata.get("disease", "N/A")
    gene = metadata.get("gene", "N/A")
    teaching_point = metadata.get("teaching_point", "")
    references = metadata.get("references", [])

    stems_html = ""
    for stem in result["stems"]:
        flag_class = stem["flag"]
        stems_html += f"""
        <tr class="{flag_class}">
            <td>{stem['id']}</td>
            <td>{stem['positions']}</td>
            <td>{stem['length']}</td>
            <td>{stem['mean_prob']:.3f}</td>
            <td><span class="flag-{flag_class}">{stem['flag'].upper()}</span></td>
        </tr>
        """

    refs_html = ""
    if references:
        refs_html = "<h3>References</h3><ul>"
        for ref in references:
            doi = ref.get("doi", "")
            pmid = ref.get("pmid", "")
            title = ref.get("title", "")
            ref_text = title if title else (doi if doi else pmid)

            if doi:
                refs_html += (
                    f'<li><a href="https://doi.org/{doi}" target="_blank">{ref_text}</a></li>'
                )
            elif pmid:
                refs_html += f'<li><a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">{ref_text}</a></li>'
            else:
                refs_html += f"<li>{ref_text}</li>"
        refs_html += "</ul>"

    structure_ascii = format_structure_ascii(result["sequence"], result["structure"], width=70)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FoldTrust Report - {result['header']}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
            color: #333;
        }}
        .container {{
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #2c3e50;
            margin-top: 30px;
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #34495e;
            margin-top: 20px;
        }}
        .metadata {{
            background: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .metadata-item {{
            margin: 8px 0;
        }}
        .metadata-label {{
            font-weight: bold;
            color: #2c3e50;
        }}
        .verdict {{
            background: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            font-size: 1.1em;
        }}
        .verdict.trust {{
            background: #d4edda;
            border-left-color: #28a745;
        }}
        .verdict.redesign {{
            background: #f8d7da;
            border-left-color: #dc3545;
        }}
        .verdict.probing {{
            background: #fff3cd;
            border-left-color: #ffc107;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: #34495e;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{
            background: #f5f5f5;
        }}
        .flag-firm {{
            background: #28a745;
            color: white;
            padding: 3px 8px;
            border-radius: 3px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        .flag-soft {{
            background: #ffc107;
            color: #333;
            padding: 3px 8px;
            border-radius: 3px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        .flag-floppy {{
            background: #dc3545;
            color: white;
            padding: 3px 8px;
            border-radius: 3px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        .structure {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            overflow-x: auto;
            white-space: pre;
            line-height: 1.4;
        }}
        .heatmap {{
            text-align: center;
            margin: 20px 0;
        }}
        .heatmap img {{
            max-width: 100%;
            border-radius: 5px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .methods {{
            background: #fff9e6;
            border-left: 4px solid #ff9800;
            padding: 15px;
            margin: 20px 0;
            font-style: italic;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 2px solid #ecf0f1;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
        .stem-counts {{
            display: flex;
            gap: 20px;
            margin: 20px 0;
        }}
        .stem-count {{
            flex: 1;
            padding: 15px;
            border-radius: 5px;
            text-align: center;
        }}
        .stem-count.firm {{
            background: #d4edda;
        }}
        .stem-count.soft {{
            background: #fff3cd;
        }}
        .stem-count.floppy {{
            background: #f8d7da;
        }}
        .stem-count-value {{
            font-size: 2em;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .stem-count-label {{
            font-size: 0.9em;
            text-transform: uppercase;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>FoldTrust Structure Reliability Report</h1>

        <div class="metadata">
            <div class="metadata-item"><span class="metadata-label">Sequence:</span> {result['header']}</div>
            <div class="metadata-item"><span class="metadata-label">Length:</span> {result['length']} nt</div>
            <div class="metadata-item"><span class="metadata-label">GC Content:</span> {result['gc_content']:.1f}%</div>
            <div class="metadata-item"><span class="metadata-label">MFE:</span> {result['mfe_energy']:.2f} kcal/mol</div>
            {f'<div class="metadata-item"><span class="metadata-label">Disease:</span> {disease}</div>' if disease != "N/A" else ''}
            {f'<div class="metadata-item"><span class="metadata-label">Gene:</span> {gene}</div>' if gene != "N/A" else ''}
        </div>

        <div class="verdict {'trust' if 'TRUST' in result['verdict'] else 'redesign' if 'REDESIGN' in result['verdict'] else 'probing'}">
            <strong>Verdict:</strong> {result['verdict']}
        </div>

        {f'<div class="metadata"><strong>Teaching Point:</strong> {teaching_point}</div>' if teaching_point else ''}

        <h2>Stem Reliability Summary</h2>

        <div class="stem-counts">
            <div class="stem-count firm">
                <div class="stem-count-value">{result['stem_counts']['firm']}</div>
                <div class="stem-count-label">Firm Stems</div>
            </div>
            <div class="stem-count soft">
                <div class="stem-count-value">{result['stem_counts']['soft']}</div>
                <div class="stem-count-label">Soft Stems</div>
            </div>
            <div class="stem-count floppy">
                <div class="stem-count-value">{result['stem_counts']['floppy']}</div>
                <div class="stem-count-label">Floppy Stems</div>
            </div>
        </div>

        <h2>Stem Analysis</h2>

        <table>
            <thead>
                <tr>
                    <th>Stem ID</th>
                    <th>Positions</th>
                    <th>Length (bp)</th>
                    <th>Mean P(pair)</th>
                    <th>Flag</th>
                </tr>
            </thead>
            <tbody>
                {stems_html}
            </tbody>
        </table>

        <div class="methods">
            <strong>Legend:</strong>
            <ul style="margin: 10px 0;">
                <li><strong>FIRM</strong> (mean P ≥ 0.85): High confidence — stem is well-supported by ensemble</li>
                <li><strong>SOFT</strong> (0.5 ≤ mean P < 0.85): Moderate confidence — alternative structures possible</li>
                <li><strong>FLOPPY</strong> (mean P < 0.5): Low confidence — MFE stem not reliable in ensemble</li>
            </ul>
        </div>

        <h2>Base-Pair Probability Matrix</h2>

        <div class="heatmap">
            <img src="pair_probabilities.png" alt="Base-pair probability heatmap">
            <p><em>Upper triangle shows probability of base-pairing between positions i and j</em></p>
        </div>

        <h2>MFE Structure</h2>

        <div class="structure">{structure_ascii}</div>

        <h2>Methods</h2>

        <div class="methods">
            <p><strong>Key principle:</strong> Ask for base-pair probabilities (or the full ensemble),
            not only the minimum free energy structure.</p>

            <p>FoldTrust uses ViennaRNA's partition function to compute the Boltzmann ensemble
            of all possible secondary structures and their base-pairing probabilities. The MFE
            structure is parsed into stems, and each stem's reliability is assessed by averaging
            the pair probabilities of its constituent base pairs.</p>

            <p><strong>Limitations:</strong> Thermodynamic models operate under simplified assumptions
            (37°C in 1M NaCl). In-cell conditions, co-transcriptional folding, RNA-binding proteins,
            and post-transcriptional modifications can alter structure. Experimental probing (SHAPE, DMS)
            and functional assays remain essential for validation.</p>
        </div>

        {refs_html}

        <div class="footer">
            <p>Generated by <strong>FoldTrust v0.1.0</strong></p>
            <p>Personal portfolio project by Kishore Anekalla</p>
            <p>Energy model: ViennaRNA default parameters (Turner 2004)</p>
        </div>
    </div>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html)


def generate_markdown_report(result: Dict, output_path: Path):
    """Generate Markdown report from analysis results."""

    metadata = result.get("metadata", {})
    disease = metadata.get("disease", "N/A")
    gene = metadata.get("gene", "N/A")
    teaching_point = metadata.get("teaching_point", "")
    references = metadata.get("references", [])

    stems_table = "| Stem ID | Positions | Length (bp) | Mean P(pair) | Flag |\n"
    stems_table += "|---------|-----------|-------------|--------------|------|\n"

    for stem in result["stems"]:
        stems_table += f"| {stem['id']} | {stem['positions']} | {stem['length']} | {stem['mean_prob']:.3f} | **{stem['flag'].upper()}** |\n"

    refs_section = ""
    if references:
        refs_section = "\n## References\n\n"
        for ref in references:
            doi = ref.get("doi", "")
            pmid = ref.get("pmid", "")
            title = ref.get("title", "")
            ref_text = title if title else (doi if doi else pmid)

            if doi:
                refs_section += f"- [{ref_text}](https://doi.org/{doi})\n"
            elif pmid:
                refs_section += f"- [{ref_text}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)\n"
            else:
                refs_section += f"- {ref_text}\n"

    structure_ascii = format_structure_ascii(result["sequence"], result["structure"], width=70)

    teaching_section = ""
    if teaching_point:
        teaching_section = f"\n### Teaching Point\n\n{teaching_point}\n"

    disease_line = f"- **Disease:** {disease}\n" if disease != "N/A" else ""
    gene_line = f"- **Gene:** {gene}\n" if gene != "N/A" else ""

    markdown = f"""# FoldTrust Structure Reliability Report

## Sequence Information

- **Sequence:** {result['header']}
- **Length:** {result['length']} nt
- **GC Content:** {result['gc_content']:.1f}%
- **MFE:** {result['mfe_energy']:.2f} kcal/mol
{disease_line}{gene_line}
## Verdict

**{result['verdict']}**
{teaching_section}
## Stem Reliability Summary

- **Firm stems:** {result['stem_counts']['firm']} (mean P ≥ 0.85)
- **Soft stems:** {result['stem_counts']['soft']} (0.5 ≤ mean P < 0.85)
- **Floppy stems:** {result['stem_counts']['floppy']} (mean P < 0.5)

## Stem Analysis

{stems_table}

**Legend:**
- **FIRM** (mean P ≥ 0.85): High confidence — stem is well-supported by ensemble
- **SOFT** (0.5 ≤ mean P < 0.85): Moderate confidence — alternative structures possible
- **FLOPPY** (mean P < 0.5): Low confidence — MFE stem not reliable in ensemble

## Base-Pair Probability Matrix

![Base-pair probability heatmap](pair_probabilities.png)

*Upper triangle shows probability of base-pairing between positions i and j*

## MFE Structure

```
{structure_ascii}
```

## Methods

**Key principle:** Ask for base-pair probabilities (or the full ensemble), not only the minimum free energy structure.

FoldTrust uses ViennaRNA's partition function to compute the Boltzmann ensemble of all possible secondary structures and their base-pairing probabilities. The MFE structure is parsed into stems, and each stem's reliability is assessed by averaging the pair probabilities of its constituent base pairs.

**Limitations:** Thermodynamic models operate under simplified assumptions (37°C in 1M NaCl). In-cell conditions, co-transcriptional folding, RNA-binding proteins, and post-transcriptional modifications can alter structure. Experimental probing (SHAPE, DMS) and functional assays remain essential for validation.

{refs_section}

---

*Generated by FoldTrust v0.1.0 — Personal portfolio project by Kishore Anekalla*

*Energy model: ViennaRNA default parameters (Turner 2004)*
"""

    with open(output_path, "w") as f:
        f.write(markdown)
