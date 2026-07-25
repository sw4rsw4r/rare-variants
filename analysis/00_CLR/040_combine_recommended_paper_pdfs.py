"""Combine the current recommended manuscript component PDFs into stable handoff files."""

from __future__ import annotations

import sys
from pathlib import Path

from pypdf import PdfReader, PdfWriter


def combine(output: Path, inputs: list[str]) -> None:
    writer = PdfWriter()
    for name in inputs:
        source = output.parent / name
        if not source.is_file():
            raise FileNotFoundError(f"Missing PDF component: {source}")
        writer.append(PdfReader(source))

    temporary = output.with_suffix(".tmp.pdf")
    with temporary.open("wb") as handle:
        writer.write(handle)
    temporary.replace(output)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 040_combine_recommended_paper_pdfs.py <package-directory>")

    root = Path(sys.argv[1]).resolve()
    combine(
        root / "Main_Figures_and_Table_PC1to5.pdf",
        [
            "Figure1_real_data_controls.pdf",
            "Figure2_simulation_operating_characteristics.pdf",
            "Figure3_LCT_common_vs_rare_controls.pdf",
            "Table1_primary_results.pdf",
        ],
    )
    combine(
        root / "Supplementary_Figures_PC1to5.pdf",
        [
            "Supplementary_FigureS1_simulation_QQ.pdf",
            "Supplementary_FigureS2_simulation_balance.pdf",
            "Supplementary_FigureS2b_ratio_sensitivity_expanded.pdf",
            "Supplementary_FigureS3_permutation_vs_no_permutation_GIF.pdf",
            "Supplementary_FigureS4_burden_SKAT_SKATO.pdf",
            "Supplementary_FigureS5_variant_firth_forest.pdf",
        ],
    )


if __name__ == "__main__":
    main()
