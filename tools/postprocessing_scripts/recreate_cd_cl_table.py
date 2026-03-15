#!/usr/bin/env python3
"""Recreate the cylinder drag/lift comparison table PDF.

The original PDF metadata indicates it was produced with Matplotlib. This
script rebuilds the same content using pandas for tabular data handling and
Matplotlib tables for layout.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REFERENCE_DRAG = 6.185330
REFERENCE_LIFT = 0.009401

HEADER_COLORS = {
    "FBM-DEF": "#f7e9e6",
    "FBM-FAC": "#e7edf8",
    "CC": "#e8f3ea",
    "CC-iso": "#f5eedf",
}

COLUMN_LABELS = [
    "Level",
    "NEL",
    "NEQ",
    "Drag",
    "Lift",
    "Rel. Error Drag",
    "Rel. Error Lift",
]


@dataclass(frozen=True)
class Section:
    name: str
    rows: list[dict[str, object]]


SECTIONS = [
    Section(
        name="FBM-DEF",
        rows=[
            {"Level": 2, "NEL": 5120, "NEQ": "-", "Drag": 5.7840, "Lift": 0.12876},
            {"Level": 3, "NEL": 40960, "NEQ": "-", "Drag": 6.0034, "Lift": 0.083599},
            {"Level": 4, "NEL": 262144, "NEQ": "-", "Drag": 6.1384, "Lift": 0.011709},
        ],
    ),
    Section(
        name="FBM-FAC",
        rows=[
            {"Level": 2, "NEL": 5120, "NEQ": "-", "Drag": 5.443000, "Lift": -0.605210},
            {"Level": 3, "NEL": 40960, "NEQ": "-", "Drag": 5.989500, "Lift": 0.063517},
            {"Level": 4, "NEL": 262144, "NEQ": "-", "Drag": 6.031400, "Lift": 0.042072},
            {"Level": 5, "NEL": 2097152, "NEQ": "-", "Drag": 6.136700, "Lift": 0.005310},
        ],
    ),
    Section(
        name="CC",
        rows=[
            {"Level": 2, "NEL": 3072, "NEQ": 95112, "Drag": 6.140670, "Lift": 0.009671},
            {"Level": 3, "NEL": 24576, "NEQ": 723984, "Drag": 6.174367, "Lift": 0.009389},
            {"Level": 4, "NEL": 196608, "NEQ": 5647392, "Drag": 6.182615, "Lift": 0.009387},
            {"Level": 5, "NEL": 1572864, "NEQ": "45E+6", "Drag": 6.184652, "Lift": 0.009396},
        ],
    ),
    Section(
        name="CC-iso",
        rows=[
            {"Level": 2, "NEL": 3072, "NEQ": 95112, "Drag": 6.180395, "Lift": 0.009881},
            {"Level": 3, "NEL": 24576, "NEQ": 723984, "Drag": 6.184546, "Lift": 0.009464},
            {"Level": 4, "NEL": 196608, "NEQ": 5647392, "Drag": 6.185213, "Lift": 0.009407},
            {"Level": 5, "NEL": 1572864, "NEQ": "45E+6", "Drag": 6.185309, "Lift": 0.009402},
        ],
    ),
]


def build_dataframe(section: Section) -> pd.DataFrame:
    frame = pd.DataFrame(section.rows)
    frame["Rel. Error Drag"] = (
        (frame["Drag"] - REFERENCE_DRAG).abs() / REFERENCE_DRAG * 100.0
    )
    frame["Rel. Error Lift"] = (
        (frame["Lift"] - REFERENCE_LIFT).abs() / REFERENCE_LIFT * 100.0
    )
    return frame


def format_integer_like(value: object) -> str:
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def format_frame(frame: pd.DataFrame) -> pd.DataFrame:
    formatted = pd.DataFrame()
    formatted["Level"] = frame["Level"].map("{:d}".format)
    formatted["NEL"] = frame["NEL"].map(format_integer_like)
    formatted["NEQ"] = frame["NEQ"].map(format_integer_like)
    formatted["Drag"] = frame["Drag"].map("{:.6f}".format)
    formatted["Lift"] = frame["Lift"].map("{:.6f}".format)
    formatted["Rel. Error Drag"] = frame["Rel. Error Drag"].map("{:.2f}%".format)
    formatted["Rel. Error Lift"] = frame["Rel. Error Lift"].map("{:.2f}%".format)
    return formatted


def add_section_table(
    ax: plt.Axes,
    *,
    section_name: str,
    formatted_frame: pd.DataFrame,
    label_y: float,
    table_y: float,
    table_h: float,
) -> None:
    ax.text(
        0.001,
        label_y,
        section_name,
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    table = ax.table(
        cellText=formatted_frame[COLUMN_LABELS].values,
        colLabels=COLUMN_LABELS,
        cellLoc="center",
        colLoc="center",
        bbox=[0.001, table_y, 0.989, table_h],
    )

    table.auto_set_font_size(False)
    table.set_fontsize(13.5)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#202020")
        cell.set_linewidth(1.1)
        if row == 0:
            cell.set_facecolor(HEADER_COLORS[section_name])
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("white")


def create_figure(output_path: Path) -> None:
    fig = plt.figure(figsize=(12.4, 9.1))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    fig.suptitle(
        (
            "Cylinder coefficients comparison "
            f"(reference: C_d={REFERENCE_DRAG:.6f}, C_l={REFERENCE_LIFT:.6f})"
        ),
        fontsize=18,
        y=0.982,
    )

    label_start = 0.885
    label_step = 0.210
    table_h = 0.145
    table_top_offset = 0.019

    for index, section in enumerate(SECTIONS):
        frame = build_dataframe(section)
        formatted = format_frame(frame)
        label_y = label_start - index * label_step
        table_y = label_y - table_top_offset - table_h
        add_section_table(
            ax,
            section_name=section.name,
            formatted_frame=formatted,
            label_y=label_y,
            table_y=table_y,
            table_h=table_h,
        )

    fig.savefig(output_path, format=output_path.suffix.lstrip("."))
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("cd_cl_table_recreated.pdf"),
        help="Output file path. The extension controls the format.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    create_figure(args.output)


if __name__ == "__main__":
    main()
