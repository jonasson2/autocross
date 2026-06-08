#!/usr/bin/env python3
"""Plot annual SST and air temperature for stations in data/station-list.txt."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, TextArea


def read_station_list(path: Path) -> list[tuple[str, str]]:
    stations: list[tuple[str, str]] = []
    with path.open(encoding="utf-8-sig") as file:
        for line_no, line in enumerate(file, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(maxsplit=1)

            if line_no == 1 and not parts[0][0].isdigit():
                continue
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_no}: expected station code and name")
            stations.append((parts[0], parts[1].strip()))
    return stations


def read_station_data(path: Path) -> tuple[list[int], list[float], list[float]]:
    years: list[int] = []
    sst: list[float] = []
    temp: list[float] = []

    with path.open(encoding="utf-8") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            years.append(int(row["year"]))
            sst.append(float(row["SST"]))
            temp.append(float(row["T"]))

    return years, sst, temp


def expand_to_year_grid(
    years: list[int],
    sst: list[float],
    temp: list[float],
    xmin: int,
    xmax: int,
) -> tuple[list[int], list[float], list[float]]:
    by_year = {year: (sst_value, temp_value) for year, sst_value, temp_value in zip(years, sst, temp)}
    full_years = list(range(xmin, xmax + 1))
    full_sst = []
    full_temp = []

    for year in full_years:
        values = by_year.get(year)
        if values is None:
            full_sst.append(math.nan)
            full_temp.append(math.nan)
        else:
            full_sst.append(values[0])
            full_temp.append(values[1])

    return full_years, full_sst, full_temp


def add_station_label(ax, name: str, xmax: int, ymin: int) -> None:
    text = TextArea(name, textprops={"fontsize": 9, "fontweight": "bold"})
    label = AnnotationBbox(
        text,
        (xmax, ymin),
        xybox=(-2, 2),
        xycoords="data",
        boxcoords="offset points",
        box_alignment=(1, 0),
        frameon=True,
        pad=0.18,
        bboxprops={"facecolor": "white", "edgecolor": "none", "alpha": 1},
    )
    ax.add_artist(label)


def display_station_code(station: str) -> str:
    return f"{int(station):03d}" if station.isdigit() else station


def main() -> int:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=script_dir / "data")
    parser.add_argument("--stations", type=Path, default=script_dir / "data" / "station-list.txt")
    parser.add_argument("--output", type=Path, default=script_dir / "sst-panel.png")
    parser.add_argument("--show", action="store_true", help="show the figure interactively")
    args = parser.parse_args()

    stations = read_station_list(args.stations)
    if len(stations) != 12:
        raise ValueError(f"expected 12 stations, found {len(stations)} in {args.stations}")

    fig, axes = plt.subplots(
        6,
        2,
        figsize=(17 / 2.54, 17 / 2.54),
        sharex=True,
        sharey=False,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )

    xmin = 1865
    xmax = 2025
    ylimits_by_station = {
        "Hraun á Skaga": (-1, 7),
        "Grímsey": (-1, 7),
        "Strandir": (-1, 7),
        "Raufarhöfn": (-1, 7),
        "Suðureyri": (-1, 7),
        "Þorvaldsstaðir": (-1, 7),
        "Stykkishólmur": (0, 8),
        "Teigarhorn": (0, 8),
        "Reykjavík": (0, 8),
        "Papey": (0, 8),
        "Grindavík": (2, 9),
        "Stórhöfði": (2, 9),
    }
    xticks = list(range(1880, 2021, 20))

    for index, ((station, name), ax) in enumerate(zip(stations, axes.flat)):
        row = index // 2
        col = index % 2
        ymin, ymax = ylimits_by_station[name]
        yticks = list(range(ymin, ymax + 1))

        data_path = args.data_dir / f"{name}.txt"
        years, sst, temp = read_station_data(data_path)
        years, sst, temp = expand_to_year_grid(years, sst, temp, xmin, xmax)

        ax.plot(years, sst, color="blue", linewidth=1.5, label="SST")
        ax.plot(years, temp, color="red", linewidth=1.5, label="T")
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        ax.set_yticklabels([""] + [str(value) for value in yticks[1:-1]] + [""])
        ax.grid(True, which="major", color="0.82", linewidth=0.7)
        add_station_label(ax, f"{name} ({display_station_code(station)})", xmax, ymin)

        ax.tick_params(
            bottom=False,
            left=False,
            top=False,
            right=False,
            labelbottom=row == 5,
            labelleft=col == 0,
            labelsize=8,
            pad=1.5,
        )

    fig.supylabel("Annual mean temperature (°C)")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.09, top=0.99, wspace=0, hspace=0)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    print(f"Wrote {args.output}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
