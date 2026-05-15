#!/usr/bin/env python3
"""Extract annual SST and air-temperature summaries from the SST workbook."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile


XLSX_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
RELS_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_RELS_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"x": XLSX_NS, "r": RELS_NS}


def default_workbook_path(script_dir: Path) -> Path:
    candidates = [
        script_dir / "data" / "sjavarhitamedaltol_manudir.xlsx",
        script_dir / "sjavarhitamedaltol_manudir.xlsx",
        script_dir.parent / "grímsey" / "sjavarhitamedaltol_manudir.xlsx",
    ]
    for path in candidates:
        if path.exists():
            return path
    return candidates[0]


def default_station_list_path(script_dir: Path) -> Path:
    return script_dir / "data" / "station-list.txt"


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

            if line_no == 1 and not re.match(r"^\d", parts[0]):
                continue
            if not re.match(r"^\d", parts[0]):
                raise ValueError(f"{path}:{line_no}: expected station code")

            station = str(int(parts[0])) if parts[0].isdigit() else parts[0]
            name = parts[1].strip() if len(parts) > 1 else ""
            stations.append((station, name))

    return stations


def column_number(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    number = 0
    for letter in letters:
        number = number * 26 + ord(letter.upper()) - ord("A") + 1
    return number


def shared_strings(xlsx: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in xlsx.namelist():
        return []

    root = ET.fromstring(xlsx.read("xl/sharedStrings.xml"))
    strings = []
    for item in root.findall("x:si", NS):
        strings.append("".join((t.text or "") for t in item.iter(f"{{{XLSX_NS}}}t")))
    return strings


def sheet_paths(xlsx: ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(xlsx.read("xl/workbook.xml"))
    rels = ET.fromstring(xlsx.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall(f"{{{PKG_RELS_NS}}}Relationship")
    }

    paths = {}
    for sheet in workbook.find("x:sheets", NS):
        rel_id = sheet.attrib[f"{{{RELS_NS}}}id"]
        target = targets[rel_id]
        if not target.startswith("xl/"):
            target = f"xl/{target}"
        paths[sheet.attrib["name"]] = target
    return paths


def cell_value(cell: ET.Element, strings: list[str]) -> str | None:
    value = cell.find("x:v", NS)
    if value is None:
        return None

    text = value.text
    if text is None:
        return None
    if cell.attrib.get("t") == "s":
        return strings[int(text)]
    return text


def normalized_int_text(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not number.is_integer():
        return None
    return str(int(number))


def numeric_text(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    if not text or text == "#":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return f"{number:.2f}"


def extract_sheet_rows(xlsx: ZipFile, sheet_path: str, strings: list[str]) -> list[tuple[str, str, str, str]]:
    sheet = ET.fromstring(xlsx.read(sheet_path))
    rows = []

    for row in sheet.find("x:sheetData", NS).findall("x:row", NS):
        values: dict[int, str | None] = {}
        for cell in row.findall("x:c", NS):
            col = column_number(cell.attrib["r"])
            if col in (1, 2, 15, 29):
                values[col] = cell_value(cell, strings)

        station = normalized_int_text(values.get(1))
        year = normalized_int_text(values.get(2))
        sst = numeric_text(values.get(15))
        temp = numeric_text(values.get(29))

        if station is None or year is None or sst is None or temp is None:
            continue

        rows.append((station, year, sst, temp))

    return rows


def station_sheet_name(station: str, paths: dict[str, str]) -> str | None:
    for wanted in (f"s{station}",):
        if wanted in paths:
            return wanted

    if not station.isdigit():
        return None

    wanted = f"s{int(station):03d}"
    if wanted in paths:
        return wanted

    pattern = re.compile(rf"^s0*{int(station)}(?:\b|[^0-9])")
    matches = [name for name in paths if pattern.match(name)]
    if len(matches) == 1:
        return matches[0]

    return None


def write_station_file(path: Path, rows: list[tuple[str, str, str, str]]) -> None:
    with path.open("w", encoding="utf-8") as file:
        file.write("station-number\tyear\tSST\tT\n")
        for station, year, sst, temp in rows:
            file.write(f"{station}\t{year}\t{sst}\t{temp}\n")


def output_file_name(station_name: str) -> str:
    name = station_name.strip()
    if not name:
        raise ValueError("station name is empty")
    return f"{name}.txt"


def main() -> int:
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stations", type=Path, default=default_station_list_path(script_dir))
    parser.add_argument("--workbook", type=Path, default=default_workbook_path(script_dir))
    parser.add_argument("--output-dir", type=Path, default=script_dir / "data")
    args = parser.parse_args()

    stations = read_station_list(args.stations)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with ZipFile(args.workbook) as xlsx:
        strings = shared_strings(xlsx)
        paths = sheet_paths(xlsx)

        missing = []
        written = 0
        for station, name in stations:
            sheet_name = station_sheet_name(station, paths)
            output_path = args.output_dir / output_file_name(name)

            if sheet_name is None:
                missing.append(station)
                write_station_file(output_path, [])
                continue

            rows = extract_sheet_rows(xlsx, paths[sheet_name], strings)
            write_station_file(output_path, rows)
            written += 1
            print(f"{station}: {len(rows)} rows from sheet {sheet_name} -> {output_path}")

    if missing:
        print("Missing sheets, wrote header-only files for:", ", ".join(missing))
    print(f"Done. Wrote {written} station files to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
