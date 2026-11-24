#!/usr/bin/env python3
"""Rebuild data/egfr_pic50.csv from ChEMBL.

Queries ChEMBL for IC50 records of EGFR (CHEMBL203), keeps records in nM with
an exact relation, converts them to pIC50 = 9 - log10(IC50 [nM]) and writes a
two-column CSV (smiles, pIC50). Duplicate molecules are merged by the median.

Only the standard library is used, so the script runs anywhere.

Examples
--------
    python scripts/fetch_chembl_egfr.py --out data/egfr_pic50.csv
    python scripts/fetch_chembl_egfr.py --max-records 200 --out /tmp/egfr_small.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data"
EGFR = "CHEMBL203"


def get_json(url: str, retries: int = 4) -> dict:
    last_error = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as handle:
                return json.load(handle)
        except Exception as exc:  # network hiccup or 5xx
            last_error = exc
            time.sleep(2 ** attempt)
    raise SystemExit(f"Request failed after {retries} attempts: {url}\n{last_error}")


def fetch_activities(target: str, page_size: int, max_records: int | None) -> list[dict]:
    activities, offset = [], 0
    while True:
        query = urllib.parse.urlencode({
            "target_chembl_id": target,
            "standard_type": "IC50",
            "standard_units": "nM",
            "limit": page_size,
            "offset": offset,
        })
        payload = get_json(f"{BASE}/activity.json?{query}")
        page = payload.get("activities", [])
        if not page:
            break
        activities.extend(page)
        total = payload.get("page_meta", {}).get("total_count")
        offset += len(page)
        print(f"  activities: {len(activities)}/{total if total else '?'}")
        if max_records and len(activities) >= max_records:
            return activities[:max_records]
        if total is not None and offset >= total:
            break
    return activities


def fetch_smiles(molecule_ids: list[str], batch: int = 50) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for start in range(0, len(molecule_ids), batch):
        chunk = molecule_ids[start:start + batch]
        query = urllib.parse.urlencode({
            "molecule_chembl_id__in": ",".join(chunk),
            "limit": len(chunk),
        })
        payload = get_json(f"{BASE}/molecule.json?{query}")
        for molecule in payload.get("molecules", []):
            chembl_id = molecule.get("molecule_chembl_id")
            smiles = (molecule.get("molecule_structures") or {}).get("canonical_smiles")
            if chembl_id and smiles:
                mapping[chembl_id] = smiles
        print(f"  structures: {min(start + batch, len(molecule_ids))}/{len(molecule_ids)}")
    return mapping


def to_pic50(value_nm: float) -> float:
    return 9.0 - math.log10(value_nm)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--target", default=EGFR, help="ChEMBL target id (default: EGFR)")
    parser.add_argument("--out", default="data/egfr_pic50.csv", type=Path)
    parser.add_argument("--page-size", default=1000, type=int)
    parser.add_argument("--max-records", default=None, type=int,
                        help="stop after N activity records (useful for a quick test)")
    args = parser.parse_args()

    print(f"Fetching IC50 records for {args.target} ...")
    activities = fetch_activities(args.target, args.page_size, args.max_records)

    # Keep exact, positive nM measurements only.
    measured = []
    for record in activities:
        value = record.get("standard_value")
        relation = record.get("standard_relation")
        chembl_id = record.get("molecule_chembl_id")
        units = record.get("standard_units")
        if not value or not chembl_id:
            continue
        if relation not in ("=", None) or units != "nM":
            continue
        try:
            value = float(value)
        except (TypeError, ValueError):
            continue
        if value <= 0:
            continue
        measured.append((chembl_id, value))

    print(f"  usable measurements: {len(measured)}")
    ids = sorted({chembl_id for chembl_id, _ in measured})
    print(f"  unique molecules: {len(ids)}")

    print("Fetching SMILES ...")
    smiles_by_id = fetch_smiles(ids)

    values_by_smiles: dict[str, list[float]] = defaultdict(list)
    for chembl_id, value in measured:
        smiles = smiles_by_id.get(chembl_id)
        if smiles:
            values_by_smiles[smiles].append(to_pic50(value))

    rows = []
    for smiles, values in values_by_smiles.items():
        values.sort()
        n = len(values)
        median = values[n // 2] if n % 2 else 0.5 * (values[n // 2 - 1] + values[n // 2])
        rows.append((smiles, round(median, 4)))
    rows.sort(key=lambda item: -item[1])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["smiles", "pIC50"])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} molecules -> {args.out}")
    if rows:
        print(f"  pIC50 range: {rows[-1][1]:.2f} - {rows[0][1]:.2f}")


if __name__ == "__main__":
    main()
