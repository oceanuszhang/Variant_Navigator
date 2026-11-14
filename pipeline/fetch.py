"""Fetch wildtype protein sequence + annotated binding-site residues from UniProt.

Usage:
    from pipeline.fetch import fetch_protein
    p = fetch_protein("EGFR")          # gene name or UniProt accession
    p.sequence, p.binding_sites
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass, field, asdict
from typing import Any

UNIPROT = "https://rest.uniprot.org/uniprotkb"
TIMEOUT = 30

# Gene name -> reviewed UniProt accession, so demo runs are deterministic
# instead of depending on search ranking.
KNOWN_ACCESSIONS = {
    "ABL1": "P00519",
    "EGFR": "P00533",        # Homo sapiens epidermal growth factor receptor
    "RPOB": "P9WGY9",        # M. tuberculosis H37Rv RNA polymerase beta subunit
    # Zymoseptoria tritici sterol 14-alpha-demethylase. UniProt holds several
    # CYP51 alleles for this species (A4ULI5 515 aa, A4ULJ1 517 aa, B2L2K3
    # 364 aa); A4ULI5 is pinned because it carries Tyr137 and is full length.
    "CYP51": "A4ULI5",
    "ERG11": "A4ULI5",       # fungal CYP51 synonym
}

# Residue-level features worth keeping as "where a drug could bind".
BINDING_FEATURE_TYPES = {
    "Binding site",
    "Active site",
    "Site",
    "Metal binding",
    "Nucleotide binding",
}


@dataclass
class Protein:
    accession: str
    gene: str
    name: str
    organism: str
    sequence: str
    binding_sites: list[dict[str, Any]] = field(default_factory=list)
    regions: list[dict[str, Any]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def residue(self, position: int) -> str:
        """1-indexed residue, as used in variant notation like T790M."""
        if not 1 <= position <= self.length:
            raise IndexError(
                f"position {position} outside {self.accession} (1-{self.length})"
            )
        return self.sequence[position - 1]

    def binding_positions(self) -> list[int]:
        out: set[int] = set()
        for site in self.binding_sites:
            out.update(range(site["start"], site["end"] + 1))
        return sorted(out)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _get(url: str) -> dict[str, Any]:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
        return json.load(resp)


def resolve_accession(query: str, organism_id: int | None = None) -> str:
    """Turn a gene name into a UniProt accession. Accessions pass through."""
    q = query.strip()
    if q.upper() in KNOWN_ACCESSIONS and organism_id is None:
        return KNOWN_ACCESSIONS[q.upper()]
    if _looks_like_accession(q):
        return q.upper()

    terms = [f"(gene:{q} OR protein_name:{q})", "reviewed:true"]
    if organism_id:
        terms.append(f"organism_id:{organism_id}")
    url = (
        f"{UNIPROT}/search?query={urllib.parse.quote(' AND '.join(terms))}"
        "&fields=accession&size=1"
    )
    results = _get(url).get("results", [])
    if not results:
        # Retry unreviewed: fungal/plant-pathogen entries often aren't in Swiss-Prot.
        terms = [t for t in terms if t != "reviewed:true"]
        url = (
            f"{UNIPROT}/search?query={urllib.parse.quote(' AND '.join(terms))}"
            "&fields=accession&size=1"
        )
        results = _get(url).get("results", [])
    if not results:
        raise LookupError(f"no UniProt entry for {query!r}")
    return results[0]["primaryAccession"]


def _looks_like_accession(q: str) -> bool:
    u = q.upper()
    return (
        6 <= len(u) <= 10
        and u[0].isalpha()
        and u[1].isdigit()
        and any(c.isdigit() for c in u[-1])
    )


def fetch_protein(query: str, organism_id: int | None = None) -> Protein:
    """Fetch sequence + binding annotations for a gene name or accession."""
    accession = resolve_accession(query, organism_id)
    entry = _get(f"{UNIPROT}/{accession}.json")

    sequence = entry["sequence"]["value"]
    binding_sites: list[dict[str, Any]] = []
    regions: list[dict[str, Any]] = []

    for feat in entry.get("features", []):
        ftype = feat.get("type", "")
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        record = {
            "type": ftype,
            "start": start,
            "end": end,
            "description": feat.get("description", ""),
            "ligand": (feat.get("ligand") or {}).get("name", ""),
        }
        if ftype in BINDING_FEATURE_TYPES:
            binding_sites.append(record)
        elif ftype in {"Domain", "Region", "Transmembrane"}:
            regions.append(record)

    genes = entry.get("genes") or [{}]
    return Protein(
        accession=accession,
        gene=(genes[0].get("geneName") or {}).get("value", query),
        name=(
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        ),
        organism=entry.get("organism", {}).get("scientificName", ""),
        sequence=sequence,
        binding_sites=binding_sites,
        regions=regions,
    )


if __name__ == "__main__":
    import sys

    target = sys.argv[1] if len(sys.argv) > 1 else "EGFR"
    org = int(sys.argv[2]) if len(sys.argv) > 2 else None
    p = fetch_protein(target, org)
    print(f"{p.accession}  {p.gene}  {p.organism}")
    print(f"{p.name}")
    print(f"length={p.length}  binding_features={len(p.binding_sites)}")
    for s in p.binding_sites[:12]:
        span = f"{s['start']}" if s["start"] == s["end"] else f"{s['start']}-{s['end']}"
        print(f"  {s['type']:16} {span:>10}  {s['ligand'] or s['description']}")
    print(f"binding residues: {p.binding_positions()[:30]}")
