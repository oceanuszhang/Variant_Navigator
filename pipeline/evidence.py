"""Paperclip literature retrieval — the evidence channel, Step 5.

DESIGN RULE, and the reason this module is separate from the scoring modules:

    Evidence never feeds the prediction.

VariantNavigator's claim is that it answers "which drugs still work?" for a
variant with *no database entry and no trial data*. If literature retrieval fed
the prediction, the positive controls would stop being controls: T790M
resistance is in thousands of papers, so a pipeline allowed to read them is
looking up the answer, not deriving it. The demo would prove nothing.

So the pipeline runs two tracks:

    Track A — prediction   ESM2 + Boltz + Ptarmigan-1, blind to literature
    Track B — evidence     this module, run AFTER a prediction is fixed

Then `concordance()` compares them. That yields three outcomes, all useful:

    agree    prediction matches published resistance data -> validation
    disagree prediction contradicts it -> a bug, or a finding; either way, look
    novel    no literature at all -> exactly the gap the project exists to fill

The third is the pitch: run a variant Paperclip has never heard of, get zero
hits, and the pipeline still returns a ranked treatment list.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

Source = Literal["pmc", "biorxiv", "medrxiv", "arxiv", "abstracts", "fda", "trials/us", "proteins"]
Signal = Literal["resistance", "efficacy", "mixed", "none"]
Concordance = Literal["agree", "disagree", "novel", "unchecked"]

TIMEOUT = 60

# Paperclip ships as a CLI whose shebang is `/usr/bin/env python3`; whichever
# python that resolves to must have `click` installed. On this machine it is the
# anaconda interpreter, not the system one, so we find a working PATH once.
# A system-wide anaconda install lives in /opt, not under $HOME, and that is the
# case on this machine: click is in /opt/anaconda3 (3.13) while `python3`
# resolves to the python.org 3.11 framework build, which has no click. Searching
# only $HOME left paperclip permanently "unavailable" with a bare
# ModuleNotFoundError. /opt paths are listed first because a system install is
# the likelier place for a shared dependency.
_CANDIDATE_DIRS = [
    Path("/opt/anaconda3/bin"),
    Path("/opt/miniforge3/bin"),
    Path("/opt/miniconda3/bin"),
    Path("/opt/homebrew/bin"),
    Path.home() / "anaconda3" / "bin",
    Path.home() / "miniforge3" / "bin",
    Path.home() / "miniconda3" / "bin",
]

RESISTANCE_TERMS = re.compile(
    r"\b(resist\w*|refractor\w*|non-?suscept\w*|insensitiv\w*|escape|"
    r"reduced suscept\w*|treatment failure|breakthrough)\b",
    re.I,
)
EFFICACY_TERMS = re.compile(
    r"\b(sensitiv\w*|suscept\w*|respon\w*|effectiv\w*|efficac\w*|"
    r"overcome\w*|retain\w* activity|active against)\b",
    re.I,
)


class PaperclipUnavailable(RuntimeError):
    """Paperclip CLI missing, not signed in, or not importable."""


@dataclass
class Citation:
    title: str
    url: str
    date: str = ""
    authors: str = ""
    abstract: str = ""
    source: str = ""

    @property
    def ident(self) -> str:
        """PMC/DOI-ish id parsed from the URL.

        The CSV export's own `id` column is unreliable — in paperclip 0.8.17 it
        carries the match context ("Matched: Abstract") rather than an id, so the
        URL is the trustworthy source.
        """
        if m := re.search(r"(PMC\d+)", self.url):
            return m.group(1)
        if m := re.search(r"(10\.\d{4,}/\S+)", self.url):
            return m.group(1).rstrip("/")
        return self.url.rstrip("/").rsplit("/", 1)[-1]

    def signal(self) -> Signal:
        text = f"{self.title} {self.abstract}"
        res = bool(RESISTANCE_TERMS.search(text))
        eff = bool(EFFICACY_TERMS.search(text))
        if res and eff:
            return "mixed"
        if res:
            return "resistance"
        if eff:
            return "efficacy"
        return "none"


@dataclass
class EvidenceBundle:
    query: str
    source: str
    citations: list[Citation] = field(default_factory=list)
    search_id: str = ""
    error: str = ""

    @property
    def found(self) -> int:
        return len(self.citations)

    def signal(self) -> Signal:
        """Majority signal across the retrieved set."""
        if not self.citations:
            return "none"
        counts: dict[str, int] = {}
        for c in self.citations:
            s = c.signal()
            if s != "none":
                counts[s] = counts.get(s, 0) + 1
        if not counts:
            return "none"
        if len(counts) > 1 and "resistance" in counts and "efficacy" in counts:
            return "mixed"
        return max(counts, key=lambda k: counts[k])  # type: ignore[return-value]

    def to_dict(self) -> dict[str, Any]:
        return {
            "query": self.query,
            "source": self.source,
            "search_id": self.search_id,
            "found": self.found,
            "signal": self.signal(),
            "error": self.error,
            "citations": [
                {"id": c.ident, "title": c.title, "url": c.url, "date": c.date,
                 "signal": c.signal()}
                for c in self.citations
            ],
        }


@lru_cache(maxsize=1)
def _env() -> dict[str, str]:
    """An environment in which the paperclip CLI can import its deps."""
    if not shutil.which("paperclip"):
        raise PaperclipUnavailable("paperclip not on PATH")
    base = dict(os.environ)
    for d in _CANDIDATE_DIRS:
        if not (d / "python3").exists():
            continue
        env = dict(base, PATH=f"{d}:{base.get('PATH', '')}")
        probe = subprocess.run(
            ["paperclip", "--version"], env=env, capture_output=True, text=True, timeout=TIMEOUT
        )
        if probe.returncode == 0:
            return env
    probe = subprocess.run(
        ["paperclip", "--version"], env=base, capture_output=True, text=True, timeout=TIMEOUT
    )
    if probe.returncode == 0:
        return base
    raise PaperclipUnavailable(
        f"paperclip present but not runnable: {probe.stderr.strip()[:200]}"
    )


def available() -> bool:
    try:
        _env()
        return True
    except (PaperclipUnavailable, subprocess.SubprocessError, OSError):
        return False


def search(query: str, *, source: Source = "pmc", n: int = 5) -> EvidenceBundle:
    """Run one Paperclip search and return parsed citations."""
    bundle = EvidenceBundle(query=query, source=source)
    try:
        env = _env()
    except (PaperclipUnavailable, subprocess.SubprocessError, OSError) as exc:
        bundle.error = str(exc)
        return bundle

    proc = subprocess.run(
        ["paperclip", "search", "-s", source, query, "-n", str(n)],
        env=env, capture_output=True, text=True, timeout=TIMEOUT,
    )
    if proc.returncode != 0:
        bundle.error = proc.stderr.strip()[:300] or "paperclip search failed"
        return bundle

    if m := re.search(r"\[(s_[0-9a-f]+)\]", proc.stdout):
        bundle.search_id = m.group(1)
    if not bundle.search_id:
        bundle.error = "no search id in paperclip output"
        return bundle

    # CSV export is the machine-readable path; the TTY output is for humans.
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "r.csv"
        dump = subprocess.run(
            ["paperclip", "results", bundle.search_id, "--save", str(out)],
            env=env, capture_output=True, text=True, timeout=TIMEOUT,
        )
        if not out.exists():
            bundle.error = dump.stderr.strip()[:300] or "csv export failed"
            return bundle
        with out.open() as fh:
            for row in csv.DictReader(fh):
                bundle.citations.append(
                    Citation(
                        title=(row.get("title") or "").strip(),
                        url=(row.get("url") or "").strip(),
                        date=(row.get("date") or "").strip(),
                        authors=(row.get("authors") or "").strip(),
                        abstract=(row.get("abstract") or "").strip(),
                        source=source,
                    )
                )
    return bundle


def drug_evidence(
    gene: str, variant: str, drug: str, *, source: Source = "pmc", n: int = 5
) -> EvidenceBundle:
    """Literature on one variant/drug pair."""
    return search(f"{gene} {variant} {drug} resistance susceptibility", source=source, n=n)


def variant_novelty(gene: str, variant: str, *, n: int = 5) -> EvidenceBundle:
    """How well characterised is this variant at all? Zero hits = novel."""
    return search(f"{gene} {variant} mutation", source="pmc", n=n)


def collect(
    gene: str,
    variant: str,
    drugs: list[str],
    *,
    source: Source = "pmc",
    n: int = 3,
) -> dict[str, EvidenceBundle]:
    """Evidence per drug, plus a '_novelty' entry for the variant itself."""
    out: dict[str, EvidenceBundle] = {"_novelty": variant_novelty(gene, variant, n=n)}
    for drug in drugs:
        out[drug] = drug_evidence(gene, variant, drug, source=source, n=n)
    return out


EXTRACT_SCHEMA = {
    "type": "object",
    "properties": {
        "drug": {"type": "string"},
        "verdict": {"type": "string", "enum": ["resistant", "sensitive", "unclear"]},
        "fold_change": {
            "type": ["number", "null"],
            "description": "Reported fold change in IC50/GR50 vs reference, if stated",
        },
        "evidence": {"type": "string", "description": "Verbatim quote supporting the verdict"},
    },
    "required": ["drug", "verdict", "evidence"],
}


@dataclass
class Extraction:
    paper_id: str
    title: str
    verdict: str
    evidence: str
    fold_change: float | None = None
    drug: str = ""


def map_extract(
    search_id: str,
    question: str,
    *,
    schema: dict[str, Any] | None = None,
    limit: int = 5,
) -> list[Extraction]:
    """Run `paperclip map` over a search and parse the per-paper JSON answers.

    This replaces keyword matching with per-paper reading. Keyword signals are
    unreliable here because a paper reporting resistance uses the vocabulary of
    efficacy throughout — every T790M paper mentions both, so the naive signal
    collapses to "mixed". An LLM read against a fixed schema does not.
    """
    import json as _json

    try:
        env = _env()
    except (PaperclipUnavailable, subprocess.SubprocessError, OSError):
        return []

    proc = subprocess.run(
        [
            "paperclip", "map", "--from", search_id, "-n", str(limit),
            "--output-schema", _json.dumps(schema or EXTRACT_SCHEMA), question,
        ],
        env=env, capture_output=True, text=True, timeout=TIMEOUT * 4,
    )
    if proc.returncode != 0:
        return []

    out: list[Extraction] = []
    # Output is a human-readable run log with one JSON object per paper, each
    # preceded by a "<id> · <ms>" line. Walk it and brace-match the objects.
    lines = proc.stdout.splitlines()
    titles: dict[int, str] = {}
    for i, line in enumerate(lines):
        if line.strip().startswith("✓"):
            titles[i] = line.strip().lstrip("✓ ").strip()

    i = 0
    last_title, last_id = "", ""
    while i < len(lines):
        line = lines[i]
        if i in titles:
            last_title = titles[i]
        if m := re.match(r"\s*(PMC\d+|10\.\S+)\s*·", line):
            last_id = m.group(1)
        if line.strip().startswith("{"):
            depth, buf = 0, []
            while i < len(lines):
                buf.append(lines[i])
                depth += lines[i].count("{") - lines[i].count("}")
                if depth <= 0:
                    break
                i += 1
            try:
                obj = _json.loads("\n".join(buf))
                out.append(
                    Extraction(
                        paper_id=last_id,
                        title=last_title,
                        verdict=obj.get("verdict", "unclear"),
                        evidence=obj.get("evidence", ""),
                        fold_change=obj.get("fold_change"),
                        drug=obj.get("drug", ""),
                    )
                )
            except _json.JSONDecodeError:
                pass
        i += 1
    return out


def verdict_from_extractions(extractions: list[Extraction]) -> Signal:
    """Majority verdict across per-paper reads."""
    if not extractions:
        return "none"
    counts: dict[str, int] = {}
    for e in extractions:
        if e.verdict in ("resistant", "sensitive"):
            counts[e.verdict] = counts.get(e.verdict, 0) + 1
    if not counts:
        return "none"
    if len(counts) > 1 and abs(counts.get("resistant", 0) - counts.get("sensitive", 0)) == 0:
        return "mixed"
    top = max(counts, key=lambda k: counts[k])
    return "resistance" if top == "resistant" else "efficacy"


def concordance(predicted: str, bundle: EvidenceBundle) -> Concordance:
    """Compare a Track-A verdict against Track-B literature.

    `predicted` is a drug_screen verdict: work / partial / fail / n/a.
    """
    if bundle.error or bundle.found == 0:
        return "novel"
    sig = bundle.signal()
    if sig in ("none", "mixed"):
        return "unchecked"
    if predicted == "fail":
        return "agree" if sig == "resistance" else "disagree"
    if predicted == "work":
        return "agree" if sig == "efficacy" else "disagree"
    return "unchecked"


def report_section(gene: str, variant: str, evidence: dict[str, EvidenceBundle],
                   verdicts: dict[str, str] | None = None) -> str:
    """Markdown block for agent_report, with citations the judges can click."""
    verdicts = verdicts or {}
    novelty = evidence.get("_novelty")
    lines = ["## Literature evidence (Paperclip)", ""]
    if novelty is not None:
        if novelty.found == 0:
            lines += [
                f"**No published characterisation of {gene} {variant} was found.** "
                "The prediction above is unvalidated by prior literature — which is "
                "the case this pipeline is built for.",
                "",
            ]
        else:
            lines += [f"{gene} {variant}: {novelty.found} papers found.", ""]

    rows = [(d, b) for d, b in evidence.items() if d != "_novelty"]
    if rows:
        lines += ["| Drug | Predicted | Literature | Agreement | Top citation |",
                  "| --- | --- | --- | --- | --- |"]
        for drug, b in rows:
            pred = verdicts.get(drug, "n/a")
            cite = (
                f"[{b.citations[0].ident}]({b.citations[0].url})" if b.citations else "—"
            )
            lines.append(
                f"| {drug} | {pred} | {b.signal()} | {concordance(pred, b)} | {cite} |"
            )
        lines.append("")

    lines += [
        "_Evidence was retrieved after the prediction was fixed and did not inform "
        "it. Signals are keyword-derived from titles and abstracts, so treat "
        "'agree' as corroboration, not proof._",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    import json
    import sys

    gene = sys.argv[1] if len(sys.argv) > 1 else "EGFR"
    variant = sys.argv[2] if len(sys.argv) > 2 else "T790M"
    drugs = sys.argv[3].split(",") if len(sys.argv) > 3 else ["erlotinib", "osimertinib"]
    ev = collect(gene, variant, drugs)
    print(json.dumps({k: v.to_dict() for k, v in ev.items()}, indent=2)[:3000])
    print()
    print(report_section(gene, variant, ev, {"erlotinib": "fail", "osimertinib": "work"}))
