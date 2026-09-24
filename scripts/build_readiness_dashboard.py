"""Build a single-file interactive readiness dashboard.

Merges the outputs of the governance tools into one HTML page that answers
"what still blocks submission?" Sections cover the headline blockers, the
prospective-pool gate, corpus redistribution terms, external outreach, and
archive integrity.

The page is self-contained: no network access, no build step, no dependencies.
Open it directly from disk.
"""

from __future__ import annotations

import argparse
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_OUTPUT = Path("artifacts/dashboard/readiness-dashboard.html")

SOURCES = {
    "archive": Path("artifacts/release-archive/software-implementation-v6-manifest.json"),
    "pool": Path("artifacts/prospective-pool/dashboard.json"),
    "licenses": Path("artifacts/license-audit/development-corpus-license-audit.json"),
    "requests": Path("artifacts/external-requests/status.json"),
}

MANUSCRIPT_DIR = Path("manuscript")
DECLARED_DELIVERABLES = (
    {
        "id": "author-summary",
        "label": "Author Summary (150-200 words)",
        "globs": ("*Author_Summary*", "*AuthorSummary*", "*author_summary*"),
    },
    {
        "id": "cover-letter",
        "label": "Cover letter",
        "globs": ("*Cover_Letter*", "*CoverLetter*", "*cover_letter*"),
    },
)

PLACEHOLDER_PATTERNS = (
    re.compile(r"\[corresponding author\]", re.IGNORECASE),
    re.compile(r"\[repository URL\]", re.IGNORECASE),
    re.compile(r"\[AUTHOR", re.IGNORECASE),
)

REDACTION = "<redacted-local-path>"


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def check_declared_deliverables(root: Path) -> list[dict[str, Any]]:
    results = []
    for item in DECLARED_DELIVERABLES:
        found: list[str] = []
        for pattern in item["globs"]:
            found.extend(path.name for path in (root / MANUSCRIPT_DIR).glob(pattern))
        results.append(
            {
                "id": item["id"],
                "label": item["label"],
                "present": bool(found),
                "files": sorted(found),
            }
        )
    return results


def find_placeholders(root: Path) -> list[dict[str, Any]]:
    findings: list[dict[str, Any]] = []
    for path in sorted((root / MANUSCRIPT_DIR).glob("*.md")):
        text = path.read_text(encoding="utf-8")
        counts: dict[str, int] = {}
        for pattern in PLACEHOLDER_PATTERNS:
            hits = pattern.findall(text)
            if hits:
                counts[pattern.pattern.strip(r"\[\]\\")] = len(hits)
        if counts:
            findings.append(
                {"file": path.relative_to(root).as_posix(), "placeholders": counts}
            )
    return findings


def redact(text: str) -> str:
    """Replace absolute local paths so a shared dashboard cannot leak identity."""
    text = re.sub(r"[A-Za-z]:\\\\?[^\"'\s,;)]{3,}", REDACTION, text)
    text = re.sub(r"/(?:home|Users)/[A-Za-z0-9._-]+[^\"'\s,;)]*", REDACTION, text)
    return text


def collect_blockers(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Every unresolved item that stands between the project and a deposit."""
    blockers: list[dict[str, Any]] = []

    licenses = data.get("licenses") or {}
    if licenses and not licenses.get("valid", False):
        blockers.append(
            {
                "id": "license-audit",
                "severity": "blocking",
                "title": "Corpus license declarations are inconsistent",
                "detail": (
                    f"{len(licenses.get('issues', []))} issue(s) reported by the audit. "
                    "A data availability statement cannot be asserted while the declarations disagree."
                ),
                "remedy": "python scripts/audit_development_corpus_license.py --strict",
            }
        )

    archive = data.get("archive") or {}
    local_paths = archive.get("local_path_findings", [])
    if local_paths:
        affected = sorted({item["path"] for item in local_paths})
        blockers.append(
            {
                "id": "local-paths",
                "severity": "blocking",
                "title": f"{len(affected)} archived file(s) embed absolute local paths",
                "detail": (
                    "A double-blind deposit must not carry the authors' filesystem layout. "
                    "The Git author account and the upstream repository account share a name, "
                    "so a home-directory path is a direct de-anonymization vector."
                ),
                "remedy": "python scripts/build_release_archive.py --fail-on-local-paths",
                "items": affected,
            }
        )

    pool = data.get("pool") or {}
    unmet = [gate for gate in pool.get("gates", []) if not gate.get("satisfied")]
    if unmet:
        blockers.append(
            {
                "id": "pool-gate",
                "severity": "blocking",
                "title": "Prospective pool is below its frozen gate",
                "detail": (
                    "Until the gate passes, no independent validation claim can be made and "
                    "predictions are not permitted on the pool."
                ),
                "remedy": "python scripts/prospective_pool_dashboard.py --markdown "
                "artifacts/prospective-pool/dashboard.md",
                "items": [
                    f"{gate['gate']}: {gate['observed']}/{gate['required']} "
                    f"(short {gate['shortfall']})"
                    for gate in unmet
                ],
            }
        )

    requests = data.get("requests") or {}
    if requests and not requests.get("sent_count", 0):
        blockers.append(
            {
                "id": "no-outreach",
                "severity": "blocking",
                "title": "No external data request has been sent",
                "detail": (
                    "Every acquisition route depends on outreach that has not happened, so "
                    "the blocked workstreams cannot advance no matter how long they wait."
                ),
                "remedy": "python scripts/external_request_status.py --markdown "
                "artifacts/external-requests/status.md",
                "items": [
                    f"{entry['draft']} ({'no queue entry' if not entry['tracked_in_queue'] else entry.get('queue_status')})"
                    for entry in requests.get("entries", [])
                ],
            }
        )

    for item in data.get("deliverables", []):
        if not item["present"]:
            blockers.append(
                {
                    "id": f"missing-{item['id']}",
                    "severity": "blocking",
                    "title": f"{item['label']} is not drafted",
                    "detail": (
                        "The target journal requires this component. Its content depends on "
                        "the unresolved decision about the manuscript's primary contribution."
                    ),
                    "remedy": f"add {MANUSCRIPT_DIR.as_posix()}/{item['id']}.*",
                }
            )

    if data.get("placeholders"):
        total = sum(sum(item["placeholders"].values()) for item in data["placeholders"])
        blockers.append(
            {
                "id": "placeholders",
                "severity": "blocking",
                "title": f"{total} placeholder(s) remain in the manuscripts",
                "detail": "Submission metadata must be real before the files are uploaded.",
                "remedy": "replace the bracketed fields with the final author and repository details",
                "items": [
                    f"{item['file']}: {', '.join(f'{k} x{v}' for k, v in item['placeholders'].items())}"
                    for item in data["placeholders"]
                ],
            }
        )

    return blockers


def build_payload(root: Path, redact_paths: bool) -> dict[str, Any]:
    data: dict[str, Any] = {
        "sources_present": {},
    }
    for key, relative in SOURCES.items():
        payload = load_json(root / relative)
        data[key] = payload
        data["sources_present"][key] = payload is not None

    data["deliverables"] = check_declared_deliverables(root)
    data["placeholders"] = find_placeholders(root)
    data["blockers"] = collect_blockers(data)

    data["generated_on"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    data["release_id"] = (data.get("archive") or {}).get("release_id", "unknown")
    data["commit"] = (data.get("archive") or {}).get("built_from_commit")

    if redact_paths:
        data = json.loads(redact(json.dumps(data)))
        data["redacted"] = True
    else:
        data["redacted"] = False

    return data


def render(payload: dict[str, Any]) -> str:
    embedded = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    # A closing script tag inside the payload would terminate the block early.
    embedded = embedded.replace("</", "<\\/")
    title = "BioCandidateRanker submission readiness"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{
  --bg: #f6f7f9; --panel: #ffffff; --ink: #14181f; --muted: #5c6675;
  --line: #e2e6ec; --ok: #1f7a4d; --warn: #a8620a; --bad: #b3261e;
  --ok-bg: #e8f5ee; --warn-bg: #fdf3e3; --bad-bg: #fdecea; --accent: #2c6fbb;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--bg); color: var(--ink);
  font: 14px/1.55 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }}
header {{ background: var(--panel); border-bottom: 1px solid var(--line); padding: 20px 28px; }}
h1 {{ margin: 0 0 4px; font-size: 20px; letter-spacing: -0.01em; }}
.meta {{ color: var(--muted); font-size: 12.5px; }}
.verdict {{ display: inline-flex; align-items: center; gap: 8px; margin-top: 14px;
  padding: 8px 14px; border-radius: 8px; font-weight: 600; font-size: 13.5px; }}
.verdict.bad {{ background: var(--bad-bg); color: var(--bad); }}
.verdict.ok {{ background: var(--ok-bg); color: var(--ok); }}
nav {{ display: flex; gap: 2px; padding: 0 28px; background: var(--panel);
  border-bottom: 1px solid var(--line); overflow-x: auto; }}
nav button {{ border: 0; background: none; padding: 12px 16px; cursor: pointer;
  font: inherit; font-weight: 500; color: var(--muted); border-bottom: 2px solid transparent; white-space: nowrap; }}
nav button:hover {{ color: var(--ink); }}
nav button[aria-selected="true"] {{ color: var(--accent); border-bottom-color: var(--accent); }}
main {{ padding: 24px 28px 60px; max-width: 1180px; }}
section[hidden] {{ display: none; }}
h2 {{ font-size: 16px; margin: 0 0 14px; }}
h3 {{ font-size: 13.5px; margin: 22px 0 10px; color: var(--muted);
  text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 12px; }}
.card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 10px; padding: 14px 16px; }}
.card .k {{ color: var(--muted); font-size: 12px; }}
.card .v {{ font-size: 24px; font-weight: 650; margin-top: 4px; letter-spacing: -0.02em; }}
.card .s {{ color: var(--muted); font-size: 12px; margin-top: 2px; }}
.card.bad .v {{ color: var(--bad); }}
.card.ok .v {{ color: var(--ok); }}
.card.warn .v {{ color: var(--warn); }}
.blocker {{ background: var(--panel); border: 1px solid var(--line); border-left: 3px solid var(--bad);
  border-radius: 8px; margin-bottom: 10px; overflow: hidden; }}
.blocker > button {{ display: block; width: 100%; text-align: left; border: 0; background: none;
  padding: 14px 16px; cursor: pointer; font: inherit; }}
.blocker .t {{ font-weight: 600; }}
.blocker .d {{ color: var(--muted); font-size: 12.5px; margin-top: 3px; }}
.blocker .body {{ padding: 0 16px 16px; border-top: 1px solid var(--line); margin-top: 2px; }}
.blocker .body[hidden] {{ display: none; }}
.blocker code {{ background: #f1f3f7; padding: 2px 6px; border-radius: 4px;
  font-size: 12px; display: inline-block; margin-top: 8px; }}
.blocker ul {{ margin: 10px 0 0; padding-left: 18px; color: var(--muted); font-size: 12.5px; }}
.blocker ul li {{ margin-bottom: 2px; }}
.bar {{ background: #eceff4; border-radius: 999px; height: 10px; overflow: hidden; margin-top: 6px; }}
.bar > span {{ display: block; height: 100%; background: var(--accent); transition: width .5s ease; }}
.bar.pass > span {{ background: var(--ok); }}
.gate {{ margin-bottom: 16px; }}
.gate .row {{ display: flex; justify-content: space-between; font-size: 13px; }}
.gate .row .n {{ font-variant-numeric: tabular-nums; color: var(--muted); }}
table {{ width: 100%; border-collapse: collapse; background: var(--panel);
  border: 1px solid var(--line); border-radius: 10px; overflow: hidden; }}
th, td {{ text-align: left; padding: 9px 12px; border-bottom: 1px solid var(--line); font-size: 13px; }}
th {{ background: #fafbfc; font-weight: 600; color: var(--muted); font-size: 12px;
  text-transform: uppercase; letter-spacing: 0.04em; cursor: pointer; user-select: none; }}
th[data-sort]:hover {{ color: var(--ink); }}
tbody tr:last-child td {{ border-bottom: 0; }}
tbody tr:hover {{ background: #fafbfc; }}
.pill {{ display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11.5px; font-weight: 600; }}
.pill.ok {{ background: var(--ok-bg); color: var(--ok); }}
.pill.warn {{ background: var(--warn-bg); color: var(--warn); }}
.pill.bad {{ background: var(--bad-bg); color: var(--bad); }}
.pill.dim {{ background: #eef1f5; color: var(--muted); }}
.controls {{ display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 12px; align-items: center; }}
.controls input[type=search] {{ flex: 1 1 240px; max-width: 340px; padding: 8px 12px;
  border: 1px solid var(--line); border-radius: 8px; font: inherit; background: var(--panel); }}
.controls label {{ display: inline-flex; align-items: center; gap: 6px; font-size: 12.5px;
  color: var(--muted); cursor: pointer; }}
.famgrid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(74px, 1fr)); gap: 6px; }}
.fam {{ border: 1px solid var(--line); border-radius: 6px; padding: 6px; text-align: center;
  font-size: 11px; background: var(--panel); }}
.fam .n {{ font-weight: 650; font-size: 13px; font-variant-numeric: tabular-nums; }}
.fam.cap {{ border-color: var(--warn); background: var(--warn-bg); }}
.fam .l {{ color: var(--muted); font-size: 9.5px; }}
.mono {{ font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; font-size: 11.5px;
  color: var(--muted); word-break: break-all; }}
.note {{ color: var(--muted); font-size: 12.5px; margin-top: 10px; }}
.empty {{ color: var(--muted); font-style: italic; }}
.count {{ color: var(--muted); font-size: 12.5px; margin-left: auto; }}
</style>
</head>
<body>
<header>
  <h1>{html.escape(title)}</h1>
  <div class="meta" id="meta"></div>
  <div id="verdict"></div>
</header>
<nav id="tabs" role="tablist"></nav>
<main id="main"></main>
<script id="payload" type="application/json">{embedded}</script>
<script>
const DATA = JSON.parse(document.getElementById('payload').textContent);
const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c =>
  ({{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}})[c]);
const fmt = (n) => typeof n === 'number' ? n.toLocaleString() : n;

/* ---------- overview ---------- */
function overview() {{
  const a = DATA.archive || {{}}, p = DATA.pool || {{}}, r = DATA.requests || {{}}, l = DATA.licenses || {{}};
  const unmet = (p.gates || []).filter(g => !g.satisfied).length;
  const cards = [
    {{ k: 'Blocking items', v: DATA.blockers.length, cls: DATA.blockers.length ? 'bad' : 'ok' }},
    {{ k: 'Unmet pool gates', v: `${{unmet}} / ${{(p.gates||[]).length}}`, cls: unmet ? 'bad' : 'ok' }},
    {{ k: 'Requests sent', v: `${{r.sent_count ?? 0}} / ${{r.draft_count ?? 0}}`, cls: (r.sent_count ?? 0) ? 'ok' : 'bad' }},
    {{ k: 'Corpora shippable', v: `${{(l.shippable_corpora||[]).length}} / ${{l.corpus_count ?? 0}}`, cls: 'warn' }},
    {{ k: 'Archived files', v: fmt(a.file_count), s: a.total_size_bytes ? `${{(a.total_size_bytes/1048576).toFixed(1)}} MB uncompressed` : '' }},
    {{ k: 'Frozen evidence verified', v: `${{(a.frozen_evidence_verified||[]).length}} / ${{(a.frozen_evidence_verified||[]).length + (a.missing_frozen_evidence||[]).length}}`, cls: (a.missing_frozen_evidence||[]).length ? 'bad' : 'ok' }},
    {{ k: 'Local-path findings', v: (a.local_path_findings||[]).length, cls: (a.local_path_findings||[]).length ? 'bad' : 'ok' }},
    {{ k: 'Placeholder fields', v: (DATA.placeholders||[]).reduce((n,p)=>n+Object.values(p.placeholders).reduce((x,y)=>x+y,0),0), cls: (DATA.placeholders||[]).length ? 'bad' : 'ok' }},
  ];
  return `
    <div class="cards">${{cards.map(c => `
      <div class="card ${{c.cls||''}}">
        <div class="k">${{esc(c.k)}}</div>
        <div class="v">${{esc(c.v)}}</div>
        ${{c.s ? `<div class="s">${{esc(c.s)}}</div>` : ''}}
      </div>`).join('')}}</div>
    <h3>What still blocks submission</h3>
    ${{DATA.blockers.length ? DATA.blockers.map((b,i) => `
      <div class="blocker">
        <button onclick="toggle('b${{i}}')" aria-expanded="false">
          <div class="t">${{esc(b.title)}}</div>
          <div class="d">${{esc(b.detail)}}</div>
        </button>
        <div class="body" id="b${{i}}" hidden>
          ${{(b.items||[]).length ? `<ul>${{b.items.map(x => `<li>${{esc(x)}}</li>`).join('')}}</ul>` : ''}}
          <code>${{esc(b.remedy)}}</code>
        </div>
      </div>`).join('') : '<p class="empty">No blocking items recorded.</p>'}}`;
}}

/* ---------- pool ---------- */
function pool() {{
  const p = DATA.pool || {{}};
  if (!p.gates) return '<p class="empty">Pool dashboard output not found. Run scripts/prospective_pool_dashboard.py --json-out artifacts/prospective-pool/dashboard.json</p>';
  const gates = p.gates.map(g => `
    <div class="gate">
      <div class="row">
        <span>${{esc(g.gate)}}</span>
        <span class="n">${{fmt(g.observed)}} / ${{fmt(g.required)}} ${{g.satisfied ? '' : `· short ${{g.shortfall}}`}}</span>
      </div>
      <div class="bar ${{g.satisfied ? 'pass' : ''}}"><span style="width:${{Math.min(100, 100*g.observed/g.required).toFixed(1)}}%"></span></div>
    </div>`).join('');
  const fams = (p.families||[]).map(f => `
    <div class="fam ${{f.at_cap ? 'cap' : ''}}" title="${{esc(f.family)}} — ${{f.records}} records, headroom ${{f.headroom}}">
      <div class="n">${{f.records}}</div>
      <div class="l">${{esc(f.family.replace('family-','#'))}}</div>
    </div>`).join('');
  const proj = p.projection || {{}};
  const opts = (p.options||[]).map(o => `
    <tr><td><strong>${{esc(o.label)}}</strong><div class="d" style="color:var(--muted);font-size:12px">${{esc(o.requires)}}</div></td>
    <td>${{esc(o.cost)}}</td></tr>`).join('');
  return `
    <h2>Prospective pool gate</h2>
    ${{gates}}
    <p class="note">Gate passes: <strong>${{p.readiness_gate_passes}}</strong> · predictions permitted:
      <strong>${{p.predictions_permitted}}</strong> · family cap ${{p.family_cap}}</p>
    <h3>Acquisition projection</h3>
    <div class="cards">
      <div class="card"><div class="k">Records needed</div><div class="v">${{proj.records_needed}}</div></div>
      <div class="card"><div class="k">Families needed</div><div class="v">${{proj.families_needed}}</div></div>
      <div class="card"><div class="k">Records per source</div><div class="v">${{proj.records_per_accepted_source}}</div></div>
      <div class="card"><div class="k">Sources needed</div><div class="v">${{proj.sources_if_records_drive_acquisition}}</div></div>
    </div>
    <h3>Resolution options</h3>
    <table><thead><tr><th>Option</th><th>Cost</th></tr></thead><tbody>${{opts}}</tbody></table>
    <h3>Per-family saturation <span class="count">${{(p.saturation||{{}}).families_at_cap?.length || 0}} at the ${{p.family_cap}}-record cap</span></h3>
    <div class="famgrid">${{fams}}</div>
    <p class="note">${{esc((p.saturation||{{}}).headroom_note || '')}}</p>`;
}}

/* ---------- licenses ---------- */
const LICENSE_CLASS = {{
  permitted_with_attribution: 'ok', permitted_unconditionally: 'ok',
  not_verified: 'bad', blocked_pending_permission: 'bad',
  not_applicable_reference_only: 'dim'
}};
let licFilter = 'all', licQuery = '', licSort = null, licDesc = true;
function licenses() {{
  const l = DATA.licenses;
  if (!l || !l.corpora) return '<p class="empty">License audit output not found. Run scripts/audit_development_corpus_license.py --json-out artifacts/license-audit/development-corpus-license-audit.json</p>';
  return `
    <h2>Corpus redistribution terms</h2>
    <p class="note">Consistent: <strong>${{l.valid}}</strong> · ${{l.blocking_corpora.length}} of ${{l.corpus_count}} corpora may not ship.
      A permissive code license does not grant redistribution rights over compiled data.</p>
    <div class="controls">
      <input type="search" placeholder="Search corpus, license, or role…" oninput="licQuery=this.value.toLowerCase();paintLicenses()">
      <label><input type="radio" name="licf" value="all" checked onchange="licFilter='all';paintLicenses()"> all</label>
      <label><input type="radio" name="licf" value="ship" onchange="licFilter='ship';paintLicenses()"> ships</label>
      <label><input type="radio" name="licf" value="blocked" onchange="licFilter='blocked';paintLicenses()"> blocked</label>
      <span class="count" id="liccount"></span>
    </div>
    <table>
      <thead><tr>
        <th data-sort="id" onclick="sortLic('id')">Corpus</th>
        <th data-sort="role" onclick="sortLic('role')">Role</th>
        <th data-sort="license_name" onclick="sortLic('license_name')">License</th>
        <th data-sort="redistribution" onclick="sortLic('redistribution')">Redistribution</th>
        <th data-sort="redistribute_in_release" onclick="sortLic('redistribute_in_release')">In release</th>
      </tr></thead>
      <tbody id="licbody"></tbody>
    </table>
    <div id="licdetail"></div>`;
}}
function licRows() {{
  let rows = (DATA.licenses.corpora || []).slice();
  if (licFilter === 'ship') rows = rows.filter(x => x.redistribute_in_release);
  if (licFilter === 'blocked') rows = rows.filter(x => !x.redistribute_in_release);
  if (licQuery) rows = rows.filter(x => [x.id,x.role,x.license_name,x.redistribution].join(' ').toLowerCase().includes(licQuery));
  if (licSort) rows.sort((a,b) => {{
    const x = a[licSort], y = b[licSort];
    const c = typeof x === 'boolean' ? (x===y?0:x?-1:1) : String(x).localeCompare(String(y));
    return licDesc ? -c : c;
  }});
  return rows;
}}
function paintLicenses() {{
  const rows = licRows();
  document.getElementById('licbody').innerHTML = rows.map(x => `
    <tr onclick="showLic('${{esc(x.id)}}')" style="cursor:pointer">
      <td><strong>${{esc(x.id)}}</strong></td>
      <td>${{esc(x.role)}}</td>
      <td>${{esc(x.license_name)}}</td>
      <td><span class="pill ${{LICENSE_CLASS[x.redistribution]||'dim'}}">${{esc(x.redistribution)}}</span></td>
      <td>${{x.redistribute_in_release ? 'yes' : 'no'}}</td>
    </tr>`).join('') || '<tr><td colspan="5" class="empty">No corpora match.</td></tr>';
  document.getElementById('liccount').textContent = `${{rows.length}} of ${{(DATA.licenses.corpora||[]).length}}`;
}}
function sortLic(key) {{ licDesc = (licSort === key) ? !licDesc : true; licSort = key; paintLicenses(); }}
function showLic(id) {{
  const x = (DATA.licenses.corpora || []).find(c => c.id === id);
  if (!x) return;
  document.getElementById('licdetail').innerHTML = `
    <div class="blocker" style="border-left-color:var(--accent);margin-top:14px">
      <div style="padding:14px 16px">
        <div class="t">${{esc(x.id)}}</div>
        <div class="d" style="color:var(--muted);font-size:12.5px;margin-top:4px">${{esc(x.license_name)}} · ${{esc(x.role)}}</div>
        ${{x.blocking_reason ? `<div class="body" style="border-top:0;margin-top:8px;padding:0"><strong>Why it cannot ship:</strong> ${{esc(x.blocking_reason)}}</div>` : ''}}
        <div class="note">${{x.evidence_count}} evidence path(s)${{(x.missing_evidence||[]).length ? ` · <span style="color:var(--bad)">missing: ${{esc(x.missing_evidence.join(', '))}}</span>` : ''}}</div>
      </div>
    </div>`;
}}

/* ---------- requests ---------- */
function requests() {{
  const r = DATA.requests;
  if (!r || !r.entries) return '<p class="empty">Request status output not found. Run scripts/external_request_status.py --json-out artifacts/external-requests/status.json</p>';
  const rows = r.entries.map(e => `
    <tr>
      <td class="mono">${{esc(e.draft.replace('docs/external-requests/',''))}}</td>
      <td>${{esc(e.recipient || '—')}} ${{e.recipient_is_placeholder ? '<span class="pill bad">placeholder</span>' : ''}}</td>
      <td>${{e.tracked_in_queue ? `<span class="pill dim">${{esc(e.queue_status)}}</span>` : '<span class="pill bad">untracked</span>'}}</td>
      <td>${{e.sent_at ? esc(e.sent_at) : '<span class="pill bad">no</span>'}}</td>
    </tr>`).join('');
  return `
    <h2>External data requests</h2>
    <p class="note">${{esc(r.headline)}}</p>
    <div class="cards">
      <div class="card ${{r.sent_count ? 'ok':'bad'}}"><div class="k">Sent</div><div class="v">${{r.sent_count}} / ${{r.draft_count}}</div></div>
      <div class="card warn"><div class="k">Ready to send</div><div class="v">${{r.ready_to_send_count}}</div></div>
      <div class="card"><div class="k">Blocked on recipient</div><div class="v">${{r.blocked_count}}</div></div>
      <div class="card ${{r.untracked_drafts.length ? 'bad':'ok'}}"><div class="k">Untracked drafts</div><div class="v">${{r.untracked_drafts.length}}</div></div>
    </div>
    <h3>Drafts</h3>
    <table><thead><tr><th>Draft</th><th>Recipient</th><th>Queue state</th><th>Sent</th></tr></thead>
      <tbody>${{rows}}</tbody></table>
    ${{r.issues.length ? `<h3>Reconciliation issues</h3><div class="blocker" style="border-left-color:var(--warn)">
      <div style="padding:14px 16px"><ul style="margin:0;padding-left:18px;color:var(--muted);font-size:12.5px">
      ${{r.issues.map(i => `<li>${{esc(i)}}</li>`).join('')}}</ul></div></div>` : ''}}`;
}}

/* ---------- archive ---------- */
function archive() {{
  const a = DATA.archive;
  if (!a) return '<p class="empty">Archive manifest not found. Run scripts/build_release_archive.py</p>';
  const bad = (a.frozen_evidence_verified||[]).filter(v => !v.sha256_matches_manifest || !v.size_matches_manifest);
  const lp = a.local_path_findings || [];
  const ownership = Object.entries((a.license_ownership||{{}}).by_corpus || {{}})
    .sort((x,y) => y[1].length - x[1].length);
  return `
    <h2>Release archive integrity</h2>
    <div class="cards">
      <div class="card"><div class="k">Release</div><div class="v" style="font-size:16px">${{esc(a.release_id)}}</div>
        <div class="s">version ${{esc(a.package_version)}}</div></div>
      <div class="card"><div class="k">Files</div><div class="v">${{fmt(a.file_count)}}</div>
        <div class="s">${{(a.total_size_bytes/1048576).toFixed(1)}} MB total</div></div>
      <div class="card ${{bad.length ? 'bad':'ok'}}"><div class="k">Frozen evidence</div>
        <div class="v">${{(a.frozen_evidence_verified||[]).length - bad.length}} / ${{(a.frozen_evidence_verified||[]).length}}</div>
        <div class="s">match the freeze</div></div>
      <div class="card ${{lp.length ? 'bad':'ok'}}"><div class="k">Local-path findings</div><div class="v">${{lp.length}}</div></div>
    </div>
    <h3>Identity</h3>
    <p class="mono">commit ${{esc(a.built_from_commit)}}<br>built ${{esc(a.built_on)}}<br>manifest ${{esc(a.source_manifest_sha256)}}</p>
    <h3>Corpus ownership of archived artifacts</h3>
    <table><thead><tr><th>Corpus</th><th>Files</th></tr></thead><tbody>
      ${{ownership.map(([k,v]) => `<tr><td>${{esc(k)}}</td><td>${{v.length}}</td></tr>`).join('')}}
    </tbody></table>
    <h3>Local-path findings <span class="count">redact before a double-blind deposit</span></h3>
    ${{lp.length ? `<table><thead><tr><th>File</th><th>Sample match</th></tr></thead><tbody>
      ${{lp.map(f => `<tr><td class="mono">${{esc(f.path)}}</td><td class="mono">${{esc(f.matches[0])}}</td></tr>`).join('')}}
    </tbody></table>` : '<p class="empty">None found.</p>'}}
    <p class="note">${{esc(a.claim_boundary || '')}}</p>`;
}}

/* ---------- render ---------- */
const SECTIONS = [
  ['overview', 'Overview', overview],
  ['pool', 'Prospective pool', pool],
  ['licenses', 'Corpus licenses', licenses],
  ['requests', 'External requests', requests],
  ['archive', 'Archive', archive],
];
function toggle(id) {{
  const el = document.getElementById(id);
  const btn = el.previousElementSibling;
  el.hidden = !el.hidden;
  btn.setAttribute('aria-expanded', String(!el.hidden));
}}
function paint() {{
  const name = location.hash.replace('#','') || 'overview';
  const section = SECTIONS.find(s => s[0] === name) || SECTIONS[0];
  document.getElementById('main').innerHTML = section[2]();
  if (section[0] === 'licenses') paintLicenses();
  [...document.querySelectorAll('#tabs button')].forEach(b =>
    b.setAttribute('aria-selected', String(b.dataset.id === section[0])));
}}
DATA.blockers = DATA.blockers || [];
document.getElementById('tabs').innerHTML = SECTIONS.map(s =>
  `<button role="tab" data-id="${{s[0]}}" onclick="location.hash='${{s[0]}}'">${{esc(s[1])}}</button>`).join('');
document.getElementById('meta').innerHTML =
  `release <strong>${{esc(DATA.release_id)}}</strong> · commit ${{esc((DATA.commit||'').slice(0,7))}} · generated ${{esc(DATA.generated_on)}}`
  + (DATA.redacted ? ' · <span class="pill warn">paths redacted</span>' : '');
const n = DATA.blockers.length;
document.getElementById('verdict').innerHTML = n
  ? `<div class="verdict bad">Not submission ready — ${{n}} blocking item${{n===1?'':'s'}}</div>`
  : `<div class="verdict ok">No blocking items recorded</div>`;
addEventListener('hashchange', paint);
paint();
</script>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--redact-paths",
        action="store_true",
        help="replace absolute local paths before embedding, for a shareable copy",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    payload = build_payload(root, args.redact_paths)
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(payload), encoding="utf-8")

    missing = sorted(key for key, present in payload["sources_present"].items() if not present)
    print(
        json.dumps(
            {
                "dashboard": output.relative_to(root).as_posix(),
                "size_bytes": output.stat().st_size,
                "blockers": len(payload["blockers"]),
                "blocker_ids": [item["id"] for item in payload["blockers"]],
                "missing_sources": missing,
                "redacted": payload["redacted"],
                "open_with": f"start {output.name}   (or drag the file into a browser)",
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
