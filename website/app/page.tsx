"use client";

import { useEffect, useMemo, useState } from "react";
import {
  Activity, BookOpen, ChevronLeft, ChevronRight, CircleDot, Database,
  Dna, ExternalLink, FlaskConical, Layers3, LoaderCircle, MapPin,
  Search, SlidersHorizontal, Stethoscope,
} from "lucide-react";

type Protein = {
  accession: string; symbol: string; name: string; membraneClass: "A" | "B" | "C";
  membraneMode: string; evidenceTier: "E1" | "E2" | "E3"; sequenceLength: number;
  pairCount: number; negativePairCount: number; evidenceCount: number; bestPairTier: string; diseaseCount: number;
  siteCount: number; defaultWeb: boolean;
};
type Disease = { relationId: string; id: string; name: string; sourceIds: string; sources: string[]; channels: string[]; level: string; otScore: number | null; pubmed: string[]; mapping: string };
type Localization = { term: string; go: string; role: string; reliability: string; source: string; version: string; semantics: string };
type Classification = { family?: string; molecularFunction?: string; biologicalProcess?: string; membraneRole?: string; specialist?: string; reactome?: string; status?: string };
type PairTuple = [string, string, string, string, string, number | null, string, string, string, number, number, string, number, number, number, string, number, number, number, string];
type Pair = { pairId: string; compoundId: string; compoundName: string; inchikey: string; formula: string; mw: number | null; status: string; flags: string; tier: string; evidenceCount: number; databaseCount: number; sources: string[]; lineageCount: number; structureCount: number; modalityCount: number; conflict: string; defaultWeb: boolean; siteCount: number; coordinateSiteCount: number; siteSides: string };
type NegativePairTuple = [string, string, string, number, number];
type NegativePair = { compoundId: string; compoundName: string; inchikey: string; evidenceCount: number; conflict: boolean };
type NegativeSummary = { pairCount: number; evidenceCount: number; conflictCount: number; pageSize: number; pageCount: number; source: string };
type Detail = Protein & { reviewedStatus: string; secondaryModes: string[]; decision: string; decisionEvidence: string; formSpecificity: string; gene: { hgnc: string; ensembl: string; ncbi: string }; sequence: string; classification: Classification; localizations: Localization[]; diseases: Disease[]; pairs: PairTuple[]; negativePairs: NegativeSummary };
type IndexData = { release: string; source: string; statistics: { formalProteins: number; formalPairs: number; formalEvidence: number; negativePairs: number; negativeEvidence: number; interactionLinkedCompounds: number; diseases: number; diseaseRelations: number; sites: number; classCounts: Record<string, number>; tierCounts: Record<string, number> }; proteins: Protein[] };

const PAGE_SIZE = 20;
const n = new Intl.NumberFormat("en-US");
const classInfo = {
  A: { title: "Membrane-spanning", color: "#7b95c6", tone: "bg-[#e8edf7] text-[#38548c]" },
  B: { title: "Lipid-anchored", color: "#49c2d9", tone: "bg-[#e5f7fa] text-[#176977]" },
  C: { title: "Peripheral / associated", color: "#67a583", tone: "bg-[#eaf3ed] text-[#3e7054]" },
};

function decodePair(p: PairTuple): Pair {
  return { pairId:p[0], compoundId:p[1], compoundName:p[2], inchikey:p[3], formula:p[4], mw:p[5], status:p[6], flags:p[7], tier:p[8], evidenceCount:p[9], databaseCount:p[10], sources:p[11] ? p[11].split(";") : [], lineageCount:p[12], structureCount:p[13], modalityCount:p[14], conflict:p[15], defaultWeb:Boolean(p[16]), siteCount:p[17], coordinateSiteCount:p[18], siteSides:p[19] };
}

function decodeNegativePair(p: NegativePairTuple): NegativePair {
  return { compoundId:p[0], compoundName:p[1], inchikey:p[2], evidenceCount:p[3], conflict:Boolean(p[4]) };
}

function Badge({ children, tone = "blue" }: { children: React.ReactNode; tone?: "blue" | "teal" | "gold" | "rose" | "gray" }) {
  const tones = { blue:"bg-[#e8edf7] text-[#38548c]", teal:"bg-[#e4f5f3] text-[#176b68]", gold:"bg-[#fff4d5] text-[#836322]", rose:"bg-[#fde9e5] text-[#9c4e45]", gray:"bg-[#edf1f3] text-[#5d707a]" };
  return <span className={`inline-flex rounded-full px-2.5 py-1 text-[11px] font-semibold ${tones[tone]}`}>{children}</span>;
}

function Metric({ value, label }: { value: number | string; label: string }) {
  return <div><p className="text-xl font-semibold tracking-tight text-[#173b6c]">{typeof value === "number" ? n.format(value) : value}</p><p className="mt-0.5 text-[11px] font-medium text-[#6c7f89]">{label}</p></div>;
}

function Pager({ page, pages, setPage }: { page: number; pages: number; setPage: (p:number)=>void }) {
  if (pages <= 1) return null;
  return <div className="flex items-center justify-between border-t border-[#e2eaed] px-4 py-3 text-xs text-[#657984]">
    <span>Page {page + 1} of {pages}</span><div className="flex gap-2"><button aria-label="Previous page" disabled={page === 0} onClick={() => setPage(page - 1)} className="page-button"><ChevronLeft className="h-4 w-4" /></button><button aria-label="Next page" disabled={page >= pages - 1} onClick={() => setPage(page + 1)} className="page-button"><ChevronRight className="h-4 w-4" /></button></div>
  </div>;
}

export default function Home() {
  const [data, setData] = useState<IndexData | null>(null);
  const [query, setQuery] = useState("");
  const [activeClass, setActiveClass] = useState("All");
  const [activeTier, setActiveTier] = useState("All");
  const [selected, setSelected] = useState<Protein | null>(null);
  const [detail, setDetail] = useState<Detail | null>(null);
  const [loadingDetail, setLoadingDetail] = useState(false);
  const [pairQuery, setPairQuery] = useState("");
  const [pairTier, setPairTier] = useState("All");
  const [pairPage, setPairPage] = useState(0);
  const [interactionMode, setInteractionMode] = useState<"positive" | "negative">("positive");
  const [negativePage, setNegativePage] = useState(0);
  const [negativeRows, setNegativeRows] = useState<NegativePair[]>([]);
  const [loadingNegative, setLoadingNegative] = useState(false);
  const [diseaseQuery, setDiseaseQuery] = useState("");
  const [diseasePage, setDiseasePage] = useState(0);

  useEffect(() => {
    fetch("/data/mempro-index.json").then(r => r.json()).then((index: IndexData) => {
      setData(index);
      const requested = new URLSearchParams(window.location.search).get("protein");
      const initial = index.proteins.find(p => p.accession === requested || p.symbol === requested) ?? index.proteins.find(p => p.symbol === "DRD2") ?? index.proteins[0];
      setSelected(initial);
    });
  }, []);

  useEffect(() => {
    if (!selected) return;
    setLoadingDetail(true); setDetail(null); setPairPage(0); setNegativePage(0); setInteractionMode("positive"); setNegativeRows([]); setDiseasePage(0); setPairQuery(""); setDiseaseQuery("");
    fetch(`/data/proteins/${selected.accession}.json`).then(r => { if (!r.ok) throw new Error("detail unavailable"); return r.json(); }).then((d: Detail) => setDetail(d)).finally(() => setLoadingDetail(false));
    const url = new URL(window.location.href); url.searchParams.set("protein", selected.accession); window.history.replaceState({}, "", url);
  }, [selected]);

  useEffect(() => {
    if (!selected || !detail || interactionMode !== "negative" || detail.negativePairs.pageCount === 0) return;
    setLoadingNegative(true);
    fetch(`/data/negative-pairs/${selected.accession}/${negativePage + 1}.json`)
      .then(r => { if (!r.ok) throw new Error("negative pair page unavailable"); return r.json(); })
      .then((payload: { rows: NegativePairTuple[] }) => setNegativeRows(payload.rows.map(decodeNegativePair)))
      .finally(() => setLoadingNegative(false));
  }, [selected, detail, interactionMode, negativePage]);

  const results = useMemo(() => {
    if (!data) return [];
    const q = query.trim().toLowerCase();
    return data.proteins.filter(p => activeClass === "All" || p.membraneClass === activeClass)
      .filter(p => activeTier === "All" || p.evidenceTier === activeTier)
      .filter(p => !q || [p.symbol,p.accession,p.name].some(v => v.toLowerCase().includes(q)));
  }, [data, query, activeClass, activeTier]);

  const pairs = useMemo(() => (detail?.pairs ?? []).map(decodePair).filter(p => pairTier === "All" || p.tier === pairTier).filter(p => { const q = pairQuery.trim().toLowerCase(); return !q || [p.compoundName,p.compoundId,p.inchikey].some(v => v.toLowerCase().includes(q)); }), [detail, pairQuery, pairTier]);
  const pairPages = Math.max(1, Math.ceil(pairs.length / PAGE_SIZE));
  const pairRows = pairs.slice(pairPage * PAGE_SIZE, (pairPage + 1) * PAGE_SIZE);
  const diseases = useMemo(() => (detail?.diseases ?? []).filter(d => { const q=diseaseQuery.trim().toLowerCase(); return !q || [d.name,d.id,d.sourceIds].some(v => v.toLowerCase().includes(q)); }), [detail,diseaseQuery]);
  const diseasePages = Math.max(1, Math.ceil(diseases.length / PAGE_SIZE));
  const diseaseRows = diseases.slice(diseasePage * PAGE_SIZE, (diseasePage + 1) * PAGE_SIZE);

  if (!data) return <main className="grid min-h-screen place-items-center bg-[#f5f8fa]"><LoaderCircle className="h-7 w-7 animate-spin text-[#168c8c]" /></main>;

  return <main className="min-h-screen bg-[#f5f8fa] text-[#243b47]">
    <header className="sticky top-0 z-40 border-b border-[#d9e3e7] bg-white/95 backdrop-blur-xl"><div className="mx-auto flex h-15 max-w-[1600px] items-center gap-6 px-5 lg:px-8"><a href="#overview" className="flex items-center gap-2.5"><span className="grid h-8 w-8 place-items-center rounded-lg bg-[#173b6c]"><Layers3 className="h-4.5 w-4.5 text-white" /></span><span className="text-lg font-bold tracking-tight text-[#173b6c]">MemPro</span></a><div className="hidden h-5 w-px bg-[#d9e3e7] md:block" /><p className="hidden text-xs font-medium text-[#647985] md:block">Membrane protein–small molecule knowledgebase</p><div className="ml-auto flex items-center gap-3"><a href="/data-access" className="hidden text-xs font-semibold text-[#365f72] hover:text-[#168c8c] sm:inline">Data downloads</a><span className="hidden text-xs text-[#71838d] md:inline">Formal release</span><Badge tone="teal">7,800 proteins</Badge></div></div></header>

    <div className="mx-auto grid max-w-[1600px] lg:grid-cols-[330px_minmax(0,1fr)]">
      <aside className="border-r border-[#dce5e8] bg-white lg:sticky lg:top-15 lg:h-[calc(100vh-60px)] lg:overflow-hidden">
        <div className="border-b border-[#e3eaed] p-5"><p className="section-kicker">Protein directory</p><div className="relative mt-3"><Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-[#7f919b]" /><input value={query} onChange={e=>setQuery(e.target.value)} placeholder="Gene, UniProt or protein name" className="search-input pl-9" /></div><div className="mt-3 grid grid-cols-2 gap-2"><select aria-label="Membrane class" value={activeClass} onChange={e=>setActiveClass(e.target.value)} className="select-input"><option value="All">All classes</option><option value="A">Class A</option><option value="B">Class B</option><option value="C">Class C</option></select><select aria-label="Evidence tier" value={activeTier} onChange={e=>setActiveTier(e.target.value)} className="select-input"><option value="All">All evidence</option><option>E1</option><option>E2</option><option>E3</option></select></div><p className="mt-3 text-xs text-[#70838d]">{n.format(results.length)} matching proteins</p></div>
        <div className="protein-list lg:h-[calc(100vh-224px)] lg:overflow-y-auto">{results.slice(0,80).map(p => <button key={p.accession} onClick={()=>setSelected(p)} className={`protein-item ${selected?.accession===p.accession ? "active" : ""}`}><span className="min-w-0"><span className="block font-semibold text-[#173b6c]">{p.symbol} <span className="font-normal text-[#778993]">· {p.accession}</span></span><span className="mt-0.5 block truncate text-xs text-[#667b86]">{p.name}</span><span className="mt-2 flex gap-1.5"><span className={`rounded px-1.5 py-0.5 text-[10px] font-bold ${classInfo[p.membraneClass].tone}`}>{p.membraneClass}</span><span className="rounded bg-[#eef2f4] px-1.5 py-0.5 text-[10px] font-bold text-[#5d707a]">{p.evidenceTier}</span></span></span><span className="text-right text-[10px] leading-5 text-[#71848e]"><b className="block text-sm text-[#36566a]">{n.format(p.pairCount)}</b>pairs</span></button>)}</div>
      </aside>

      <article className="min-w-0">
        {selected && <>
          <section id="overview" className="record-hero"><div className="mx-auto max-w-[1180px] px-5 py-8 lg:px-9 lg:py-10"><div className="flex flex-wrap items-start justify-between gap-6"><div className="max-w-3xl"><div className="flex flex-wrap items-center gap-2"><Badge tone="blue">Class {selected.membraneClass} · {classInfo[selected.membraneClass].title}</Badge><Badge tone={selected.evidenceTier==="E1"?"teal":selected.evidenceTier==="E2"?"gold":"gray"}>{selected.evidenceTier}</Badge></div><h1 className="mt-4 text-4xl font-semibold tracking-[-0.035em] text-[#173b6c] sm:text-5xl">{selected.symbol}</h1><p className="mt-2 text-lg text-[#536d7b]">{selected.name}</p><p className="mt-3 font-mono text-sm text-[#70838d]">UniProt {selected.accession}</p></div><div className="grid grid-cols-2 gap-x-8 gap-y-4 rounded-2xl border border-[#d7e3e6] bg-white/75 px-5 py-4 shadow-sm sm:grid-cols-5"><Metric value={selected.sequenceLength} label="amino acids" /><Metric value={selected.pairCount} label="positive pairs" /><Metric value={selected.negativePairCount} label="negative pairs" /><Metric value={selected.diseaseCount} label="disease relations" /><Metric value={selected.siteCount} label="site assertions" /></div></div></div></section>

          <nav className="record-nav"><div className="mx-auto flex max-w-[1180px] gap-1 overflow-x-auto px-5 lg:px-9">{[["overview","Overview"],["interactions","Interactions"],["diseases","Diseases"],["localization","Localization"],["sequence","Sequence"]].map(([id,label])=><a key={id} href={`#${id}`} className="record-nav-link">{label}</a>)}</div></nav>

          {loadingDetail && <div className="grid min-h-[420px] place-items-center"><div className="text-center"><LoaderCircle className="mx-auto h-7 w-7 animate-spin text-[#168c8c]" /><p className="mt-3 text-sm text-[#6c808a]">Loading the complete record…</p></div></div>}
          {detail && <div className="mx-auto max-w-[1180px] space-y-8 px-5 py-8 lg:px-9">
            <section className="content-section reveal"><div className="section-heading"><div><p className="section-kicker">Identity and classification</p><h2>Protein overview</h2></div><Dna className="section-icon" /></div><div className="grid gap-6 lg:grid-cols-[1.15fr_0.85fr]"><div className="definition-grid"><div><dt>Membrane mechanism</dt><dd>{detail.membraneMode.replaceAll("_"," ")}</dd></div><div><dt>Membrane evidence</dt><dd>{detail.evidenceTier} · {detail.decisionEvidence || "source-traceable classification"}</dd></div><div><dt>Structural family</dt><dd>{detail.classification.family || "unresolved"}</dd></div><div><dt>Membrane role</dt><dd>{detail.classification.membraneRole || "unresolved"}</dd></div><div><dt>Molecular function</dt><dd>{detail.classification.molecularFunction || "unresolved"}</dd></div><div><dt>Biological process</dt><dd>{detail.classification.biologicalProcess || "unresolved"}</dd></div></div><div className="rounded-2xl bg-[#f3f7f8] p-5"><p className="text-xs font-bold uppercase tracking-[0.1em] text-[#6d808a]">Identifiers</p><dl className="mt-4 space-y-3 text-sm"><div className="flex justify-between gap-4"><dt>HGNC</dt><dd className="text-right font-medium text-[#173b6c]">{detail.gene.hgnc || "—"}</dd></div><div className="flex justify-between gap-4"><dt>Ensembl</dt><dd className="max-w-[65%] text-right font-medium text-[#173b6c]">{detail.gene.ensembl || "—"}</dd></div><div className="flex justify-between gap-4"><dt>NCBI Gene</dt><dd className="text-right font-medium text-[#173b6c]">{detail.gene.ncbi || "—"}</dd></div><div className="flex justify-between gap-4"><dt>Form scope</dt><dd className="text-right font-medium text-[#173b6c]">{detail.formSpecificity.replaceAll("_"," ") || "—"}</dd></div></dl></div></div></section>

            <section id="interactions" className="content-section reveal scroll-mt-28">
              <div className="section-heading"><div><p className="section-kicker">Protein–compound evidence</p><h2>Interactions</h2></div><FlaskConical className="section-icon" /></div>
              <div className="evidence-profile" aria-label="Positive and negative pair counts">
                <button className={interactionMode === "positive" ? "active" : ""} onClick={() => setInteractionMode("positive")}><span>Positive pairs</span><b>{n.format(detail.pairCount)}</b><i style={{width:`${Math.max(3, detail.pairCount / Math.max(detail.pairCount, detail.negativePairCount) * 100)}%`}} /></button>
                <button className={interactionMode === "negative" ? "active negative" : "negative"} onClick={() => setInteractionMode("negative")}><span>Negative pairs</span><b>{n.format(detail.negativePairCount)}</b><i style={{width:`${Math.max(3, detail.negativePairCount / Math.max(detail.pairCount, detail.negativePairCount) * 100)}%`}} /></button>
              </div>

              {interactionMode === "positive" ? <>
                <div className="table-toolbar"><div className="relative flex-1"><Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-[#7d909a]" /><input value={pairQuery} onChange={e=>{setPairQuery(e.target.value);setPairPage(0)}} placeholder="Search compound name, MemPro ID or InChIKey" className="search-input pl-9" /></div><select value={pairTier} onChange={e=>{setPairTier(e.target.value);setPairPage(0)}} className="select-input w-32"><option value="All">All BE tiers</option><option>BE1</option><option>BE2</option><option>BE3</option></select><span className="whitespace-nowrap text-xs text-[#6e818b]">{n.format(pairs.length)} pairs</span></div>
                <div className="data-table"><div className="pair-head"><span>Compound</span><span>Evidence</span><span>Sources</span><span>Sites</span></div>{pairRows.map(p=><details key={p.pairId} className="data-row"><summary className="pair-grid"><span className="min-w-0"><b className="block truncate text-[#173b6c]">{p.compoundName}</b><small className="block truncate">{p.compoundId}{p.formula ? ` · ${p.formula}` : ""}</small></span><span><Badge tone={p.tier==="BE1"?"teal":p.tier==="BE2"?"blue":"gold"}>{p.tier}</Badge><small className="mt-1 block">{p.evidenceCount} records</small></span><span><b>{p.databaseCount}</b><small className="block truncate">{p.sources.join(" · ") || "source unavailable"}</small></span><span><b>{p.siteCount}</b><small className="block">{p.coordinateSiteCount} coordinate-complete</small></span></summary><div className="row-detail"><div><span>InChIKey</span><b>{p.inchikey || "—"}</b></div><div><span>Putative experiment lineages</span><b>{p.lineageCount}</b></div><div><span>Independent structures</span><b>{p.structureCount}</b></div><div><span>Evidence modalities</span><b>{p.modalityCount}</b></div><div><span>Membrane-side site labels</span><b>{p.siteSides || "not annotated"}</b></div><div><span>Release scope</span><b>{p.defaultWeb ? "primary" : "extended"}</b></div></div></details>)}</div><Pager page={pairPage} pages={pairPages} setPage={setPairPage} />
              </> : <>
                <div className="negative-note"><CircleDot className="h-4 w-4" /><span>Inactive or negative under the reported assay conditions. Positive evidence, when present, remains visible as a conflict.</span></div>
                {loadingNegative ? <div className="empty-state"><LoaderCircle className="mx-auto h-5 w-5 animate-spin" /></div> : detail.negativePairs.pairCount ? <div className="data-table"><div className="negative-head"><span>Compound</span><span>Negative evidence</span><span>Source</span><span>Status</span></div>{negativeRows.map(p=><details key={p.compoundId} className="data-row"><summary className="negative-grid"><span className="min-w-0"><b className="block truncate text-[#173b6c]">{p.compoundName}</b><small className="block truncate">{p.compoundId}</small></span><span><b>{p.evidenceCount}</b><small className="block">records</small></span><span><b>{detail.negativePairs.source}</b></span><span>{p.conflict ? <Badge tone="rose">positive–negative conflict</Badge> : <Badge tone="gray">negative</Badge>}</span></summary><div className="row-detail"><div><span>InChIKey</span><b>{p.inchikey || "—"}</b></div><div><span>Interpretation</span><b>{p.conflict ? "Context-dependent evidence; both signs retained" : "Negative in the reported assay context"}</b></div></div></details>)}</div> : <div className="empty-state">No formal negative pair for this protein.</div>}
                <Pager page={negativePage} pages={Math.max(1, detail.negativePairs.pageCount)} setPage={setNegativePage} />
              </>}
            </section>

            <section id="diseases" className="content-section reveal scroll-mt-28"><div className="section-heading"><div><p className="section-kicker">Canonical disease context</p><h2>Disease associations</h2><p>MONDO identity, original source IDs, evidence channels and source databases are preserved together.</p></div><Stethoscope className="section-icon" /></div><div className="table-toolbar"><div className="relative flex-1"><Search className="absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-[#7d909a]" /><input value={diseaseQuery} onChange={e=>{setDiseaseQuery(e.target.value);setDiseasePage(0)}} placeholder="Search disease name, MONDO or source ID" className="search-input pl-9" /></div><span className="text-xs text-[#6e818b]">{n.format(diseases.length)} relations</span></div><div className="data-table"><div className="disease-head"><span>Disease</span><span>Level</span><span>Evidence channels</span><span>Sources</span></div>{diseaseRows.length ? diseaseRows.map(d=><details key={d.relationId} className="data-row"><summary className="disease-grid"><span><b className="block text-[#173b6c]">{d.name}</b><small>{d.id || "source-only disease"}</small></span><span><Badge tone={d.level==="very high"?"teal":d.level==="high"?"blue":"gold"}>{d.level}</Badge></span><span className="tag-list">{d.channels.slice(0,3).map(c=><i key={c}>{c.replaceAll("_"," ")}</i>)}</span><span><b>{d.sources.join(" · ") || "—"}</b><small className="block">{d.sourceIds}</small></span></summary><div className="row-detail"><div><span>Canonical mapping</span><b>{d.mapping || "unresolved"}</b></div><div><span>Open Targets score</span><b>{d.otScore ?? "—"}</b></div><div><span>PubMed IDs</span><b>{d.pubmed.join(" · ") || "—"}</b></div><div><span>Relation ID</span><b>{d.relationId}</b></div></div></details>) : <div className="empty-state">No included disease association for this protein.</div>}</div><Pager page={diseasePage} pages={diseasePages} setPage={setDiseasePage} /></section>

            <section id="localization" className="content-section reveal scroll-mt-28"><div className="section-heading"><div><p className="section-kicker">Observed cellular context</p><h2>Subcellular localization</h2></div><MapPin className="section-icon" /></div>{detail.localizations.length ? <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-3">{detail.localizations.map((l,i)=><div key={`${l.term}-${i}`} className="localization-card"><CircleDot className="h-4 w-4 text-[#168c8c]" /><div><b>{l.term}</b><p>{[l.go,l.role,l.reliability].filter(Boolean).join(" · ")}</p><small>{l.source}{l.version ? ` · ${l.version}` : ""}</small></div></div>)}</div> : <div className="empty-state">No mapped HPA subcellular localization record.</div>}</section>

            <section id="sequence" className="content-section reveal scroll-mt-28"><div className="section-heading"><div><p className="section-kicker">Canonical UniProt sequence</p><h2>Protein sequence</h2><p>{n.format(detail.sequence.length)} amino acids</p></div><BookOpen className="section-icon" /></div><div className="sequence-view"><div className="sequence-ruler"><span>1</span><span>{detail.sequence.length}</span></div><div className="sequence-text">{detail.sequence.match(/.{1,10}/g)?.map((chunk,i)=><span key={i} title={`Residues ${i*10+1}–${Math.min((i+1)*10,detail.sequence.length)}`}>{chunk}</span>)}</div></div></section>
          </div>}
        </>}
      </article>
    </div>

    <footer className="border-t border-[#dce5e8] bg-white"><div className="mx-auto flex max-w-[1600px] flex-wrap justify-between gap-3 px-5 py-5 text-xs text-[#70838d] lg:px-8"><span>{data.release}</span><span>{n.format(data.statistics.formalPairs)} pairs · {n.format(data.statistics.diseaseRelations)} disease relations · {n.format(data.statistics.sites)} site assertions</span></div></footer>
  </main>;
}
