import { useCallback, useEffect, useState } from 'react'
import { api, ApiError, API_BASE } from './api'
import type { Header, ManifestRow, MemorySnapshot, RollupNode, Series } from './api'
import type { DiffResponse } from './api'
import { RunPicker } from './components/RunPicker'
import { Badges } from './components/Badges'
import { Icicle } from './components/Icicle'
import { NodeDetails } from './components/NodeDetails'
import { TimeSeries } from './components/TimeSeries'
import { MemoryPanel } from './components/MemoryPanel'
import { DiffTable } from './components/DiffTable'

type Tab = 'rollup' | 'series' | 'memory' | 'diff'

export default function App() {
  const [runs, setRuns] = useState<ManifestRow[]>([])
  const [fatal, setFatal] = useState<string | null>(null)
  const [schema, setSchema] = useState<string>('')

  const [primary, setPrimary] = useState<string | null>(null)
  const [base, setBase] = useState<string | null>(null)
  const [tab, setTab] = useState<Tab>('rollup')

  const [header, setHeader] = useState<Header | null>(null)
  const [rollup, setRollup] = useState<RollupNode | null>(null)
  const [series, setSeries] = useState<Series | null>(null)
  const [memory, setMemory] = useState<MemorySnapshot | null>(null)
  const [selected, setSelected] = useState<RollupNode | null>(null)
  const [diff, setDiff] = useState<DiffResponse | null>(null)

  // -- load the run list (and backend health) ---------------------------------
  const loadRuns = useCallback(async () => {
    setFatal(null)
    try {
      const h = await api.health()
      setSchema(h.schema_version)
      const rows = await api.runs()
      setRuns(rows)
      setPrimary((cur) => cur ?? (rows.length ? rows[0].run : null))
    } catch (e) {
      setFatal(e instanceof ApiError ? e.message : String(e))
    }
  }, [])

  useEffect(() => {
    void loadRuns()
  }, [loadRuns])

  // -- load everything for the primary run ------------------------------------
  useEffect(() => {
    if (!primary) return
    let cancelled = false
    setSelected(null)
    void (async () => {
      try {
        const [hdr, roll, mem] = await Promise.all([
          api.header(primary),
          api.rollup(primary),
          api.memory(primary),
        ])
        if (cancelled) return
        setHeader(hdr)
        setRollup(roll)
        setSelected(roll)
        setMemory(mem)
      } catch (e) {
        if (!cancelled) setFatal(e instanceof ApiError ? e.message : String(e))
      }
      // series is optional (404 when the run had no -perStep buffer)
      try {
        const s = await api.series(primary)
        if (!cancelled) setSeries(s)
      } catch {
        if (!cancelled) setSeries(null)
      }
    })()
    return () => {
      cancelled = true
    }
  }, [primary])

  // -- load the diff when both runs are chosen and the diff tab is active -----
  useEffect(() => {
    if (tab !== 'diff' || !primary || !base || base === primary) {
      setDiff(null)
      return
    }
    let cancelled = false
    void (async () => {
      try {
        const d = await api.diff(base, primary)
        if (!cancelled) setDiff(d)
      } catch (e) {
        if (!cancelled) setFatal(e instanceof ApiError ? e.message : String(e))
      }
    })()
    return () => {
      cancelled = true
    }
  }, [tab, primary, base])

  const tabs: Tab[] = ['rollup', 'series', 'memory', 'diff']

  return (
    <div className="app">
      <header className="topbar">
        <div className="brand">
          Ladruno <span className="brand-sub">profiler</span>
        </div>
        <div className="topbar-meta">
          <span className="mono dim">{API_BASE}</span>
          {schema && <span className="dim"> · schema {schema}</span>}
          <button className="btn" onClick={() => void loadRuns()}>↻ runs</button>
        </div>
      </header>

      {fatal && (
        <div className="banner error">
          <strong>Backend error:</strong> {fatal}
          <div className="banner-hint">
            Start it with <code>python profiler_api.py --file profile.h5</code> (or set{' '}
            <code>VITE_API_BASE</code>), then click <em>↻ runs</em>.
          </div>
        </div>
      )}

      <section className="picker-wrap">
        <RunPicker rows={runs} primary={primary} base={base} onPrimary={setPrimary} onBase={setBase} />
      </section>

      {header && <Badges header={header} />}

      <nav className="tabs">
        {tabs.map((t) => (
          <button
            key={t}
            className={`tab ${t === tab ? 'active' : ''}`}
            onClick={() => setTab(t)}
            disabled={t === 'diff' && (!base || base === primary)}
            title={t === 'diff' && (!base || base === primary) ? 'pick a "vs" run to compare' : undefined}
          >
            {t}
          </button>
        ))}
      </nav>

      <main className="content">
        {tab === 'rollup' && rollup && (
          <div className="rollup-split">
            <Icicle root={rollup} selectedPath={selected?.path ?? null} onSelect={setSelected} />
            <NodeDetails node={selected} />
          </div>
        )}

        {tab === 'series' &&
          (series ? (
            <TimeSeries series={series} />
          ) : (
            <div className="muted pad">
              No per-step series for this run — re-run with <code>profiler('start','-perStep')</code>.
            </div>
          ))}

        {tab === 'memory' && memory && <MemoryPanel mem={memory} />}

        {tab === 'diff' &&
          (diff ? (
            <DiffTable diff={diff} />
          ) : (
            <div className="muted pad">Pick a “vs” run in the table above to compare against the viewed run.</div>
          ))}
      </main>
    </div>
  )
}
