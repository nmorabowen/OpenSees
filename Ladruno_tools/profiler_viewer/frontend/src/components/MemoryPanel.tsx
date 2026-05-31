import type { MemorySnapshot } from '../api'
import { fmtBytes, fmtInt } from '../format'

interface Props {
  mem: MemorySnapshot
  // optional baseline run (the "vs" selection) for a leak check via census diff
  base?: MemorySnapshot | null
  baseLabel?: string
}

// End-of-run memory: per-type live byte counters + peak, and the TaggedObject
// live-component census (classTag -> count) as a sorted bar list. When a baseline
// run is given, a census diff flags classTags that grew (a leak signal: more
// objects live than the comparison run).
export function MemoryPanel({ mem, base, baseLabel }: Props) {
  const census = [...mem.components_live].sort((a, b) => b.count - a.count)
  const maxCount = census.reduce((m, c) => Math.max(m, c.count), 0) || 1

  const baseCount = new Map((base?.components_live ?? []).map((c) => [c.classTag, c.count]))
  const delta = (tag: number, count: number) => (base ? count - (baseCount.get(tag) ?? 0) : 0)
  const grew = base ? census.filter((c) => delta(c.classTag, c.count) > 0) : []

  return (
    <div className="memory">
      <div className="mem-cards">
        <div className="mem-card"><span className="mem-k">matrix_live</span><span className="mem-v">{fmtBytes(mem.matrix_live)}</span></div>
        <div className="mem-card"><span className="mem-k">vector_live</span><span className="mem-v">{fmtBytes(mem.vector_live)}</span></div>
        <div className="mem-card"><span className="mem-k">id_live</span><span className="mem-v">{fmtBytes(mem.id_live)}</span></div>
        <div className="mem-card peak"><span className="mem-k">peak</span><span className="mem-v">{fmtBytes(mem.peak_bytes)}</span></div>
      </div>

      {base && (
        <div className={`leak-badge ${grew.length ? 'warn' : 'good'}`}>
          {grew.length
            ? `⚠ leak check vs ${baseLabel ?? 'base'}: ${grew.length} class(es) have more live objects `
              + `(${grew.map((c) => `${c.classTag}:+${delta(c.classTag, c.count)}`).join(', ')})`
            : `✓ leak check vs ${baseLabel ?? 'base'}: no class has more live objects`}
        </div>
      )}

      <div className="details-sub">live-component census (by classTag)</div>
      {census.length === 0 ? (
        <div className="muted">No census recorded (run without <code>-memory</code>, or armed after the model was built).</div>
      ) : (
        <table className="censustab">
          <thead>
            <tr>
              <th>classTag</th>
              <th className="num">live count</th>
              {base && <th className="num">Δ vs base</th>}
              <th></th>
            </tr>
          </thead>
          <tbody>
            {census.map((c) => {
              const d = delta(c.classTag, c.count)
              return (
                <tr key={c.classTag}>
                  <td className="mono">{c.classTag}</td>
                  <td className="num">{fmtInt(c.count)}</td>
                  {base && (
                    <td className={`num ${d > 0 ? 'slower' : d < 0 ? 'faster' : 'flat'}`}>
                      {d > 0 ? '+' : ''}{d === 0 ? '—' : fmtInt(d)}
                    </td>
                  )}
                  <td className="barcell">
                    <span className="bar" style={{ width: `${(c.count / maxCount) * 100}%` }} />
                  </td>
                </tr>
              )
            })}
          </tbody>
        </table>
      )}
    </div>
  )
}
