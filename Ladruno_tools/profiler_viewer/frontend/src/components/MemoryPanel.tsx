import type { MemorySnapshot } from '../api'
import { fmtBytes, fmtInt } from '../format'

interface Props {
  mem: MemorySnapshot
}

// End-of-run memory: per-type live byte counters + peak, and the TaggedObject
// live-component census (classTag -> count) as a sorted bar list.
export function MemoryPanel({ mem }: Props) {
  const census = [...mem.components_live].sort((a, b) => b.count - a.count)
  const maxCount = census.reduce((m, c) => Math.max(m, c.count), 0) || 1

  return (
    <div className="memory">
      <div className="mem-cards">
        <div className="mem-card"><span className="mem-k">matrix_live</span><span className="mem-v">{fmtBytes(mem.matrix_live)}</span></div>
        <div className="mem-card"><span className="mem-k">vector_live</span><span className="mem-v">{fmtBytes(mem.vector_live)}</span></div>
        <div className="mem-card"><span className="mem-k">id_live</span><span className="mem-v">{fmtBytes(mem.id_live)}</span></div>
        <div className="mem-card peak"><span className="mem-k">peak</span><span className="mem-v">{fmtBytes(mem.peak_bytes)}</span></div>
      </div>

      <div className="details-sub">live-component census (by classTag)</div>
      {census.length === 0 ? (
        <div className="muted">No census recorded (run without <code>-memory</code>, or armed after the model was built).</div>
      ) : (
        <table className="censustab">
          <thead>
            <tr><th>classTag</th><th className="num">live count</th><th></th></tr>
          </thead>
          <tbody>
            {census.map((c) => (
              <tr key={c.classTag}>
                <td className="mono">{c.classTag}</td>
                <td className="num">{fmtInt(c.count)}</td>
                <td className="barcell">
                  <span className="bar" style={{ width: `${(c.count / maxCount) * 100}%` }} />
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </div>
  )
}
