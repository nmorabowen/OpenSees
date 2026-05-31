import type { DiffResponse } from '../api'
import { fmtMs } from '../format'

interface Props {
  diff: DiffResponse
}

function depthOf(path: string): number {
  return path.split('/').length - 1
}

// Run-to-run diff: per-node wall delta along the stable rollup path (P0#6).
// Negative d_abs (candidate faster) is green; positive (slower) is red.
export function DiffTable({ diff }: Props) {
  return (
    <div className="diff">
      <div className="diff-head">
        <span className="mono">{diff.base}</span> → <span className="mono">{diff.cand}</span>
        <span className="muted"> (base → candidate; green = faster)</span>
      </div>
      <table className="difftab">
        <thead>
          <tr>
            <th>node</th>
            <th className="num">base</th>
            <th className="num">cand</th>
            <th className="num">Δ</th>
            <th className="num">Δ%</th>
            <th className="num">Δshare</th>
            <th>status</th>
          </tr>
        </thead>
        <tbody>
          {diff.rows.map((r) => {
            const tone = r.d_abs_ms < -1e-9 ? 'faster' : r.d_abs_ms > 1e-9 ? 'slower' : 'flat'
            return (
              <tr key={r.path}>
                <td className="mono" style={{ paddingLeft: `${4 + depthOf(r.path) * 14}px` }}>{r.name}</td>
                <td className="num">{fmtMs(r.base_ms)}</td>
                <td className="num">{fmtMs(r.cand_ms)}</td>
                <td className={`num ${tone}`}>{r.d_abs_ms >= 0 ? '+' : ''}{fmtMs(r.d_abs_ms)}</td>
                <td className={`num ${tone}`}>{r.d_pct === null ? '—' : `${r.d_pct >= 0 ? '+' : ''}${r.d_pct.toFixed(1)}%`}</td>
                <td className="num">{r.d_share >= 0 ? '+' : ''}{(r.d_share * 100).toFixed(1)}%</td>
                <td><span className={`chip ${r.status}`}>{r.status}</span></td>
              </tr>
            )
          })}
        </tbody>
      </table>
    </div>
  )
}
