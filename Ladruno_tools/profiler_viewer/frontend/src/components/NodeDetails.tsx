import type { RollupNode } from '../api'
import { fmtInt, fmtMs, fmtPct } from '../format'

interface Props {
  node: RollupNode | null
}

// Detail panel for the icicle-selected node: scalar metrics + the per-element
// -class breakdown (deep mode, P0#1) when present.
export function NodeDetails({ node }: Props) {
  if (!node) return <div className="details muted">Select a frame to inspect it.</div>

  const rows = node.elem_by_type
    ? [...node.elem_by_type].sort((a, b) => b.wall_ms - a.wall_ms)
    : []

  return (
    <div className="details">
      <div className="details-path mono">{node.path}</div>
      <dl className="metrics">
        <div><dt>wall</dt><dd>{fmtMs(node.wall_ms)}</dd></div>
        <div><dt>share</dt><dd>{fmtPct(node.share)}</dd></div>
        <div><dt>calls</dt><dd>{fmtInt(node.calls)}</dd></div>
        <div><dt>min / max</dt><dd>{fmtMs(node.wall_ms_min)} / {fmtMs(node.wall_ms_max)}</dd></div>
        <div><dt>/ step</dt><dd>{fmtMs(node.wall_ms_per_step)}</dd></div>
        <div><dt>/ DOF</dt><dd>{fmtMs(node.wall_ms_per_dof)}</dd></div>
      </dl>

      {rows.length > 0 && (
        <>
          <div className="details-sub">per-element-class (deep)</div>
          <table className="elemtab">
            <thead>
              <tr>
                <th>classTag</th>
                <th className="num">evals</th>
                <th className="num">wall</th>
                <th className="num">ns / eval</th>
                <th>form</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((e) => (
                <tr key={e.classTag}>
                  <td className="mono">{e.classTag}</td>
                  <td className="num">{fmtInt(e.count)}</td>
                  <td className="num">{fmtMs(e.wall_ms)}</td>
                  <td className="num">{e.wall_ns_per_elem.toFixed(0)}</td>
                  <td>
                    <span className={`chip ${e.fb_coupled ? 'fb' : 'disp'}`}>
                      {e.fb_coupled ? 'force-based' : 'disp-based'}
                    </span>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </div>
  )
}
