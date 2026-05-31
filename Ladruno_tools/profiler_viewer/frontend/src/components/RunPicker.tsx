import type { ManifestRow } from '../api'
import { fmtInt, fmtMs } from '../format'

interface Props {
  rows: ManifestRow[]
  primary: string | null
  base: string | null // the "compare against" run for the diff view
  onPrimary: (run: string) => void
  onBase: (run: string | null) => void
}

// Run picker: one row per run. The "view" column selects the primary run shown
// in every panel; the "vs" column selects an optional baseline for the diff.
export function RunPicker({ rows, primary, base, onPrimary, onBase }: Props) {
  return (
    <table className="picker">
      <thead>
        <tr>
          <th>view</th>
          <th>vs</th>
          <th>run</th>
          <th>integrator / algorithm</th>
          <th className="num">nDOF</th>
          <th className="num">nElem</th>
          <th className="num">nSteps</th>
          <th className="num">wall</th>
        </tr>
      </thead>
      <tbody>
        {rows.map((r) => (
          <tr key={r.run} className={r.run === primary ? 'sel' : undefined}>
            <td>
              <input
                type="radio"
                name="primary"
                checked={r.run === primary}
                onChange={() => onPrimary(r.run)}
                aria-label={`view ${r.run}`}
              />
            </td>
            <td>
              <input
                type="radio"
                name="base"
                checked={r.run === base}
                // click an already-selected baseline to clear it
                onClick={() => onBase(r.run === base ? null : r.run)}
                onChange={() => onBase(r.run)}
                aria-label={`compare against ${r.run}`}
              />
            </td>
            <td className="mono">{r.run}</td>
            <td>
              {r.integrator || '—'} / {r.algorithm || '—'}
            </td>
            <td className="num">{fmtInt(r.nDOF)}</td>
            <td className="num">{fmtInt(r.nElem)}</td>
            <td className="num">{fmtInt(r.nSteps)}</td>
            <td className="num">{fmtMs(r.wall_ms_total)}</td>
          </tr>
        ))}
      </tbody>
    </table>
  )
}
