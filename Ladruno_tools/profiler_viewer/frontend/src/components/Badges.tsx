import type { Header } from '../api'
import { fmtInt, fmtSci } from '../format'

interface Props {
  header: Header
}

type Tone = 'info' | 'warn' | 'good'

function Badge({ label, value, tone = 'info', title }: {
  label: string
  value: string
  tone?: Tone
  title?: string
}) {
  return (
    <span className={`badge ${tone}`} title={title}>
      <span className="badge-label">{label}</span>
      <span className="badge-value">{value}</span>
    </span>
  )
}

// Header badges: the at-a-glance run shape + the explicit-dynamics levers
// (dt vs dt_cr, oversample) the FEM panel flagged as the highest-value signals.
export function Badges({ header }: Props) {
  const { integrator, algorithm, solver, nDOF, nElem, nSteps, threads,
          dt_min, dt_cr, oversample_ratio } = header

  // oversample_ratio = dt_cr / dt: >1 means the step is finer than critical, so
  // dt could be raised (a cheap explicit speed-up). Flag it.
  const over = oversample_ratio
  const overTone: Tone = over > 1.5 ? 'warn' : over > 0 ? 'good' : 'info'
  const overTitle =
    over > 1.5
      ? `Running ~${over.toFixed(1)}× finer than the critical step — dt could likely be raised.`
      : over > 0
        ? 'Step size is near the critical step (well-tuned).'
        : 'No critical-step estimate reported for this run.'

  return (
    <div className="badges">
      {(integrator || algorithm || solver) && (
        <Badge label="scheme" value={`${integrator || '—'} · ${algorithm || '—'} · ${solver || '—'}`} />
      )}
      <Badge label="nDOF" value={fmtInt(nDOF)} />
      {nElem > 0 && <Badge label="nElem" value={fmtInt(nElem)} />}
      {nSteps > 0 && <Badge label="nSteps" value={fmtInt(nSteps)} />}
      {threads > 0 && <Badge label="threads" value={fmtInt(threads)} />}
      {dt_cr > 0 && (
        <Badge label="dt / dt_cr" value={`${fmtSci(dt_min)} / ${fmtSci(dt_cr)}`}
               tone={overTone} title={overTitle} />
      )}
      {over > 0 && (
        <Badge label="oversample" value={`${over.toFixed(2)}×`} tone={overTone} title={overTitle} />
      )}
    </div>
  )
}
