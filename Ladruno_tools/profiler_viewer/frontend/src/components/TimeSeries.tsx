import type { Series } from '../api'
import { fmtMs, fmtSci } from '../format'
import { LineChart, type LineSeries } from './LineChart'

interface Props {
  series: Series
  // optional baseline run (the "vs" selection) overlaid as dashed/faded lines
  base?: Series | null
  baseLabel?: string
}

const PALETTE = ['#5aa9e6', '#ffb454', '#7ed957', '#ff6b6b', '#c792ea', '#4dd0e1', '#f06292']

// Per-step time-history panel: one multi-line chart for per-phase wall time,
// plus iterations-per-step and dt-per-step. When a baseline run is given, its
// curves are overlaid dashed/faded for run-to-run comparison.
export function TimeSeries({ series, base, baseLabel }: Props) {
  const x = series.step
  const suffix = baseLabel ? ` (${baseLabel})` : ' (base)'

  const phaseSeries: LineSeries[] = series.phases.map((p, j) => ({
    label: p,
    color: PALETTE[j % PALETTE.length],
    y: series.wall_ms.map((row) => row[j] ?? 0),
  }))
  // overlay only the phases the baseline shares (matched by name → same color)
  if (base) {
    base.phases.forEach((p, bj) => {
      const j = series.phases.indexOf(p)
      if (j >= 0) {
        phaseSeries.push({
          label: p + suffix,
          color: PALETTE[j % PALETTE.length],
          x: base.step,
          y: base.wall_ms.map((row) => row[bj] ?? 0),
          dashed: true,
        })
      }
    })
  }

  const iterSeries: LineSeries[] = [{ label: 'iters', color: '#ffb454', y: series.iters }]
  if (base?.iters.some((v) => v > 0))
    iterSeries.push({ label: 'iters' + suffix, color: '#ffb454', x: base.step, y: base.iters, dashed: true })

  const dtSeries: LineSeries[] = [{ label: 'dt', color: '#7ed957', y: series.dt }]
  if (base?.dt.some((v) => v > 0))
    dtSeries.push({ label: 'dt' + suffix, color: '#7ed957', x: base.step, y: base.dt, dashed: true })

  return (
    <div className="series">
      {phaseSeries.length > 0 && (
        <LineChart title="per-phase wall time / step" x={x} series={phaseSeries}
                   yFormat={fmtMs} xLabel="step" height={240} />
      )}
      {series.iters.some((v) => v > 0) && (
        <LineChart title="iterations / step" x={x} series={iterSeries}
                   yFormat={(v) => v.toFixed(0)} xLabel="step" height={160} />
      )}
      {series.dt.some((v) => v > 0) && (
        <LineChart title="dt / step" x={x} series={dtSeries}
                   yFormat={fmtSci} xLabel="step" height={160} />
      )}
    </div>
  )
}
