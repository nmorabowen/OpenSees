import type { Series } from '../api'
import { fmtMs, fmtSci } from '../format'
import { LineChart, type LineSeries } from './LineChart'

interface Props {
  series: Series
}

const PALETTE = ['#5aa9e6', '#ffb454', '#7ed957', '#ff6b6b', '#c792ea', '#4dd0e1', '#f06292']

// Per-step time-history panel: one multi-line chart for per-phase wall time,
// plus iterations-per-step and dt-per-step.
export function TimeSeries({ series }: Props) {
  const x = series.step

  const phaseSeries: LineSeries[] = series.phases.map((p, j) => ({
    label: p,
    color: PALETTE[j % PALETTE.length],
    y: series.wall_ms.map((row) => row[j] ?? 0),
  }))

  return (
    <div className="series">
      {phaseSeries.length > 0 && (
        <LineChart title="per-phase wall time / step" x={x} series={phaseSeries}
                   yFormat={fmtMs} xLabel="step" height={240} />
      )}
      {series.iters.some((v) => v > 0) && (
        <LineChart title="iterations / step" x={x}
                   series={[{ label: 'iters', color: '#ffb454', y: series.iters }]}
                   yFormat={(v) => v.toFixed(0)} xLabel="step" height={160} />
      )}
      {series.dt.some((v) => v > 0) && (
        <LineChart title="dt / step" x={x}
                   series={[{ label: 'dt', color: '#7ed957', y: series.dt }]}
                   yFormat={fmtSci} xLabel="step" height={160} />
      )}
    </div>
  )
}
