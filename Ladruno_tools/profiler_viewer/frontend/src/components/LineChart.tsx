import { useEffect, useRef, useState } from 'react'

export interface LineSeries {
  label: string
  color: string
  y: number[]
}

interface Props {
  title: string
  x: number[]
  series: LineSeries[]
  height?: number
  yFormat?: (v: number) => string
  xLabel?: string
}

const M = { top: 8, right: 12, bottom: 26, left: 56 }
const MAX_POINTS = 1500 // downsample long series so the SVG path stays light

function niceTicks(min: number, max: number, count: number): number[] {
  if (min === max) return [min]
  const span = max - min
  const step0 = span / count
  const mag = Math.pow(10, Math.floor(Math.log10(step0)))
  const norm = step0 / mag
  const step = (norm >= 5 ? 5 : norm >= 2 ? 2 : 1) * mag
  const start = Math.ceil(min / step) * step
  const out: number[] = []
  for (let v = start; v <= max + step * 0.5; v += step) out.push(v)
  return out
}

// Generic responsive SVG multi-line chart with simple axes + legend. Pure SVG,
// no charting dependency.
export function LineChart({ title, x, series, height = 200, yFormat, xLabel }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(760)

  useEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver((e) => {
      const w = e[0]?.contentRect.width
      if (w && w > 0) setWidth(w)
    })
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  const n = x.length
  const stride = n > MAX_POINTS ? Math.ceil(n / MAX_POINTS) : 1

  const xMin = x.length ? x[0] : 0
  const xMax = x.length ? x[n - 1] : 1
  let yMin = 0
  let yMax = 0
  for (const s of series) for (const v of s.y) if (v > yMax) yMax = v
  if (yMax === yMin) yMax = 1

  const pw = Math.max(10, width - M.left - M.right)
  const ph = Math.max(10, height - M.top - M.bottom)
  const sx = (v: number) => M.left + (xMax === xMin ? 0 : (v - xMin) / (xMax - xMin)) * pw
  const sy = (v: number) => M.top + ph - (v - yMin) / (yMax - yMin) * ph

  const fmtY = yFormat ?? ((v: number) => v.toPrecision(3))
  const yTicks = niceTicks(yMin, yMax, 4)
  const xTicks = niceTicks(xMin, xMax, 6)

  const pathFor = (s: LineSeries): string => {
    let d = ''
    for (let i = 0; i < n; i += stride) {
      d += `${d ? 'L' : 'M'}${sx(x[i]).toFixed(1)} ${sy(s.y[i] ?? 0).toFixed(1)}`
    }
    return d
  }

  return (
    <div className="chart" ref={wrapRef}>
      <div className="chart-title">{title}</div>
      <svg width={width} height={height} role="img" aria-label={title}>
        {/* y grid + labels */}
        {yTicks.map((t) => (
          <g key={`y${t}`}>
            <line x1={M.left} y1={sy(t)} x2={width - M.right} y2={sy(t)} className="grid" />
            <text x={M.left - 6} y={sy(t) + 3} className="axis-label num-anchor">{fmtY(t)}</text>
          </g>
        ))}
        {/* x labels */}
        {xTicks.map((t) => (
          <text key={`x${t}`} x={sx(t)} y={height - 8} className="axis-label mid-anchor">
            {t.toLocaleString('en-US')}
          </text>
        ))}
        {xLabel && <text x={(M.left + width - M.right) / 2} y={height} className="axis-cap mid-anchor">{xLabel}</text>}
        {/* series */}
        {series.map((s) => (
          <path key={s.label} d={pathFor(s)} fill="none" stroke={s.color} strokeWidth={1.4} />
        ))}
      </svg>
      <div className="legend">
        {series.map((s) => (
          <span key={s.label} className="legend-item">
            <span className="swatch" style={{ background: s.color }} />
            {s.label}
          </span>
        ))}
      </div>
    </div>
  )
}
