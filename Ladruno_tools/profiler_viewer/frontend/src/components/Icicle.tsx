import { useEffect, useRef, useState } from 'react'
import type { RollupNode } from '../api'
import { colorFor, fmtMs, fmtPct } from '../format'

const ROW = 26 // px per depth level

interface Rect {
  node: RollupNode
  x: number
  w: number
  depth: number
}

// Flatten the rollup tree into positioned rectangles. Each child's width is
// proportional to its wall_ms against its PARENT's wall_ms, so a parent whose
// children don't account for all its time shows a trailing gap (= self time).
function layout(root: RollupNode, totalWidth: number): { rects: Rect[]; depth: number } {
  const rects: Rect[] = []
  let maxDepth = 0
  const walk = (node: RollupNode, x: number, w: number, depth: number) => {
    rects.push({ node, x, w, depth })
    if (depth > maxDepth) maxDepth = depth
    if (w <= 0 || node.wall_ms <= 0) return
    let cx = x
    for (const c of node.children) {
      const cw = (c.wall_ms / node.wall_ms) * w
      walk(c, cx, cw, depth + 1)
      cx += cw
    }
  }
  walk(root, 0, totalWidth, 0)
  return { rects, depth: maxDepth }
}

interface Props {
  root: RollupNode
  selectedPath: string | null
  onSelect: (node: RollupNode) => void
}

// SVG icicle (flame graph). Click a frame to select it (drives the detail panel).
export function Icicle({ root, selectedPath, onSelect }: Props) {
  const wrapRef = useRef<HTMLDivElement>(null)
  const [width, setWidth] = useState(900)

  useEffect(() => {
    const el = wrapRef.current
    if (!el) return
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width
      if (w && w > 0) setWidth(w)
    })
    ro.observe(el)
    return () => ro.disconnect()
  }, [])

  const { rects, depth } = layout(root, width)
  const height = (depth + 1) * ROW

  return (
    <div className="icicle" ref={wrapRef}>
      <svg width={width} height={height} role="img" aria-label="rollup flame graph">
        {rects.map((r) => {
          const sel = r.node.path === selectedPath
          const labelFits = r.w > 44
          const chars = Math.max(0, Math.floor(r.w / 6.6) - 1)
          const label = r.node.name.length > chars ? r.node.name.slice(0, chars) + '…' : r.node.name
          return (
            <g
              key={r.node.path}
              transform={`translate(${r.x.toFixed(2)},${r.depth * ROW})`}
              className="frame"
              onClick={() => onSelect(r.node)}
            >
              <title>
                {`${r.node.path}\n${fmtMs(r.node.wall_ms)}  ·  ${fmtPct(r.node.share)} of root  ·  ${r.node.calls} call(s)`}
              </title>
              <rect
                width={Math.max(0.5, r.w - 1)}
                height={ROW - 1}
                rx={2}
                fill={colorFor(r.node.name, r.depth)}
                stroke={sel ? '#ffd24a' : 'rgba(0,0,0,0.35)'}
                strokeWidth={sel ? 2 : 0.5}
              />
              {labelFits && (
                <text x={5} y={ROW / 2 + 4} className="frame-label">
                  {label}
                </text>
              )}
            </g>
          )
        })}
      </svg>
    </div>
  )
}
