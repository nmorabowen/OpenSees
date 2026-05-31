// Small number/label formatters shared across panels.

export function fmtMs(ms: number): string {
  if (ms === 0) return '0'
  if (ms < 0.001) return `${(ms * 1e6).toFixed(0)} ns`
  if (ms < 1) return `${(ms * 1e3).toFixed(1)} µs`
  if (ms < 1000) return `${ms.toFixed(2)} ms`
  return `${(ms / 1000).toFixed(3)} s`
}

export function fmtPct(frac: number): string {
  return `${(frac * 100).toFixed(1)}%`
}

export function fmtBytes(b: number): string {
  if (b < 1024) return `${b} B`
  if (b < 1024 * 1024) return `${(b / 1024).toFixed(1)} KiB`
  if (b < 1024 * 1024 * 1024) return `${(b / (1024 * 1024)).toFixed(2)} MiB`
  return `${(b / (1024 * 1024 * 1024)).toFixed(2)} GiB`
}

export function fmtInt(n: number): string {
  return n.toLocaleString('en-US')
}

// Compact scientific for dt-like small values.
export function fmtSci(x: number): string {
  if (x === 0) return '0'
  if (x >= 1e-3 && x < 1e4) return x.toPrecision(3)
  return x.toExponential(2)
}

// Deterministic color from a string (used to tint icicle frames by node name).
export function colorFor(name: string, depth: number): string {
  let h = 0
  for (let i = 0; i < name.length; i++) h = (h * 31 + name.charCodeAt(i)) >>> 0
  const hue = (h % 360 + depth * 47) % 360
  return `hsl(${hue} 55% ${42 + (depth % 3) * 5}%)`
}
