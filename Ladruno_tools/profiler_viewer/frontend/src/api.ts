// Typed client for the profiler FastAPI backend (P7, profiler_api.py).
// The types mirror the JSON wire format documented in profiler_results.py, so
// the contract is checked at compile time. The frontend never reads HDF5.

// Backend base URL. Defaults to the backend's own default port; override with a
// VITE_API_BASE env at build/dev time. The backend sets permissive CORS, so a
// cross-origin fetch from the Vite dev server works without a proxy.
export const API_BASE: string =
  (import.meta.env.VITE_API_BASE as string | undefined) ?? 'http://127.0.0.1:8000'

export interface Health {
  status: string
  file: string
  schema_version: string
  runs: string[]
}

// One picker row (GET /runs).
export interface ManifestRow {
  run: string
  model: string
  integrator: string
  algorithm: string
  solver: string
  nDOF: number
  nElem: number
  nSteps: number
  threads: number
  dt_min: number
  dt_cr: number
  oversample_ratio: number
  wall_ms_total: number
  timestamp: string
}

// Full run header (GET /runs/{run}/header).
export interface Header {
  run: string
  model: string
  engine_sha: string
  integrator: string
  algorithm: string
  solver: string
  units: string
  timestamp: string
  nSteps: number
  nDOF: number
  nElem: number
  nNode: number
  nnz: number
  threads: number
  dt_min: number
  dt_max: number
  dt_cr: number
  oversample_ratio: number
  wall_ms_total: number
  cpu_ms_total: number
}

// Per-element-class breakdown row (deep mode, P0#1).
export interface ElemRow {
  classTag: number
  count: number
  wall_ms: number
  wall_ns_per_elem: number
  fb_coupled: boolean
}

// One rollup node (GET /runs/{run}/rollup); recursive.
export interface RollupNode {
  name: string
  path: string
  calls: number
  wall_ms: number
  wall_ms_min: number
  wall_ms_max: number
  cpu_ms: number
  share: number
  wall_ms_per_step: number
  wall_ms_per_dof: number
  elem_by_type?: ElemRow[]
  children: RollupNode[]
}

// Per-step time history (GET /runs/{run}/series); null when no -perStep buffer.
export interface Series {
  step: number[]
  t: number[]
  dt: number[]
  iters: number[]
  mem_live_bytes?: number[]
  mem_peak_bytes?: number[]
  phases: string[]
  wall_ms: number[][] // [nSteps][nPhase]
}

export interface ComponentLive {
  classTag: number
  count: number
}

// End-of-run memory snapshot (GET /runs/{run}/memory).
export interface MemorySnapshot {
  matrix_live: number
  vector_live: number
  id_live: number
  peak_bytes: number
  components_live: ComponentLive[]
}

export interface DiffRow {
  path: string
  name: string
  base_ms: number
  cand_ms: number
  d_abs_ms: number
  d_pct: number | null
  d_share: number
  status: 'added' | 'removed' | 'changed'
}

export interface DiffResponse {
  base: string
  cand: string
  rows: DiffRow[]
}

// HTTP error carrying the backend's status + detail (404 missing run/series, 503 no file).
export class ApiError extends Error {
  status: number
  constructor(status: number, message: string) {
    super(message)
    this.status = status
    this.name = 'ApiError'
  }
}

async function get<T>(path: string): Promise<T> {
  let resp: Response
  try {
    resp = await fetch(`${API_BASE}${path}`)
  } catch {
    throw new ApiError(0, `cannot reach backend at ${API_BASE} (is profiler_api.py running?)`)
  }
  if (!resp.ok) {
    let detail = `${resp.status} ${resp.statusText}`
    try {
      const body = (await resp.json()) as { detail?: string }
      if (body.detail) detail = body.detail
    } catch {
      /* non-JSON body; keep the status text */
    }
    throw new ApiError(resp.status, detail)
  }
  return (await resp.json()) as T
}

export const api = {
  health: () => get<Health>('/health'),
  runs: () => get<ManifestRow[]>('/runs'),
  header: (run: string) => get<Header>(`/runs/${encodeURIComponent(run)}/header`),
  rollup: (run: string) => get<RollupNode>(`/runs/${encodeURIComponent(run)}/rollup`),
  series: (run: string) => get<Series>(`/runs/${encodeURIComponent(run)}/series`),
  memory: (run: string) => get<MemorySnapshot>(`/runs/${encodeURIComponent(run)}/memory`),
  diff: (base: string, cand: string) =>
    get<DiffResponse>(`/diff?base=${encodeURIComponent(base)}&cand=${encodeURIComponent(cand)}`),
}
