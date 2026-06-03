"""relax_to_static — drive a model to static equilibrium with LadrunoDynamicRelaxation.

The stock OpenSees transient driver has no convergence early-exit, so static-rest
detection is script-owned (ADR-21 decision 4). This helper runs the transient
analysis in chunks and watches the **velocity and displacement decay** at a set of
monitored nodes — at static rest both go to zero. It is the documented entry point
for the integrator.

Required model recipe (mirrors CentralDifferenceLadruno):
    ops.constraints('Transformation')   # or Plain
    ops.numberer('Plain')
    ops.system('Diagonal')              # makes the M*-solve a diagonal apply
    ops.test('NormUnbalance', 1e-6, 1)
    ops.algorithm('Linear')             # exactly one solve / step
    ops.integrator('LadrunoDynamicRelaxation')   # Gershgorin M*, kinetic damping
    # then a CONSTANT load pattern, then:
    relax_to_static(ops, nodes=[2], ndf=2)

Usage:
    from relax_to_static import relax_to_static
    ok, steps, why = relax_to_static(ops, nodes=[2], ndf=2)
"""


def relax_to_static(ops, nodes, ndf, chunk=50, dt=1.0, max_iter=200000,
                    tol_vel=1.0e-6, tol_du=1.0e-8, verbose=False):
    """Run LadrunoDynamicRelaxation to static rest.

    Parameters
    ----------
    ops      : the opensees / openseespy module (already holding a built model
               with `integrator('LadrunoDynamicRelaxation', ...)` set).
    nodes    : iterable of node tags to monitor for rest.
    ndf      : number of dofs per monitored node.
    chunk    : pseudo-time steps per analyze() call.
    dt       : pseudo-time step (match the integrator's -dt; default 1.0).
    tol_vel  : max |nodal velocity| at rest.
    tol_du   : max |displacement change over a chunk| at rest.

    Returns (converged: bool, steps_run: int, reason: str).

    Rest needs BOTH gates: kinetic damping zeros velocities at each KE peak, so a
    momentarily-small velocity mid-convergence is filtered by the displacement-
    change gate (which is still large until the structure actually settles).
    """
    ops.analysis('Transient')

    def _disp_snapshot():
        return {n: [ops.nodeDisp(n, d + 1) for d in range(ndf)] for n in nodes}

    prev = _disp_snapshot()
    done = 0
    while done < max_iter:
        if ops.analyze(chunk, dt) != 0:
            return False, done, 'analyze failed (divergence / dt too large)'
        done += chunk

        vmax = 0.0
        dumax = 0.0
        for n in nodes:
            for d in range(ndf):
                vmax = max(vmax, abs(ops.nodeVel(n, d + 1)))
                u = ops.nodeDisp(n, d + 1)
                dumax = max(dumax, abs(u - prev[n][d]))
                prev[n][d] = u

        if verbose:
            print(f"  relax step {done:>7d}: max|v|={vmax:.3e}  max|du|={dumax:.3e}")

        if vmax < tol_vel and dumax < tol_du:
            return True, done, 'converged'

    return False, done, 'reached max_iter without rest'
