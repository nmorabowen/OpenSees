# WP-163/165 multi-rank gate -- model (OpenSeesMP, one Tcl interpreter per rank).
# Checker: mp_stage_check.py. Run under a real launcher, e.g. on Esmeralda:
#   srun --mpi=pmix_v3 openseesmp.sh mp_stage_model.tcl <out_dir>      (2 ranks)
#
# Each rank builds its OWN small model (no domain partitioning: this exercises the
# recorders' launcher-environment path, which is what openseesmp uses):
#   rank 0: truss 1-2 (nodes 1, 2)      rank 1: truss 11-12 (nodes 11, 12)
# and every rank defines region 1 = nodes {1, 2}, so rank 1 owns NONE of it.
#
# Gates (see the checker):
#   * per-rank part files from the launcher env (SLURM_PROCID inside an srun step)
#   * INFO RUN_ID identical on both ranks, RUN_ID_SCOPE = launcher       (WP-165 MP-9)
#   * a mid-run pattern with an imposed displacement keeps ONE MODEL_STAGE (WP-165 R6)
#   * the -R 1 recorder on rank 1 writes an EMPTY_PARTITION stage        (WP-165 MP-8)
#   * the Monitor sink is per rank: mon.part-<rank>.h5                    (WP-163 MP-6)

set out [lindex $argv 0]
set pid [getPID]
set np [getNP]

wipe
model basic -ndm 2 -ndf 2
uniaxialMaterial Elastic 1 1000.0
if {$pid == 0} {
    node 1 0.0 0.0; fix 1 1 1
    node 2 1.0 0.0; fix 2 0 1
    element truss 1 1 2 1.0 1
    set free 2
} else {
    node 11 0.0 5.0; fix 11 1 1
    node 12 1.0 5.0; fix 12 0 1
    element truss 11 11 12 1.0 1
    set free 12
}
region 1 -node 1 2

recorder ladruno $out/stage.ladruno -N displacement -G energy
recorder ladruno $out/region.ladruno -R 1 -N displacement
recorder Monitor -node $free -dof 1 -sink $out/mon.h5

timeSeries Linear 1
pattern Plain 1 1 { load $free 1.0 0.0 }
constraints Transformation
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-12 10
algorithm Newton
integrator LoadControl 0.25
analysis Static
analyze 2

# R6: a pattern holding an SP moves the domain-change stamp but not the topology
pattern Plain 2 1 { sp $free 1 0.001 }
analyze 2

wipe
puts "MP_STAGE_RANK $pid OF $np DONE"
