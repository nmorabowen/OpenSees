# WP-163 M4 gate -- a SEQUENTIAL OpenSees run inside an `sbatch --ntasks=N` job
# WITHOUT srun. The batch shell exports SLURM_NTASKS=N and SLURM_PROCID=0 but no
# real job step; the recorder must keep the declared filename (not part-0 of N).
# Checker: mp_stage_check.py (expects <out_dir>/m4_seq.ladruno, PARTITIONED = 0).
#   OpenSees m4_sequential_model.tcl <out_dir>      (inside the batch script)
set out [lindex $argv 0]
wipe
model basic -ndm 2 -ndf 2
uniaxialMaterial Elastic 1 1000.0
node 1 0.0 0.0; fix 1 1 1
node 2 1.0 0.0; fix 2 0 1
element truss 1 1 2 1.0 1
recorder ladruno $out/m4_seq.ladruno -N displacement
timeSeries Linear 1
pattern Plain 1 1 { load 2 1.0 0.0 }
system FullGeneral
integrator LoadControl 1.0
algorithm Linear
analysis Static
analyze 1
wipe
puts "M4_SEQ DONE"
