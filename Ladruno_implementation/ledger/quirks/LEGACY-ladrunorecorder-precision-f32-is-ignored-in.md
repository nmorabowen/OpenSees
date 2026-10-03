---
wp: LEGACY
title: "LadrunoRecorder -precision f32 is ignored in -envelope mode — STORED_PRECISION now honest (FIXED)"
legacy_seq: 65
---
### LadrunoRecorder `-precision f32` is ignored in `-envelope` mode — STORED_PRECISION now honest (FIXED)
- **Bites:** `-precision f32` only changes the dtype of the streaming per-step DATA
  datasets (`StreamingSink::createTimeSeries3d`, `Ladruno_Sinks.cpp` — `H5T_IEEE_F32LE`).
  In `-envelope` mode there are no streaming DATA datasets; the only result datasets are
  the EnvelopeSink MIN/MAX/ABSMAX, which are **always f64**. But `initialize()` stamped
  `INFO/STORED_PRECISION` purely from the `store_data_f32` flag → an `-envelope -precision f32`
  file claimed `f32` while every dataset in it was f64. A reader trusting the attribute to
  pick its diff tolerance would be misled.
- **Fix:** `STORED_PRECISION` is now `f32` only when `store_data_f32 && !envelope_mode`
  (it must describe what is actually on disk); a one-time warning is emitted if `-precision f32`
  is combined with `-envelope`. (Honoring f32 *inside* the envelope datasets is a separate,
  judgment-dependent enhancement — not done; the label-honesty fix is unambiguous.) 2026-06-03.
