# ReMASTER Patch Notes (TIPB)
Date: 2026-02-04

This repo uses **ReMASTER** (as a BEAST2 package) for SC simulation. We found a mismatch between simulated typed trees and the hazard model assumed by **TISCTreePrior**.

## Summary
Root cause: `remaster.reactionboxes.ContinuousCoalescentReactionBox.getNextReactionTime()` uses an invalid update when sampling waiting times for **piecewise-constant** `Reaction` propensities. The previous implementation updates a uniform `u` by subtracting interval CDF mass, which can drive `u <= 0` and produce NaN reaction times. In practice, NaN times bias/suppress reactions (notably migration), yielding unrealistic typed trees (often at the tip-type lower bound).

Fix: replace the sampler with a correct inverse-CDF implementation:
- sample `w ~ Exp(1)` once (`w = -log(u)`)
- walk intervals and subtract hazards `w -= prop * Δt` until `w` is within the current interval.

The patch is recorded as:
- `patches/0001-fix-piecewise-reaction-waiting-time-sampler.patch`

## Evidence / Diagnostics
TIPB-side diagnostics (generated per replicate oracle run):
- `paper_validation/scripts/benchmarks/oracle_diagnostics.py`
- Key CSV: `diagnostics/remaster_waittime_sanity.csv`

Project write-up:
- `docs/plans/2026-02-04-remaster_treeprior_alignment.md`
- Beginner-friendly explainer (propensity / waiting-time sampling / when it triggers):
  - `docs/plans/2026-02-04-remaster_waiting_time_sampler_bug_explainer.md`

## How To Rebuild (vendored reference copy)
The reference source lives under `paper_validation/deps/remaster/` but is generally not committed as a full dependency snapshot.

Build a zip containing `lib/remaster.v2.7.4.TIPB.jar`:

```bash
cd paper_validation/deps/remaster
export JAVA_HOME=/opt/homebrew/opt/openjdk@17/libexec/openjdk.jdk/Contents/Home
export PATH="$JAVA_HOME/bin:$PATH"
ant clean build
```

The output zip is:
- `paper_validation/deps/remaster/dist/remaster.v2.7.4.TIPB.zip`

## Installing Into BEAST (local machine)
On macOS, BEAST packages are stored under:
- `~/Library/Application Support/BEAST/2.7/remaster/`

From the built zip, update:
- `lib/remaster.v2.7.4.TIPB.jar`

and (optional) the source jar:
- `remaster.v2.7.4.TIPB.src.jar`

Keep a backup of the original jars before replacing.
