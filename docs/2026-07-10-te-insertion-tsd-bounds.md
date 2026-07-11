# TE insertion TSD bounds

Date: 2026-07-10 America/Los_Angeles

## Purpose

Expose TEvarSim target-site duplication length controls through the
`simulate-data te-insertion` wrapper so mPing simulations can use biologically
appropriate short TSDs by default.

## Status

- Added `--tsd-min` and `--tsd-max` to `te-insertion`.
- Defaults are now `--tsd-min 3` and `--tsd-max 5`.
- Both values are passed to `tevarsim Simulate`.
- Bounds are validated before launching TEvarSim.
- README and unit tests were updated.

## Commands

Planned verification:

```bash
pixi run pytest tests/test_te_insertion.py -v
pixi run simulate-data te-insertion --help
```

## Logic

TEvarSim defaults to a wider TSD range. mPing insertions are expected to have
short target-site duplications, so the wrapper now defaults to the mPing-focused
range while still allowing callers to override it for other TE families.

## Next Steps

Use explicit values in production SLURM scripts when a run should record its
TSD assumptions in the command line:

```bash
--tsd-min 3 \
--tsd-max 5
```
