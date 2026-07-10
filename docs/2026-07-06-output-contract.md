# Stable output contract for simulated genomes

Date: 2026-07-06 15:56 America/Los_Angeles

## Purpose

Make TE/SV genome-modification outputs predictable for downstream short-read
and long-read simulation workflows.

## Status

- `te-insertion` writes `<output>/final_genome.fa` and `<output>/truth_te.vcf`.
- `sv-placement` writes `<output>/final_genome.fa` and `<output>/truth_sv.vcf`.
- Tool-native intermediate outputs are preserved in the output directory.
- README examples and SLURM scripts now use `final_genome.fa` for downstream
  references.

## Commands

Planned verification from repository root:

```bash
pixi run format
pixi run lint
pixi run test
```

## Logic

Downstream read simulators and combined TE/SV workflows should not need to know
tool-specific filenames such as `Sim_Chr1.fa` or
`sv_output_modified_genome.fa`.

## Next Steps

- Verify against mocked unit tests.
- Smoke-test with real external tools when TEvarSim, SURVIVOR, ART, and PBSIM3
  are available on the target HPC environment.
