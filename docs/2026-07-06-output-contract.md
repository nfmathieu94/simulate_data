# Stable output contract for simulated genomes

Date: 2026-07-06 15:56 America/Los_Angeles

Last verified: 2026-09-02 21:21 PDT (America/Los_Angeles)

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

Current verification from repository root:

```bash
pixi run format
pixi run lint
pixi run test
```

## Logic

Downstream read simulators and combined TE/SV workflows should not need to know
tool-specific filenames such as `Sim_Chr1.fa` or
`sv_output_modified_genome.fa`.

## Current outcome

The stable output names are covered by the test suite. TEvarSim, SURVIVOR,
ART, PBSIM3, and CCS are installed in the Pixi environment, and the linked
operational workspace contains outputs generated with these tools. No output
contract work remains open here.
