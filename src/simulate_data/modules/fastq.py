"""Simulate random FASTQ reads."""

import argparse
import random
import sys


def register_parser(parser):
    """Define the ``fastq`` subcommand's CLI arguments."""
    parser.add_argument(
        "-n",
        "--num-reads",
        type=int,
        default=10,
        help="Number of reads to generate (default: 10)",
    )
    parser.add_argument(
        "-l",
        "--read-length",
        type=int,
        default=150,
        help="Length of each read (default: 150)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output FASTQ file (default: stdout)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility",
    )
    return parser


def main(args):
    """Generate and write simulated FASTQ reads."""
    rng = random.Random(args.seed)
    bases = "ACGT"
    qual_chars = "".join(chr(33 + q) for q in range(41))

    out = args.output
    for i in range(args.num_reads):
        seq = "".join(rng.choices(bases, k=args.read_length))
        quals = "".join(rng.choices(qual_chars, k=args.read_length))
        out.write(f"@simulated_read_{i + 1}\n{seq}\n+\n{quals}\n")

    if out is not sys.stdout:
        out.close()
