# Variant Config Setup

`nova` uses a JSON configuration file to specify the types and numbers of variants to simulate.

There are three types of variant generators:
- `random`: Generate random sequences of a given length and GC content
  - `n`: Number of sequences to generate
  - `length`: Length of the sequences
  - `gc_content`: GC content of the sequences
- `simple`: Generate sequences with a simple repeat pattern
  - `n`: Number of sequences to generate
  - `repeat`: The repeat pattern
  - `units`: The number of units in the repeat
- `predefined`: Use predefined sequences from a FASTA file
  - `fasta`: Path to the FASTA file containing the sequences
  - `spec`: A dictionary specifying the number of each sequence to generate

For all variant generators, you can specify a `mutations` dictionary to apply mutations to the sequences.
- `mutations`:
      - `substitution_rate`: The rate of single base substitutions to apply to the sequences

## Example
```json
{
  "random": [
    {
      "n": 50,
      "length": 100,
      "gc_content": 0.4
    },
    {
      "n": 50,
      "length": 300
    }
  ],
  "simple": [
    {
      "n": 100,
      "repeat": "CAG",
      "units": 20
    }
  ],
  "predefined": [
    {
      "Alu": {
        "fasta": "predefined/dfam_AluY_homininae.fasta",
        "spec": {
          "AluYa5": 30,
          "AluYb8": 20
        }
      }
    },
    {
      "L1": {
        "fasta": "predefined/dfam_L1HS_homininae.fasta",
        "spec": {
          "L1HS_full": 2
        },
        "mutations": {
          "substitution_rate": 0.03
        }
      }
    }
  ]
}
```

This example generates:
- 50 random sequences of 100bp with 40% GC content
- 50 random sequences of 300bp with default GC content
- 100 CAG repeat sequences (20 units each)
- 30 AluYa5 + 20 AluYb8 sequences
- 2 L1HS full-length sequences with 3% divergence (single base substitutions) from the reference sequence

## Validation

You can validate your configs using the `validate-config` command:

```bash
nova validate-config examples/variant_config.json # valid config
nova validate-config tests/test_data/test_config_broken.json # broken example
```

`nova simulate` also validates configs before running, so this feature is mostly for convenience.