# CLAUDE.md - Nova Project Navigation

Quick reference and navigation guide for Claude Code. For detailed information, see linked documents.

## Project Overview

`nova` is a de novo variant insertion simulator for benchmarking structural variant detection in ultra-low frequency scenarios using long-read sequencing data.

→ **Full details**: [claude_docs/project_overview.md](claude_docs/project_overview.md)

## Quick Start

```bash
# Activate environment
conda activate nova

# Run simulation
nova simulate input.bam config.json -o output/

# Validate configuration
nova validate-config config.json
```

## General Guidelines

- Think critically
- Avoid "yes-man" behavior: raise concerns and potential pitfalls
- Favor code cleanliness and conciseness: avoid over-engineering, consult with user when unsure

→ **Research findings**: [claude_docs/window_limit_experiments.md](claude_docs/window_limit_experiments.md)

## Documentation Map

### Core Documentation
- **Project Overview**: [claude_docs/project_overview.md](claude_docs/project_overview.md)
- **Architecture Guide**: [claude_docs/architecture.md](claude_docs/architecture.md)

### Development Guides
- **Testing Philosophy**: [claude_docs/testing_philosophy.md](claude_docs/testing_philosophy.md)
- **Testing Setup**: [tests/README.md](tests/README.md)
- **Analysis Scripts**: [scripts/README.md](scripts/README.md)

## Documentation Guidelines for Claude

### When to Create Documentation

| Content Type | Location | When to Create |
|-------------|----------|----------------|
| Quick reference | `CLAUDE.md` | Critical params/commands only |
| Experiments | `claude_docs/` | Significant findings with data |
| Architecture | `claude_docs/` | Major design changes |
| Test guides | `tests/README.md` | Testing requirements change |
| Script docs | `scripts/README.md` | New analysis tools added |

### Documentation Standards

1. **Document, don't comment**: Inline comments should be minimal, write implementation details to docstrings and .md files when needed
2. **New experiments**: Create `claude_docs/descriptive_name.md`
3. **Keep CLAUDE.md concise**: Under 200 lines, navigation focus
4. **Cross-reference**: Link between related documents
5. **Use templates**: Follow [claude_docs/documentation_style_guide.md](claude_docs/documentation_style_guide.md)

### Quick Rules
- ✅ Create docs for: experiments, major features, user requests
- ❌ Don't create docs for: minor fixes, temporary notes, code comments
- 📝 Always update CLAUDE.md when: adding new docs, changing key params

## Known Issues
**Environment**: If `which python` doesn't show the nova conda environment path, prompt the user to exit and run `conda activate nova` before proceeding.
**Snakemake output crashes Claude Code**: Running `snakemake` commands can cause BrokenPipeError crashes. Allow the user to test `snakemake` commands instead.