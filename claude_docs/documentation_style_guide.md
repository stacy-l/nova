# Documentation Style Guide

*Created: 2025-06-30*

## Purpose

This guide ensures consistent, discoverable, and maintainable documentation across the nova project. All Claude instances and contributors should follow these guidelines.

## Document Organization

### File Placement

| Content Type | Location | Example |
|-------------|----------|---------|
| Navigation & quick reference | `CLAUDE.md` | Command examples, key parameters |
| Experimental results | `claude_docs/` | `window_limit_experiments.md` |
| Architecture & design | `claude_docs/` | `architecture.md` |
| Testing information | `tests/README.md` | Test setup, requirements |
| Script documentation | `scripts/README.md` | Usage, parameters |
| API documentation | Source code docstrings | Function/class docs |

### When to Create New Documentation

#### DO Create Docs When:
✅ Completing a significant experiment with results
✅ Adding a new major feature or system
✅ User explicitly requests documentation
✅ Discovering non-obvious behavior that needs explanation
✅ Creating reusable analysis workflows

#### DON'T Create Docs When:
❌ Making minor bug fixes
❌ Doing routine maintenance
❌ Information belongs in code comments
❌ Content would duplicate existing docs

## Document Structure

### Standard Sections

1. **Title** - Clear, descriptive
2. **Executive Summary** - For docs >500 words
3. **Background/Context** - Why this exists
4. **Main Content** - Organized with clear headers
5. **Examples** - Practical usage
6. **Related Documentation** - Cross-references

### Markdown Template

```markdown
# Document Title

*Created: YYYY-MM-DD*

## Executive Summary

[2-3 sentence overview for long documents]

## Background

[Context and motivation]

## Main Content

### Section 1
[Content with examples]

### Section 2
[Content with examples]

## Related Documentation
- [Link 1](relative/path.md)
- [Link 2](relative/path.md)
```

## Writing Style

### General Guidelines

1. **Be Concise**: Get to the point quickly
2. **Use Examples**: Show, don't just tell
3. **Active Voice**: "Nova simulates..." not "Simulations are performed..."
4. **Present Tense**: For current behavior
5. **Technical Accuracy**: Verify claims with code

### Formatting Standards

#### Code Blocks
```bash
# Always include language identifier
nova simulate input.bam config.json -o output/
```

#### Command Examples
- Include `$` prefix for shell commands only when showing output
- Use realistic filenames and paths
- Show both simple and complex usage

#### Tables
Use tables for structured comparisons:
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max-reads-per-window` | 1 | Clustering limit |

#### Emphasis
- **Bold** for important concepts or warnings
- *Italics* for introducing new terms
- `Code font` for parameters, filenames, commands

## File Naming

### Conventions
- Use lowercase with underscores: `window_limit_experiments.md`
- Be specific: `memory_optimization_report.md` not `optimization.md`
- Include dates in experimental reports: `*_2024_01.md` if multiple versions

### Bad Examples
- `notes.md` - Too vague
- `TODO.md` - Use issue tracker instead
- `misc.md` - Always be specific

## Cross-References

### Internal Links
```markdown
See [Architecture Guide](architecture.md) for details.
```

### Section Links
```markdown
See [Window Limit Parameter](#window-limit-parameter) above.
```

### External Links
```markdown
Uses [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) for variant calling.
```

## CLAUDE.md Maintenance

### Keep It Focused
- Under 200 lines
- Navigation and critical info only
- Link to details rather than include them

### Update Triggers
- New major feature added
- Documentation structure changes
- Critical parameter changes
- New experimental results

### Don't Add
- Detailed explanations (link instead)
- Historical information
- Temporary notes

## Experimental Reports

### Required Elements
1. **Problem Statement**: What question are we answering?
2. **Methodology**: How did we test?
3. **Results**: What did we find?
4. **Conclusions**: What does it mean?
5. **Recommendations**: What should users do?

### Data Presentation
- Use visualizations when helpful
- Include summary statistics
- Provide raw data location
- Explain significance

## Version Control

### Commit Messages
```
docs: Add window limit experiments report

- Document clustering issue and resolution
- Include performance metrics
- Add recommendations for users
```

### Update Strategy
- Update existing docs for minor changes
- Create new docs for major revisions
- Archive old experimental reports if needed

## Quality Checklist

Before committing documentation:

- [ ] Follows file placement guidelines
- [ ] Uses standard structure
- [ ] Includes examples where appropriate
- [ ] Cross-references related docs
- [ ] Spell-checked and proofread
- [ ] Technical accuracy verified
- [ ] CLAUDE.md updated if needed

## Examples of Good Documentation

### Good: Specific and Actionable
```markdown
## Window Limit Parameter

Use `--max-reads-per-window=1` to prevent genomic clustering:

```bash
nova simulate input.bam config.json --max-reads-per-window 1 -o output/
```

This prevents multiple nova reads from being selected from the same 1MB window,
eliminating false positives from variant caller clustering.
```

### Bad: Vague and Theoretical
```markdown
## Parameters

Various parameters can be adjusted to change the behavior of the simulation.
Different values may produce different results.
```

## Maintenance

### Regular Reviews
- Quarterly documentation audit
- Remove outdated content
- Update cross-references
- Consolidate redundant docs

### User Feedback
- Track common questions
- Update docs to address gaps
- Clarify confusing sections

## Conclusion

Good documentation is:
- **Discoverable**: Easy to find when needed
- **Actionable**: Helps users accomplish tasks
- **Maintainable**: Easy to keep current
- **Consistent**: Follows predictable patterns

When in doubt, optimize for clarity and usefulness over completeness.