# AnaCoDa Versioning

AnaCoDa uses a four-component version: `MAJOR.MINOR.PATCH.BUILD`.

| Component | When to bump |
|-----------|--------------|
| MAJOR     | Breaking changes to the public R or `.rst` file API. |
| MINOR     | New user-visible features (new priors, new models, new R methods, new restart-file sections). |
| PATCH     | Bug fixes that change behavior but not API.  Examples: the legacy-format restart-MCMC fix; the numElongationMixtures persistence fix. |
| BUILD     | Internal cleanups, refactors, build-system changes, doc-only changes.  Bump even when no code behavior changes if you want a fresh "this build is different" signal. |

## Practical cadence

Bump on every PR that materially changes the package.  Update both:

- `DESCRIPTION`: `Version:` line.
- `DESCRIPTION`: `PreviousGitCommitHash:` to the commit SHA of the previous bump.
- `DESCRIPTION`: `Date:` to today (UTC).

The git commit SHA is also captured automatically at `./configure` time
and embedded in every restart file via the `>buildInfo:` block, so it is
unambiguous in retrospect even when `Version:` itself was not bumped.

## Restart-file format generations

Because `.rst` files outlive the AnaCoDa builds that wrote them, the
reader can also identify the *format generation* of a `.rst` from its
structure (presence/absence of certain sections, arity of category
rows).  This is independent of the `Version:` string and is the
authoritative compatibility key.  See `Parameter::detectRestartFileGeneration()`.

| Generation              | Recognizing signals                                         |
|-------------------------|-------------------------------------------------------------|
| MODERN_WITH_BUILDINFO   | `>buildInfo:` block at top.                                 |
| MODERN_NO_BUILDINFO     | No `>buildInfo:`, but `>numSynthesisRateCategories:` present and 4-field `>categories:` rows. |
| LEGACY_2022_RELEASE     | No `>buildInfo:`, no `>numSynthesisRateCategories:`, 2-field `>categories:` rows.  Matches the official 2022 AnaCoDa release (`0.1.4.0`). |
| UNKNOWN                 | Anything else; reader emits a warning and uses best-effort fallbacks. |

When `Version:` and the file-structure generation disagree, the
structure wins.  Tests cover all four generations.

## Automation (future)

A `tools/bump_version.sh` script could automate the DESCRIPTION updates
and the tag.  Not implemented yet; current cadence is manual.
