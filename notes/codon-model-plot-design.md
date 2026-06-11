# Codon-frequency vs log(phi) plotting: generalization design

Status: design agreed 2026-06-07. Stage 1 DONE (grouping builder + panels= arg,
commits bb17bb3/9d97f5c); Stage 2 DONE for per-AA panel (codonModelData, 3dbea5f).
Remaining: extend codonModelData to split groupings; then ggplot engine.
Branch: feat/compact-trace-and-codon-plots (RMF).

## Problem

The codon-frequency-vs-log(phi) model plots are hard-coded. Three near-parallel
code paths implement the *same* pipeline (per-gene codon proportions -> bin by phi
-> summarize -> draw + overlay model curve):

- `plotSinglePanel` -- one AA, all its codons, denominator = AA total (standard plot)
- `plotWobbleSplitPanel` -- configurable codon sub-groups, denominator = group total
- legacy `"original"` AnaCoDa path

Layouts are a `switch` in `.getSplitPanels()`. Binning is fixed (quantiles of 5%),
the summary glyph is fixed (point + sd error bar), and the model overlay is a
posterior-mean line with no predictive band. Each new plot = new switch case +
bespoke draw code.

## Direction (agreed)

End state is **ggplot2**; we are reinventing what its grammar already provides.
But we migrate in stages so there is always a stable reference to diff against:

1. **Generalize the base-R code** (this stays the "standard to compare to").
2. **Extract a tidy data-content layer** (the ggplot-ready interface).
3. (later) **ggplot2 + patchwork renderer**, validated against the stage-1/2 base-R output.
4. Retire base-R once ggplot reaches visual parity.

Dependency note: RMF already Imports dplyr, tidyr, magrittr, rlang -- the tidy
data layer adds NO new deps. Only the eventual ggplot renderer adds ggplot2
(+ patchwork for multi-panel assembly).

## Architecture: named data content + layout

### Data content (stage 2): one tidy long data frame

A single function `(genome, phi, parameter) -> data.frame` is the source of truth.
Columns (sketch):

    aa  codon  codon.group  role    phi    proportion  lower  upper  n
    S   TCA    4            data    -0.31  0.22        0.18   0.26   140   # binned point+range
    S   TCA    4            model   <grid> 0.21        0.19   0.24   NA    # curve + PPI
    S   TCA    4            gene     0.7   0.25        NA     NA      1    # per-gene (violins)

`aes(x = phi, y = proportion, color = codon)` maps directly. The three
extensibility wishes then fall out of the grammar, not custom code:

- violins   -> geom_violin over role == "gene" rows
- bin size  -> argument to the data-prep (quantile n or explicit breaks)
- 90% PPI   -> geom_ribbon(lower, upper) over role == "model" (needs trace draws)

### Layout: grouping description + AA lists

The plot is described as `(grouping, aa-list)` panel groups -- the user's framing.
The `type` field already inside `.getSplitPanels` IS the grouping description;
it is just buried in a switch. Named vocabulary (maps to existing `type`):

| grouping name        | existing type | draws                                            |
|----------------------|---------------|--------------------------------------------------|
| individual           | (single)      | each codon of the AA, denom = AA total            |
| wobble-purine        | ag            | A/G-wobble codons within each block               |
| wobble-pyrimidine    | tc            | C/T-wobble codons within each block               |
| ry-aggregate         | ry / ry-all6  | purine-block vs pyrimidine-block proportion       |
| block-4v2            | 4v2           | 4-codon block vs 2-codon block                    |
| r-block              | r-block       | R-wobble block split                              |

Stage 1: replace the `switch` with a builder that consumes
`group(grouping, aa, label)` rows. Pre-designed plots become named
(data-spec, layout-spec) pairs; users compose with the SAME builder.

    plotROCOptions(panels = list(
        group("wobble-purine", c("E","K","Q"),             "2 Codon AAs"),
        group("wobble-purine", c("A","G","P","S","T","V"), "4 Codon AAs"),
        group("block-4v2",     c("L","R"),                 "6 Codon AAs")),
      subtitle = "Purine (A/G) Wobble Codons")
    # split-ag is exactly this list under a preset name
    plotROCOptions(layout = "split-ag")
    # standard one-AA plot:
    plotROCOptions(panels = list(group("individual", AA.STANDARD)))

Wrinkle: split-6codon shows 3 groupings per single AA (4v2 / r-block / Y-wobble
for Leu, Arg). It does not fit "one grouping per row"; keep it a named preset
(specialist page) rather than forcing it into the builder grammar.

## Honest trade-offs for the eventual ggplot port

Most maps cleanly (point/line/ribbon/violin geoms, facet_wrap by AA). Harder,
because they were hand-built in base R:

- within-pair rescaling (ag/ct): easy -- a data transform at the data layer.
- codon-group row banding + right-margin "N Codon AAs" labels: needs patchwork
  assembly + strip theming (or ggh4x nested strips). Real work.
- in-panel gene-count histogram strip: secondary geom_col on a rescaled axis,
  or an aligned patchwork sub-strip. Fiddly.

Decide port-vs-simplify per feature at ggplot time.

## Regression net

S. cerevisiae snapshot baseline (vdiffr) from a reused A-RMF S288c fit
(downsampled genes + thinned trace; trace draws retained for PPI). Captured
BEFORE the stage-1 refactor; stages 1 and 2 must keep snapshots byte-identical.
The vignette (`vignettes/codon-model-plots.Rmd`) embeds pre-rendered PNGs
(house style: eval=FALSE fit chunk + committed PNGs).

## Current option API (this branch)

plotROCOptions(layout, show.gene.hist, show.date, color.codon.groups, aa.include)
- color.codon.groups (renamed from show.group.borders, default FALSE): codon-group
  color coding of frames/axes/row-labels; off by default since the split already
  separates AAs by codon count.
