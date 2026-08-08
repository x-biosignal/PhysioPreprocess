## Test environments

* Local: Ubuntu 24.04, R 4.5.2 (x86_64-pc-linux-gnu)
* win-builder: R-devel and R-release (planned)
* macbuilder: R-release on macOS (planned)
* R-hub v2 (GitHub Actions, `r-hub/actions/run-check`): linux, windows, macos (planned)

## R CMD check results

`R CMD check --as-cran` gives:

0 ERRORs | 0 WARNINGs | 1 NOTE

on the CRAN incoming pipeline, where `PhysioCore` (this package's only
non-mainstream strong dependency) will already be available. The local
`--as-cran` run reports two WARNINGs instead; both are environmental /
expected and are explained below. No user-facing package problems remain.

### CRAN incoming feasibility (New submission)

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Yusuke Matsui <mail.to.matsui@gmail.com>'

New submission

Strong dependencies not in mainstream repositories:
  PhysioCore
```

This is the first submission of PhysioPreprocess to CRAN, hence "New
submission". PhysioPreprocess depends on `PhysioCore`, a sibling package from
the same PhysioExperiment ecosystem that is submitted to CRAN immediately
before this package (see "Reverse dependencies" below). Because `PhysioCore`
is not yet on a mainstream repository, the local check elevates this item to a
WARNING; once `PhysioCore` is accepted first it is a mainstream strong
dependency and this item reduces to the standard new-submission NOTE. Both
outcomes are expected and justifiable.

The incoming-feasibility check may additionally report "possibly misspelled
words" for domain terms such as Butterworth, resampling, and ICA. These are
correct signal-processing terms; the acronym ICA is spelled out as
"independent component analysis" in the Description. It may also report "unable
to verify current time", which is an infrastructure artifact unrelated to the
package.

The URLs in DESCRIPTION (the package homepage on GitHub and its r-universe
page under the `x-biosignal` organisation) may be reported as "possibly
invalid" (HTTP 404). These public resources become live at the coordinated
public release, which happens together with CRAN acceptance; the URLs are
correct and point at the package's canonical locations.

### qpdf WARNING (local only)

```
* checking CRAN incoming feasibility ... WARNING
'qpdf' is needed for checks on size reduction of PDFs
```

This reflects `qpdf` being absent from the local check host only. The package
ships no PDFs, and this WARNING does not occur on CRAN, win-builder,
macbuilder, or R-hub, where `qpdf` is present.

## Reverse dependencies

PhysioPreprocess is part of the PhysioExperiment ecosystem. Its dependency
`PhysioCore` is submitted first; PhysioPreprocess and its siblings then follow.
The intended CRAN submission order is:

  PhysioCore -> PhysioIO / PhysioPreprocess -> PhysioAnalysis

There are no reverse dependencies currently on CRAN.

## Additional notes

* The maintainer is the sole author.
* Examples that require the optional 'fastICA' package (in Suggests) are guarded
  with `requireNamespace()` and wrapped in `\donttest{}` so they do not fail
  when the suggested package is unavailable.
