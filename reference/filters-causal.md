# Causal and stateful IIR filtering primitives

This file provides one-pass, causal, state-carrying IIR filtering that
mirrors `scipy.signal` (`lfilter`, `lfilter_zi`, `sosfilt`,
`sosfilt_zi`, `sosfiltfilt`). Unlike the zero-phase
[`butterworthFilter`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
default (which uses forward-backward `filtfilt` and is inherently
non-causal), these functions process a signal strictly forward in time
and can carry filter state across consecutive chunks. This is the
numerical core that a real-time streaming backend (e.g. a future
`PhysioStream`) reuses so that offline and online filtering produce
identical results.

## Details

All recursions use the Direct-Form-II-Transposed (DF2T) structure, which
is the same structure `scipy.signal` and
[`signal::filter`](https://rdrr.io/pkg/signal/man/filter.html) use, so
the outputs agree to machine precision.

## References

Oppenheim, A.V. & Schafer, R.W. (2010). *Discrete-Time Signal
Processing*, 3rd ed. Prentice Hall.

Virtanen, P. et al. (2020). SciPy 1.0. *Nature Methods* 17, 261-272.
(reference implementation: `scipy.signal.sosfilt` / `lfilter` /
`lfilter_zi` / `sosfilt_zi`).
