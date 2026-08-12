# Re-referencing Operations for EEG Data

Functions for changing the reference electrode in EEG recordings.
Re-referencing is a common preprocessing step that affects the spatial
distribution of the signal.

## Usage

``` r
rereference(
  x,
  ref_type = c("average", "robust", "median", "channel", "channels", "REST"),
  ref_channels = NULL,
  exclude = NULL,
  input_assay = NULL,
  output_assay = "rereferenced",
  keep_ref = TRUE,
  robust_noise_sd = 4,
  robust_max_iter = 5
)
```

## Arguments

- x:

  A PhysioExperiment object.

- ref_type:

  Type of re-referencing: "average" (common average reference), "robust"
  (PREP robust average reference: iteratively detect and exclude
  high-variance bad channels from the average; Bigdely-Shamlo et al.,
  2015), "median" (channel-wise median reference, robust to outlier
  channels), "channel" (single channel), "channels" (average of
  specified channels), or "REST" (Reference Electrode Standardization
  Technique).

- ref_channels:

  For "channel" or "channels" type, the channel name(s) or index/indices
  to use as reference.

- exclude:

  Channels to exclude from average reference calculation (e.g., non-EEG
  channels like EOG, EMG).

- input_assay:

  Input assay name. If NULL, uses default assay.

- output_assay:

  Output assay name. Default is "rereferenced".

- keep_ref:

  Logical. If TRUE, keeps the original reference channel(s) in the
  output (zeroed). If FALSE, removes them.

- robust_noise_sd:

  For "robust", the threshold (robust SDs above the median channel
  variance) for flagging a bad channel (default: 4).

- robust_max_iter:

  For "robust", the maximum number of PREP iterations (default: 5).

## Value

A `PhysioExperiment` object with re-referenced data stored in a new
assay named `output_assay`. The `reference` and `previous_reference`
metadata fields are updated. When `keep_ref = FALSE` and using
channel-based reference, the reference channel(s) are removed and a new
object with reduced channel count is returned.

## Details

Re-referencing transforms the data by subtracting a reference signal
from each channel. The choice of reference affects the spatial
distribution and interpretation of the signal.

**Average reference** ("average"): Subtracts the mean of all channels at
each time point. This is commonly used for high-density EEG and provides
a reference-independent measure, but requires good spatial sampling.

**Single channel reference** ("channel"): Subtracts the signal from a
specified electrode. Common choices include Cz, linked mastoids
(A1+A2)/2, or nose reference.

**Multi-channel reference** ("channels"): Subtracts the average of
multiple specified channels. Useful for linked mastoids or other custom
references.

## References

Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the Brain." 2nd
ed. Oxford University Press. Re-reference EEG data

Changes the reference electrode for EEG recordings. Supports common
re-referencing schemes including average reference, linked mastoids, and
single electrode reference.

Nunez, P.L. & Srinivasan, R. (2006). "Electric Fields of the Brain." 2nd
ed. Oxford University Press.

## See also

[`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
for querying the current reference,
[`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md)
for checking average reference status,
[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
for frequency-domain preprocessing.

## Examples

``` r
# Create example EEG data
set.seed(123)
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(1000), nrow = 100, ncol = 10)),
  colData = S4Vectors::DataFrame(
    label = c("Fp1", "Fp2", "F3", "F4", "C3", "C4", "P3", "P4", "O1", "O2"),
    type = rep("EEG", 10)
  ),
  samplingRate = 256
)

# Apply average reference
pe_avg <- rereference(pe, ref_type = "average")

# Re-reference to a single channel (Cz)
pe_cz <- rereference(pe, ref_type = "channel", ref_channels = "C3")

# Re-reference to linked mastoids (if available)
# pe_linked <- rereference(pe, ref_type = "channels",
#                          ref_channels = c("M1", "M2"))
```
