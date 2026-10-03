# PhysioPreprocess: preprocessing for PhysioExperiment objects

PhysioPreprocess is the signal-conditioning layer of the ecosystem. It
operates directly on
[PhysioExperiment::PhysioExperiment](https://x-biosignal.r-universe.dev/PhysioExperiment/reference/PhysioExperiment.html)
containers, reading an input assay and writing the processed signal into
a new named assay so that every step is traceable and the raw data is
never overwritten.

## Filtering

[`butterworthFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/butterworthFilter.md)
(IIR, SOS or ba form),
[`firFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/firFilter.md)
(linear-phase FIR),
[`notchFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/notchFilter.md)
(fixed power-line notch),
[`filterSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/filterSignals.md)
(moving average), and
[`detrendSignal()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignal.md)
/
[`detrendSignals()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detrendSignals.md)
/
[`removeBaseline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeBaseline.md)
for trend and offset removal.

## Causal and stateful filtering

[`sosDesign()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosDesign.md),
[`sosfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfilt.md),
[`sosfiltInit()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltInit.md),
[`sosfiltfilt()`](https://x-biosignal.github.io/PhysioPreprocess/reference/sosfiltfilt.md),
[`lfilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilter.md),
[`lfilterInit()`](https://x-biosignal.github.io/PhysioPreprocess/reference/lfilterInit.md),
and
[`StreamFilter()`](https://x-biosignal.github.io/PhysioPreprocess/reference/StreamFilter.md)
for block-wise, state-carrying filtering (the primitives used by the
streaming layer).

## Resampling

[`resample()`](https://x-biosignal.github.io/PhysioPreprocess/reference/resample.md)
(arbitrary rate),
[`decimate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/decimate.md)
(integer downsampling with anti-aliasing),
[`interpolate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolate.md)
(integer upsampling), with per-assay rate tracking via
[`assaySamplingRates()`](https://x-biosignal.github.io/PhysioPreprocess/reference/assaySamplingRates.md)
and
[`setAssaySamplingRate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/setAssaySamplingRate.md).

## Artifact handling

[`detectBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectBadChannels.md)
and
[`interpolateBadChannels()`](https://x-biosignal.github.io/PhysioPreprocess/reference/interpolateBadChannels.md),
[`detectArtifacts()`](https://x-biosignal.github.io/PhysioPreprocess/reference/detectArtifacts.md)
on continuous data,
[`rejectBadEpochs()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rejectBadEpochs.md)
and
[`baselineCorrect()`](https://x-biosignal.github.io/PhysioPreprocess/reference/baselineCorrect.md)
on epoched data, and the configurable
[`cleanData()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanData.md)
wrapper.

## ICA and ASR

[`runICA()`](https://x-biosignal.github.io/PhysioPreprocess/reference/runICA.md)
/
[`removeICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/removeICAComponents.md)
and the earlier
[`icaDecompose()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaDecompose.md)
/
[`icaRemove()`](https://x-biosignal.github.io/PhysioPreprocess/reference/icaRemove.md)
/
[`classifyICAComponents()`](https://x-biosignal.github.io/PhysioPreprocess/reference/classifyICAComponents.md)
interface for independent-component artifact removal;
[`asrCalibrate()`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrCalibrate.md),
[`asrProcess()`](https://x-biosignal.github.io/PhysioPreprocess/reference/asrProcess.md)
and
[`cleanRawdata()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanRawdata.md)
for Artifact Subspace Reconstruction.

## Line-noise removal

[`zapLine()`](https://x-biosignal.github.io/PhysioPreprocess/reference/zapLine.md)
(ZapLine / DSS) and
[`cleanLine()`](https://x-biosignal.github.io/PhysioPreprocess/reference/cleanLine.md)
(sinusoidal regression) for deep power-line attenuation without a wide
spectral notch.

## Re-referencing and pipelines

[`rereference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/rereference.md)
with
[`getCurrentReference()`](https://x-biosignal.github.io/PhysioPreprocess/reference/getCurrentReference.md)
/
[`isAverageReferenced()`](https://x-biosignal.github.io/PhysioPreprocess/reference/isAverageReferenced.md),
and
[`createPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/createPipeline.md)
/
[`applyPipeline()`](https://x-biosignal.github.io/PhysioPreprocess/reference/applyPipeline.md)
to compose reproducible, replayable preprocessing chains.

## Where to go next

The data container itself lives in PhysioExperiment. See the
`preprocessing-guide` and `ica-artifact-removal` vignettes for worked,
end-to-end examples.

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioPreprocess>

- <https://x-biosignal.r-universe.dev/PhysioPreprocess>

- Report bugs at
  <https://github.com/x-biosignal/PhysioPreprocess/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>
