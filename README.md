# Do flashes on simultaneous Heidelberg plates avoid Earth's shadow?

A pre-registered test for Project Lost Light. The Bruce astrograph at Heidelberg-Königstuhl photographed the same
field at the same moment through two lenses onto two plates. A brief flash of light lands on both plates, whereas an
emulsion defect affects only one, so a source on both plates of a twin pair is not a plate defect. If such sources
include glints of sunlight from objects in high orbit, they should be scarcer where Earth's shadow falls at
geostationary radius. This repository measures the twin pairs of 1900 to 1957, selects the sources on both plates,
and tests for that deficit.

The pre-registration is `draft/preregistration_v2.md`, rendered from its template (`draft/preregistration_v2.src.md`)
with every number filled from `work/results.json`. Every setting of the analysis, with its reason, is in
`config/analysis.json`.

## The freeze and this repository's history

`FREEZE.json` records the SHA-256 of the code, the settings, the pair lists, and the inputs manifest, with one checksum
over that list, which is posted with the pre-registration. `python3 tools/freeze.py check FREEZE.json` confirms that a
copy matches it. This repository starts from the frozen files as one commit. Its working history was left out because
it held an accidentally committed crash dump, and the commit named in `FREEZE.json` belongs to that history.

## What is here

- `heidelberg/lib`, `heidelberg/stages`: the measurement (detection, plate solutions, photometry, artificial flashes),
  the selection, the shadow geometry, and the likelihood and within-pair tests.
- `heidelberg/producers`: the calibration (toys on the record geometry: the null, power, the claim rule's size, and
  the upper limit's coverage).
- `config/analysis.json`: every setting that can change a result, each with its reason; `config/compute.json`: the
  settings that cannot (processes, memory, timeouts).
- `heidelberg/record_pairs.csv`, `heidelberg/record_control_pairs.csv`: the record set and the different-night
  controls; `heidelberg/measure_pairs_tuning.csv`: the tuning pairs on which every rule was set.
- `heidelberg/inputs.sha256`: the SHA-256 of every input file (plate tables, catalogue snapshots, ephemerides, saved
  catalogue answers).
- `Snakefile`: the whole analysis as one workflow; `tools/`: comparison, sealing, and rendering scripts; `aws/`: the
  launcher for the jobs that read plate scans.

## How to run it

Everything runs in one container image (`Dockerfile`; the package versions are pinned in `requirements.lock`):

    docker build --platform linux/amd64 -t lost-light:dev .
    ./run.sh snakemake draft -j4                 # the documents, from the results file
    ./run.sh snakemake check_sample -j4          # the record set and its geometry, rebuilt from the plate tables

The steps that read plate scans (measuring, the record run's two parts) and the calibration run on Amazon Web
Services through `aws/launch.sh`, which starts one machine for a Snakemake target and terminates it when the target
finishes or at its time cap; outputs go to S3 as the job runs. The inputs are copied from an S3 mirror and checked
against `heidelberg/inputs.sha256` before anything runs.

## The record run and its blinding

Part one measures every pair and makes the selection blind; the checksums of the transient sample and of the coverage
map are registered before any analysis, and the selection is time-locked with drand (`tools/seal.py`). Part two, on a
new machine, checks every file against the registered checksums, runs the tests, and locks its results to drand round
37,911,412 (1 April 2027, UTC); the encrypted file's checksum is posted before that round. The unlock
(`tools/unseal.py`) opens only the file with that checksum.

## After the unlock

On 1 April 2027, once the encrypted results' checksum is recorded (`tools/registration.py seal`), one target opens them
and applies the frozen decision rule:

    ./run.sh snakemake after_unlock --config run=record -j4   # unseal, verdict, rate, and the Zenodo deposit

The same chain runs on a dry run's seal (`--config run=dry`), which writes its own seal record.

## Data and licences

The plate scans are those of the Heidelberg Digitized Astronomical Plates, served by the German Astrophysical Virtual
Observatory (GAVO). The code is released under the MIT licence (`LICENSE`) and the derived data under CC BY 4.0
(`LICENSE-data`).
