<!-- The pre-registration in the call's six sections (the template of the other Lost Light proposals): the file to edit.
Every number is a placeholder filled from work/results.json (tools/results.py) by tools/render.py (Snakefile rule draft).
Built on 1 November 2026 from the previous repository's full-length draft; the comparison with Solano et al. [7],
the replication recommendations, and the limitations are stated here in brief, with every setting in config/analysis.json. Claude drafts; Donal edits, and his edits are final. Items in square brackets
are to be settled before posting. -->

# Do flashes on simultaneous Heidelberg plates avoid Earth's shadow?

## 1. Abstract

(1) Rationale and plan. Transients on the First Palomar Sky Survey (POSS-I) are reported to be scarcer within Earth's shadow, as expected if some were glints from objects in high orbit [9]. The Bruce astrograph at Heidelberg-Königstuhl photographed each field on two plates simultaneously, so a genuine flash appears on both, whereas a plate defect appears on only one. We will test whether sources found on both plates of 1,321 twin pairs from 1900 to 1957 are scarcer where Earth's shadow fell at geostationary radius. If no sources are glints, we expect to limit their share to 4% at 95% confidence.

(2) Team and facilities. Dr Donal Hill, a data scientist in industry (PhD in particle physics, Oxford; formerly at CERN), is the sole analyst. Dr Markus Demleitner, an astrophysicist at the German Astrophysical Virtual Observatory, advises on the plate archive. The scans are public, and no ethics approval is required.

(3) Timeline. After acceptance, we will measure the plates and select the sources blind in February 2027, and run the statistical tests in March. The results remain sealed until 1 April 2027.

(4) Data sharing. The code, selected sources, coverage maps, null tests, and results will be made public.

## 2. Introduction

A transient is a point source recorded on one photographic exposure, which is absent from images of the same sky taken at other times. Villarroel et al. [9] reported that transients on the First Palomar Sky Survey (POSS-I) are scarcer in Earth's shadow. Brief glints of sunlight from an object in orbit cannot occur in Earth's shadow, whereas plate defects "do not care where the Earth's shadow falls on the plate" [8]. Earth's shadow therefore offers a test of the physical nature of transients. Since the POSS-I plates predate the first artificial satellite, a deficit caused by glints would indicate reflective objects in high orbit before 1957.

Two further analyses have reported the same deficit. Bruehl et al. [2] separated transients from plate defects using machine learning, and found that the transients remain scarcer in Earth's shadow. Doherty [4] re-analysed the original authors' catalogue and also reported the deficit. Watters et al. [13] argued that plate defects and selection effects could account for these signals. Villarroel et al. [11] replied that the subset analysed by Watters et al. is too small to test the shadow deficit, and that the principal findings stand.

On a single plate, a brief flash and an emulsion defect can look alike. The two can be separated only through the appearance of each image, whether via cuts on image shape or machine learning [2]; no independent exposure of the same moment exists to test that separation. The number of defects that survive among the transients thus cannot be counted on the plates themselves. A shadow test on single plates depends on how well background is modelled, and each background candidate passing the selection dilutes the shadow deficit that genuine glints would cause. A second plate, exposed at the same moment through a separate lens, removes such defects instead of modelling them.

The Bruce astrograph at Heidelberg-Königstuhl took such plates. It carried two lenses side by side, each exposing its own blue-sensitive plate of the same field at the same time. The two lenses are so close together that an object at geostationary radius (42,164 kilometres from Earth's centre) appears in the same direction through both. A flash of light from such an object lands at the same position on both plates, whereas an emulsion defect, a speck of dust, or a scanning fault affects only one; requiring a source on both plates can thus remove plate-level defects. Plates of the same field taken hours or nights apart cannot make this distinction, because a flash on one plate is missing from the other, just as a defect is.

Our study uses two separate sets of Bruce twin pairs. The main sample, described under Participants/Sample, holds the 1,321 pairs on which the hypotheses below are tested. The 46 tuning pairs, none of which Earth's shadow overlaps, were used to fix every rule of the analysis, and are not part of the main sample. Their results are given under Pilot Data.

The Lost Light challenge asks for a primary search on photographic plates from observatories other than Palomar, exposed before the launch of Sputnik in October 1957. The Bruce twin plates meet this requirement. However, of the Bruce twin pairs exposed from 1949, when POSS-I began, Earth's shadow overlaps only one, exposed on 30 September 1957, four days before the launch of Sputnik. Every other shadow pair is from before 1949, the last of them from 1937. This gap reflects the archive itself: the Bruce astrograph took 2,647 twin pairs before 1935, but only 121 from 1935 until the launch of Sputnik.

The challenge prefers telescopes of at least 0.8 m aperture. The Bruce lenses are 400 mm across, so their plates are shallower than those of POSS-I and record only bright flashes. We accept this loss of depth in return for a stronger rejection of plate-level background.

Even with strong plate-level background rejection, some sources other than glints can still appear on both plates: variable stars and novae, asteroids, and the trails of meteors and aircraft. We remove these using catalogues, computed asteroid positions, and rules on image shape. These rules, together with the requirement of a source on both plates, form our selection.

Two types of background remain after this selection, and must be measured. The first are chance coincidences: two unrelated detections, one on each plate, that happen to lie within three arcseconds of each other. We measure this background on every pair by shifting the detections on one plate by 20 arcseconds and repeating the selection. Detections already found on both plates, which include any glint or static source, are left out of the shift, so every source that passes afterwards is a chance coincidence. We repeat this process in eight directions, and take the average count as the pair's chance estimate.

The second remaining background is from static sources, namely stars or galaxies that are genuine but missing from today's catalogues. To bound them, we pair plates from the main sample with plates of the same fields taken on other nights. No flash can appear on both plates of a different-night pair, whereas a static source can. Across the 39 tuning pairs that the selection accepted, 33 sources passed compared to 33.6 expected from chance coincidences alone. Neither static sources nor glints thus added a measurable excess. This count of fewer than one source per pair compares with an average of 167 transients per plate in the POSS-I sample of Villarroel et al. [9]. The shallower Bruce plates account for part of this difference, but only a second plate can remove plate-level defects.

In our approach, the size of the remaining background is thus measured on the plates themselves. A count alone, however, cannot single out glints, because a glint and a static source both contribute to a count in the same way. Only Earth's shadow can separate them.

Since each pair holds fewer than one selected source on average, our shadow test pools all of the pairs. The 442 pairs crossed by the shadow hold as much shadowed sky as 208 pairs wholly in shadow. At this pre-registration stage, the main sample has of course not yet been measured. Scaling the tuning pairs' chance estimate to the full dataset, we forecast around 970 selected sources, of which 150 would fall in Earth's shadow if none are glints. The Bruce plates thus give a shadow test that is fully independent of POSS-I, on a large sample with a background we can measure. If there are no glints, the test is expected to limit them to about 4% of the selected sources.

### Hypotheses

A selected source is one that appears on both plates of a pair, and passes every rule of our selection. Our test treats each selected source as either background or a glint. Background sources, such as chance coincidences and static sources, do not depend on Earth's shadow. Glints, by contrast, cannot occur in Earth's shadow, and we assume that they are otherwise spread evenly over the sky.

The common exposure of a pair is the time during which both of its plates were exposed. The shadow share, s, of a given location on a pair is the fraction of the common exposure that the location spent in Earth's shadow. The shadow share thus runs from s = 0 in full sunlight, to s = 1 in full shadow throughout the common exposure. Near its edge, Earth's shadow is partial, and s then counts the fraction of the Sun's disc that Earth hides, averaged over the common exposure. We compute s where the line of sight from Königstuhl reaches geostationary radius. Seen from Königstuhl, Earth's shadow moves across the sky during an exposure, because the observatory moves with Earth's rotation. During the common exposures of the shadow pairs, which have a median exposure time of 150 minutes, the shadow moves by a median of 3.6 degrees, which is large compared to the 6 x 6 degree area tested. We therefore evaluate the shadow at 21 instants spread evenly over each common exposure, rather than at its midpoint alone.

At a location with shadow share s, we model the rate of glints as reduced by a factor 1 − κs. If κ = 1, glints vanish in shadow, and if κ = 0, the shadow has no effect on them. We express the strength of a glint signal as the glint share, namely the fraction of the selected sources on an average pair that are glints. If η is the expected number of glints on an average pair, the glint share is η / (η + B). In this expression, B is the average chance estimate over the main sample's pairs, which does not depend on Earth's shadow. B converts a glint share into a number of glints on the plates tested, and the limit expressed as a rate of glints per square degree per hour is a limit on η itself. The background level in the fit is a separate quantity, fitted to the data as described under Planned Analyses. The main sample has not yet been measured, so our forecasts of the power and expected limit use B = 0.86 per pair, the chance estimate averaged over the 39 tuning pairs that the selection accepted. This estimate leaves out any source found on both plates, so genuine glints found on both plates do not contribute to it.

#### H1 (primary outcome): selected sources are suppressed in Earth's shadow

H1 states that κ > 0: selected sources are suppressed in Earth's shadow at geostationary radius, as seen from Königstuhl. The null hypothesis, H0, is κ = 0, under which the rate of selected sources does not depend on Earth's shadow. A deficit is claimed only if both shadow tests described under Planned Analyses reject H0 (each one-sided at the 5% level), and if none of the checks point to another cause. A claim would be evidence of a shadow deficit from an instrument, site, emulsion, and era different from those of POSS-I. However, it would not establish what the selected sources are.

If neither test rejects H0 and every check passes, the result is an upper limit on the glint share, and on the rate of glints per square degree per hour, for bright flashes recorded on the shadow pairs' plates between 1900 and 1957. Applying this limit to the POSS-I transients requires converting between the two surveys' brightness and colour, and assuming that the population of transients did not change between the two eras. Any other result is inconclusive, and no claim will be made.

#### H2 (secondary outcome): a deficit follows the parallax of geostationary radius

Because of the parallax of the orbit, Earth's shadow at geostationary radius falls on a different part of the sky as seen from Königstuhl than as seen from Earth's centre. Glints from objects at geostationary radius would be missing only where the shadow falls as seen from Königstuhl. H2 states that any deficit is stronger with the shadow computed from Königstuhl than with the shadow computed from Earth's centre. We repeat both tests with the shadow computed from Earth's centre and report their statistics beside those of the primary test. Like the other secondary analyses, this comparison is exploratory and does not change the outcome.

## 3. Methods

### Participants/Sample

The units of the study are twin pairs of Bruce plates. Before any plate was measured, we identified 2,768 Bruce pairs that were exposed together for at least 30 minutes between 1900 and the launch of Sputnik in October 1957. For each pair, the test area is the central 6 x 6 degrees of the plates.

A shadow pair is one in which Earth's shadow at geostationary radius (as seen from Königstuhl) overlapped the test area during the exposure. Most shadow pairs are only partly in shadow, so they also contribute sunlit area. We matched each shadow pair with two sunlit pairs, which Earth's shadow did not overlap. The shadow pairs were taken one at a time, in random order, and each was assigned the sunlit pair closest to it in declination, season, and year, where declination is the sky's equivalent of latitude. A second pass, in the same order, assigned each shadow pair a second sunlit pair.

The sunlit pairs provide the likelihood test with a reference: they show the rate of selected sources on similar plates, at similar positions on the sky and times of year, with no shadow effect. Without them, this reference would come only from the sunlit parts of the shadow pairs themselves, which always lie next to the shadow and are wholly missing on 42 pairs whose test area lies in full shadow for at least 95% of the exposure.

The within-pair test does not use the sunlit pairs, since it compares positions within each shadow pair. In simulated samples, the median expected upper limit on the glint share is 6% with one sunlit pair per shadow pair and 4% with two, with no further improvement in limit found using additional sunlit pairs. We therefore elect to use two sunlit pairs per shadow pair.

We remove pairs whose times are uncertain, and pairs exposed with the Sun above an altitude of −12 degrees or with the plate's centre below an altitude of 10 degrees. The main sample consists of the 442 shadow pairs and 879 sunlit pairs that remain. The scans of the main sample were copied from the GAVO archive in November 2026. A pair whose scans prove unreadable will be removed before the selected sources are registered, and reported as a shadow or sunlit pair.

A further 29 pairs, referred to as different-night pairs, each combine one plate of a main-sample pair, drawn at random, with a plate of the same field taken on a different night. No flash can appear on both plates of such a pair, so the sources passing our selection requirements on these pairs measure how many static sources the selection lets through on average across the main sample's sky.

Each plate's time comes from the Heidelberg plate catalogue, which agrees within 30 seconds with the scanned observation logbooks on the twelve plates checked. The logbooks confirm that the times were copied correctly, but not that the observatory's clocks were correct. We therefore also time plates using the trails of known asteroids, compared with their computed tracks.

On 19 plates of the tuning pairs, this asteroid clock gives a median offset from the catalogue's times of −1.3 minutes, with a standard error of 1 minute. Part one of our analysis plan measures the asteroid clock on the main sample before the coverage map is produced. If the median offset differs from zero by more than three standard errors and by more than five minutes, every plate's time is shifted by it; otherwise the catalogue's times are kept. The times used, and the decision, will be registered with the selected sources.

### Sample size and power

We estimate the power of the test using simulated samples. Each simulated sample contains the same 1,321 pairs as the main sample, with the same sky positions, exposure times, and coverage by Earth's shadow. The background of each simulated pair is drawn around the average measured on the tuning pairs, and glints are added at a chosen glint share. The power is the percentage of simulated samples in which each test alone finds a deficit, or in which both do, which is the condition for a claim (the null tests and the comparison of data quality are not simulated):

| Glint share | Likelihood test alone (κ = 1) | Within-pair test alone (κ = 1) | Claim (κ = 1) | Claim (κ = 0.5) |
|---|---|---|---|---|
| 10% | 16% | 17% | 5% | 2% |
| 15% | 53% | 25% | 18% | 6% |
| 25% | 97% | 55% | 54% | 18% |

The test is therefore likely to claim a deficit only if a quarter or more of the selected sources are glints that vanish in shadow. For smaller glint shares, including zero, the most likely outcome is an upper limit on the glint share.

In simulated samples with no glints, the median 95% upper limit is a glint share of 4%, and the limit lies between 2% and 8% in two thirds of simulated samples. With no glints and a correct background model, a deficit was claimed in none of 1,000 simulated samples, where the 5% significance level allows up to one in twenty.

### Materials and Procedures

The scans are those of the Heidelberg Digitized Astronomical Plates, served by the German Astrophysical Virtual Observatory (GAVO). We keep a compressed copy of each scan, confirmed to restore every pixel and header card of the original, such that every scan is read only once from the archive.

Each plate is measured on its own. Its plate solution, the conversion from positions on the plate to positions on the sky, comes from matching its stars to the Gaia star catalogue [5]. A candidate is a pair of detections, one on each plate, that are each other's nearest neighbour within three arcseconds. We measure the brightness on each plate at the positions of both detections, and keep the candidate only if all four measurements are at least 7 times the noise (the spread of brightness measured at random positions on the plate). Solano et al. [7] required 30 times the noise on a single plate. Since the second plate already rejects defects, a lower threshold keeps fainter flashes, and a value of 7 was set on the tuning pairs.

We follow the definition of a transient given by Solano et al. [7]: a point source, seen at one epoch only, with no counterpart in modern optical catalogues, which is not an artefact and is confirmed on a second, independent image. We therefore remove a candidate with any Gaia source within five arcseconds. We also remove a candidate with a Pan-STARRS source within five arcseconds that is at most two magnitudes fainter than the candidate.

A candidate is kept only if it appears as a point source on both plates, and if its two images are as similar to each other as two images of a star of the same brightness. A flash lasting less than a second can look sharper than a star [10], so the cut on image width removes only images broader than a star's. A selected source must also lie away from stars, catalogued objects, known asteroids, and active meteor radiants. Four further cuts remove flaws in the plate or scan: marks on the scanner, which recur at the same scan pixel on at least three scans; copies of stars at a fixed offset, caused by a double exposure or a guiding jump; sources surrounded by a dark ring, a defect of the plate's development; and locations where the plate solution is poor, with fewer than 50% of the 30 nearest bright detections matching Gaia. One function applies every cut to candidates, chance coincidences, and artificial flashes alike.

Where our selection departs from that of Solano et al. [7], the reasons are as follows. Image shape is judged against each plate's own stars, since the limits that Solano et al. set on image size (in pixels) are specific to the POSS-I scans. Bright stars and catalogue counterparts are taken from Gaia, at positions moved to each plate's epoch. Known asteroids are removed using our own ephemeris, computed for each plate's times, so that our selection depends on no external service. The Lost Light call also points to the criteria of Villarroel et al. [12]. These were set for images from the Zwicky Transient Facility, a modern survey, and give no thresholds for photographic plates. Where they do carry over, our vetoes are wider: they remove sources within 1.5 arcseconds of a catalogued object and within 10 arcseconds of a known asteroid, against our 5 and 60 arcseconds.

To measure how many genuine flashes the selection keeps, we introduce around 2,000 artificial flashes into both scans of every pair, each at the same position on the sky in both scans. Each artificial flash is built from the images of nearby stars of the Sun's colour on the same plate. The efficiency is the fraction of star-like artificial flashes that pass every selection step applied to the real candidates. We define brightness as a magnitude in the Gaia G band, within which larger numbers are fainter. On the tuning pairs, the efficiency fell to half its peak at a median brightness of G = 14.75. The artificial flashes span G = 11 to 17.3. Their faint end lies well beyond the half-efficiency point, so they trace the decline in efficiency towards faint flashes. A pair is refused if its artificial flashes fail the criteria for looking like real sources.

We count chance coincidences on each pair by shifting the detections on one plate that have no counterpart on the other plate by 20 arcseconds, and repeating every selection step. Detections found on both plates, which include any genuine flash, are left out of this step, so every source it selects is a chance coincidence. We repeat this in eight directions, and the pair's chance estimate is taken as the mean count over the shifts.

Requiring a source on both plates applies the coincidence method of Bothe and Geiger [1]. Walther Bothe, coincidentally a professor in Heidelberg from 1932 (pardon the pun), shared the 1954 Nobel Prize in Physics for the method. Coincidences between independent detectors have since become a standard method for rejecting background in particle physics, including in flagship experiments at CERN. Those experiments measure accidental coincidences by delaying one detector's signal, the equivalent of shifting one plate's detections in our analysis.

Aside from how many background sources to expect, our tests also need to know where a selected source could have been found on each pair. The control positions provide this: they are the positions of artificial flashes that are found on both plates and pass every cut on position. The control positions thus share the real candidates' losses to bright sky, plate edges, crowding, and bright stars.

This registration freezes the design, code, and settings, and records their SHA-256 checksums. The frozen files are archived on Zenodo (https://doi.org/10.5281/zenodo.24090376), and the SHA-256 checksum of the freeze is 1ed306efc7064d6c9bee5ac26c1e76326309f6834b77d9bc3ad9dd49969c094a. Every setting of the analysis is listed, with its reason, in the repository's configuration file (config/analysis.json). The analysis of the main sample will then run in two parts. Part one measures every pair and selects the sources blind. It then computes the SHA-256 checksum of the list of selected sources, and encrypts the selection's outputs (including the list) using drand timelock encryption [6], which cannot be opened before a set time. That time falls after part one ends and before part two starts, and is fixed when part one is scheduled.

Before any analysis, we publish the lock time and the checksums of every file that part one registers in a timestamped post on ResearchHub. These files include the list's checksum, the encrypted file's checksum, the plate times used, and the coverage map, which records each pair's plates, common exposure, test area, and role (shadow, sunlit, or different-night).

Part two runs on a new machine. It checks every file against the registered checksums, halts if any file differs, opens the encrypted selection, and halts unless the list of selected sources matches its registered checksum. It then runs the statistical tests, and encrypts their results on the selected sources to drand round 37,911,412, at the start of 1 April 2027 (UTC). The encrypted file's checksum is posted before that round, so every run of part two is on record. The outcomes of the checks that can prevent a claim, described under Planned Analyses, are reported before the unlock.

After the unlock, we publish the list of selected sources, and confirm that it matches the checksum registered in part one. Any change to the code or settings after this registration will be made before part one, and posted as a dated amendment with its reason.

### Planned Analyses

Each plate pair i has an expected number of background sources Bᵢ, which pass the selection on both of its plates. Bᵢ includes the chance coincidences between the two plates' separate backgrounds, such as dust on one plate and a defect on the other, so it is measured for the pair as a whole rather than for each plate. These background sources are spread over the test area by a smooth background map wᵢ, fitted to the chance coincidences. Glints add an expected number η per pair, scaled by the pair's efficiency times its common exposure relative to the average pair (Eᵢ), and reduced in shadow by the factor 1 − κs. The density of selected sources at position x on pair i is therefore

λᵢ(x) = Bᵢ wᵢ(x) + η Eᵢ [1 − κ s(x)].

The likelihood test fits η, κ, and the backgrounds to the positions of the selected sources on all plate pairs. Each pair's background is held near its chance estimate, after every chance estimate is multiplied by a common factor. This factor is the number of selected sources divided by the total of the chance estimates. The overall level of the background is thus fitted from the data, and only its variation from pair to pair comes from the chance estimates. An error common to all chance estimates therefore cannot appear as a shadow deficit.

The test statistic, q, is twice the improvement in log-likelihood over a fit with κ = 0 if the fitted κ is positive, and zero otherwise. This is equivalent to the one-sided profile likelihood ratio used in searches for new particles at CERN's Large Hadron Collider [3]. Since each pair holds fewer than one selected source on average, we calibrate it with simulated samples rather than with its asymptotic distribution. We compute its p-value from simulated samples at each glint share on a grid, and take the largest p-value over the grid.

The within-pair test, adapted from Bruehl et al. [2], compares the selected sources only with the control positions on the same pair. Without a shadow effect, the shadow shares of a pair's n selected sources sum on average to n times the mean shadow share of its control positions, each weighted by the background map. The test statistic, D, sums the shortfall from this value over all pairs.

We require both tests because they rest on different comparisons. The likelihood test draws much of its power from the comparison between shadow and sunlit pairs, which can be shifted by anything that varies with distance from the point opposite the Sun. The within-pair test compares sources only with control positions on the same pair, so differences between pairs in the level of the background or in the efficiency cannot shift it. It depends on the other pairs only through the shape of the background map, which is fitted to all pairs together.

The result is one of three outcomes, decided by the code frozen in this registration. A deficit is claimed if both tests give a p-value below 0.05, both null tests pass, and shadow and sunlit pairs are comparable in terms of data quality. The result is null if neither test gives a p-value below 0.05 and every check passes, in which case the upper limit is the result. In every other case the result is inconclusive. No claim is then made, the reason is stated, and the upper limit is reported as resting on a background model that the failed condition calls into question.

Deviations from the planned sample are handled without redrawing it. Pairs whose scans cannot be read, whose plates cannot be solved, or whose artificial flashes fail the criteria are left out of both tests, reported by role, and not replaced.

For glints that vanish in shadow (κ = 1), the upper limit is the largest glint share that is not ruled out at 95% confidence. We set this limit with a Neyman construction on simulated samples. We also express the limit as a rate of glints per square degree per hour, both for all flashes from G = 11 to 17.3 and for flashes brighter than G = 14.75, where the efficiency falls to half its peak. We always report κ, η, both p-values, the common factor applied to the chance estimates, and the upper limit.

Two types of check can prevent a claim. Firstly, we rerun both tests on two populations that cannot contain glints:

- Chance coincidences found with one of the eight shifts. These are left out when fitting the background map, so that the test on them is independent of the fit.
- Artificial flashes that pass the selection, at most five per pair, chosen at random. They are introduced without regard to Earth's shadow, so a deficit among them would mean that the selection itself keeps fewer sources in shadow, for example if the shadowed parts of the plates happen to be of poorer quality.

A null test fails if the likelihood test gives a p-value below 0.01, or if the within-pair test gives one below 0.005 in either direction, so an excess fails as well as a deficit. A failure on either population prevents a claim.

Secondly, we compare the shadow and sunlit pairs across five measures of data quality: how often both plates were solved, how often the selection accepted the pair, the depth (the brightness at which half of the Gaia stars are detected), the noise in the brightness measurements, and the efficiency. For each measure, we take the difference between shadow and sunlit pairs in units of their spread. A difference larger than 0.5 in either direction on any measure prevents a claim.

The following secondary analyses are also reported alongside the result. They help to interpret the main result, but do not change the outcome:

- A positive control, in which each artificial flash is removed at random with a probability equal to its shadow share. This creates a known deficit in shadow, which both tests should find.
- Both tests with the shadow seen from Earth's centre. This tests whether a deficit follows the parallax expected at geostationary radius.
- The test of Bruehl et al. [2] as published, which compares selected sources in shadow with random positions pooled over pairs, for a direct comparison with that work.
- The asteroid clock's offsets, plate by plate, and the upper limit on static sources from the different-night pairs.

The design follows the recommendations of Villarroel [8] for replications. Each transient is judged by its shape against the plate's own stars and must be point-like, and the catalogue vetoes are optical. Whole-plate scans are used, with the analysis confined to the central 6 x 6 degrees of each pair, and the plate times are checked against the observation logbooks and the asteroid clock. The upper limit comes from simulated samples that carry the uncertainties of the background and its map. The control positions are artificial flashes, which are neither stars nor galaxies and are not built from the observed transients.

The design departs from three of the recommendations. Firstly, the sample is chosen by Earth's shadow rather than for homogeneous coverage of the sky, but each source is compared with control positions on the same plates, and the sunlit pairs are matched in declination, season, and year.

Secondly, Villarroel [8] recommends normalising by the total counts of transients and controls over all plates. Our within-pair test instead compares each pair's selected sources only with that pair's control positions. It thus gives up the comparison between pairs, which the likelihood test already makes. For direct comparison, we also report the test of Bruehl et al. [2], which normalises by the total counts.

Thirdly, Villarroel [8] suggests judging the power against the POSS-I result of 7.6 standard deviations [9]. A significance depends on the size of the sample it was measured on, so it does not carry over to plates of a different depth, era, and number. We therefore state the power against the glint share, which is a property of the Heidelberg sources themselves.

Below we summarise the study's main limitations, and what guards against each:

- Photographic emulsion records a very brief flash more efficiently than steady light of the same energy. Converting the upper limit into a rate of glints of a given energy therefore requires a correction, which we carry as a range of 0.5 to 2.0 magnitudes. The correction affects only this conversion, not the shadow tests or the glint share.
- The background is modelled by a single smooth map across the test area, fitted to all pairs together. The null test on chance coincidences checks this map, and its failure prevents a claim. We also simulated samples in which the true background differs from the fitted map, being either uniform or half as steep. A deficit was then claimed in 6% and 2% of these samples respectively, against a significance level of 5%. These simulations leave out the null test on chance coincidences, which can only prevent claims, so the rates are upper bounds.
- Static sources - genuine stars or galaxies missing from modern catalogues - do not depend on Earth's shadow, but are likely to be more common in crowded fields. More static sources appearing on sunlit pairs compared to shadow pairs could mimic a deficit in the likelihood test, but not in the within-pair test, which a claim also requires. Five of the six different-night pairs measured alongside the tuning pairs passed our selection; on these pairs, static sources exceed the chance estimate by at most 1 per pair at 95% confidence. The main analysis will measure this rate again on 29 different-night pairs. Static sources are not counted by the chance estimates, so we checked whether they could produce a false deficit. In simulated samples with no glints, but with more static sources on every pair than the bound above, the likelihood test alone found a deficit in 4% of samples. At the 5% significance level, the test is expected to do so in 5% of samples with no signal at all, so these static sources do not raise the rate of false deficits.
- The sunlit pairs lie further from the point opposite the Sun than the shadow pairs, which the likelihood test allows for through each pair's own efficiency and chance estimate. Uncatalogued asteroids, the main remaining difference, could hide part of a genuine deficit but not create one, and near that point they move far enough during an exposure to fail the cut on image shape. The within-pair test compares positions within each pair, so it is not shifted by such differences.
- Observers sometimes marked objects on the plates in ink, and a mark drawn at the same object on both plates of a pair would appear on both. Such marks were made without regard to Earth's shadow, so, like static sources, they can add to the background but cannot produce a deficit.

### Ethics and Data Management

The study uses public archival scans and catalogues, and involves no human or animal subjects; with no personal data, no confidentiality measures are required. GAVO's operators agreed to the volume of scans we will read. We acknowledge the Heidelberg Digitized Astronomical Plates project and GAVO at their request, and will acknowledge Project Lost Light in the manuscript. We will publish the code, frozen snapshots of the catalogues, the selected sources with their coverage map, every null test, and the results. The measurements of each plate are not deposited, since the code regenerates them exactly from the public scans. The code will be hosted on GitHub under the MIT licence, and archived on Zenodo with a DOI, together with the data. Our compressed copies of the scans are kept on Amazon Web Services. The manuscript will be deposited on arXiv and submitted to an open-access journal.

AI assistance (Anthropic Claude Max personal subscription) was used to build the code and checks behind this pre-registration, and to draft and edit the text. The analyst set the design, made every decision, and checked the results. No AI models were used to judge transient candidates: every cut is a deterministic computation with logged inputs and outputs. The same division of work, and the same exclusion of AI models from judging candidates, will apply to all work undertaken after this pre-registration.

There are no conflicts of interest in this work.

## 4. Pilot Data

We set every analysis rule on the 46 tuning pairs, none of which Earth's shadow overlaps. On these pairs, the selected sources are consistent with chance coincidences alone: 33 sources were selected, against a chance estimate of 33.6. Of the tuning pairs, 44 had both plate solutions accepted. The selection accepted 39 of these, and refused the remainder because their artificial flashes failed the criteria. On the different-night pairs, the 95% upper limit on static sources, above the chance estimate, is around 1 per pair.

End-to-end analysis dry runs have been executed on the tuning pairs, covering all steps from measuring the plates to generating time-locked results. These runs confirm that the data pipeline runs correctly, and produces the same results each time.

## 5. Budget

All costs are for computing and storage on Amazon Web Services. We measured the cost of each step, using a dry run of the full analysis on 52 pairs (the 46 tuning pairs and six different-night pairs). The estimates below scale these costs up to the 1,350 pairs of the main analysis, at on-demand prices, with an equal allowance covering a rerun of any step after a technical failure:

| Item | Cost |
|---|---|
| Copying the 2,644 scans from GAVO, compressed (2,583 GB) | $4 |
| Measuring every plate (2.2 hours on a 192-core machine) | $29 |
| Part one: the blind selection | $9 |
| Part two: the statistical tests | $15 |
| Calibration on the frozen code | $4 |
| Two further dry runs on the tuning pairs | $5 |
| Storing the compressed scans and the measurements for twelve months | $96 |
| Estimated cost | $163 |
| Allowance for reruns after technical failures | $163 |
| Total requested | $326 |

The analyst's time is not charged, no publication charge is requested, and the computing already spent on the pilot, calibrations, and dry runs (about $133) is not included.

## 6. References

1. Bothe, W. & Geiger, H. Über das Wesen des Comptoneffekts. *Z. Phys.* **32**, 639 (1925).
2. Bruehl, S., Doherty, B., Streblyanska, A. & Villarroel, B. Use of machine learning to enhance detection of transient astronomical phenomena in historical observatory images. *Sci. Rep.* **16**, 31890 (2026).
3. Cowan, G., Cranmer, K., Gross, E. & Vitells, O. Asymptotic formulae for likelihood-based tests of new physics. *Eur. Phys. J. C* **71**, 1554 (2011).
4. Doherty, B. Independent replication of nuclear test-transient correlations and Earth shadow deficit in POSS-I photographic plates. Preprint at arXiv:2604.00056 (2026).
5. Gaia Collaboration, Vallenari, A. et al. Gaia Data Release 3. *Astron. Astrophys.* **674**, A1 (2023).
6. Gailly, N., Melissaris, K. & Romailler, Y. tlock: practical timelock encryption from threshold BLS. Preprint at IACR Cryptology ePrint Archive 2023/189 (2023).
7. Solano, E., Villarroel, B. & Rodrigo, C. Discovering vanishing objects in POSS I red images using the Virtual Observatory. *Mon. Not. R. Astron. Soc.* **515**, 1380 (2022).
8. Villarroel, B. Replications, validations and exciting developments (section "Recommendations for future replication efforts"). The VASCO network, vasconsite.wordpress.com (11 July 2026).
9. Villarroel, B. et al. Aligned, multiple-transient events in the First Palomar Sky Survey. *Publ. Astron. Soc. Pac.* **137**, 104504 (2025).
10. Villarroel, B., Solano, E. & Marcy, G. W. On the image profiles of transients in the Palomar Sky Survey. Preprint at arXiv:2507.15896 (2025).
11. Villarroel, B., Streblyanska, A., Bruehl, S. & Geier, S. A response to "Critical evaluation of studies alleging evidence for technosignatures in the POSS1-E photographic plates" by Watters et al. (2026). Preprint at arXiv:2602.15171 (2026).
12. Villarroel, B., Watters, W. A., Streblyanska, A., Solano, E., Geier, S. & Mattsson, L. A cost-effective search for extraterrestrial probes in the Solar system. *Mon. Not. R. Astron. Soc.* **546**, staf1158 (2026).
13. Watters, W. A. et al. Critical evaluation of studies alleging evidence for technosignatures in the POSS1-E photographic plates. Preprint at arXiv:2601.21946 (2026).
