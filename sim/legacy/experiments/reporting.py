"""Descriptions tied to the scope of each experiment and its computed metrics.

Keep these in one place so rebuilding a report cannot restore an outdated claim.
The numerical metrics remain the experiment's own output.
"""
from pathlib import Path


def refine_summary(summary):
    t = summary["id"]
    m = summary.get("metrics", {})
    captions = {}
    if t == "t01_static_pyramid":
        summary.update(
            title="A stationary pyramid produces depth-like patterns",
            finding=f"The reconstructed pipeline produces layered patterns although the simulated scene has no motion or subsurface. Mean best window-fit scores are {m['mean_best_adjR2_pyramid']:.3f} inside the ground footprint and {m['mean_best_adjR2_desert']:.3f} outside it.",
            limitations="One stylized scene and one random seed. The footprint is a ground-plane mask, not a scatterer attribution map: layover mixes locations. A high fit score is not a probability of a buried structure.",
        )
        captions["t01_transect_rg"] = "<b>Through the near face and desert.</b> Both focusing variants produce depth-like patterns along this range transect. The image coordinate mixes height and ground range through layover."
        captions["t01_maps"] = "<b>Where the best fits occur.</b> The input amplitude, maximum window-fit score, and winning model depth on a 1 m target grid. High scores also occur away from the pyramid. The winning depth is a processing output, not a surveyed location."
        captions["t01_trajectories"] = "<b>Three static trajectories.</b> Reference/offset patch shifts at a near-face, side-face, and desert location. None of the simulated scatterers moves. The curves show that changing the sub-aperture can change the measured patch displacement."
    elif t == "t02_desert_null":
        summary.update(
            title="Empty desert produces a similar score distribution",
            finding=f"The mean best score is {m['mean_best_adjR2_pyramid']:.3f} on the pyramid footprint and {m['mean_best_adjR2_desert_scene']:.3f} in the desert-only scene. Their empirical distribution distance is {m['ks_statistic']:.4f}. The absence of buried objects does not prevent structured output.",
            limitations="One realization per scene; the 32-pixel patches overlap on a 4-pixel grid, so the samples are dependent. The saved conventional KS p-value is descriptive only here. Similar distributions do not establish statistical equivalence or validate a detector.",
        )
        captions["t02_histograms"] = "<b>Overlapping distributions.</b> Best window-fit scores and their winning model depths for the pyramid and desert controls. Compare their effect sizes; the many overlapping patches are not independent experimental replicates."
        captions["t02_maps"] = "<b>Same scale, two scenes.</b> Best-fit score maps for a pyramid scene and a desert-only scene. Both inputs contain no buried objects. This comparison tests this score's specificity, rather than proving the two scenes identical."
    elif t == "t03_mechanism":
        summary.update(
            title="Surface spacing predicts the apparent depth",
            finding=f"In the two-scatterer examples, peak depth tracks azimuth separation at approximately {m['depth_per_metre_separation']:.2f} model-metres per metre. Farther scatterers can contribute outside the patch, and the sampled depth axis folds high spatial frequencies back by aliasing.",
            limitations="Ideal point scatterers and a prescribed spectral model. The wall's comb fit depends on reflectivity: the ideal-wall fitted spacing differs from the prediction, while the randomized variants are closer. These examples identify a mechanism in this reconstruction, not the cause of a particular published feature.",
        )
        captions["t03_wall"] = "<b>A periodic wall.</b> Dotted lines mark the predicted course harmonics. The fitted fundamental is about 6.9 m for the two randomized variants versus a 7.0 m prediction; the ideal-wall fit is about 9.8 m. The prediction is not equally accurate in every scene."
        captions["t03_heading"] = "<b>Rotate the same wall.</b> Fitted spacings in these three simulations track the projected-spacing prediction. This is a controlled geometry test; a comparison of real acquisitions must also account for visibility, resolution, and registration."
    elif t == "t04_parameters":
        summary.update(
            title="Processing choices set the axis and change the patterns",
            finding="Patch size changes spatial smoothness. The declared wavelength scales the model-depth axis, and the sub-aperture width changes its sampling and the resulting pattern. The wavelength panel illustrates an exact algebraic rescaling; a separate end-to-end test verifies that relation.",
            limitations="One scene per sweep, with one parameter changed at a time. Sensitivity alone does not invalidate an inverse method; the relevant question is whether calibrated, independently validated targets remain correctly localized. No such validation is supplied by this experiment.",
            method="Desert and pyramid transects with patch widths 16, 32, and 64 pixels and sub-aperture fractions 30%, 50%, and 70%. The wavelength panel relabels a single computed result at 0.24, 0.48, and 0.96 m. Band edges are snapped to FFT bins.",
        )
        captions["t04_lambda"] = "<b>An exact change of scale.</b> A single computed score image with its depth labels multiplied by the wavelength ratio. This illustrates the formula; it is not three independently simulated detections. Without an independently justified wavelength, the depth scale is unconstrained."
        captions["t04_bandwidth"] = "<b>Different sub-aperture widths.</b> The same transect processed at 30%, 50%, and 70% bandwidth. Each panel uses its own Nyquist depth extent. A common-axis comparison and known targets would be needed to assess localization stability."
    elif t == "t05_motion":
        best = max(max(v) for v in m["snr_vs_static"].values())
        summary.update(
            title="Motion changes the output; frequency becomes model depth",
            finding=f"At 1 cm imposed amplitude, the response reaches {best:.2f} times the static trajectory rms. The focusing step places vibration near {m['depth_per_hz_m']:.2f} model-metres per hertz; the 2 Hz case peaks at a second harmonic. This is a response measurement, not a calibrated detection limit.",
            limitations="One sinusoidal patch and static background clutter. Subtracting the exactly known motion-free baseline is an advantage available in simulation. No sensor-noise distribution, false-alarm calibration, or cavity-to-surface model is included. The result does not establish a general limit on radar vibrometry.",
        )
        captions["t05_sensitivity"] = f"<b>Response relative to a known baseline.</b> The change in sampled patch trajectories with and without imposed vibration. The dashed line is static trajectory rms. The maximum ratio is {best:.2f}; this ratio is not a detection probability. Shaded amplitudes are illustrative reference values."
        captions["t05_depth_is_frequency"] = "<b>The isolated motion response.</b> Surface oscillations at 1, 2, and 5 Hz are placed along a model-depth axis. Dotted lines mark the fundamental-frequency prediction; the 2 Hz case is dominated by its second harmonic. These are moving surface scatterers with no buried target."
    elif t == "t06_ordering":
        shape = m["ellipse_shape_gate_pass"]
        full = m["ellipse_full_gate_pass"]
        summary.update(
            title="The full ellipse gate rejects the sampled static pixels",
            finding=f"The corrected per-window shape gate accepts {shape['pyramid pixels']:.1%} of pyramid and {shape['empty desert']:.1%} of desert pixels. After the 0.005 px minor-axis floor, the full gate accepts {full['pyramid pixels']:.1%} and {full['empty desert']:.1%}. The earlier shuffle comparison is retained as a diagnostic only.",
            limitations="One scene per class and 1,500 dependent pixels. Arbitrary look permutations change correlations between overlapping sub-apertures, so the shuffle statistic is not a calibrated physical null test. The new independent-scene controls replace it as the primary null comparison. Shape-only acceptance is not a false positive of the full gate; sensitivity must also be tested on physical positives.",
            method="For each of 26 contiguous W25 windows, choose the best ellipse mode from 1–6, then apply adjusted R² ≥ 0.25, axis ratio ≥ 0.1 and minor axis ≥ 0.005 px. Pixel acceptance requires at least one fully accepted window. Legacy shuffle results use 60 arbitrary order permutations and are descriptive only.",
        )
        captions["t06_permutation"] = "<b>Check all parts of a gate.</b> The right panel excludes the amplitude floor; including it rejects all sampled static pixels. The two shuffle panels are historical diagnostics: permutations do not preserve overlapping-look covariance, so their nominal 5% reference is not a validated false-alarm rate. Per-window selection has been corrected and the scene results rerun."
    elif t == "t07_wells":
        summary.update(
            role="illustration",
            tag="illustration",
            eyebrow="Appendix · Constructed rendering",
            title="Illustration: surface points and rendering choices",
            question="What can an intentionally constructed surface scene look like after depth smoothing and chosen display settings?",
            finding=f"In this constructed scene, dominant point scatterers produce dark columns ({m['power_ratio_point_wells_vs_desert']:.2f} times desert power), while mixed face edges produce bright bands ({m['power_ratio_bright_bands_vs_desert']:.2f} times desert power). The displayed volume ends at {m['depth_axis_max_m_paper']:.1f} model-metres, or {m['depth_axis_max_m_slides']:.1f} m after a chosen wavelength rescaling.",
            limitations="Eight bright points are placed by hand. The 2.2 m wavelength, smoothing and display thresholds are chosen for illustration. This example does not apply T6's acceptance gate. Neither the point count nor the roughly 630 m display floor is an independent prediction. Visual resemblance does not establish the origin of the Khafre features, validate the surface targets' realism, or reproduce the original processing. This illustration neither reproduces nor measures spiral ramps.",
        )
        captions["t07_compare"] = "<b>Published and synthetic images.</b> Left: a Khafre Research Project slide, reproduced at reduced size for critique; copyright its authors. Right: a synthetic slice with no buried objects, using a similar display style. Resemblance is an illustration, not a reproduction of the original data or proof of its origin."
        captions["t07_plan"] = "<b>The constructed scene and its power map.</b> Eight added point scatterers are circled in the input. Dominant points produce dark patches; mixed face edges produce bright bands. The dashed line marks a transect. Locations and display settings were chosen for this illustration."
        captions["t07_slices"] = "<b>Sections through the synthetic volume.</b> Bright bands and isolated point scatterers produce different column-like patterns. The left axes use a chosen wavelength of 2.2 m; the right axes retain 0.48 m. All panels stop at the computed depth limit. This rescaling does not establish the calibration used for the published slides."
        captions["khafre_cad_ramps"] = "<b>Published interpretation.</b> A CAD illustration of cylinders and spiral ramps beneath a wireframe pyramid. The illustration is distinct from a measured radar image. Reproduced at reduced size for critique; copyright the Khafre Research Project."
        captions["khafre_cad_coils"] = "<b>Another published CAD view.</b> Cylinders, coils, and basal blocks depict the project's interpretation. Their shapes and depths are not independently validated here. Reproduced at reduced size for critique; copyright the Khafre Research Project."
        captions["t07_isosurface"] = "<b>A rendering of the synthetic volume.</b> Low-power regions (blue) and high-power regions (orange) appear as columns or walls after smoothing along depth. Their common lower boundary is the end of this computed volume. Isosurface shape depends on the chosen threshold and smoothing."
        captions["t07_rungs"] = "<b>Spectra in selected regions of the constructed scene.</b> The comb fits differ across regions and do not cleanly follow one course spacing. This comparison does not measure or explain a published spiral structure."
    for figure in summary.get("figures", []):
        stem = Path(figure["file"]).stem
        if stem in captions:
            figure["caption"] = captions[stem]
    return summary
