"""Catalogue of the published figures: script, data tier and the plotted-value tables each script writes.

Tier "repo": drawn from the de-identified repository data (data/derived).
Tier "restricted": needs the restricted project data (Beijing-3 products, station polygons, unit coordinates or the
Landsat annual composites); the code is included and was verified on the project data.
"""
FIGURES = {
    'Fig1': ('fig01_study_area.py', 'restricted', 'study area, unit and cross-section', ['fig1_counts.json']),
    'Fig2': ('fig02_workflow.py', 'repo', 'workflow of the study', []),
    'Fig3': ('fig03_event_study.py', 'repo', 'mixed-pixel event study', []),
    'Fig4': ('fig04_time_axis.py', 'repo', 'one time axis, two records, three surfaces', ['fig3_cohorts.csv', 'fig3_lanes.csv', 'fig3_longitudinal.csv']),
    'Fig5': ('fig05_baselines.py', 'restricted', 'matched and unmatched baselines', ['fig4_summary.json', 'fig4e_designs.csv']),
    'Fig6': ('fig06_spatial_pattern.py', 'restricted', 'spatial pattern and station means', ['fig5_summary.json']),
    'Fig7': ('fig07_spectral_evidence.py', 'repo', 'spectral, red-edge, distance and clearance evidence',
             ['fig6a_relative_reflectance.csv', 'fig6b_red_edge.csv', 'fig6c_distance_08.csv', 'fig6d_clearance.csv']),
    'Fig8': ('fig08_imaging_share.py', 'repo', 'share reproducible by imaging alone', ['fig7_scenario_shares.csv', 'fig7_summary.csv', 'fig7_clearance.csv']),
    'FigS1': ('figS01_calibration.py', 'restricted', 'radiometric traceability and calibration', ['figS1fg_scene_ndvi_metrics.csv', 'figS1h_band_cv_rmse.csv']),
    'FigS2': ('figS02_unmixing_tests.py', 'repo', 'panel end-member and unmixing tests', ['figS6c_fe_offset.csv', 'figS6_bias_rmse_150m.csv']),
    'FigS3': ('figS03_installation_year.py', 'restricted', 'installation-year detection', []),
    'FigS4': ('figS04_unit_distributions.py', 'repo', 'unit-level distributions', []),
    'FigS5': ('figS05_common_footprints.py', 'repo', 'common-footprint changes', ['supp_diag/figS8_common_by_clearance.csv']),
    'FigS6': ('figS06_mixed_minus_interrow.py', 'repo', 'mixed minus inter-row contrast by row-band fraction', ['supp_diag/figS23_mixed_minus_interrow.csv']),
    'FigS7': ('figS07_single_year_cohorts.py', 'repo', 'single-year cohorts', ['fig3_cohorts.csv', 'fig3_lanes.csv', 'fig3_longitudinal.csv']),
    'FigS8': ('figS08_coverage.py', 'restricted', 'coverage decomposition', []),
    'FigS9': ('figS09_spectral_contrasts.py', 'repo', 'visible and NIR spectral contrasts', ['supp_diag/figS21_spectral.csv']),
    'FigS10': ('figS10_distance_bands.py', 'repo', 'distance bands', ['supp_diag/figS14_distance.csv']),
    'FigS11': ('figS11_imaging_boundary.py', 'repo', 'imaging and boundary sensitivity', ['supp_diag/figS22_optics_boundary.csv']),
    'FigS12': ('figS12_rainfall.py', 'repo', 'mixed-pixel rainfall diagnostic', []),
}
# Fig. 6 and Fig. S8 write the summaries that the numbers read (fig5_summary.json, figS27_coverage.csv); in the
# project layout they run before scripts/31_numbers.py, and Fig. 2 prints two numbers of numbers.json.
ORDER = ['Fig6', 'FigS8'] + [n for n in FIGURES if n not in ('Fig6', 'FigS8', 'Fig2')] + ['Fig2']
