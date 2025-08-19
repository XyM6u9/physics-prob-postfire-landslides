Codes for Probabilistic Post-Fire Shallow Landslide Susceptibility Modeling Considering Spatiotemporal Land Cover Uncertainties: A Case of January 2025 Palisades Wildfire in Southern California

Code and workflows for probabilistic, physics-based modeling of post-fire shallow-landslide susceptibility and its uncertainty.
Core pieces:
  RFHydrorealization.py — generates random fields of post-fire hydraulic conductivity
  RFroot.py — generates random fields of post-fire root cohesion
  Main_SF_para_prob_eff.m — MATLAB model for shallow-landslide susceptibility (factor of safety and failure probability)

Inputs:
  Place required geospatial inputs in RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx derive from ArcGIS Pro

Quick start:
  1) Generate hydraulic conductivity random fields
    RFHydrorealization.py produces ensembles of post-fire hydraulic conductivity (or multipliers) over the study grid.
  2) Generate root cohesion random fields
    RFroot.py produces ensembles of post-fire root cohesion (or multipliers) that evolve over time.
3) Run the physical model in MATLAB
  Open Main_SF_para_prob_eff.m, and set the paths.

Outputs:
  Fs maps per time step and ensemble member
