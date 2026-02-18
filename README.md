## **Update After First-Round Peer Review**

The revised manuscript includes additional modeling components for **Topographic Wetness Index (TWI)** and **Soil Water Repellency (SWR)**. As a result, the main simulation script has been updated.

Please use **`Main_SF_para_prob_eff_revisedpaper.m`** instead of **`Main_SF_para_prob_eff.m`** to reproduce the results reported in the revised manuscript.

---

## **First Submission**

Code for *Probabilistic Post-Fire Shallow Landslide Susceptibility Modeling Considering Spatiotemporal Land Cover Uncertainties: A Case of the January 2025 Palisades Wildfire in Southern California*.

Probabilistic, physics-based modeling of post-fire shallow-landslide susceptibility and its uncertainty.

### Core scripts
- **`RFHydrorealization.py`** — generates random fields of post-fire hydraulic conductivity  
- **`RFroot.py`** — generates random fields of post-fire root cohesion  
- **`Main_SF_para_prob_eff.m`** — MATLAB model for shallow-landslide susceptibility (factor of safety and failure probability)

### Inputs
- Place the required geospatial inputs in **`RasterT_Palisad4_SpatialJoin6_TableToExcel.xlsx`** (derived from ArcGIS Pro).

### Quick start
1. Generate hydraulic conductivity random fields  
   - Run **`RFHydrorealization.py`** to produce ensembles of post-fire hydraulic conductivity (or multipliers) over the study grid.
2. Generate root cohesion random fields  
   - Run **`RFroot.py`** to produce ensembles of post-fire root cohesion (or multipliers) that evolve over time.
3. Run the physical model in MATLAB  
   - Open **`Main_SF_para_prob_eff.m`** and set the required paths.

### Outputs
- FS maps for each time step and ensemble member

