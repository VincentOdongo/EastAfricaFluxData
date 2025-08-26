# EA Dryland Carbon & Water

This directory holds everything you need to reproduce the analyses in **“Contrasting carbon and water dynamics in two East African dryland ecosystems”**.

## Subfolders

### code/  
All MATLAB scripts for processing, analysis, and plotting:
- `KE_KAP_flux_rangeland.m`  
  • Reads and gap‐fills Kapiti (rangeland) flux data  
- `KE_AUS_flux_cropland.m`  
  • Reads and gap‐fills Ausquest (cropland) flux data  
- `KE_KAP_rangeland_biomet.m` / `KE_AUS_cropland_biomet.m`  
  • Imports biomet sensor files (soil moisture, PPFD, etc.)  
- `wue_analysis.m`  
  • Calculates and bins WUE vs. VPD, generates Figure 8  
- `cue_comparison.m`  
  • Calculates CUE barplots and t-tests (Figures 6 & 7)  
- (etc. – list any other custom .m files you have)

### data/  
Interim raw files used in these scripts:
- `KE_KAP_gapfilledNEE_LE.txt`  
- `KE_AUS_gapfilledNEE_LE.txt`  
- `KE_AUS_cropland_biomet.csv`  
- `KE_KAP_rangeland_biomet.xlsx`

> **Note:** final FLUXNET releases (with DOIs) will be linked from the top‐level README once available.

---

**Usage**  
1. Run any script in `code/` from MATLAB (scripts assume their working directory is this folder).  
2. All file paths point to `../data/…`; no further edits should be needed.  
3. Figures will pop up in MATLAB’s figure window; you can export using `export_fig` or the built-in “Save As.”

---

*Last updated: July 2025 by Vincent Odongo*  
