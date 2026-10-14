# IODP U1418 public data: attribution and extraction

Zindorf, M., März, C., and Schnetger, B. (2020). “Data report:
wavelength-dispersive X-ray fluorescence–based geochemical data, Site U1418,
IODP Expedition 341, Gulf of Alaska”. DOI: https://doi.org/10.2204/iodp.proc.341.206.2020
Official report: https://publications.iodp.org/proceedings/341/206/341_206.html
Official workbook: https://publications.iodp.org/proceedings/341/SUPP_MAT/206/XRF/U1418%20XRF%20Data.xlsx

The authors' linked XRF data are redistributed under **CC BY 4.0**:
https://creativecommons.org/licenses/by/4.0/ . Attribution does not imply
endorsement. This data licence is separate from the package's GPL-3 code.

The bundled iodp_u1418_major_trace.csv is an unchanged copy of the independently
cleared existing extraction, not a fresh supplement download. SHA-256:
06b31a651aba7df1dc8b3c3f379018c4b4793fe778ce4920a4d039831eb062a3
Original MD5: 14f5a0bdb8e19c18165986fb88ea9271
LF-normalized MD5: 8817ffe358c53b939c81e4e4ca00718a

Extraction selected workbook sheet 1, A:E metadata, F:O ten major oxides,
Q:AF sixteen trace elements and AG reported Sum Of Conc. P and AJ:AS derived
elemental wt% columns were omitted. Nonempty Label identifier rows were retained.
A:E became depth_m, sample_type, sample_length_cm, sample_id, text_id; AG became
oxide_total_wt_pct. Numerical cells and blanks were preserved. There are 32
fields and 284 identified rows/unique sample IDs, whereas the article reports
285 analyzed samples. This count difference is unresolved.

Major oxides are dry-sediment wt%; trace elements are separate ppm; depth_m is
CCSF-A depth. Trace zeros/blanks and reported totals remain unchanged. The example
closes only SiO2, TiO2, Al2O3, Fe2O3, MnO, MgO, CaO, Na2O, K2O, P2O5 to a selected
subcomposition. It never changes the CSV, adds a remainder, imputes censoring or
claims whole-rock mass closure. All 2,840 selected oxide values are positive.

Manual refresh recipe (not run by package/example/check): download the official
workbook and first verify historical SHA-256
f862fcec7b3772519eb4124ad051a19a9a49a3390026f7c84e48f2ab623eaa0d.
If it differs, stop and review the changed source rather than accepting a refresh.
Select/rename the sheet columns and identified rows as above. Compare ordered
schema, sample IDs and every numerical/blank cell with this pinned CSV. CSV
serialization may change byte identity even when cells agree; report both cell
comparison and new byte hash if reserializing. The workbook hash and previous
cell comparison are historical provenance evidence: a fresh official workbook
retrieval/comparison was unavailable at clearance, and equivalence has not been
newly demonstrated. No automatic network retrieval is included.
