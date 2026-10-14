# Mineralogy–assay reference catalogue

Metadata only; evidence assessed 14 October 2026. These six existing extractions
are references for measurement and consistency checks, **not fitted calibration
populations**. Long-table row counts count measurements, not specimens. No raw
six-source data, source documents, specimen keys or extraction scripts are bundled.
A shared material key does not establish a common analytical mass or spatial basis.

| Citation / source | Extract rows / material units | Assay and modal basis; matching | Rights status and scientific limits |
| --- | --- | --- | --- |
| Austin (2024), [Cloncurry METAL dataset](https://doi.org/10.60623/82trleue), Cannington subset; [methods paper](https://doi.org/10.5194/essd-16-5027-2024) | 32,436 / 153 exact specimen matches | TIMA **vol%** versus polished-surface pXRF elemental **ppm**, separate 1SE; exact specimen-key join | Dataset CC BY 4.0 verified at [AuScope](https://repository.data.auscope.org.au/dataset/cloncurry-metal-multimodal-integrated-mineral-system-characterisation-data). Holes/Unclassified retained; spot chemistry is not bulk mass chemistry. Pb **143 numeric + 10 nondetects**; Zn153 numeric; Ag11 numeric +142 nondetects. |
| Tibane, Mamba and Nkosi (2026), [Zandrivierspoort Formation BIF](https://doi.org/10.1016/j.dib.2026.112831) | 440 / 20 sample pairs | Reported QXRD mineral **wt%**, XRF oxide/total FeO **wt%**; sample-ID join | Article CC BY 4.0 verified in [primary deposited article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13213779/). Source-recorded ideal Fe lower bounds exceed assay in all20 pairs; anomalous totals near199%, negative LOI and nonclosing XRD. No fresh forward calculation or supplement-byte verification claimed. |
| Belachew, Belachew and Tamene (2026), [Ashashire gold deposit](https://doi.org/10.11648/j.ajche.20261404.15) | 156 / 6 representative cores | Au fire assay/AAS **g/t**, mixed elemental assays; QEMSCAN labelled **wt%**, image-mode-to-mass conversion undocumented; core-key join | Article CC BY 4.0 verified at [publisher](https://www.sciencepg.com/article/10.11648/j.ajche.20261404.15). Triplicates are not independent cores; no quantified per-core gold-bearing mineral masses, censored/trace modes and one subtotal above100%. |
| [Moholwa manganese data (2022)](https://doi.org/10.1016/j.dib.2022.108566) | 40 / 3 ore-type aggregates | Average ICP-OES mixed element/oxide **wt%**, natural-ore QXRD **wt%**; aggregate ore-type matching | CC BY 4.0 is provenance-recorded, not newly primary-verified. Partial sums, absent/censored phases and missing phase chemistry; three averages cannot be treated as independent specimen populations. |
| Zhang et al. (2026), [Multi-Scale Process Mineralogy of Cd and Ag in a Pb-Zn Ore](https://doi.org/10.3390/min16030316) | 29 / 1 homogenized bulk batch | 19 chemical results **wt%/ppm**, 10 reported TIMA mineral **mass%** (96.12% subtotal); subsamples from one bulk material | Article-specific CC BY 4.0 claim and full measurement description remain unresolved after inaccessible publisher full text/PDF. Source-recorded ideal listed sulfides exceed assay. Do not normalize away the contradiction or count29 measurements as specimens. |
| Yang, Xiao, Yang and Zhang (2016), [Processing Mineralogy Study on Lead and Zinc Oxide Ore in Sichuan](https://doi.org/10.3390/met6040093) | 27 / 1 bulk composite | 12 chemistry results **wt%**, Ag **g/t**, original Total_Fe label; 15 microscope/MLA reported mineral **weight%**; common raw-ore material | Article-specific rights remain unresolved. Elemental phase allocations excluded. Missing supported phase chemistry yields partial lower bounds; no demonstrated Zhang-style incompatibility. |

Nondetects are qualifiers, not zero or structural absence. For Cannington the
153 Pb records all match numeric Galena modes, but only143 assays are numeric:
the earlier all153-numeric provenance statement was incorrect. No failed join
or changed selection is required to explain the difference.

## Extraction identities and QA

Source-relative reference filenames below identify existing extraction evidence
at recovered source commit `989577aaa449d28ac554d806ef57daa47122069c`.
Hashes are independently inspected CSV byte identities; they do **not** prove
scientific validity, fresh workbook/PDF extraction or CSV/Parquet parity.
Each extract has a corresponding source provenance record and extraction script;
neither is copied into this package. Units and qualifiers must be preserved in
any independently cleared reuse, without normalization of inconsistent totals.

| Extraction reference | SHA-256 |
| --- | --- |
| cloncurry_metal_2024_cannington.csv | `d0b26aded9d71715a2786e6a04dbc000f1ddf887f427876a4cda9fd318f55db3` |
| zandrivierspoort_bif.csv | `51c010fb07576f6fc5dcc8bef681ee051b9a51e77dc0cb34d52a705f46a3d671` |
| ashashire_gold.csv | `586f34c5434469ab046e0e85c5b0341e1b338006a977b80b09d437e060c53dea` |
| moholwa_2022.csv | `88e553b1967bd6f3c03642fafc9ce457b1c6c7bbf4c31b7fa95a9388aac36c72` |
| zhang_2026_pbzn.csv | `e4bbcd7ee3ed239fcdfc62fa878a66e86f8a4e2601603b213103c843f87bcd3a` |
| sichuan_2016_pbzn.csv | `eefc7d28cf9d059b51e822a92faaaf7ec2905e5eebc7f846f1bc72c4c0dc37eb` |

Registry `source_registry.csv` SHA-256:
`c77ce7b376ccd0b6b9ee22f253450d43045a1884d5b199118cf993f3ba6e6c1b`.
Cannington extraction reference `scripts/extract_cloncurry_metal.py` SHA-256:
`f0162c186da21dd8563f3a31b116b74093345a4149fbca6e09e1d93a1f0906c0`.
Recorded pair-check reference `analysis/extracted_pair_checks.R` SHA-256:
`83b9b10d2b4405c372eedf3844cef6c85f0119443be897bfbc7c522e4848a602`.
Those scripts were inspected, not executed in this catalogue audit.

## Additional source triage

| Lead | Recorded access/rights | Suitable interpretation |
| --- | --- | --- |
| [Rocklea Dome C3DMM](https://doi.org/10.25919/5ed83bf55be6a), [methods](https://doi.org/10.5194/essd-13-1371-2021) | Dataset CC BY 4.0 recorded; no new licence audit here | Public RC XRF/LOI and hyperspectral indices. Fe oxide/hematite–goethite indices are spectral products, not quantitative mineral mass%. Earlier XRD discussed as validation does not supply the public numeric join. |
| [Nkout iron ore](https://doi.org/10.1016/j.oregeorev.2014.02.015) | Complete numeric pairing inaccessible; separately licensed data unresolved | Reports51 samples but accessible summaries do not establish a complete specimen-level chemistry/modal join. No extraction or redistribution clearance. |
| [Peterson ferromanganese ore](https://doi.org/10.1080/25726641.2023.2222909) | Article CC BY-NC-ND 4.0 recorded | Abstract assemblage/separation evidence without verified open numeric modal join; no quantitative extraction admitted. |
| USGS iron/mafic assays; Amaruq; Wiwirowno laterite and other source-recorded leads | Chemistry-only or unaudited leads; no new rights/pairing clearance | Context or potential statistical comparisons only. Do not promote to paired mineral calibration or validation populations. |

Original `SOURCE_TRIAGE.md` SHA-256:
`d619e5aa82fb4455f6e66afa8dbf9c89de2a46a96d6603169303f2ea2d3750e6`.
This bounded catalogue does not expand that audit. **Only the separately attributed
IODP CSV is bundled numerical data.** Its dry-sediment oxide/trace basis and
extraction limitations are in [IODP-LICENSE-PROVENANCE.md](IODP-LICENSE-PROVENANCE.md).
