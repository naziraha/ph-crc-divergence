# Commit map: project history → public repository

The analysis documents (plan, addenda, deviation log) and script headers cite commits of the private project repository, in which this repository was developed. The public history was produced from it with `git filter-repo`, keeping only the files of this repository: commit dates, authors and messages are unchanged, commits that touched only excluded files (manuscript drafts, internal notes) were dropped, and every remaining commit has a new hash. This table maps each original hash to its public equivalent ("—" = commit dropped because it contained only excluded files).

| Original | Public | Date (UTC+8) | Message |
|---|---|---|---|
| `4f1f32f` | `66d1c4f` | 2026-10-04 11:35 | Pre-specification locked before step 4 |
| `913635a` | `994e95c` | 2026-10-04 11:37 | Step 4: persistent homology profiles per stage (descriptive) |
| `06db10a` | `a149d04` | 2026-10-04 11:54 | Step 5: primary permutation test (999 permutations) |
| `0fbc7bf` | `82a369a` | 2026-10-04 12:47 | Step 5: sensitivity analyses S1-S9 |
| `6ea31ea` | `a3bb307` | 2026-10-04 13:54 | Pre-specification addendum: GSE20916 and GSE37364, locked before download |
| `aae528a` | `7cf500a` | 2026-10-04 13:56 | Validation: GSE20916 step 1 data preparation |
| `d73e5d1` | `179a53f` | 2026-10-04 14:02 | Addendum 2 (stage-aware adjustment, GSE20916 micro rules) and deviations log, before step 2 on validation cohorts |
| `d244f36` | `b6c4056` | 2026-10-04 14:05 | Validation: GSE20916 step 2 identity screen and composition (primary + stage-aware) |
| `30d0dab` | `de4bd2c` | 2026-10-04 14:09 | Validation: GSE20916 step 3 embeddings and micro pooling decision |
| `ac2a591` | `d8bba05` | 2026-10-04 14:12 | Deviation D4: GSE20916 microdissected step 4 (all-sample diagrams, crypt spread), before step 4 |
| `e817a3c` | `43af563` | 2026-10-04 14:26 | Validation: GSE20916 steps 4-5 (bulk test, micro descriptive), D5 |
| `d38db92` | `b436fcf` | 2026-10-04 18:12 | Addendum 3 (secondary DL meta-analysis); GSE37364 steps 1-2 |
| `3bcfc5d` | `e5d7df2` | 2026-10-04 19:12 | GSE37364 steps 3-5; exploratory D3 (GSE41258 stage-aware); secondary meta-analysis |
| `a74a5cd` | `—` | 2026-10-04 19:23 | Results summary steps 1-5 for Prof Peng |
| `26e16a4` | `—` | 2026-10-04 19:47 | Word version of results summary for Prof Peng (title page, embedded figures, Word tables) |
| `1611030` | `a2eda43` | 2026-10-06 10:28 | Deviation D6: Ripser cross-check of GSE41258 step 4 diagrams (cross-check only) |
| `d456d5c` | `831def6` | 2026-10-06 10:29 | D6: GUDHI vs Ripser cross-check of GSE41258 step 4 diagrams - agreement to float precision |
| `15cbf68` | `d1d5b5a` | 2026-10-06 10:37 | Addendum 4: standardised effect sizes with bootstrap CIs and pooling (descriptive), before computation |
| `bddafc0` | `da56bc5` | 2026-10-06 10:57 | Addendum 4 outputs kept, marked invalid (bootstrap duplicates); point estimates descriptive only; negative I2 shown as 0 |
| `a446a34` | `12f6911` | 2026-10-06 10:58 | Deviation D7: empirical power by relabelling dilution of the raw Normal->Adenoma benchmark; classification rule (before computation) |
| `72e4d85` | `061b2c5` | 2026-10-06 11:57 | D7 power results: raw interpretable; shape-only withdrawn (relabelling mixture artefact) |
| `b98407b` | `1c51972` | 2026-10-06 11:57 | D7 status: raw power interpretable, shape-only withdrawn |
| `9d814f9` | `d989051` | 2026-10-06 13:11 | Deviation D8: sample-size sensitivity of the shape-only test (before computation); GEO metadata fetch script |
| `0761c37` | `b7d8f1d` | 2026-10-06 13:18 | D8 results: shape-only detection vs sample size |
| `a96b652` | `—` | 2026-10-06 13:51 | Revised summary sent to Prof Peng |
| `ed90115` | `—` | 2026-10-06 14:09 | Ripser cross-check report for Prof Peng (docx + pdf) |
| `713e6ca` | `—` | 2026-10-06 14:33 | Technical report v1 for Prof Peng (TDA review): data, exact settings, pre-registration, deviations, systematic-error table, proposed positive control |
| `d5673a7` | `6a540d9` | 2026-10-06 14:41 | Addendum 5: validation re-run, scan-date batch check, positive/negative controls with acceptance criteria, H0-only/H1-only tests (before running) |
| `32cc053` | `9172e21` | 2026-10-06 14:47 | Addendum 5 part 1: reference-code validation re-run (report + console log) |
| `11aadde` | `443e515` | 2026-10-06 15:31 | Addendum 5 part 2 as-specified run (Excel date swap discovered) + deviation D9 parsing correction, before re-run |
| `1ec4d3b` | `67eaaaa` | 2026-10-06 15:32 | Addendum 5 part 2: scan dates with D9 parsing correction (Normal->Adenoma confounded with scan day) |
| `0a3a86a` | `6e7bf8e` | 2026-10-06 15:35 | Addendum 5 part 2b appended (batch-only normals comparison + size-matched N->A reference), before running |
| `5cae936` | `9abc2cf` | 2026-10-06 15:47 | Addendum 5 part 2b results (batch-only normals vs size-matched N->A); part 3 script |
| `fe9dc84` | `016272a` | 2026-10-06 18:36 | Addendum 5 part 3: positive/negative control results |
| `5c5791a` | `90c24a9` | 2026-10-06 19:33 | Addendum 5 part 3b appended: level-2 diagnostic with bug check (before running) |
| `1d684b4` | `d1dca66` | 2026-10-06 19:35 | Addendum 5 part 3b: level-2 diagnostic - insertion present; lost mainly by 10-PC truncation and gene selection |
| `d09a883` | `fa68746` | 2026-10-06 20:06 | Addendum 5 part 4: H0-only and H1-only tests (secondary) |
| `c35db60` | `—` | 2026-10-06 20:42 | Technical report v2 for Prof Peng: executive summary, addendum 5 checks, updated interpretation, next steps (v1 unchanged) |
| `0bf9847` | `—` | 2026-10-06 20:49 | Technical report v2: executive summary revised per PI (wording); summary kept to one page (10-pt formatting only) |
| `e6e7dd2` | `—` | 2026-10-06 20:57 | Technical report v2: real Word bullets (converter --fix-lists), summary kept on page 2 (10 pt, 1-pt paragraph spacing); v1 not rebuilt |
| `93613bd` | `6765159` | 2026-10-06 21:10 | Addendum 6: GSE117606 as fourth cohort (before any expression data downloaded or opened) |
| `abde6eb` | `76a5d2a` | 2026-10-06 21:21 | GSE117606: CEL header table (no intensities; archive kept outside OneDrive) and addendum 6 section 7.1 batch check |
| `8c9b93f` | `5d1b5df` | 2026-10-06 21:27 | GSE117606 step 1 (addendum 6): series matrix, Entrez->symbol mapping, main set 65/59/71, array-quality flags |
| `11474e8` | `f82c571` | 2026-10-06 21:30 | Deviation D10: GSE117606 step 3 will report PC correlations with array quality (descriptive, recorded before step 3) |
| `845dcfe` | `87e6ae2` | 2026-10-06 21:32 | GSE117606 step 2 (addendum 6): identity screen (15 flags, none excluded), composition scores, stage-blind and stage-aware adjustment |
| `a155f15` | `b646e7d` | 2026-10-06 21:35 | GSE117606: PI identity decision (keep all 195) recorded; deviation D11 (sensitivity without GSM3304812), before step 3 |
| `5002344` | `9bef682` | 2026-10-06 21:35 | GSE117606 step 2 summary: PI identity decision recorded (keep all 195) |
| `b734dde` | `558bb62` | 2026-10-06 21:38 | GSE117606 step 3 (addendum 6): embeddings (unadjusted, adjusted, stage-aware), D10 quality correlations, SSA/set-aside projection, sensitivity embeddings |
| `60e161e` | `1e0f4ab` | 2026-10-06 21:42 | Deviation D12: GSE117606 quality-covariate adjustment sensitivity (declared before step 4) |
| `120a7ba` | `310ed72` | 2026-10-06 21:43 | GSE117606 step 4 (addendum 6): descriptive persistence profiles, unadjusted and adjusted |
| `b78121b` | `c854097` | 2026-10-06 22:46 | GSE117606 step 5 primary (addendum 6): unadjusted + adjusted, raw + shape-only, 999 permutations, Holm |
| `8a8bd7f` | `da4199a` | 2026-10-07 05:11 | GSE117606 step 5 secondary analyses (stage-aware, paired, noQAflag, D11, D12, PC/gene variants) and 4-cohort meta-analyses (addendum 6) |
| `c165512` | `cf4dd3d` | 2026-10-07 08:19 | Addendum 7: Step 13 benchmark pre-registration (PERMDISP, PERMANOVA, stage separation, MST/H0) - before any computation |
| `6fac06a` | `9b1c923` | 2026-10-07 11:18 | Step 13: addendum 7 section 4 simulation validation - all acceptance criteria passed (no real data) |
| `2f8cf74` | `a6e8d80` | 2026-10-07 11:31 | Step 13 benchmark (addendum 7): PERMDISP, PERMANOVA, stage separation, MST/H0, size-matched; 4 cohorts |
| `8e8df5b` | `99fed42` | 2026-10-08 09:06 | Addendum 8: gradual vs switch-like change pre-registration (cross-fitted gene trajectories) - before any computation |
| `248d295` | `edb2610` | 2026-10-08 09:08 | D13: addendum 8 sigmoid-fitting implementation (profile grid + refinement), declared before validation |
| `cca77b1` | `f2d502b` | 2026-10-08 11:45 | Addendum 8 section 6 validation: 13/14 criteria passed; sigmoid power 0.62 < 0.80 (D14) - real data on hold |
| `0275b7a` | `87a9ccf` | 2026-10-08 11:48 | Addendum 8 Amendment 1: fixed-window resolvability rule (D14), committed before validation rerun |
| `1f15bf9` | `c494b28` | 2026-10-08 13:52 | Addendum 8 validation rerun (Amendment 1): 13/14 passed; gap-scenario unresolvable 0.793 < 0.80 (D15) - real data on hold, PI to decide |
| `a3af4f9` | `03bed32` | 2026-10-08 14:57 | Addendum 8 Amendment 2: conditional interpretation of switch-gene calls (D15), before any real data |
| `e73638d` | `024d9f0` | 2026-10-08 15:40 | Addendum 8 real-data analysis (4 cohorts): N->A support for gradual change; A->C inconclusive (gene level gradual, population precision limited); no replicated switch or non-monotonic genes |
| `5a9b03a` | `62eba6d` | 2026-10-08 15:58 | Addendum 9: carcinoma heterogeneity (variability, CMS, dispersion decomposition) pre-registration - before any computation |
| `9590f65` | `c3daff2` | 2026-10-08 16:06 | D16: addendum 9 CMScaller implementation details (explicit centring; example-data validation definition), before validation |
| `4f1ce27` | `8692b68` | 2026-10-08 16:09 | Addendum 9 section 6 validation: 12/17 passed; V1 power, LOO over-correction, CMScaller 0.889 (D17) - real data on hold |
| `12dfc3c` | `315fcfb` | 2026-10-09 09:33 | Addendum 9 Amendment 1: bias-adjusted PERMDISP distances, V1 edge 4 SD, CMScaller accepted at 0.889 (D17) - before rerun |
| `3844e5e` | `def4248` | 2026-10-09 09:35 | Addendum 9 section 6 validation rerun (Amendment 1): all 17 criteria passed (no real data) |
| `699796d` | `83eb288` | 2026-10-09 09:45 | Addendum 9 real data: primary conclusion inconclusive (2 within-subtype, 1 subtype diversity; robust to confident-only); 117 replicated more-variable genes, CMS-template enriched |
| `547017a` | `—` | 2026-10-10 12:53 | Add index of real-data outputs |
| `bca0b9e` | `—` | 2026-10-10 14:09 | Plain-language project status summary for Prof Peng |
| `296eaed` | `c55c232` | 2026-10-10 14:25 | Addendum 10: TCGA-COAD carcinoma stages I-IV with the frozen pipeline - pre-registration before any TCGA computation |
| `438d686` | `24fffe6` | 2026-10-10 14:30 | D18: TCGA-COAD standard (non-centred) Xena HiSeqV2 used instead of PANCAN file; 25% filter in step 3 as for GEO |
| `35f3517` | `b821bec` | 2026-10-10 14:31 | Addendum 10 step 1: TCGA-COAD preparation (274 tumours: I 45, II 110, III 80, IV 39; 41 normals) |
| `c825eff` | `fa1cd2a` | 2026-10-10 14:35 | Addendum 10 step 2: TCGA-COAD composition scores and adjustment; identity screen flags 46/274 tumours (decision pending) |
| `cdf3bcf` | `326c337` | 2026-10-10 14:38 | D19: TCGA-COAD keep all 274 tumours; sensitivity without 18 low-colon-identity tumours (before step 3) |
| `cd2cf7a` | `410ff91` | 2026-10-10 14:40 | Addendum 10 step 3: TCGA-COAD embeddings (main, noLowColon, MSS, withnormal; sensitivity PCs/genes) |
| `ec07b9d` | `2c5ec7c` | 2026-10-10 14:43 | Addendum 10 step 4: TCGA-COAD persistence profiles per stage (descriptive) |
| `92f9de4` | `7c3d076` | 2026-10-10 14:45 | D20: TCGA positive-control conditions and acceptance criteria, declared before step 5 |
| `250314e` | `6c2baab` | 2026-10-10 15:37 | Addendum 10 step 5: TCGA-COAD permutation test - no transition significant (any analysis); controls 5/6 pass (SP shape-only 4/20 > 3/20) |
| `e9789ba` | `ec91fd5` | 2026-10-10 16:21 | D21: TCGA control outcome accepted with limitation; earlier TCGA pilot superseded (pseudo-replication) |
| `eb33a2f` | `5772926` | 2026-10-10 16:30 | Addendum 10 section 4.1: TCGA spread (PERMDISP 0/6 sig) and H0/H1-only (1/24 sig: H1-only adjusted shape III->IV) |
| `a3e5aeb` | `622c7f4` | 2026-10-10 16:37 | Addendum 10 section 4.2: TCGA CMS by stage (coverage 95.5%; CMS4 rises, CMS3 falls with stage) |
| `99017b1` | `9db44c4` | 2026-10-10 16:39 | Addendum 10 CMS summary: caution that CMS4 calls on unadjusted expression partly reflect stromal content |
| `69d178e` | `73ec95b` | 2026-10-10 16:43 | D22: TCGA normal -> Stage I main test (one per patient) and paired normal vs own tumour (any stage) replacing infeasible paired Stage I analysis |
| `207b048` | `2fa521c` | 2026-10-10 16:47 | Addendum 10 section 4.3: TCGA normal vs Stage I (shape and raw significant, both versions) and paired normal vs own tumour |
| `fc9a8d2` | `05b3bfa` | 2026-10-10 16:51 | Addendum 10 section 4.4: pilot method p depends on subsample count (0.95 -> 0.001); patient-level p = 0.885 raw / 0.941 shape; D21 updated |
| `37a269a` | `—` | 2026-10-10 16:55 | Objectives vs findings summary (what was found, not found, and why) |
| `3cf85e9` | `ad30480` | 2026-10-10 17:03 | Addendum 9: descriptive annotation of replicated variable genes (CMS template class, median mean shift) |
| `baf4bc4` | `1735880` | 2026-10-10 17:14 | Addendum 11: chromosome 20 and differentiation axis pre-registration (A5 rule i: top 5 arms) - before any computation |
| `ca18e5b` | `8fe3ef7` | 2026-10-10 17:16 | D23: chromosome-arm annotation from MSigDB C1 positional sets (platform files lack locations); GISTIC2 download checksum, 272/274 matched |
| `384c81e` | `fbc24ef` | 2026-10-10 17:18 | Addendum 11 Amendment 1 (Part C residual composition; arm check of 8 genes) before running; D24 simulation W3 0.090 > 0.089, Part B on hold |
| `6523fd8` | `0e24bd4` | 2026-10-10 17:21 | Addendum 11 Part A: 20q CONFIRMED (rank top 5 4/4, BF 4/4, TCGA rho 0.85); comparators 20p/13q/8q confirmed; Part C: no residual composition signal |
| `1051818` | `62fecf1` | 2026-10-10 18:08 | D24 accepted with restriction (D_resid p descriptive); addendum 11 Amendment 2 (chromosome-arm f) before computing; Part A wording |
| `9373c3f` | `aa54dc1` | 2026-10-10 18:13 | Addendum 11 Part B + Amendment 2: differentiation axis not confirmed (1/4; full sets 3/4); all dispersion shares minor; f_CMS reproduces addendum 9 |
| `d2ed028` | `—` | 2026-10-10 18:20 | objectives_vs_findings: updated wording after addendum 11 (subtype, chromosome arms, differentiation hypothesis-generating) and limitations |
| `693974d` | `—` | 2026-10-10 18:21 | One-page plain-language summary of findings (addenda 1-11) for Prof Peng |
| `d0e85ab` | `—` | 2026-10-10 18:26 | Revised paper outline (draft 1): biology-focused framing, TDA as pre-registered negative finding and benchmark |
| `bef3266` | `—` | 2026-10-10 21:45 | Novelty literature review (PubMed): TDA transitions, DNB in CRC, TDA in CRC omics, multi-cohort negative studies; revised novelty |
| `32b80eb` | `—` | 2026-10-10 21:52 | Literature review round 2: broader search (Consensus/Semantic Scholar/OpenAlex, arXiv, web); 7 new verified records; novelty conclusion and qualified wording |
| `73b4e71` | `2d5d4ca` | 2026-10-10 22:31 | Addendum 12: quantitative PH outputs (landscapes, barcode summaries, dimension split) and standard figures - pre-registration before any computation |
| `4018ea9` | `4350e26` | 2026-10-10 22:34 | D25: addendum 12 validation criteria interpretation and implementation details, before running |
| `5751d5f` | `5af82c2` | 2026-10-10 22:40 | Addendum 12 validation: 14/16 passed; H0 landscape underpowered (SP raw 10/20, CL 12/20) (D26) - real data on hold |
| `d5fe9e6` | `08a6011` | 2026-10-10 22:44 | Addendum 12 Amendment 1: landscape tests H1 only, H0 landscapes figures only (D26), before real data |
| `c531104` | `e6fae01` | 2026-10-10 23:15 | Addendum 12 real data: H1 landscape concordance 83% (supportive); 23/336 barcode summaries significant (mostly H0); dimension split all cohorts (H0 reproduces 94%; one isolated H1-only); loops rule not met; 60 figures |
| `cb6859b` | `—` | 2026-10-10 23:20 | Results outline v1 (hybrid framing; exact numbers with source files and commits; figure plan) |
| `c619448` | `—` | 2026-10-10 23:57 | Paper package: original and paper figures, pertinent tables, captions, plain-language reading guide |
| `99b4b35` | `c0d7b75` | 2026-10-11 16:20 | D27: plan steps/cohorts not run or changed; outline and guide: 5-PC sensitivity result and global-structure scope limitation |
| `7b740eb` | `—` | 2026-10-11 16:21 | Correct 5-PC statement: S5 already run in GSE20916, GSE37364, TCGA (rule not met); GSE117606 S5-S8 pre-specified, not run |
| `611b239` | `082b945` | 2026-10-11 16:28 | Sensitivity table for all GEO cohorts (GSE117606 S5-S8 already run, 80/80); R1 wording on sensitivity to analysis choices; fix stale GSE117606 step05 summary |
| `be79af3` | `9bcd90c` | 2026-10-11 17:13 | Source data tables (GEO): PH diagram values, primary validation summary with rule_met, centroid distances, arm scores, per-sample tables; all 19 reproduction checks passed |
| `78642b0` | `e1a85f9` | 2026-10-11 17:24 | Fig 2: GSE41258 H1 diagrams per stage (adjusted, raw and shape-only) and H0 death boxes, from source data; sources and caption |
| `5d718a2` | `7627009` | 2026-10-11 17:30 | paper/final structure; Fig 2 moved to main_figures; journal caption and plain-language reading guide (.md + .docx); INDEX.md; captions_to_docx helper |
| `5018359` | `626bdd0` | 2026-10-11 17:46 | Fig 2 v2 (standard H0+H1 diagrams, representative subsample, PC insets) with caption, guide (.docx) and sources; v1 renamed Fig2_v1 and moved to supplementary; INDEX: both 'draft, comparing' |
| `672e0aa` | `cf50541` | 2026-10-11 17:51 | Fig 2 v2: rows a and b show the same samples (raw-mode representative subsample); shape-only diagrams recomputed for display, raw recomputation checked against source data; caption, guide, sources updated |
| `8fa03ee` | `48b04e4` | 2026-10-11 17:58 | Fig 2 v2 approved; Fig 3 (GSE41258 barcodes, Betti curves, H1 landscapes) with caption, reading guide (.docx) and sources; INDEX row |
| `7bbd066` | `2ed8b9f` | 2026-10-11 18:06 | Fig 3 approved; main Tables 1 (cohorts) and 2 (GEO primary validation) as editable Word + CSV, with captions, reading guides (.docx), sources and INDEX rows |
| `fa87662` | `3f71e8b` | 2026-10-11 20:45 | Table 2: shape-only H1 landscape columns only, patients in both groups, footnote on shared-patient handling (from code); Table 1: GSE117606 75 vs 70 footnote, GSE117607 checked; Supplementary Table S_H1_landscape (4 columns) |
| `fb12175` | `6dbb326` | 2026-10-11 21:38 | D28 description correction: 'label-permutation test (one sample per patient per stage)' replaces 'patient-level' in paper/final, results outline and guide; Fig 2 v1/v2 and Fig 3 re-rendered (text only); conservative-effect sentence in reading guides |
| `bed83c4` | `—` | 2026-10-11 21:55 | paper_package: CAPTIONS wording (D28); README marks the package superseded by paper/final/ (not for submission) |
| `9d26533` | `7b600e3` | 2026-10-12 06:42 | Table 1: sample-set footnote (main/bulk, from step 1 code); Table 2: Dimension detail removed from Word table (kept in CSV), footnote adjusted; Table 2 guide: single-version shape differences paragraph; Tables 1-2 approved |
| `6622f2c` | `ef548e0` | 2026-10-12 06:54 | Fig 4 heterogeneity (adjusted; per-sample centroid distances + PERMDISP D) and unadjusted candidate supplementary figure, with captions, reading guides (.docx), sources and INDEX rows |
| `12f3bda` | `bc58a66` | 2026-10-12 06:56 | Fig 4 guide: nearest-neighbour vs centre-distance medians (descriptive); explanation corrected to match numbers; values added to Fig4_sources.md |
| `e98f81e` | `—` | 2026-10-12 06:57 | Fig4_sources.md: nearest-neighbour and centre-distance medians and ratios (descriptive check) |
| `bacc031` | `e8b4716` | 2026-10-12 07:09 | Fig 4: 'between-sample heterogeneity' titles and wording; pattern description (A->C 3/4, N->A 2/4, no C->M increase); composition check (descriptive) in guide and sources; Fig 4 approved |
| `c62afba` | `6635f45` | 2026-10-12 07:16 | Supplementary Table S_cohort_differences (four GEO cohorts side by side, existing numbers only) with caption, reading guide, sources and INDEX row |
| `052f546` | `77544aa` | 2026-10-12 07:27 | TableS_cohort_differences approved; Fig 5 CMS subtypes and carcinoma heterogeneity with caption, reading guide (.docx), sources and INDEX row |
| `b45b70f` | `—` | 2026-10-12 09:14 | Fig 5 caption/guide: why GSE41258 A->C value differs from Fig 4 (bias adjustment, single uncorrected test vs Holm) |
| `dea68c0` | `e8c515c` | 2026-10-12 09:22 | Fig 5 panel a: carcinoma distance to all-carcinoma vs own-CMS centre (10 PCs; checked vs addendum 9 D_within); former PC1/PC2 panel moved to supplementary FigS_CMS_PC1_PC2; captions, guides, sources, INDEX |
| `462255d` | `e558513` | 2026-10-12 09:32 | Fig 5 approved; Fig 6 (117 variable genes: stage-mean and within-stage SD heatmaps, four example genes by fixed rule) with caption, reading guide (.docx), sources and INDEX row |
| `a3928f3` | `d7c9516` | 2026-10-12 09:53 | Fig 6: exact addendum 9 replication rule and per-gene cohort counts; grey normal/adenoma points; fixed heatmap label set; Supplementary Table S_117_genes with caption, guide, sources, INDEX |
| `7e39a42` | `e97c099` | 2026-10-12 12:27 | Fig 6 approved; Fig 7 chromosome arms (arm variability, 20q scores, TCGA copy-number validation, dispersion shares) with caption, reading guide (.docx), sources and INDEX row |
| `e4bcfbb` | `5037e39` | 2026-10-12 12:59 | Supplementary Table S_power_controls (addendum 5 controls, addendum 10 TCGA control, step 5 subsample sizes) with caption, reading guide, sources and INDEX row |
| `570535f` | `f2b46b6` | 2026-10-12 13:08 | TableS_power_controls: part D (D7 raw-mode power, D8 shape-only sample-size sensitivity); caption, guide, sources, INDEX |
| `ca57ac0` | `c9f0fa1` | 2026-10-12 13:20 | Fig 7d: 20q and CMS only; four-arm results moved to Supplementary Table S_four_arms (negative f, D24 caution); GSE41258 20q variability described as borderline; caption, guide, sources, INDEX |
| `a08a7f8` | `—` | 2026-10-12 13:24 | Fig 7: panel a wording (variability up across most arms in 3 cohorts, near zero in GSE117606; 20q in top four with 20p, 13q, 8q); Fig 7 approved |
| `3530a10` | `54e63c5` | 2026-10-12 13:29 | Fig 1 workflow schematic (study design, pipeline, decision path) with caption, reading guide (.docx), sources and INDEX row |
| `94aa69e` | `757aa69` | 2026-10-12 13:51 | Supplementary Table S_TCGA_stages (step 5, PERMDISP, CMS by stage, normal vs Stage I secondary); Fig 1 TCGA label widened and re-rendered; captions, guides, sources, INDEX |
| `a3e22da` | `100108f` | 2026-10-12 14:25 | Fig 1: 'microarray' on GEO platform lines; (Figs 2-3) and (Table 2) in pipeline boxes; caption updated; Fig 1 approved |
| `1145aee` | `f3c0abd` | 2026-10-12 14:36 | Supplementary Table S_sensitivity: step 5 sensitivity analyses (PC/gene settings, sample sets and adjustments, TCGA) |
| `3760f00` | `5ee2ba8` | 2026-10-12 14:40 | TableS_sensitivity: analogue rows not bold, labelled 'analogue only'; approved |
| `52a81cb` | `—` | 2026-10-12 14:44 | Fig 7 guide 20q sentence; supplementary statuses ready for co-author review; README for Prof Peng |
| `1309d94` | `—` | 2026-10-13 21:52 | Methods draft (written from code) with provenance table and ambiguity list |
| `ad12fb8` | `b6f352d` | 2026-10-13 22:17 | Methods review: precise W1 description, data provenance, PERMDISP and symbol checks, Supplementary Methods, prespecified wording |
| `7142899` | `338e63d` | 2026-10-14 06:10 | docs: wording only - 'pre-registered/pre-registration' replaced by 'prespecified before analysis' (analysis plan file and its history unchanged) |
| `188b844` | `e12669a` | 2026-10-14 06:27 | Cross-cohort analyses with harmonised HGNC symbols (122 genes; S7 deviation); precise Wasserstein wording; S8 duplicate check; untrack personal draft and superseded figure script |
| `566df08` | `—` | 2026-10-14 06:40 | Fig 6 and Fig 7 approved after symbol harmonisation (formatting pass pending) |

