# Validation of the H2 antagonist run

Rows: 6 (one per drug x ATC code x route) for 5 drugs.
Receptor: Inactive human histamine H2 receptor with famotidine bound (cryo-EM, 3.0 A) (7UL3).

| check                                                                          | result   | detail                                                                                                                                                 |
|:-------------------------------------------------------------------------------|:---------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
| SMILES parse and InChIKey agree                                                | pass     | 6 rows, mismatches: none                                                                                                                               |
| one structure per drug                                                         | pass     | 5 distinct structures for 5 drugs                                                                                                                      |
| RDKit molecular weight reproducible                                            | pass     | largest difference 0.004 Da                                                                                                                            |
| QED reproducible from the SMILES                                               | pass     | largest difference 0.0005                                                                                                                              |
| QED close to ChEMBL's own value                                                | pass     | largest difference 0.005 (cimetidine); ChEMBL rounds to two decimals and uses its own parent structure                                                 |
| Crippen log P and PubChem XLogP3 agree within 1.5 log units                    | pass     | largest gap 1.16 (ranitidine) over 5 drugs                                                                                                             |
| redocking recovers the crystallographic famotidine pose (< 2 A)                | pass     | selected pose 1.15 A (rank 3 of 20); closest pose 0.99 A at rank 6                                                                                     |
| all poses inside the orthosteric site (centroid < 5 A from the crystal ligand) | pass     | {"cimetidine": 3.85, "ranitidine": 0.25, "famotidine": 0.76, "nizatidine": 0.45, "roxatidine": 4.44}                                                   |
| poses contact Asp98, the conserved aminergic anchor                            | pass     | contacting: ['cimetidine', 'ranitidine', 'famotidine', 'nizatidine', 'roxatidine']; not: []                                                            |
| at most two non-intravenous routes per ATC code (rules 3 and 4)                | pass     | {('cimetidine', 'A02BA01'): 1, ('famotidine', 'A02BA03'): 1, ('nizatidine', 'A02BA04'): 1, ('ranitidine', 'A02BA02'): 1, ('roxatidine', 'A02BA06'): 1} |
| one row per ATC code and route (rule 5)                                        | pass     | 6 rows for 5 drugs                                                                                                                                     |
| intravenous bioavailability is 1 and oral is below it                          | pass     | IV rows [1.0], oral rows [0.575, 0.5, 0.45, 0.7]                                                                                                       |
| every value has a source and every gap has a reason                            | pass     | none                                                                                                                                                   |
| therapeutic index equals log(toxicity) - log(1/activity)                       | pass     | checked 5 rows, off: none                                                                                                                              |

## Docking

| drug       | atc     |   vina_score_kcal_per_mol |   selected_pose_rank |   poses_in_orthosteric_site |   site_offset_A |   n_poses |   best_minus_second | pocket_contacts                                                                                                      |   n_hydrogen_bonds |   redock_rmsd_to_crystal_A |   redock_best_pose_rmsd_A |   redock_best_pose_rank | receptor                                             | box                                                    |
|:-----------|:--------|--------------------------:|---------------------:|----------------------------:|----------------:|----------:|--------------------:|:---------------------------------------------------------------------------------------------------------------------|-------------------:|---------------------------:|--------------------------:|------------------------:|:-----------------------------------------------------|:-------------------------------------------------------|
| cimetidine | A02BA01 |                    -5.668 |                    2 |                           7 |            3.85 |        20 |                0    | ASP98; VAL99; CYS102; THR103; LEU149; VAL176; ASP186; THR190; TRP247; TYR250; PHE251; PHE254                         |                  4 |                     nan    |                    nan    |                     nan | PDB 7UL3 chain A (human H2 receptor, UniProt P25021) | 24.0 A cube centred on the crystallographic famotidine |
| ranitidine | A02BA02 |                    -5.79  |                    1 |                          12 |            0.25 |        20 |               -0.15 | ASP98; VAL99; CYS102; THR103; CYS174; LYS175; VAL176; ASP186; THR190; TRP247; TYR250; PHE251; PHE254; LEU274; TYR278 |                  3 |                     nan    |                    nan    |                     nan | PDB 7UL3 chain A (human H2 receptor, UniProt P25021) | 24.0 A cube centred on the crystallographic famotidine |
| famotidine | A02BA03 |                    -5.682 |                    3 |                           7 |            0.76 |        20 |                0.23 | ASP98; VAL99; CYS102; THR103; VAL176; ASP186; THR190; TYR250; PHE251; ARG257; GLU270; LEU274                         |                  8 |                       1.15 |                      0.99 |                       6 | PDB 7UL3 chain A (human H2 receptor, UniProt P25021) | 24.0 A cube centred on the crystallographic famotidine |
| nizatidine | A02BA04 |                    -5.606 |                    1 |                          12 |            0.45 |        20 |               -0.02 | TYR94; ASP98; VAL99; CYS102; CYS174; LYS175; VAL176; ASP186; THR190; TRP247; TYR250; PHE251; PHE254; LEU274; TYR278  |                  3 |                     nan    |                    nan    |                     nan | PDB 7UL3 chain A (human H2 receptor, UniProt P25021) | 24.0 A cube centred on the crystallographic famotidine |
| roxatidine | A02BA06 |                    -7.132 |                    1 |                          12 |            4.44 |        20 |               -0.14 | ASP98; VAL99; CYS102; VAL176; ASP186; TRP247; TYR250; PHE251; PHE254; LEU274; TYR278                                 |                  1 |                     nan    |                    nan    |                     nan | PDB 7UL3 chain A (human H2 receptor, UniProt P25021) | 24.0 A cube centred on the crystallographic famotidine |

## What is measured and what is illustrated

- Pose coordinates, Vina scores, pocket contacts and hydrogen-bond distances are computed from the
  structure and the docking run; the binding page shows them unchanged.
- The path the ligand takes into the pocket is an illustration along the receptor's extracellular
  axis. It is not a molecular dynamics trajectory and no energy is computed along it.
- FAERS counts are spontaneous reports: a reporting signal, never an incidence rate.
- The rat oral LD50 is an acute single-dose animal value; the therapeutic index built from it and
  the human label dose is a ratio for comparing these five drugs with each other, not a clinical
  safety margin.
