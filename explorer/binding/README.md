# Ligand Binding Theatre

Drugs grouped by the receptor they act on, each one animated into its binding pocket, with the
receptor built from its own topology rather than from a cartoon.

```
python3 curation/h2_features.py                               # writes h2_features/binding_data.json
python3 explorer/binding/build_binding.py \
        --features h2_features/binding_data.json              # writes explorer/binding/data/receptors.json
python3 -m http.server 8777 --directory explorer              # open /binding/index.html
```

## Where each thing on screen comes from

| On screen | Source | Honest? |
|---|---|---|
| Ligand atoms and bonds | RDKit ETKDG conformer, MMFF-minimised, of the PubChem SMILES | real 3D structure of the free ligand |
| Helix lengths (`TM5 175-202`) | UniProt `Transmembrane` features of the target's accession | real residue ranges |
| Orange residue markers | UniProt `Site` / `Binding site` features, placed at the residue's own fractional depth in its helix | real residue, schematic position |
| Disulfide line | UniProt `Disulfide bond` | real pairing, schematic position |
| Affinity (pChEMBL, Kd) | ChEMBL median over activities at that target, only when two or more documents agree | real |
| Last beat (TM6 holds or swings) | ChEMBL `action_type` of the drug at that target | real label, schematic motion |
| Where the helices sit around the bundle | a canonical class-A arrangement written into `binding.js` | **schematic**, said so on the page |

So the contact distances printed on the dashed lines are indicative, not docking output. The page says
this in its own caveat panel and lists the PDB entries for the same accession, which is where measured
geometry should come from. Nothing here is a substitute for docking; it is a way to see which residues
the evidence says matter, and how one group of drugs compares with another on affinity.

## The four beats

1. **approach** - the ligand drifts in from the extracellular bulk.
2. **pocket mouth** - it descends between the helices that line the orthosteric site.
3. **contacts** - dashed lines to each annotated residue, with the affinity.
4. **verdict** - for an antagonist or inverse agonist TM6 (blue) stays put: the pocket is occupied and the
   receptor is not switched on. For an agonist TM6 swings out, the hallmark of class-A activation.

## Adding another group

`build_binding.py` is not specific to one receptor: run `curation/h2_features.py` with another ATC
level-3 group (`ATC_GROUP=C07AB python3 curation/h2_features.py`), point `--features` at the new
`binding_data.json`, and the page grows another target-group chip. A target with no UniProt accession,
or a non-GPCR target, still draws: the bundle arrangement falls back to an even fan, and the caveat
panel says which annotations were missing.
