Class 12: Structural Bioinformatics II
================

\#\#Preping for docking We want to produce a protein-only PDB file and a
drug only PDB file.

``` r
library(bio3d)
#download the pdb file
get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

    ## [1] "./1hsg.pdb"

``` r
pdb <- read.pdb("1hsg.pdb")
protein<- atom.select(pdb, "protein", value= TRUE)
write.pdb(protein, file= "1hsg_protein.pdb")
```

and for the ligand

``` r
ligand <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(ligand, file = "1hsg_ligand.pdb")
```

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg.pdb")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

q1 non-prtein residues ans: Non-protein/nucleic resid values: \[ HOH
(127), MK1 (1) \]

q2 yes you can see the whole of the binding site. the resolution is too
high for the hydrogen atoms to appear in crystal structures

\#\#process our docking results

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "1hsg_protein.pdb")
```

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```
