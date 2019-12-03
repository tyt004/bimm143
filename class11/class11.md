class 11: structural bioinformatics 1
================

\#\#The PDB database for biomoleculer structure data

Q1: Download a CSV file from the PDB site (accessible from “Analyze” -\>
“PDB Statistics” \> “by Experimental Method and Molecular Type”. Move
this CSV file into your RStudio project and determine the percentage of
structures solved by X-Ray and Electron Microscopy. Also can you
determine what proportion of structures are protein?

download CSV (accessible from “Analyze” -\> “PDB Statistics” \>“by
Experimental Method and Molecular Type”)

``` r
#read csv
data<- read.csv("Data Export Summary.csv")
data
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other
    ## 1               X-Ray   131278          2059               6759     8
    ## 2                 NMR    11235          1303                261     8
    ## 3 Electron Microscopy     2899            32                999     0
    ## 4               Other      280             4                  6    13
    ## 5        Multi Method      144             5                  2     1
    ##    Total
    ## 1 140104
    ## 2  12807
    ## 3   3930
    ## 4    303
    ## 5    152

total number of enteries

``` r
total<-sum(data$Total)
```

proportion from each method

``` r
(data$Total/sum(data$Total))*100
```

    ## [1] 89.0702879  8.1419744  2.4984742  0.1926305  0.0966331

propotion that are protein

``` r
round(sum(data$Proteins)/sum(data$Total)*100,2)
```

    ## [1] 92.71

Q2: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?

Direction 1.Go to PDB and Download the PDB file 2. Open VMD and go to
the main window medium(there is three: big, medium, and small) 3. Use
file broswer to open desired file 4. To change graphics go to VMD and
click graphics 5. Change drawing methods

other things you can change coloring method selecteed atoms is the atoms
you ant to change ex: water not protein and not water create rep: makes
a new file for you to change w/o affecting the other

Q3: Water molecules normally have 3 atoms. Why do we see just one atom
per water molecule in this structure? Because the water molecules are
too small for the resolution Q4: There is a conserved water molecule in
the binding site. Can you identify this water molecule? What residue
number does this water molecule have (see note below)? there are four
resid192

``` r
library(bio3d)

pdb<-read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

`{r}, eval=FALSE ligand<-`
