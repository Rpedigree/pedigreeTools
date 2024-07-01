### Example: Pedigree analyses for self-polinated species.

Here we illustrate the use of the function ```getASelfing``` in the computation of the A matrix from an example pedigree which contains information on number of selfing cycles (```nCycles```), and illustrate the differences in the A matrix before and after taking selfing into account. 


#### Section A: Creation of example ```data.frame``` pedigrees with and without selfing cycles

Below we create an example pedigree ( identical to the one in Example 1 ) before taking selfing into account.

```R
pedNoCycles <-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
                      dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
                      label=as.character(1:14))
```

| Subject  |      Sire     |  Dam |
|----------:|-------------:|------:|
| 1	| NA| 	NA
| 2	| NA| 	NA
| 3	| NA| 	NA
| 4	| NA| 	NA
| 5	| NA| 	NA
| 6	| 1| 2
| 7	| 3	| 2
| 8	| 5	| NA
| 9	| 6| 	7
| 10	| 4| 7
| 11	| 8	| NA
| 12	| 1	| 9
| 13	| 10| 9
| 14	| 8| 13



Here we create an example pedigree similar to the one used previously with an additional column called ```nCycles``` that takes selfing into account. 

Note that this pedigree must be complete and sorted ( see Example 1 )before using it with the functions in Section B.

```R
pedCycles <-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
                      dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
                      label=as.character(1:14),nCycles=c(0,0,0,0,0,0,0,5,0,0,0,0,3,0))
```

| Subject  |      Sire (Parent 1)    |  Dam (Parent 2) | nCycles |
|----------:|-------------:|------:|-------:|
| 1	| NA| 	NA | 0|
| 2	| NA| 	NA | 0|
| 3	| NA| 	NA | 0|
| 4	| NA| 	NA | 0|
| 5	| NA| 	NA | 0|
| 6	| 1| 2 | 0|
| 7	| 3	| 2 | 0|
| 8	| 5	| NA | 5|
| 9	| 6| 	7 | 0|
| 10	| 4| 7 | 0|
| 11	| 8	| NA | 0|
| 12	| 1	| 9 | 0|
| 13	| 10| 9 | 3|
| 14	| 8| 13 | 0|


#### Section B: Computation of the additive relationship matrix using the function ```getASelfing```

Here we create the matrix Aself from ```pedCycles``` using getAselfing which computes the additive relationship matrix from the pedigree given the information contained in ```nCycles```.

```R

library(pedigreeTools)

Aself <- getASelfing(ID=pedCycles$label,Par1=pedCycles$sire,Par2=pedCycles$dam,nCycles=pedCycles$nCycles,nCyclesDefault=0)

```
We also create the matrix A from the pedigree before selfing following the method outlined in example 1:

```R
pedNoCycles = pedigree(label = pedNoCycles$label, sire = pedNoCycles$sire, dam = pedNoCycles$dam)
A <- getA(pedNoCycles)
```
Visualizing the difference matrix of ```Aself``` and ```A``` shows us that only the members 8,11,13 and 14 are different, as is expected.

```R

diff <- Aself - A
image(diff)

```
<img src="https://github.com/siddharth51292/pedigreeTools/blob/patch-2/inst/examples/diffMatrix.png" width="500">

[Home](https://github.com/Rpedigree/pedigreeTools)
 

 
