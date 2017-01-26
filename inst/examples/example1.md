### Example 1: Completing & sorting a pedigree

A valid pedigree object is a pedigree that is : 
* Complete : Any individual that appears as ancestor must also appear as a row in the Pedigree. 
* Sorted :  Ancestors must preceed progeny in the Pedigree.
* Formal Class : It is an S4 object of formal class 'Pedigree'. ( See Section B )

The function ```editPed()``` can be used to complete and sort pedigrees. In the following example we use the small pedigree printed below to illustrate the use of the function ```editPed()``` in the creation of a valid pedigree object from an incomplete and scrambled pedigree.

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


#### Section A : Scrambling the above pedigree and rendering it incomplete. 

Here is some example code that : 
* Creates an incomplete version of the pedigree by omitting rows that correspond to four of the five 1st generation parents. 
* Scrambles the incomplete Pedigree. 
```R
        ## Create incomplete version of example pedigree 
        PedInc <- data.frame(sire=as.character(c(NA,1,3,5,6,4,8,1,10,8)),
                  dam= as.character(c(NA,2,2,NA,7,7,NA,9,9,13)),
                  label=as.character(5:14))
        
        ## Scramble incomplete pedigree
        PedScram <- PedInc[sample(replace=FALSE, 1:10),] 

```
#### Section B : Reconstructing the complete, sorted pedigree and converting it into a valid 'Pedigree' object

Below is some more example code that :
* Completes and sorts the incomplete pedigree using ```editPed``` and returns a DataFrame.
* Uses the function ```pedigree``` to convert the data frame representation into an S4 object of formal class 'Pedigree'. 

```R    
        library(pedigreeR)
        
        ## Complete and sort incomplete Pedigree using editPed
        PedEdit<- editPed(sire=PedScram$sire, dam= PedScram$dam, label=PedScram$label) 
        
        ## Converted the data frame PedEdit into an S4 object of formal class 'Pedigree'
        PedFinal<- with(PedEdit, pedigree(label=label, sire=sire, dam=dam))

```
The object PedFinal is a valid pedigree object that can be used for further analysis with PedigreeR and pedigreemm. This is what it should look like for this example : 

```R
PedFinal

   sire  dam
1  <NA> <NA>
2  <NA> <NA>
3  <NA> <NA>
4  <NA> <NA>
5  <NA> <NA>
7     3    2
8     5 <NA>
6     1    2
11    8 <NA>
9     6    7
10    4    7
13   10    9
12    1    9
14    8   13

str(PedFinal)

Formal class 'pedigree' [package "pedigreeR"] with 3 slots
  ..@ sire : int [1:14] NA NA NA NA NA 3 5 1 7 8 ...
  ..@ dam  : int [1:14] NA NA NA NA NA 2 NA 2 NA 6 ...
  ..@ label: chr [1:14] "1" "2" "3" "4" ...

```

[Home](https://github.com/Rpedigree/pedigreeR)
 
