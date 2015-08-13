#### Example: completing & sorting a pedigree

Before creating a pedigree object, the pedigree must be: (i) completed, that is, any individual that appears as ancestor must also appear as a row in the pedigree, and (ii) sorted, so that ancestors preceed progeny in the pedigree. The functions ```editPed()``` can be used to complete and sort pedigrees. In the following example we use the small pedigree printed below to illustrate the use of the function editPed() and the creation of a pedigree object.

| Subject  |      Sire     |  Dam |
|----------:|-------------:|------:|
| 1	| NA| 	12
| 2	| 11| 10
| 3	| 1	| 2
| 4	| 1	| 2

```R
 library(pedigreeR)
 animal=c()
 sire=c()
 dam=c()
 PED=editPed() 
 head(PED)
```

[Home](https://github.com/Rpedigree/pedigreeR)
 
