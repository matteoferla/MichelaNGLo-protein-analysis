## Probability homozygous

gnomAD variants are shown with the idea that they are phenotypically neutral at that zygosity.

Homozygously neutral is certainly true for variants over .5%,
but stuff gets murky the less common they are.

The gnomAD control database does not have every individual on Earth.
So this means that there may be one or more peoples that are homozygous for a given allele found
heterozygously only in gnomAD.
If these people showed a disease phenotype in a country with a rare disease genome sequencing program,
such as the UK, they may be in the Genomics England database.

However, doing a square of the allele frequency to get the homozygous frequency
does not account for heterogenity of the population.
If the homozygous frequency is small enough the assumption that allele frequency
is the same as heterozygous frequency is a totally sane thing to do.

The probability of finding one or more persons with homozygous alleles would
follow a negative binomial, e.g. one minus the binomial probability of getting zero successes from population sized trials with the squared frequency as probability.

Using the ethnicity split may give a better result, however, 
whereas the number of individuals in each subgroup is known within gnomAD,
for a larger population it is much trickier
as one can either assume the population has the same split as gnomAD (bogus assumption)
or one sifts through census data which has different groups (e.g. finnish and nfe vs. caucasian).

The maths to calculate that would be one minus elementwise product of 
the vector of probabilities per subpopulation group, each calculated as the binomial pmf of zero k/successes from a sample/n that is total subpopulation members
with a probability that is the squared allele frequecy for that group.

The problem then is that it assumes each subgroup is homogeneous with no consanguinity.

Finally, there is the problem that finding out if there is one or more homozygous persons
is kind of pointless...

