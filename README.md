# CRISPR-selectivity
Project to investigate targeting selectivity in the CRISPR CAS9 system using mean field models and other algorithms

Folders:

-Algorithms

-Model 1

-Model 2

-Model 3

-nt_data


Models 1,2,3 are mean field models of increasing complexity which calculate the targeting selectivity between a targeting sequence and a targeted genome (strand). Model 1 describes the genome using mean field frequency parameters with equal binding energy across all complementary base pairs. The later models incorporate more features such as dimer frequencies, non-degenerate binding energies and ultimate the 'Santa Lucia rules' of DNA hybridisation in model 3.

Several algorithms were constructed to compare the results of the mean field models with simulations using real DNA sequence data in human chromosomes. nt_data contains the sequence data.

Project Supervisor: Stefano Angioletti-Uberti
Thanks to: Enrico Fontana, Harrison Zhu
