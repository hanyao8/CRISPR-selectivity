# CRISPR-selectivity
Project to investigate targeting selectivity in the CRISPR CAS9 system using mean field models and other algorithms


Models 1,2,3 are mean field models of increasing complexity which calculate the targeting selectivity between a targeting sequence and a targeted genome (strand). Model 1 describes the genome using mean field frequency parameters with equal binding energy across all complementary base pairs. The later models incorporate more features such as dimer frequencies, non-degenerate binding energies and ultimate the 'Santa Lucia rules' of DNA hybridisation in model 3.

Several algorithms were constructed to compare the results of the mean field models with simulations using real DNA sequence data in human chromosomes. Due to the size of genome data, they are omitted from the folder 'nt_data'.

Project Supervisor: Stefano Angioletti-Uberti

With thanks to: Enrico Fontana, Harrison Zhu, Mike Chen
