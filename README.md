Code for identifying resistance mutations in Sars-Cov-2 against 3CLpro and RdRp based on the [Stanford Sars-Cov-2 database](https://covdb.stanford.edu/drms/rdrp/)  

Currently only the `create_BN_import_nanopore.R` script is relevant (for Illumina, MIK or other external HF see `create_BN_import_nanopore.R`). Also, we only create a BN import file for 3Clpro mutations. Although the code can be used to identify any RdRP resistance mutations as well.  

Possible future developments:
- Use the Stanford SQL database and extract relevant information. See for example [here](https://github.com/hivdb/covid-drdb-payload), [here](https://github.com/hivdb/covid-drdb-payload/issues/1022) and [here](https://github.com/hivdb/covid-drdb/blob/master/derived_tables/08_rdrp_resistance_mutations.sql).
- Integrate code with other routine analysis.
- Run the code "oppsett" by "oppsett" whenever we get new results.
