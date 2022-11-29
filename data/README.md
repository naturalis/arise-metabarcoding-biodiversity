This section contains the following intermediate and result files:

- [ARISE_Sample_information_logbook.xlsx](ARISE_Sample_information_logbook.xlsx) - sample information from the soil data (locations ect)
- [ASVtab_raw.csv](ASVtab_raw.csv) - csv data from the soil data
- [OTU97tab_tax.csv](OTU97tab_tax.csv) - OTU table from the soil data
- [Primers.xlsx](Primers.xlsx) - Primers used by BaseClear for NovaSeq platform

Prepare your data:

To use this pipeline you need to know the primers that has been used and then seperate them by markers. 
The data we use is seperated in 4 different kind of markers, they are presented in this table.

| Markers | #I7: 1 – 12    | #I7: 13 – 24    | #I7: 25 – 36 | #I7: 37 – 48 |
| :-----: | :---: | :---: | :---: | :---: |
| I5: 1 – 8 | ITS locatie 1 21022-375801-001 NI030   | ITS locatie 3 21022-375801-003 NI032   | CO1 locatie 1 21022-375801-001 NI057 | ITS nioo locatie 1 21022-375801-001 NI059 |
| I5: 9 – 16 | ITS locatie 2 21022-375801-002 NI031   | ITS locatie 4 21022-375801-004 NI033   | 18s locatie 1 21022-375801-001 NI058 | 18s nioo locatie 1 21022-375801-001 NI060 |
| I5: 17 – 24 | 16s locatie 1 21022-375801-001 NI034   | 16s locatie 3 21022-375801-003 NI036   |  |  |
| I5: 25 – 32 | 16s locatie 2 21022-375801-002 NI035   | 16s locatie 4 21022-375801-004 NI037   |  |  |

It is important to seprate the data before u run the pipeline. 
