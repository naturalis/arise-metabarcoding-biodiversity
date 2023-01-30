This section contains the following intermediate and result files:

- [Primers.xlsx](Primers.xlsx) - Primers used by BaseClear for NovaSeq platform
- [UNITE_database] - This directory contains UNITE database general fasta files 

### Prepare your data:

To use this pipeline you need to know the primers that has been used and then seperate them by markers. The data we use is seperated in 4 different kind of markers, they are presented in this table.

| Markers | #I7: 1 – 12    | #I7: 13 – 24    | #I7: 25 – 36 | #I7: 37 – 48 |
| :-----: | :---: | :---: | :---: | :---: |
| I5: 1 – 8 | ITS locatie 1 21022-375801-001 NI030   | ITS locatie 3 21022-375801-003 NI032   | CO1 locatie 1 21022-375801-001 NI057 | ITS nioo locatie 1 21022-375801-001 NI059 |
| I5: 9 – 16 | ITS locatie 2 21022-375801-002 NI031   | ITS locatie 4 21022-375801-004 NI033   | 18s locatie 1 21022-375801-001 NI058 | 18s nioo locatie 1 21022-375801-001 NI060 |
| I5: 17 – 24 | 16s locatie 1 21022-375801-001 NI034   | 16s locatie 3 21022-375801-003 NI036   |  |  |
| I5: 25 – 32 | 16s locatie 2 21022-375801-002 NI035   | 16s locatie 4 21022-375801-004 NI037   |  |  |

It is important to seprate the data before u run the pipeline. In this case we are ging to seperate the data by ITS markers: `NI030, NI031, NI032, NI033`
To seperate the data:

    find . -name "*NI030*" -exec mv -t ./ITS {} +

Where `"*NI030*"` is the marker that we wanted to move to `./ITS` directory.
Use this to move all needed data with the right markers to a specific directory. This will be the data that we will be using for this pipeline. 

### Set up data on the server

To use the server we need to set up the data on the MaaS server. To do so use scp

    scp *.gz ubuntu@145.136.253.38:winny.thoen/Data/

Where `*.gz` are all the files within the current directory and `ubuntu@145.136.253.38:~/winny.thoen/Data/` is the directory on the MaaS server. It is not possible to place your data from a local computer in your specific account on the server. To do so de data can be moved from `~/winny.thoen/Data` to `/home/winny.thoen/arise-metabarcoding-biodiversity/data/raw_sequencesNovaSeq`. Where `winny.thoen` is the personal account where the data needs to be stored. To move the data:
    
    sudo mv *.gz /home/winny.thoen/arise-metabarcoding-biodiversity/data/raw_sequencesNovaSeq
    
If the data is to big, `mv` doesn't work. Use find:

    sudo find . -name '*gz*' -exec mv {} /home/winny.thoen/arise-metabarcoding-biodiversity/data/raw_sequencesNovaSeq \;

    
After moving the data check the permissions before using the dada2 pipeline. If necessary change the permissions using `chmod`.


