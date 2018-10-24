#!/bin/bash
### Script to download data from EMBL GeneCore using Aspera  
ASCP="/home/tomas/.aspera/connect/bin/ascp"
/home/tomas/.aspera/connect/bin/ascp -k 1 -i /home/tomas/.aspera/connect/etc/asperaweb_id_dsa.openssh \
https://w3-09.genecore.embl.de:/aspera/faspex/received/1604 /home/tomas/CEITEC_lab/Dis3L2/Nandan_new_set/