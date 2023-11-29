STAR --runThreadN 10 --genomeDir /media/user/sdg/home/qiuguo/10X/SHP1/H9M/H9M/H9M_STAR/ \
 --readFilesCommand zcat \
 --readFilesIn  /media/user/sdg/home/qiuguo/10X/SHP1/H9M/H9M/raw/SC94/SC94_R2.fq.gz /media/user/sdg/home/qiuguo/10X/SHP1/H9M/H9M/raw/SC94/SC94_R1.fq.gz \
 --soloType CB_UMI_Simple
 --soloCBwhitelist /path/to/cell/barcode/whitelist \
 --outFileNamePrefix ./SC94