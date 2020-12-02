This directory contains applications for stand-alone use, 
built specifically for a Linux 64-bit machine.

For help on the bigBed and bigWig applications see:
http://genome.ucsc.edu/goldenPath/help/bigBed.html
http://genome.ucsc.edu/goldenPath/help/bigWig.html

View the file 'FOOTER.txt' to see the usage statement for 
each of the applications.

##############################################################################
Thank you to Bob Harris for permission to distribute a binary
version of the lastz and lastz_D programs, from:

   https://github.com/lastz/lastz

Version 1.04.00 as of April 2018:

-rwxrwxr-x 1  625283 Apr  6 11:15 lastz-1.04.00
-rwxrwxr-x 1  628835 Apr  6 11:15 lastz_D-1.04.00

$ md5sum lastz*
429e61ffdf1612b7f0f0c8c2095609a7  lastz-1.04.00
4f9a558a65c3a07d0f992cd39b3a27e1  lastz_D-1.04.00

##############################################################################
This entire directory can by copied with the rsync command
into the local directory ./

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

Or from our mirror site:

rsync -aP rsync://hgdownload-sd.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

Individual programs can by copied by adding their name, for example:

rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSize ./
