# **ComPy** 
## *Compare .bam and .vcf files using a sqlite3 database* 
 
 
### Overview
1. Install ComPy 
2. Run ComPy 
	1. Compare files 
	2. Extract data from database 
	3. Remove data from database 
 
 
### **1. Install ComPy** 
 
> ***You have to make sure to run ComPy on a Linux machine!*** 
> 
> If you don't have permission to install packages make sure to have anaconda installed 
 
#### *Download ComPy* 
 
If you have permissions to install you may easily use
 
    pip install ComPy
 
to install ComPy on your machine. Otherwise, you can install ComPy with anaconda
 
    conda install -c biconda ComPy
 
The third possibility is to manually clone this GitHub repository, but in this case you have to install all dependencies manually by hand.
 
> **I highly recommend Anaconda or something similar (like miniconda) to isntall ComPy in the most easies way.** 
 
 
#### *An important thing to note:* 
 
Compy will save all plots, the database and the logs per default to you home directory! (/home/*yourusername*/ComPy/)