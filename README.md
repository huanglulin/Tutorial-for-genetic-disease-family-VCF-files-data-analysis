
# Tutorial for genetic disease family VCF files data analysis

This is a basic tutorial for genetic disease family VCF files data analysis. There are six families for tutorial. Family 1-4 can be analyzed by trioFilter.py (https://github.com/WGLab/Trio-filter). TrioFilter can find out the candidate genes automatically with the merged family VCF file, the pedigree file (.ped) and the config file (config.ini).

## Case 1. IHA (idiopathic hemolytic anemia) sample, which is a recessive disease. 84060 is case (son), 84615   88962   92157 are unaffected brother, mother and father, respectively.
Here is the pedigree diagram of this family:

 
     VCF file: /share/datasets/GeneticDiseaseFamilyVcfs/case1/case1.vcf
     PED file: /share/datasets/GeneticDiseaseFamilyVcfs/case1/example.ped
     config.ini file: /share/datasets/GeneticDiseaseFamilyVcfs/case1/config.ini

Use this command to run the analysis:

$ trioFilter.py config.ini

     result: /share/datasets/GeneticDiseaseFamilyVcfs/case1/example.recessive.csv
     
Only one candidate gene in the result list.

## Case 2. TAF family, which is an X-linked family. LID57247 LID57244 LID57243 LID57250 are affected brother, affected brother, unaffected father, unaffected mother, respectively.

Here is the pedigree diagram of this family:
 

     VCF file: /share/datasets/GeneticDiseaseFamilyVcfs/case2/case2.vcf
     PED file: /share/datasets/GeneticDiseaseFamilyVcfs/case2/example.ped
     config.ini file: /share/datasets/GeneticDiseaseFamilyVcfs/case2/config.ini
     result: /share/datasets/GeneticDiseaseFamilyVcfs/case2/example.recessive.csv

## Case 3. De novo mutation example with muscle weakness and difficult swallowing. ll, llF, llM are affected son, father and mother respectively.
Here is the pedigree diagram of this family:
 

     VCF file: /share/datasets/GeneticDiseaseFamilyVcfs/case3/case3.vcf
     PED file: /share/datasets/GeneticDiseaseFamilyVcfs/case3/example.ped
     config.ini file: /share/datasets/GeneticDiseaseFamilyVcfs/case3/config.ini
     result: /share/datasets/GeneticDiseaseFamilyVcfs/case3/example.denovo.csv

## Case 4. Dominant mutation for Gina family. 400 401 399 402 403 are unaffected father, affected mother, affected daughter, affected son, affected daughter, respectively.
Here is the pedigree diagram of this family:
 


TNXB and NAGA are in candidate gene list.
          VCF file: /share/datasets/GeneticDiseaseFamilyVcfs/case4/case4.vcf
          PED file: /share/datasets/GeneticDiseaseFamilyVcfs/case4/example.ped
          config.ini file: /share/datasets/GeneticDiseaseFamilyVcfs/case4/config.ini
          result: /share/datasets/GeneticDiseaseFamilyVcfs/case4/example.dominant.csv

## Case 5. Compound heterozygous variants in NAGLU in both brother and sister, who are affected with progressive cognitive decline.
Here is the pedigree diagram of this family:
 


    VCF file: /share/datasets/GeneticDiseaseFamilyVcfs/case5/ruinew2.extract_consensus.vcf and yongkangnew2.extract_consensus.vcf
    Reference: https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0066-9
   
   There had some errors when merge these two files together. I did the analysis of this family with the following steps.
   
### Step 1: using annovar to annotate each vcf file:

$ table_annovar.pl ruinew2.vcf /home/lulinhuang/bin/annovar/humandb/ -buildver hg19 -out ruinew2 -remove -protocol refGene,cytoBand,gnomad_exome,gnomad_genome,dbnsfp33a -operation g,r,f,f,f -nastring . –vcfinput
$ table_annovar.pl yongkangnew2.extract_consensus.vcf /home/lulinhuang/bin/annovar/humandb/ -buildver hg19 -out yongkangnew2 -remove -protocol refGene,cytoBand,gnomad_exome,gnomad_genome,dbnsfp33a -operation g,r,f,f,f -nastring . -vcfinput

### Step 2: extract the nonsynonymous variants in patient ruinew2:

$ grep -E 'Func.refGene|splicing|nonsynonymous|stop|frameshift' ruinew2.hg19_multianno.txt>ruinew2.txt 

Step 3: extract the rare nonsynonymous variants in patient ruinew2

$ awk '$13<0.01' ruinew2.txt > MAFruinew2_new.txt

### Step 4: extract the same rare nonsynonymous variants in patients ruinew2 and yongkangnew2 using a python script

$ python mfind.py MAFruinew2_new.txt 9 yongkangnew2.hg19_multianno.txt 9 > same.txt

import os, sys, string;
from collections import defaultdict

def getKeywords(f, ind2):
   mf = open(f, 'r')
   line = mf.readline();
   kwdict = defaultdict(str)
   while line:
      line = string.strip(line);
      lsp = line.split('\t');
      kwdict[lsp[ind2]] = 1;
      line = mf.readline();

   return kwdict;

if __name__=='__main__':
  if len(sys.argv)<5:
     print 'Usage: python', sys.argv[0], 'file1 ind1 file2 ind2'
     print '\t The two indexes is 0-based'
     print '\t The second index is to find keywords'
     print '\t The output will be filtered line from the first file'
     sys.exit(1);

  kwdict = getKeywords(sys.argv[3], int(sys.argv[4]));

  mf = open(sys.argv[1], 'r'); ind1 = int(sys.argv[2])
  line = mf.readline();
  while line:
     line = string.strip(line);
     lsp = line.split('\t');

     if lsp[ind1] in kwdict: print line;

     line = mf.readline();

### Step 5: extract the genes which had autosome recessive inherent model:
#### find out the genes which mutated at least twice for compound heterozygous mutations
cut -f7 same.txt | sort | uniq -c | sort >memesame.txt
awk '$1>1' memesame.txt>mimisame.txt
#### find out the homozygous mutations:
grep –E “1/1:” same.txt>samehom.txt

### Step 6: Using phenolyzer to find out the disease genes
Go to the webside: http://phenolyzer.wglab.org/, fill in the disease/phenotype and gene list like this:	
 
Then we can get the results. The NAGLU gene is the best candidate gene.
http://phenolyzer.wglab.org/done/44416/OHx3OoeCwlK_T1em/index.html
 
http://phenolyzer.wglab.org/done/44416/OHx3OoeCwlK_T1em/index.html
 


## Case 6. KBG syndrome due to de novo mutation in ANKRD11. You can perform single-case analysis using proband only, or perform family-based analysis
Here is the pedigree diagram of this family:
 
    VCF file:/share/datasets/GeneticDiseaseFamilyVcfs/case6/Unaffected_brother.vcf  Unaffected_father.vcf  Unaffected_mother.vcf  Unaffected_sister1.vcf  Unaffected_sister2.vcf  proband.vcf
    Reference: http://molecularcasestudies.cshlp.org/content/2/6/a001131.full

### Step 1: Using BCF tools to merge all these vcf files to a single vcf file family6.vcf.
### Step 2: Using annovar to do annotate of the VCF file.
$ table_annovar.pl family6.vcf /home/lulinhuang/bin/annovar/humandb/ -buildver hg19 -out family6 -remove -protocol refGene,cytoBand,gnomad_exome,gnomad_genome,dbnsfp33a -operation g,r,f,f,f -nastring . –vcfinput
### Step 3: extract the nonsynonymous variants of this family:
grep -E 'Func.refGene|splicing|nonsynonymous|stop|frameshift' family6.vcf.hg19_multianno.txt>family6misssense.txt
### Step 4: extract the rare nonsynonymous variants of this family:
awk '$24<0.0001' family6misssense.txt > MAFfamily6.txt
awk '$12<0.0001' MAFfamily6.txt > MAFfamily6_2.txt
### Step 5: Sort the MAFfamily6_2.txt file in the excel to filter the heterozygote mutations only exists in the proband. There are 233 variants.
### Step 6: Using phenolyzer to find out the disease genes

Go to the webside: http://phenolyzer.wglab.org/, using KBG syndrome disease/phenotype and the filtered gene list.

 Here are the results. The AMKRD11 is the best candidate gene.
http://phenolyzer.wglab.org/done/44497/kvNlqRSQF7OtHJD3/index.html
http://phenolyzer.wglab.org/done/44497/kvNlqRSQF7OtHJD3/index.html
 

 
