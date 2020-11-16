##code for manuscript "Evolution of late-stage metastatic melanoma is dominated by tetraploidization and aneuploidy"##

The code in this repository provides all steps necessary to reproduce the neoantigenic mutation component done in Vergara et al, "Evolution of late-stage metastatic melanoma is dominated by tetraploidization and aneuploidy".  

Patient SK-G will be used as an example to generate neoepitope work.

1) Run CooVar (PMID: 23116482) to generate wildtype and mutated protein sequences for each lesion
-------------------------------------------------------------------------------------------------

perl /home/ivergara/tools/coovar-0.07/coovar.pl -e /home/ivergara/reference/Homo_sapiens.GRCh37.75.mainChrom.gtf -r /home/ivergara/human_g1k_v37_mainChrom.fasta -t DATA/CALLS/SK-G_SRR2159462.coovar -o SK-G_SRR2159462_OUT --circos

perl /home/ivergara/tools/coovar-0.07/coovar.pl -e /home/ivergara/reference/Homo_sapiens.GRCh37.75.mainChrom.gtf -r /home/ivergara/human_g1k_v37_mainChrom.fasta -t DATA/CALLS/SK-G_SRR2159463.coovar -o SK-G_SRR2159463_OUT --circos

perl /home/ivergara/tools/coovar-0.07/coovar.pl -e /home/ivergara/reference/Homo_sapiens.GRCh37.75.mainChrom.gtf -r /home/ivergara/human_g1k_v37_mainChrom.fasta -t DATA/CALLS/SK-G_SRR2159465.coovar -o SK-G_SRR2159465_OUT --circos


2) Create folders for patient
-----------------------------

mkdir EPITOPES_17
mkdir EPITOPES_17/SK-G
mkdir EPITOPES_REF_17
mkdir EPITOPES_REF_17/SK-G

3) Extract 17-mers from mutated and wild-type protein sequences for each lesion
-------------------------------------------------------------------------------

 perl generate-epitope-sequence.pl SK-G_SRR2159462_OUT/categorized-gvs.gvf SK-G_SRR2159462_OUT/transcripts/variant_peptides.fasta  DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_17/SK-G/SK-G_SRR2159462_epitopes.out
perl generate-epitope-sequence.pl SK-G_SRR2159462_OUT/categorized-gvs.gvf SK-G_SRR2159462_OUT/transcripts/reference_peptides.fasta DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_REF_17/SK-G/SK-G_SRR2159462_epitopes.out
 
perl generate-epitope-sequence.pl SK-G_SRR2159463_OUT/categorized-gvs.gvf SK-G_SRR2159463_OUT/transcripts/variant_peptides.fasta  DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_17/SK-G/SK-G_SRR2159463_epitopes.out
perl generate-epitope-sequence.pl SK-G_SRR2159463_OUT/categorized-gvs.gvf SK-G_SRR2159463_OUT/transcripts/reference_peptides.fasta DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_REF_17/SK-G/SK-G_SRR2159463_epitopes.out
 
perl generate-epitope-sequence.pl SK-G_SRR2159465_OUT/categorized-gvs.gvf SK-G_SRR2159465_OUT/transcripts/variant_peptides.fasta  DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_17/SK-G/SK-G_SRR2159465_epitopes.out
perl generate-epitope-sequence.pl SK-G_SRR2159465_OUT/categorized-gvs.gvf SK-G_SRR2159465_OUT/transcripts/reference_peptides.fasta DATA/TSV/NS/SK-G.all.ns.tsv 17 > EPITOPES_REF_17/SK-G/SK-G_SRR2159465_epitopes.out
 
4) Remove redundant epitopes before running netMHC4.0
-----------------------------------------------------

perl remove-redundancy.pl EPITOPES_17/SK-G
perl remove-redundancy.pl EPITOPES_REF_17/SK-G 

Expected content within folders EPITOPES_17 and EPITOPES_REF_17 can be found within OUTPUT.


5) Run netMHC on mutated and wildtype peptides
----------------------------------------------

/home/ivergara/tools/netMHC-4.0/netMHC -a HLA-A0201,HLA-A0301,HLA-B1402,HLA-B0702,HLA-C0702,HLA-C0802 -l 9 -tdir SK-G_I_17_var_tmp -f ./EPITOPES_17/SK-G/input_netMHC_noRedundancy.fsa > ./EPITOPES_17/SK-G/I_input_netMHC_noRedundancy.out

/home/ivergara/tools/netMHC-4.0/netMHC -a HLA-A0201,HLA-A0301,HLA-B1402,HLA-B0702,HLA-C0702,HLA-C0802 -l 9 -tdir SK-G_I_17_ref_tmp -f ./EPITOPES_REF_17/SK-G/input_netMHC_noRedundancy.fsa > ./EPITOPES_REF_17/SK-G/I_input_netMHC_noRedundancy.out 

6) Parse netMHC output
----------------------

perl parse-netMHC-classI-output.pl EPITOPES_17/SK-G/I_input_netMHC_noRedundancy.out >  EPITOPES_17/SK-G/I_input_netMHC_noRedundancy.tab
perl parse-netMHC-classI-output.pl EPITOPES_REF_17/SK-G/I_input_netMHC_noRedundancy.out >  EPITOPES_REF_17/SK-G/I_input_netMHC_noRedundancy.tab

7) Build neoepitope call matrix for patient
-------------------------------------------

perl build-matrix-epitopes-netmhc.pl SK-G 17 > SK-G_neoepitopes_matrix_17.txt

The expected final output file can be found within OUTPUT.


