We want to make input PacBio_amplicon files with many templates present in pan-cov libs

Input file (here, called input_seqs.csv) is a table of the insert sequences, the target names, and the insert length

First, we copied input_seqs.csv to processed_seqs.csv, and opened it up in BBEdit for initial processing.

We place a "template" PacBio amplicon (PacBio_amplicon_MERS.gb) at the end of each line:

find:

```
\r
```

replace:

```
\rLOCUS       MERS                       1121 bp    DNA        linear       12-DEC-2022\nSOURCE      Tyler Starr\n  ORGANISM  .\nFEATURES             Location/Qualifiers\n     termini5        1..33\n                     /locus_tag="termini5"\n                     /label="termini5"\n     spacer          679..1099\n                     /locus_tag="spacer"\n                     /label="spacer"\n     barcode         1100..1115\n                     /locus_tag="barcode"\n                     /label="barcode"\n     termini3        1116..1121\n                     /locus_tag="termini3"\n                     /label="termini3"\n     gene            34..678\n                     /locus_tag="gene"\n                     /label="gene"\nORIGIN\n        1 ggccgcggag gcggagggtc ggctagccat atgCAAGCGG AAGGGGTCGA GTGCGATTTT\n       61 AGTCCACTGT TATCCGGCAC ACCGCCTCAA GTCTATAATT TTAAGAGACT GGTTTTCACA\n      121 AACTGCAATT ATAATCTAAC GAAATTGTTA AGCTTATTTA GCGTCAATGA TTTCACCTGC\n      181 TCTCAAATCT CTCCCGCCGC CATAGCAAGC AACTGTTACA GCAGTTTGAT TCTTGATTAC\n      241 TTCTCTTACC CATTATCTAT GAAATCTGAC CTGTCTGTAA GCTCAGCAGG GCCGATCTCC\n      301 CAATTTAACT ATAAGCAATC CTTTTCAAAT CCAACCTGCC TTATTCTGGC CACCGTGCCT\n      361 CACAATCTAA CCACGATAAC GAAGCCCTTA AAATATAGCT ACATCAATAA ATGCTCCCGT\n      421 CTTCTTTCTG ACGACAGGAC CGAGGTGCCT CAGTTAGTTA ATGCCAACCA GTATTCCCCA\n      481 TGCGTATCCA TCGTTCCGAG TACGGTATGG GAGGATGGTG ATTACTATAG AAAGCAGTTA\n      541 AGCCCGTTGG AAGGCGGCGG GTGGTTGGTT GCCAGCGGTT CAACAGTTGC TATGACAGAG\n      601 CAGCTTCAAA TGGGCTTCGG AATTACAGTC CAGTATGGCA CTGATACCAA TTCTGTTTGC\n      661 CCAAAACTTG AGTTCGCTct cgaggggggc ggttccgaac aaaagcttat ttctgaagag\n      721 gacttgtaat agagatctga taacaacagt gtagatgtaa caaaatcgac tttgttccca\n      781 ctgtactttt agctcgtaca aaatacaata tacttttcat ttctccgtaa acaacatgtt\n      841 ttcccatgta atatcctttt ctatttttcg ttccgttacc aactttacac atactttata\n      901 tagctattca cttctataca ctaaaaaact aagacaattt taattttgct gcctgccata\n      961 tttcaatttg ttataaattc ctataattta tcctattagt agctaaaaaa agatgaatgt\n     1021 gaatcgaatc ctaagagaat taatgatacg gcgaccaccg agatctacac tctttcccta\n     1081 cacgacgctc ttccgatctN NNNNNNNNNN NNNNNgcggc c \n//\r
```

We now have our file with the "replacement" name, sequence, and sequence length ready to be pulled with wildcard matches into the same template replacement genbank formatted file.


Next, replace the template's gene sequence with the new insert sequecne

find:
```
^([^,\r]+\,\d+)\,(\w+)(\rLOCUS       MERS                       1121 bp    DNA        linear       12-DEC-2022\nSOURCE      Tyler Starr\n  ORGANISM  .\nFEATURES             Location/Qualifiers\n     termini5        1..33\n                     /locus_tag="termini5"\n                     /label="termini5"\n     spacer          679..1099\n                     /locus_tag="spacer"\n                     /label="spacer"\n     barcode         1100..1115\n                     /locus_tag="barcode"\n                     /label="barcode"\n     termini3        1116..1121\n                     /locus_tag="termini3"\n                     /label="termini3"\n     gene            34..678\n                     /locus_tag="gene"\n                     /label="gene"\nORIGIN\n        1 ggccgcggag gcggagggtc ggctagccat atg)CAAGCGG AAGGGGTCGA GTGCGATTTT\n       61 AGTCCACTGT TATCCGGCAC ACCGCCTCAA GTCTATAATT TTAAGAGACT GGTTTTCACA\n      121 AACTGCAATT ATAATCTAAC GAAATTGTTA AGCTTATTTA GCGTCAATGA TTTCACCTGC\n      181 TCTCAAATCT CTCCCGCCGC CATAGCAAGC AACTGTTACA GCAGTTTGAT TCTTGATTAC\n      241 TTCTCTTACC CATTATCTAT GAAATCTGAC CTGTCTGTAA GCTCAGCAGG GCCGATCTCC\n      301 CAATTTAACT ATAAGCAATC CTTTTCAAAT CCAACCTGCC TTATTCTGGC CACCGTGCCT\n      361 CACAATCTAA CCACGATAAC GAAGCCCTTA AAATATAGCT ACATCAATAA ATGCTCCCGT\n      421 CTTCTTTCTG ACGACAGGAC CGAGGTGCCT CAGTTAGTTA ATGCCAACCA GTATTCCCCA\n      481 TGCGTATCCA TCGTTCCGAG TACGGTATGG GAGGATGGTG ATTACTATAG AAAGCAGTTA\n      541 AGCCCGTTGG AAGGCGGCGG GTGGTTGGTT GCCAGCGGTT CAACAGTTGC TATGACAGAG\n      601 CAGCTTCAAA TGGGCTTCGG AATTACAGTC CAGTATGGCA CTGATACCAA TTCTGTTTGC\n      661 CCAAAACTTG AGTTCGCT
```

replace:
```
\1\3\2
```


Next, let's wrap the gb formatted text into a third element of a true "CSV" for further scripting:

find:

```
\r^LOCUS
```

replace:


```
,"LOCUS
```


find:

```
//
```

replace:


```
//"
```

Then, we replaced the length column to be not the length of the newly inserted sequence, but rather the difference in length between the template MERS sequence (645 bp) and the replacement length, i.e. the number that we need to subtract from higher numbers in the map to reset the number indexing after we replaced the insert sequence.

Next, with the lovely help of ChatGPT, we generated Python code to parse the genbank format column element of the csv while replacing numbers and names as needed. This code is in generate_amplicons_from_csv.ipynb. The output of this is processed_genbank_maps.gb, which we copied into the main `../data` folder to use as our input PacBio amplicons.
