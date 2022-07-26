# NASA_GeneLab_Summer2022
Work related to NASA GeneLab summer 2022 internship with Mike Lee on amplicon and metagenomics sequencing.

https://github.com/AstrobioMike/GL-2022-summer-internship/wiki

---

## HackMD Pages

### Week 2
- First practice page - https://hackmd.io/sH_CnHPbQRm6kt_vA3KvSg
- Amplicon Processing Workflow Explained - https://hackmd.io/_fPaZACSRA2z9_pkvy75IA
- Metagenomics Proccesing Workflow Explained - https://hackmd.io/yLMrIK4URyOtiS1apGMzxQ

### Week 3
- Amplicon/Metagenomics Processing Tutorial Notes/Qs - https://hackmd.io/GFIAAk19TNKSZPchx3sDpQ
- Snakemake Command-line Code for Test Dataset - https://hackmd.io/P_AZynOaR06bjm14A04noA
- GeneLab GLDS-286 Metagenomic Processing - https://hackmd.io/-rMzenqzTPOTr17dIqgWQg

### Week 4 
- Current Metagenomics V&V Code Notes - https://hackmd.io/xm-gBK5eRS2DOE_jJsLJyA

---

## **Amplicon/Metagenomics Code Instructions**

### **Installation**

- Install conda

--> If on Mac: 
```
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
--> If on Windows: (in WSL environment)
```
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
--> Then, install using bash command
```
bash Miniconda3-latest-*.sh
```
(For more conda information, visit Mike's [conda intro page](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda))

- Install [dp_tools](https://github.com/AstrobioMike/GL-2022-summer-internship/wiki/Working-towards-Jonathan's-validation-structure)
```
curl -L -o dp_tools-condaEnv.yaml https://raw.githubusercontent.com/J-81/dp_tools/main/condaEnv.yaml
conda env create -f dp_tools-condaEnv.yaml
```

- Download test data to run on

--> Metagenomics test data:
```
curl -L -o GL-metagenomics-output-for-validation-testing.tar.gz https://figshare.com/ndownloader/files/36039989
tar -xzvf GL-metagenomics-output-for-validation-testing.tar.gz
cd GL-metagenomics-output-for-validation-testing/
```
--> Amplicon test data:
```
curl -L -o GL-amplicon-output-for-validation-testing.tar.gz https://figshare.com/ndownloader/files/36252156
tar -xzvf GL-amplicon-output-for-validation-testing.tar.gz
cd GL-amplicon-output-for-validation-testing/
```

### **To Run Metagenomics Code** - follow [Rose Carion's instructions](https://github.com/rosecarion/GL_2022_Internship) (also re-written below)
- Download "metagenomics_user_input_validation.py"

- Make sure dp_tools is activated and you've cd'ed into appropriate directory above
```
conda activate dp_tools
```

- When running "metagenomics_user_input_validation.py", user inputs are required:

1. GLDS ID (e.g. 'GLDS-276')
2. Path to the sample names file (e.g. '/home/inikabhatia/GeneLab-Internship/GL-metagenomics-output-for-validation-testing/unique-sample-IDs.txt')
3. Expected additional filename prefix that was added to the files that describe multiple samples, default is a back slash (Press 'return' key if none)
4. Specify whether the test data is single-ended ('y' for 'yes', 'n' for 'no')

### **To Run Amplicon Code** - follow [Rose Carion's instructions](https://github.com/rosecarion/GL_2022_Internship) (also re-written below)
- Download "amplicon_user_input_validation.py"

- Make sure dp_tools is activated and you've cd'ed into appropriate directory above
```
conda activate dp_tools
```

- When running "amplicon_user_input_validation.py", user inputs are required:

1. GLDS ID (e.g. 'GLDS-276')
2. Path to the sample names file (e.g. '/home/inikabhatia/GeneLab-Internship/GL-metagenomics-output-for-validation-testing/unique-sample-IDs.txt')
3. Output file prefix if there is one (Press 'return' key if none)
4. Specify whether primers trimmed prior to GeneLab processing ('y' for 'yes', 'n' for 'no')
5. Specify whether the test data is single-ended ('y' for 'yes', 'n' for 'no')
