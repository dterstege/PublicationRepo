Analysis of WholeBrain Outputs
==============================

WholeBrainAnalysis.m was developed for the compilation and analysis of datasets collected using the Whole Brain software suite developed by Daniel Fürth and modified using [scripts](https://github.com/dterstege/PublicationRepo/tree/main/Terstege2022A/WholeBrain/R) and [plugins](https://github.com/dterstege/CavalieriPointMask) developed by Dylan Terstege.

Whole Brain is a neuroanatomical information system.  Neuroanatomical data obtained from microscope images is encoded and stored in stereotactic coordinate form within the Allen Mouse Brain Atlas. The software suite can be accessed at http://www.wholebrainsoftware.org

## Table of Contents

| Section  | Contents |
| ------------- | ------------- |
| [**1. Inputs**](#in) | How to organize your data for analysis |
| [**2. Operation**](#op) | How to run the analysis |
| [**3. Outputs**](#out) | How the analysis will present your data |
| [**4. Requirements**](#req) | What is required to run the analysis |

<a name="in"/>

## 1. Inputs

**Input Files**

After running the [*R* script](https://github.com/dterstege/PublicationRepo/tree/main/Terstege2022A/WholeBrain/R), users will have two output csv files generated from each image.  One of these files will have the suffix "cells.csv" while the other has the suffix "grids.csv".  

**File Organization**

Users should organize these output files so that each mouse has its own folder (referred to in the analysis as the "parent folder" for that mouse).  Within this folder, there should be two subfolders titled "cell" and "grid" (caps sensitive folder naming).  *Cell* should contain all "cells.csv" files for that animal, while *grid* contains all "grids.csv" files for that animal.

**Variables in MATLAB Script**

<a name="op"/>

## 2. Operation

<a name="out"/>

## 3. Outputs

<a name="req"/>

## 4. Requirements

Please ensure that the following packages and files are located in a MATLAB-accessible folder.  It is reccommended that the user copies these files to the MATLAB Path directory.

**Required Packages**:

- [fdr_bh.m](https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh)
- [degrees_und.m](https://sites.google.com/site/bctnet/)
- [targeteddeletioncurve.m](https://github.com/dterstege/TargetedNodeDeletionToolbox)
- [mcl.m](https://github.com/AndrasHartmann/MMCL)
- [deduce_mcl_clusters.m](https://github.com/AndrasHartmann/MMCL)
- [smallworldassessment.m](https://github.com/dterstege/SmallWorldAssessment)

**Required Files** ([included in this repository](https://github.com/dterstege/PublicationRepo/tree/main/Terstege2022A/WholeBrain/MATLAB/Atlases)):

- WB_Atlas.mat
- a user-defined atlas (ex: User_Atlas_98reg.mat)
