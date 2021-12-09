Analysis of WholeBrain Outputs
==============================

WholeBrainAnalysis.m was developed for the compilation and analysis of datasets collected using the Whole Brain software suite developed by Daniel Fürth and modified using [scripts](https://github.com/dterstege/PublicationRepo/tree/main/Terstege2022A/WholeBrain/R) and [plugins](https://github.com/dterstege/CavalieriPointMask) developed by Dylan Terstege.

Whole Brain is a neuroanatomical information system.  Neuroanatomical data obtained from microscope images is encoded and stored in stereotactic coordinate form within the Allen Mouse Brain Atlas. The software suite can be accessed at http://www.wholebrainsoftware.org

Current Version: V 1.0.1



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

*Please note that the organization of the User Atlas is very important. When analyzing data using pairwise correlations, including a region which is missing from one of the subjects will have severe consequences of the resulting correlation. Atlas organization should be suited to only analyze regions present across all subjects in the analysis*

**File Organization**

Users should organize these output files so that each mouse has its own folder (referred to in the analysis as the "parent folder" for that mouse).  Within this folder, there should be two subfolders titled "cell" and "grid" (caps sensitive folder naming).  *Cell* should contain all "cells.csv" files for that animal, while *grid* contains all "grids.csv" files for that animal.

**Variables in MATLAB Script**

Users should read the MATLAB documentation in full prior to beginning the analysis.  At the beginning of each block of code in the script, there is a blurb outlining the required packages, giving insight into the analysis parameters, and providing general instructions.  This blurb also outlines the user input variables in each section.  These variables ask the user to input things from group identifiers, the number of mice per group, and individual mouse IDs, to technical details such as network thresholds.  Adjust these values to suit your experiment.

<a name="op"/>

## 2. Operation

It is recommended that users read the MATLAB documentation in full prior to beginning the analysis.  It is also recommended that users run each section independently, which allows users to utlize the *Save* and *Load* commands built into the last two sections of the code (*Sections X and Y*).  The loading of experimental data can be time consuming, depending on group sizes.  Therefore, it can be useful to run this section and then save the resulting structure array, in case that users would like to re-run the analysis using different parameters.

<a name="out"/>

## 3. Outputs

For group analyses, correlation matrices, adjacency matrices, and circle plots are all automatically generated while running *Section 03. Basic Network Analyses*.  All other data can be found nested in the *WB* structure element in the MATLAB workspace.  All "per animal" data will be organized in the same order as it was input.  This order is stored under the user input ID structure elements for later reference.

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

**Version Requirements**

This analysis was developped using MATLAB R2020a
