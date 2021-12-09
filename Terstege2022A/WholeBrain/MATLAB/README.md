Analysis of WholeBrain Outputs
==============================

WholeBrainAnalysis.m was developed for the compilation and analysis of datasets collected using the Whole Brain software suite developed by Daniel Furth (http://www.wholebrainsoftware.org) and incorporates a modification using scripts and plugins (https://github.com/dterstege/CavalieriPointMask) developed by Dylan Terstege.

Whole Brain is a neuroanatomical information system.  Neuroanatomical data obtained from microscope images is encoded and stored in stereotactic coordinate form within the Allen Mouse Brain Atlas. The software suite can be accessed at http://www.wholebrainsoftware.org

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
