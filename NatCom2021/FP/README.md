## Fiber Photometry Analysis

The attached scripts were used in the fiber photometry analyses conducted during the described experiments.

Non-experimental sample data has been provided to demonstrate the required data organization.

### Analysis Information

Inputs: 
- Neurophotometrics Photometry Output (.csv)
- ANYmaze Outputs (.csv)

Analyses (data stored in the MATLAB structure element "FP"):
- dFF = delta F over F
- peaks = information about metrics related to peaks in the photometry traces

Operational Definitions:
- dFF: Isosbestic data is fit to a biexponential decay which is linearly scaled to the raw 470 signal data.
- Baseline Correction: Minimum dFF value from the test segment is set to the mean value from the user-defined baseline segment. All other values are scaled acordingly.
- Peaks: Peaks are defined as any point at which the dFF trace makes a downward deflection after rising to at least two standard deviations above the median dFF value 

### Notes to consider when running the analysis:

- Analyses conducted using MATLAB R2020a on macOS 10.15.7
- Please read the entire documentation prior to running the script and adjust lines as needed
- Run each section of the analysis scripts in order
- After running 'FP_Part1.m', organize data from all groups into a single excel workbook titled 'test_summary.xlsx'. See attached file for an exmple of how to organize this file 

## Contact Us

- **Dylan Terstege** (developped analyses) - ![twitter-icon_16x16](https://user-images.githubusercontent.com/44174532/113163958-e3d3e400-91fd-11eb-8d79-17906d8d3f25.png)[@dterstege](https://twitter.com/dterstege) - ![Mail](https://user-images.githubusercontent.com/44174532/113164412-50e77980-91fe-11eb-9282-dd83852578ce.png)
<dylan.terstege@ucalgary.ca>


Principal Investigator:
- Jonathan Epp (supervisor) - https://epplab.com
