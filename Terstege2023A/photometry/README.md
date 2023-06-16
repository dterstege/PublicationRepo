## Fiber Photometry Analysis

The attached scripts were used in the fiber photometry analyses conducted during the described experiments.

### Analysis Information

Inputs: 
- Neurophotometrics Photometry Output (.csv)
- ANYmaze Outputs (.csv)

Analyses (data stored in the MATLAB structure element "FP"):
- dFF = delta F over F
- AUC = area under the curve

Operational Definitions:
- dFF: Isosbestic data is fit to a biexponential decay which is linearly scaled to the raw 470 signal data.

### Notes to consider when running the analysis:

- Analyses conducted using MATLAB R2020a on macOS 10.15.7
- Please read the entire documentation prior to running the script and adjust lines as needed
- Run each section of the analysis scripts in order

## Contact Us

- **Dylan Terstege** (Developed Analyses) - ![twitter-icon_16x16](https://user-images.githubusercontent.com/44174532/113163958-e3d3e400-91fd-11eb-8d79-17906d8d3f25.png)[@dterstege](https://twitter.com/dterstege) - ![Mail](https://user-images.githubusercontent.com/44174532/113164412-50e77980-91fe-11eb-9282-dd83852578ce.png)
<dylan.terstege@ucalgary.ca>

- Jonathan Epp (Principal Investigator / Corresponding Author) - https://epplab.com
