Photometry
----------

Analyses used with the fiber photometry experiments

It is recommended that users read through the documentation in the first section of each script in full prior to running the analyses.  These sections have been copied and pasted below.

From SarginFP_ANYmaze_v111.m:

~~~
%                       *GENERAL INFORMATION*
%
%   Behavioural Paradigm: Align photometry data with ANYmaze binary outputs
%
%   Run block-wise
%   Jump to 2nd last section (Section X) to save data structure at any time
%   Load previously saved data structures using last section (Section Y)
%
%                        *ANALYSIS INFORMATION*
%
%   Inputs:
%       Doric Fiber Photometry Output (.csv)
%       ANYmaze Outputs (.csv)
%
%   Analyses:
%   (Data stored in MATLAB structure "FP" under "outputs")
%       dFF = delta F over F
%           Full trial info
%           Also separated into behavioural epochs
%       peaks = information about peak frequency and mean peak height
%           Full trial info
%           Also separated into behavioural epochs
%       Analyses are sorted into "behaviourdependent" and
%       "behaviourindependent", with independent being over the course of
%       the entire trial and dependent being linked to the SimBA outputs
%       
%       Option to output data in z-score format.
%           Choosing this option will make all "dFF" variables in the
%           output structure a z-scored dF/F.  The z-score is calculated
%           from non-baseline corrected dF/F, using the mean and standard
%           deviation from the entire trace
%
%   Operational Definitions:
%       dFF:
%           Isosbestic data is fit to a biexponential decay
%           This biexponential decay is linearly scaled to the raw 470 data
%           (raw = sig, fitted = sigfit)
%           dFF = (sig-sigfit)./sigfit
%           Baseline correction:
%               Trial is split into test and baseline segments
%               Mean dFF value during middle portion of baseline is
%               calculated (X minutes)
%               Minimum dFF value from the test segment is calculated
%               The difference between these two values is then added to
%               every dFF value during the test segment
%               This effectively raises the minimum dFF value from the test
%               segment to the mean baseline segment value
%       peaks:
%           May require the addition of a lowpass filter depending on the
%           dataset
%           Peaks defined as any point at which the dFF trace makes a
%           downward deflection after rising
%           Peak detection restricted to peaks exceeding a height of two
%           standard deviations above the median dFF value
%           Mean peak height is relative to 0
%           Behaviour-specific peak frequencies are normalized by duration
%           of behavioural epoch
%
%                           *ASSUMPTIONS*
%                      *IMPORTANT - PLEASE READ*
%
%   Like most analyses, this code was developed with a number of
%   assumptions about the organization of data files. The intention was to
%   make this as flexible as possible in terms of the number of groups,
%   group sizes, and input. However, it makes the assumption that all
%   photometry inputs are of the same length. If this is not the case,
%   please manually adjust the total number of frames in the recording by
%   trimming the few frames from the end which account for this.
%   There are a number of instances in which the number of frames for
%   different recordings (photometry and behaviour) and hard-written into
%   the code. Each of these instances are documented and can be found by
%   searching (control+F) this file with the search term ''. Be sure
%   to adjust these values as needed to suit the analysis.
%
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.1
%
%   V1.0.1 - Basic Behavioural Alignment. Assumes consistent framerate
%   V1.1.1 - Integration of a pickup window for excluding experimenter-introduced noise;
%
%   Developed for the Sargin Lab, University of Calgary
%   (https://sarginlab.com)
%   Datasets used during script development were collected by Matt Dawson
%
%   Created 08/18/2021 Dylan Terstege (https://github.com/dterstege)
%   Epp Lab, University of Calgary (https://epplab.com)
%   Contact: dylan.terstege@ucalgary.ca
~~~

From SarginFP_SimBA_v163.m
~~~
%                       *GENERAL INFORMATION*
%
%   Behavioural Paradigm: Social Interaction - Familiar vs Intruder
%       Design: photometry recording from resident mouse during a
%       "baseline" period. A familiar mouse is then introduced. During this
%       time, anogenital and head/torso interactions are recorded using
%       SimBA (https://github.com/sgoldenlab/simba). After a period of time
%       the familiar mouse is removed and the resident mouse is allowed
%       time to settle in the apparatus. Minutes later, an unfamiliar
%       intruder mouse is introduced. Once again, anogenital and head/torso
%       interactions are recorded during this time. The intruder mouse is
%       then removed and the resident mouse is allowed to settle once more
%       before the conclusion of the test
%
%   Run block-wise
%   Jump to 2nd last section (Section X) to save data structure at any time
%   Load previously saved data structures using last section (Section Y)
%
%                        *ANALYSIS INFORMATION*
%
%   Inputs:
%       Doric Fiber Photometry Output (.csv)
%       SimBA Outputs:
%           Mouse interacting with familiar mouse (.csv)
%           Mouse interacting with intruder mouse (.csv)
%               Input the number of frames trimmed from the beginning of
%               the associated video files
%
%   Analyses:
%   (Data stored in MATLAB structure "FP" under "outputs")
%       dFF = delta F over F
%           Full trial info
%           Also separated into behavioural epochs
%       peaks = information about peak frequency and mean peak height
%           Full trial info
%           Also separated into behavioural epochs
%       Analyses are sorted into "behaviourdependent" and
%       "behaviourindependent", with independent being over the course of
%       the entire trial and dependent being linked to the SimBA outputs
%       
%       Option to output data in z-score format.
%           Choosing this option will make all "dFF" variables in the
%           output structure a z-scored dF/F.  The z-score is calculated
%           from non-baseline corrected dF/F, using the mean and standard
%           deviation from the entire trace
%
%   Operational Definitions:
%       dFF:
%           Isosbestic data is fit to a biexponential decay
%           This biexponential decay is linearly scaled to the raw 470 data
%           (raw = sig, fitted = sigfit)
%           dFF = (sig-sigfit)./sigfit
%           Baseline correction:
%               Trial is split into test and baseline segments
%               Mean dFF value during middle portion of baseline is
%               calculated (X minutes)
%               Minimum dFF value from the test segment is calculated
%               The difference between these two values is then added to
%               every dFF value during the test segment
%               This effectively raises the minimum dFF value from the test
%               segment to the mean baseline segment value
%       peaks:
%           May require the addition of a lowpass filter depending on the
%           dataset
%           Peaks defined as any point at which the dFF trace makes a
%           downward deflection after rising
%           Peak detection restricted to peaks exceeding a height of two
%           standard deviations above the median dFF value
%           Mean peak height is relative to 0
%           Behaviour-specific peak frequencies are normalized by duration
%           of behavioural epoch
%
%                           *ASSUMPTIONS*
%                      *IMPORTANT - PLEASE READ*
%
%   Like most analyses, this code was developed with a number of
%   assumptions about the organization of data files. The intention was to
%   make this as flexible as possible in terms of the number of groups,
%   group sizes, and input. However, it makes the assumption that all
%   photometry inputs are of the same length. If this is not the case,
%   please manually adjust the total number of frames in the recording by
%   trimming the few frames from the end which account for this.
%   There are a number of instances in which the number of frames for
%   different recordings (photometry and behaviour) and hard-written into
%   the code. Each of these instances are documented and can be found by
%   searching (control+F) this file with the search term '#frames'. Be sure
%   to adjust these values as needed to suit the analysis.
%
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.6.3
%
%   V1.1 - Basic Behavioural Alignment. Assumes consistent framerate
%   1.1.1 - Z-score integration; bout-by-bout information;
%   behaviour-dependent dF/F saved to variable; behaviour stretch has been
%   fixed
%   1.1.2 - mean peak amplitude defaults to the baseline-corrected value;
%   %dF/F output option added
%   1.1.3 - timestamp of peak amplitude during each segment in which peaks
%   are detected; option to define time ranges of interest added to
%   beahviour-independent section (Z)
%
%   V1.2 - Utilizes MatLab Image Processing Toolkit
%   (https://www.mathworks.com/products/image.html)
%   V1.2.1 - alignment using nearest neighbour interpolation; no longer asking user
%   for number of frames trimmed -- number of total frames is assumed by
%   fixed framerate and the number of dummy frames to be added to the
%   beginning is calculated based on the discrepency between the total
%   number of frames imported and the expected number of frames
%   
%   V1.3 - Back to previous re-alignment style; total number of frames
%   still assumed by fixed framerate and the number dummy frames added to
%   the beginning is calculated based on the discrepency between the total
%   number of frames imported and the expected number of frames; fixed
%   "extra data" issues
%
%   V1.4 - Rolling Window Peak Detection
%
%   V1.5 - Trace Analyses during Epoch, independent of behaviour
%   For simplicity, this is nested within the "Behaviour Dependent" section
%   of the analysis (simply because we've already defined the epochs during
%   this section)
%
%   Developed for the Sargin Lab, University of Calgary
%   (https://sarginlab.com)
%   Datasets used during script development were collected by Matt Dawson
%
%   Created 06/03/2021 Dylan Terstege (https://github.com/dterstege)
%   Epp Lab, University of Calgary (https://epplab.com)
%   Contact: dylan.terstege@ucalgary.ca
~~~
