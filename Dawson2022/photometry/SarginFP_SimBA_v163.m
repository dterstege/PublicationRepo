%%                          SarginFP_SimBA.m
%
%   Fiber Photometry Analysis
%   
%   It is recommended that users read through the documentation in full
%   prior to using this analysis. This applies to each section in the code
%
%   All outputs are accessed through the structure element "FP" in the
%   MATLAB workspace
%
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

disp('Running SarginFP_SimBA');

%%  01. Initialization - Group Information
%   Default information for each input is listed below
%
%   Group Information:
%       Input number of groups
%           default: groupnum = 2
%       Input group identifiers
%           default: groups = ["control" "manipulation"]
%       Input number of animals per group
%           default: pergroup = [10 10]
%
%   Analysis Information:
%       Z-Score Use:
%           0 if not necessary
%           1 if necessary
%           default: zscore = 0
%       %dF/F Use:
%           0 if not necessary
%           1 if necessary
%           default: percentdFF = 0
%       NOTE: ONLY SELECT ONE OF THESE OPTIONS

%clear structure to start fresh
clear FP

%input number of groups
FP.info.groupnum = 1;
%input group identifiers
FP.info.groups = ["GH"]; %manually adjust
%input number of animals per group
FP.info.pergroup = [10]; %manually adjust
%indicate whether z scores are wanted
FP.info.zscore = 0;
%indicate whether %dF/F is wanted
FP.info.percentdFF = 0;


%peak detection rolling window time in seconds
FP.info.windowtime = 30;

%process successful
disp('Process 1. Initialization - Group Data Complete');

%%  02. Initialization - Subject Information
%   Somewhat time intensive
%   Requires manual selection of files
%   Default information for each input is listed below
%
%   Variables to be adjusted:
%       Number of frames for each section (see descriptions)
%       Individual Information:
%           Identify animal ID
%           Identify Doric output file
%           Identify SimBA output files
%           Number of frames in the dataset

%internal initialization
for ii = 1:FP.info.groupnum(1)
    %create structure arrays for raw doric timestamp inputs
    %adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
    FP.inputs.doric.time.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw doric calcium signal inputs
    %adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
    FP.inputs.doric.sig.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw doric isosbestic inputs
    %adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
    FP.inputs.doric.iso.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw simba - familiar  inputs
    %adjust the number of frames (default = 6000, aka 5mins @20fps) as needed
    FP.inputs.familiar.headtorso.(FP.info.groups(ii))=zeros(6000,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw simba - intruder  inputs
    %adjust the number of frames (default = 6000, aka 5mins @20fps) as needed
    FP.inputs.intruder.headtorso.(FP.info.groups(ii))=zeros(6000,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw simba - familiar  inputs
    %adjust the number of frames (default = 6000, aka 5mins @20fps) as needed
    FP.inputs.familiar.anogenital.(FP.info.groups(ii))=zeros(6000,FP.info.pergroup(ii)); %#frames
    %create structure arrays for raw simba - intruder  inputs
    %adjust the number of frames (default = 6000, aka 5mins @20fps) as needed
    FP.inputs.intruder.anogenital.(FP.info.groups(ii))=zeros(6000,FP.info.pergroup(ii)); %#frames
    %create info table to indicate whether behaviour occurs in the order of
    %intruder-familiar or familiar-intruder
    FP.info.order.(FP.info.groups(ii))=zeros(FP.info.pergroup(ii),1);
end

%import and organize raw data for first group
for ii=1:FP.info.groupnum(1)
    for iii=1:FP.info.pergroup(ii)
        %import raw data
        prompt=strcat("Select Doric file for mouse #",num2str(iii)," in ",FP.info.groups(ii)," group:");
        disp(prompt);
        [doricfile,pathname]=uigetfile(('*.csv'),prompt); %load doric file
        doricdata=readmatrix(strcat(pathname,doricfile));
        %trim if exceeds maximum trial duration (default = 180705, aka 25mins
        %@120fps)
        if size(doricdata,1)>180705 %#frames
            doricdata=doricdata(1:180705,:); %#frames
        else
            %do nothing
        end
        prompt=strcat("Select SimBA - Familiar file for mouse #",num2str(iii)," in ",FP.info.groups(ii)," group:");
        disp(prompt);
        [familiarfile,pathname]=uigetfile(('*.csv'), prompt); %load simba familiar file
        familiardata=readmatrix(strcat(pathname,familiarfile));
        prompt=strcat("Select SimBA - Intruder file for mouse #",num2str(iii)," in ",FP.info.groups(ii)," group:");
        disp(prompt);
        [intruderfile,pathname]=uigetfile(('*.csv'), prompt); %load simba intruder file    
        intruderdata=readmatrix(strcat(pathname,intruderfile));
        prompt="Was the first mouse familiar (1) or an intruder (0)?: ";
        firstmouse=input(prompt);
        if firstmouse==1
            FP.info.order.(FP.info.groups(ii))(iii)=1;
        else
            %do nothing
        end
        
        famtrim=6000-size(familiardata,1); %#frames
        inttrim=6000-size(intruderdata,1); %#frames        
        
        %organize raw data
        %adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
        FP.inputs.doric.time.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,1); %timestamps from doric output
        FP.inputs.doric.sig.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,2); %calcium-dependent signal from doric output
        FP.inputs.doric.iso.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,4); %isosbestic sigal from doric output
        %adjust the number of frames (default = 6000, aka 5mins @20fps) as needed
        FP.inputs.familiar.headtorso.(FP.info.groups(ii))((famtrim+1):6000,iii)=familiardata(1:(6000-famtrim),501); %resident_headtorso column %#frames
        FP.inputs.familiar.anogenital.(FP.info.groups(ii))((famtrim+1):6000,iii)=familiardata(1:(6000-famtrim),499); %resident_anogenital column %#frames
        FP.inputs.intruder.headtorso.(FP.info.groups(ii))((inttrim+1):6000,iii)=intruderdata(1:(6000-inttrim),501); %resident_headtorso column %#frames
        FP.inputs.intruder.anogenital.(FP.info.groups(ii))((inttrim+1):6000,iii)=intruderdata(1:(6000-inttrim),499); %resident_anogenital column %#frames
        
    end
end


%clear remaining unnecessary variables
clear ii iii prompt pathname doricfile familiarfile intruderfile
clear doricdata familiardata intruderdata famtrim inttrim firstmouse

%process successful
disp('Process 2. Initialization - Subject Information Complete');


%%  3. Photometry Processing - Behaviour Independent
%   The generation of delta F/F data from the photometry inputs
%   
%   User-inputted baseline period
%       default: baseline = [60 180]
%           first value is the number of seconds to trim from the start of
%           the trial, to minimize the influence of handling stress on the
%           baseline value
%           second value is the duration of the baseline period in seconds
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

%initialization
FP.info.baseline = [60 180];  %[start duration]

%batch processing
for ii=1:FP.info.groupnum(1)
            
        %create structure arrays for photometry outputs
        %adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
        FP.outputs.behaviourindependent.traces.dFF_notcorrected.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourindependent.analyses.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourindependent.analyses.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourindependent.analyses.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourindependent.analyses.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourindependent.analyses.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        
    for iii=1:FP.info.pergroup(ii)
        %full trace
        time=FP.inputs.doric.time.(FP.info.groups(ii))(:,iii);
        sig=FP.inputs.doric.sig.(FP.info.groups(ii))(:,iii);
        iso=FP.inputs.doric.iso.(FP.info.groups(ii))(:,iii);
        temp_x=1:length(iso);
        temp_x=temp_x'; %transpose data frame
        isofit=fit(temp_x,iso,'exp2'); %fit isosbestic data with a biexponential decay
        sigfit=robustfit(isofit(temp_x),sig); %linearly scale biexponential decay to the raw 470 data using robust fit
        sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);
        dFF=(sig-sigfit)./sigfit;   %delta F/F
        
        %add framerate to info
        FP.info.samplingrate=time(2)-time(1);
        
        %rolling window frame number
        FP.info.windowframes=round(FP.info.windowtime/FP.info.samplingrate);
        
        %baseline adjustment
        [~,idx1]=min(abs(time-(FP.info.baseline(1)))); %identify beginning of baseline
        [~,idx2]=min(abs(time-((FP.info.baseline(2)+FP.info.baseline(1))))); %identify end of baseline
        baseline=dFF(idx1:idx2); %baseline period
        test=dFF(idx2:end); %test period
        meanbaseline=mean(baseline); %calculate mean baseline dFF
        mintest=min(test); %identify minimum dFF during the test segment
        adjust=diff([mintest,meanbaseline]); %calculate the difference between baseline mean and test minimum
        cordFF=dFF+abs(adjust); %add the absolute value of difference to every point in the trace
        
        %z-score
        %if selected, z-score will be applied to the analysis
        %variables cordFF and dFF will no longer represent raw dF/F traces
        %and will instead be zscores without baseline correction
        if FP.info.zscore(1) == 1
            cordFF=zscore(dFF);
            dFF=zscore(dFF);
        else
            %proceed with non-z-scored data
        end
        
        % %dF/F
        %if selected, analyses will be conducted using %dF/F
        %variable dFF will represent the raw traces while cordFF will
        %represent the %dF/F output
        if FP.info.percentdFF(1) == 1
            cordFF=((dFF-meanbaseline)./meanbaseline)*100; %percent change in dFF from baseline
        else
            %proceed with non-%dF/F data
        end
        
        %area under the curve calculation
        auc=trapz(cordFF); %area under corrected dFF curve
        
        %peak detection and analysis
        minpk=movmedian(dFF,FP.info.windowframes)+(2*movstd(dFF,FP.info.windowframes)); %minimum height to for peak detection
        %consideration is set to two standard deviations above the median 
        %dFF value for the test segment of the recording
        tracedif=ge(dFF,minpk);
        pks=dFF(tracedif);
        %since minimum peak height was calculated using the un-adjusted
        %test portion, we must use the un-corrected dFF trace to find
        %peaks. Their indicies can be found using position 2 of the
        %findpeaks function
        pkfreq=size(pks,1)/(FP.info.samplingrate(1)*size(time,1)); %total peak frequency
        meanpk=mean(cordFF(tracedif)); %mean peak height
        maxpkidx=find(dFF==max(dFF),1,'first'); %find index of max peak. If multiple, it will select first
        maxpk=max(cordFF); %max dFF
        maxpktime=time(maxpkidx); %time of max

        %input data into output structures
        FP.outputs.behaviourindependent.traces.dFF_notcorrected.(FP.info.groups(ii))(:,iii)=dFF;
        FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))(:,iii)=cordFF;
        FP.outputs.behaviourindependent.analyses.AUC.(FP.info.groups(ii))(iii)=auc;
        FP.outputs.behaviourindependent.analyses.PeakFrequency.(FP.info.groups(ii))(iii)=pkfreq;
        FP.outputs.behaviourindependent.analyses.MeanPeakHeight.(FP.info.groups(ii))(iii)=meanpk;
        FP.outputs.behaviourindependent.analyses.MaxPeakHeight.(FP.info.groups(ii))(iii)=maxpk;
        FP.outputs.behaviourindependent.analyses.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=maxpktime;
        
    end
end


%clear remaining unnecessary variables
clear ii iii sig iso temp_x isofit sigfit dFF mintest adjust test cordFF
clear baseline idx1 idx2 meanbaseline time minpeak pks auc pkfreq meanpk
clear minpk locs maxpkidx maxpk maxpktime tracedif

%process successful
disp('Process 3. Photometry Processing - Behaviour Independent Complete');

%%  4. Photometry Processing - Behaviour Epoch Dependent
%   
%   Integration of behavioural information into photometry analyses
%   Must run all previous processes before running this process
%
%   User-defined onset and duration of events
%       first = 300 (the onset time in seconds of the first epoch
%       second = 900 (the onset time in seconds of the second epoch
%       duration.familiar = 300 (the duration in seconds of the familiar epoch
%       duration.intruder = 300 (the duration in seconds of the intruder
%       epoch
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
%   NaN used as a placeholder for behaviour-specific dF/F traces. The trace
%   is limited to any values prior to the first NaN occurence

%initialization
FP.info.epoch.first = 300; %start
FP.info.epoch.second = 900; %start
FP.info.epoch.duration.familiar = 300; %duration
FP.info.epoch.duration.intruder = 300; %duration

%batch processing
for ii=1:FP.info.groupnum(1)

        %create structure arrays for photometry outputs
        FP.outputs.behaviourdependent.behaviour.anogenital.familiar.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourdependent.behaviour.anogenital.intruder.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourdependent.behaviour.headtorso.familiar.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourdependent.behaviour.headtorso.intruder.(FP.info.groups(ii))=zeros(180705,FP.info.pergroup(ii)); %#frames
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.dFF.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar/FP.info.samplingrate(1))),FP.info.pergroup(ii)); %placeholder array
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.dFF.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar/FP.info.samplingrate(1))),FP.info.pergroup(ii)); %placeholder array
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.dFF.(FP.info.groups(ii))=NaN(round((FP.info.epoch.duration.intruder/FP.info.samplingrate(1))),FP.info.pergroup(ii)); %placeholder array
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.dFF.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder/FP.info.samplingrate(1))),FP.info.pergroup(ii)); %placeholder array
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.familiarepoch.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.familiarepoch.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.familiarepoch.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.familiarepoch.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.familiarepoch.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.intruderepoch.AUC.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.intruderepoch.PeakFrequency.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.intruderepoch.MeanPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.intruderepoch.MaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        FP.outputs.behaviourdependent.analyses.intruderepoch.TimeMaxPeakHeight.(FP.info.groups(ii))=zeros(1,FP.info.pergroup(ii));
        
    for iii=1:FP.info.pergroup(ii)
            
        %define epochs
        if FP.info.order.(FP.info.groups(ii))(iii)==1
            FP.info.epoch.familiar = [FP.info.epoch.first,FP.info.epoch.duration.familiar]; %familiar first
            FP.info.epoch.intruder = [FP.info.epoch.second,FP.info.epoch.duration.intruder]; %intruder second
        elseif FP.info.order.(FP.info.groups(ii))(iii)==0
            FP.info.epoch.familiar = [FP.info.epoch.second,FP.info.epoch.duration.familiar]; %intruder first
            FP.info.epoch.intruder = [FP.info.epoch.first,FP.info.epoch.duration.intruder]; %familiar second
        end 
        
        %load data
        dFF=FP.outputs.behaviourindependent.traces.dFF_notcorrected.(FP.info.groups(ii))(:,iii); %load not corrected dFF
        time=FP.inputs.doric.time.(FP.info.groups(ii))(:,iii);
        cordFF=FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))(:,iii); %load corrected dFF
        htf=FP.inputs.familiar.headtorso.(FP.info.groups(ii))(:,iii); %load familiar resident_headtorso column
        agf=FP.inputs.familiar.anogenital.(FP.info.groups(ii))(:,iii); %load familiar resident_anogenital column
        hti=FP.inputs.intruder.headtorso.(FP.info.groups(ii))(:,iii); %load intruder resident_headtorso column
        agi=FP.inputs.intruder.anogenital.(FP.info.groups(ii))(:,iii); %load intruder resident_anogenital column
        
        %identify time segments of behavioural epochs
        [~,idxf1]=min(abs(time-(FP.info.epoch.familiar(1)))); %identify beginning of familiar
        [~,idxf2]=min(abs(time-((FP.info.epoch.familiar(2)+FP.info.epoch.familiar(1))))); %identify end of familiar
        [~,idxi1]=min(abs(time-(FP.info.epoch.intruder(1)))); %identify beginning of intruder
        [~,idxi2]=min(abs(time-((FP.info.epoch.intruder(2)+FP.info.epoch.intruder(1))))); %identify end of intruder
        fcordFF=cordFF(idxf1:idxf2); %corrected dFF during the familiar epoch
        icordFF=cordFF(idxi1:idxi2); %corrected dFF during the intruder epoch
        fdFF=dFF(idxf1:idxf2); %corrected dFF during the familiar epoch
        idFF=dFF(idxi1:idxi2); %corrected dFF during the intruder epoch
        ftime=time(idxf1:idxf2); %photometry time during familiar
        itime=time(idxi1:idxi2); %photometry time during familiar
        
        %metrics during these epochs
        %familiar
        FP.outputs.behaviourdependent.analyses.familiarepoch.AUC.(FP.info.groups(ii))(iii)=trapz(fcordFF);
        minpk=movmedian(fdFF,FP.info.windowframes)+(2*movstd(fdFF,FP.info.windowframes));
        tracedif=ge(fdFF,minpk);
        pks=fdFF(tracedif);
        FP.outputs.behaviourdependent.analyses.familiarepoch.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(FP.info.samplingrate(1)*size(ftime,1));
        FP.outputs.behaviourdependent.analyses.familiarepoch.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(fcordFF(tracedif));
        maxpkidx=find(fdFF==max(fdFF),1,'first');
        FP.outputs.behaviourdependent.analyses.familiarepoch.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(fcordFF);
        FP.outputs.behaviourdependent.analyses.familiarepoch.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=ftime(maxpkidx);
        %intruder
        FP.outputs.behaviourdependent.analyses.intruderepoch.AUC.(FP.info.groups(ii))(iii)=trapz(icordFF);
        minpk=movmedian(idFF,FP.info.windowframes)+(2*movstd(idFF,FP.info.windowframes));
        tracedif=ge(idFF,minpk);
        pks=idFF(tracedif);
        FP.outputs.behaviourdependent.analyses.intruderepoch.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(FP.info.samplingrate(1)*size(itime,1));
        FP.outputs.behaviourdependent.analyses.intruderepoch.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(icordFF(tracedif));
        maxpkidx=find(idFF==max(idFF),1,'first');
        FP.outputs.behaviourdependent.analyses.intruderepoch.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(icordFF);
        FP.outputs.behaviourdependent.analyses.intruderepoch.TimeMaxPeakHeight.(FP.info.groups(ii))=itime(maxpkidx);
        
        %find indices of behaviours during epochs
        %stretch and align datasets
        %sets arrays to the same length and gap fills
        % *** NOTE: WILL NOT WORK PROPERLY IF THE FRAMERATE OF THE
        % BEHAVIOURAL RECORDING EXCEEDS THAT OF THE PHOTOMETRY RECORDING
        % ***
        %familiar
        famtime=(FP.info.epoch.familiar(1):(FP.info.epoch.familiar(2)/size(htf,1)):(FP.info.epoch.familiar(1)+FP.info.epoch.familiar(2))); %dummy time series
        famtime=famtime(2:end)'; %trim starting zero
        htfstretch=zeros(size(ftime,1),1); %generate sink to contain the new dataframes
        agfstretch=zeros(size(ftime,1),1); %generate sink to contain the new dataframes
        for iv=1:size(ftime,1) %batch one row at a time
            temp_time=ftime(iv);   
            [~,ix]=min(abs(famtime-temp_time)); %find index to align closest timestamps
            htfstretch(iv)=htf(ix,:); %apply index to head/torso
            agfstretch(iv)=agf(ix,:); %apply index to anogenital
        end
        %familiar
        inttime=(FP.info.epoch.intruder(1):(FP.info.epoch.intruder(2)/size(hti,1)):(FP.info.epoch.intruder(1)+FP.info.epoch.intruder(2))); %dummy time series
        inttime=inttime(2:end)'; %trim starting zero
        htistretch=zeros(size(itime,1),1); %generate sink to contain the new dataframes
        agistretch=zeros(size(itime,1),1); %generate sink to contain the new dataframes
        for iv=1:size(itime,1) %batch one row at a time
            temp_time=itime(iv);   
            [~,ix]=min(abs(inttime-temp_time)); %find index to align closest timestamps
            htistretch(iv)=hti(ix,:); %apply index to head/torso
            agistretch(iv)=agi(ix,:); %apply index to anogenital
        end

        %find indices of bouts of behaviour of interest
        htfstretch=(htfstretch==1); %turns this into a logical array
        agfstretch=(agfstretch==1); %turns this into a logical array
        htistretch=(htistretch==1); %turns this into a logical array
        agistretch=(agistretch==1); %turns this into a logical array
        
        %save stretched binaries for visualization
        FP.outputs.behaviourdependent.behaviour.anogenital.familiar.(FP.info.groups(ii))((idxf1:idxf2),iii)=agfstretch;
        FP.outputs.behaviourdependent.behaviour.anogenital.intruder.(FP.info.groups(ii))((idxi1:idxi2),iii)=agistretch;
        FP.outputs.behaviourdependent.behaviour.headtorso.familiar.(FP.info.groups(ii))((idxf1:idxf2),iii)=htfstretch;
        FP.outputs.behaviourdependent.behaviour.headtorso.intruder.(FP.info.groups(ii))((idxi1:idxi2),iii)=htistretch;
        
        %use incides of behaviours to isolate behaviour-specific photometry
        %signal
        %corrected
        htf=fcordFF(htfstretch);
        agf=fcordFF(agfstretch);
        hti=icordFF(htistretch);
        agi=icordFF(agistretch);
        %not corrected
        nchtf=fdFF(htfstretch);
        ncagf=fdFF(agfstretch);
        nchti=idFF(htistretch);
        ncagi=idFF(agistretch);
        
        %write not corrected to structure for visualiaztion
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.dFF.(FP.info.groups(ii))(1:size(nchtf,1),iii)=nchtf;
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.dFF.(FP.info.groups(ii))(1:size(ncagf,1),iii)=ncagf;
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.dFF.(FP.info.groups(ii))(1:size(nchti,1),iii)=nchti;
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.dFF.(FP.info.groups(ii))(1:size(ncagi,1),iii)=ncagi;
        
        htfminpk=movmedian(nchtf,FP.info.windowframes)+(2*movstd(nchtf,FP.info.windowframes));
        agfminpk=movmedian(ncagf,FP.info.windowframes)+(2*movstd(ncagf,FP.info.windowframes));
        htiminpk=movmedian(nchti,FP.info.windowframes)+(2*movstd(nchti,FP.info.windowframes));
        agiminpk=movmedian(ncagi,FP.info.windowframes)+(2*movstd(ncagi,FP.info.windowframes));
        
        %analyze and input data
        %familiar
        %head/torso
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.AUC.(FP.info.groups(ii))(iii)=trapz(htf)/(size(htf,1)*(time(2)-time(1)));
        tracedif=ge(nchtf,htfminpk);
        pks=dFF(tracedif);
        maxpkidx=find(cordFF==max(htf),1,'first'); %find index of max peak. If multiple, it will select first
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(htf); %max dFF
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=time(maxpkidx); %time of max
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(size(nchtf,1)*(time(2)-time(1)));
        FP.outputs.behaviourdependent.analyses.headtorso.familiar.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(cordFF(tracedif));
        %anogenital
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.AUC.(FP.info.groups(ii))(iii)=trapz(agf)/(size(agf,1)*(time(2)-time(1)));
        tracedif=ge(ncagf,agfminpk);
        pks=dFF(tracedif);
        maxpkidx=find(cordFF==max(agf),1,'first'); %find index of max peak. If multiple, it will select first
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(agf); %max dFF
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=time(maxpkidx); %time of max
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(size(ncagf,1)*(time(2)-time(1)));
        FP.outputs.behaviourdependent.analyses.anogenital.familiar.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(cordFF(tracedif));
        
        %intruder
        %head/torso
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.AUC.(FP.info.groups(ii))(iii)=trapz(hti)/(size(hti,1)*(time(2)-time(1)));
        tracedif=ge(nchti,htiminpk);
        pks=dFF(tracedif);
        maxpkidx=find(cordFF==max(hti),1,'first'); %find index of max peak. If multiple, it will select first
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(hti); %max dFF
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=time(maxpkidx); %time of max
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(size(nchti,1)*(time(2)-time(1)));
        FP.outputs.behaviourdependent.analyses.headtorso.intruder.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(cordFF(tracedif));
        %anogenital
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.AUC.(FP.info.groups(ii))(iii)=trapz(agi)/(size(agi,1)*(time(2)-time(1)));
        tracedif=ge(ncagi,agiminpk);
        pks=dFF(tracedif);
        maxpkidx=find(cordFF==max(agi),1,'first'); %find index of max peak. If multiple, it will select first
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.MaxPeakHeight.(FP.info.groups(ii))(iii)=max(agi); %max dFF
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))(iii)=time(maxpkidx); %time of max
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.PeakFrequency.(FP.info.groups(ii))(iii)=size(pks,1)/(size(ncagi,1)*(time(2)-time(1)));
        FP.outputs.behaviourdependent.analyses.anogenital.intruder.MeanPeakHeight.(FP.info.groups(ii))(iii)=mean(cordFF(tracedif));
        
    end
end


%clear remaining unnecessary variables
clear ii iii iv ix fdFF idFF
clear time famtime inttime dFF cordFF fcordFF icordFF
clear htf agf hti agi nchtf ncagf nchti ncagi htfstretch agfstretch htistretch agistretch
clear idxf1 idxf2 idxi1 idxi2 idx2
clear temp_time test minpk pks ftime itime locs maxpkidx tracedif htfminpk htiminpk agfminpk agiminpk

%process successful
disp('Process 4. Photometry Processing - Behaviour Dependent Complete');

%%  5. Photometry Processing - Bout-by-Bout Information

%   The following section will detail photometry outputs for each bout of
%   anogenital or head/torso interactions.
%
%   Each column represents a different mouse, with bout information being
%   presented in chronological order down the column. Ex: if the AUC during
%   the first bout for mouse 1 is 3, the cell at the intersection of the
%   first column and row will be 3. NaN represents a placeholder value
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

%batch processing
for ii=1:FP.info.groupnum(1)
            
        %create structure arrays for photometry outputs %NaN placeholders
        FP.outputs.boutbybout.analyses.anogenital.familiar.AUC.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.familiar.AUC.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.intruder.AUC.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.intruder.AUC.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.familiar.PeakFrequency.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.familiar.PeakFrequency.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.intruder.PeakFrequency.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.intruder.PeakFrequency.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.familiar.MeanPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.familiar.MeanPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.intruder.MeanPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.intruder.MeanPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.familiar.MaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.familiar.MaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.intruder.MaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.intruder.MaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.familiar(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.anogenital.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        FP.outputs.boutbybout.analyses.headtorso.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))=NaN((round(FP.info.epoch.duration.intruder(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        
    for iii=1:FP.info.pergroup(ii)
        
        %define epochs
        if FP.info.order.(FP.info.groups(ii))(iii)==1
            FP.info.epoch.familiar = [FP.info.epoch.first,FP.info.epoch.duration.familiar]; %familiar first
            FP.info.epoch.intruder = [FP.info.epoch.second,FP.info.epoch.duration.intruder]; %intruder second
        else
            FP.info.epoch.familiar = [FP.info.epoch.second,FP.info.epoch.duration.familiar]; %intruder first
            FP.info.epoch.intruder = [FP.info.epoch.first,FP.info.epoch.duration.intruder]; %familiar second
        end 
        
        %load data
        time=FP.inputs.doric.time.(FP.info.groups(ii))(:,iii);
        cordFF=FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))(:,iii); %load corrected dFF
        htf=FP.inputs.familiar.headtorso.(FP.info.groups(ii))(:,iii); %load familiar resident_headtorso column
        agf=FP.inputs.familiar.anogenital.(FP.info.groups(ii))(:,iii); %load familiar resident_anogenital column
        hti=FP.inputs.intruder.headtorso.(FP.info.groups(ii))(:,iii); %load intruder resident_headtorso column
        agi=FP.inputs.intruder.anogenital.(FP.info.groups(ii))(:,iii); %load intruder resident_anogenital column
        
        %identify time segments of behavioural epochs
        [~,idxf1]=min(abs(time-(FP.info.epoch.familiar(1)))); %identify beginning of familiar
        [~,idxf2]=min(abs(time-((FP.info.epoch.familiar(2)+FP.info.epoch.familiar(1))))); %identify end of familiar
        [~,idxi1]=min(abs(time-(FP.info.epoch.intruder(1)))); %identify beginning of intruder
        [~,idxi2]=min(abs(time-((FP.info.epoch.intruder(2)+FP.info.epoch.intruder(1))))); %identify end of intruder
        fcordFF=cordFF(idxf1:idxf2); %corrected dFF during the familiar epoch
        icordFF=cordFF(idxi1:idxi2); %corrected dFF during the intruder epoch
        fdFF=cordFF(idxf1:idxf2); %dFF during the familiar epoch
        fminpk=movmedian(fdFF,FP.info.windowframes)+(2*movstd(fdFF,FP.info.windowframes)); %minimum height to for peak detection
        idFF=cordFF(idxi1:idxi2); %dFF during the intruder epoch
        iminpk=movmedian(idFF,FP.info.windowframes)+(2*movstd(idFF,FP.info.windowframes)); %minimum height to for peak detection
        ftime=time(idxf1:idxf2); %photometry time during familiar
        itime=time(idxi1:idxi2); %photometry time during familiar
        
        %find indices of behaviours during epochs
        %stretch and align datasets
        %sets arrays to the same length and gap fills
        % *** NOTE: WILL NOT WORK PROPERLY IF THE FRAMERATE OF THE
        % BEHAVIOURAL RECORDING EXCEEDS THAT OF THE PHOTOMETRY RECORDING
        % ***
        %familiar
        famtime=(FP.info.epoch.familiar(1):(FP.info.epoch.familiar(2)/size(htf,1)):(FP.info.epoch.familiar(1)+FP.info.epoch.familiar(2))); %dummy time series
        famtime=famtime(2:end)'; %trim starting zero
        htfstretch=zeros(size(ftime,1),1); %generate sink to contain the new dataframes
        agfstretch=zeros(size(ftime,1),1); %generate sink to contain the new dataframes
        for iv=1:size(ftime,1) %batch one row at a time
            temp_time=ftime(iv);   
            [~,ix]=min(abs(famtime-temp_time)); %find index to align closest timestamps
            htfstretch(iv)=htf(ix,:); %apply index to head/torso
            agfstretch(iv)=agf(ix,:); %apply index to anogenital
        end
        %familiar
        inttime=(FP.info.epoch.intruder(1):(FP.info.epoch.intruder(2)/size(hti,1)):(FP.info.epoch.intruder(1)+FP.info.epoch.intruder(2))); %dummy time series
        inttime=inttime(2:end)'; %trim starting zero
        htistretch=zeros(size(itime,1),1); %generate sink to contain the new dataframes
        agistretch=zeros(size(itime,1),1); %generate sink to contain the new dataframes
        for iv=1:size(itime,1) %batch one row at a time
            temp_time=itime(iv);   
            [~,ix]=min(abs(inttime-temp_time)); %find index to align closest timestamps
            htistretch(iv)=hti(ix,:); %apply index to head/torso
            agistretch(iv)=agi(ix,:); %apply index to anogenital
        end

        %find indices of bouts of behaviour of interest
        htfstretch=(htfstretch==1); %turns this into a logical array
        agfstretch=(agfstretch==1); %turns this into a logical array
        htistretch=(htistretch==1); %turns this into a logical array
        agistretch=(agistretch==1); %turns this into a logical array
        
        %find indices of when behaviours change
        %assumes that the first frame will be a 0, meaning that there is no
        %interaction occuring
        
        %anogenital familiar
        %initialization
        temp1=find(diff(agfstretch)==0);
        behavchange=setdiff(min(temp1):max(temp1),temp1); %find indices of changes in behavioural state
        behavlist=[1,behavchange,size(agfstretch,1)];
        behavstart=behavlist(1:2:end); %start of a behavioural sequence
        behavend=behavlist(2:2:end);  %end of a behavioural sequence
        if size(behavstart,2)>size(behavend,2) %safety check to ensure that arrays are of the same length
            behavstart=behavstart(1:size(behavend,2));
        else
            behavend=behavend(1:size(behavstart,2));
        end
        
        %batch through bouts 
        for iv=1:(size(behavstart,2)-1)
            tempdFF=fcordFF(behavend(iv)+1:behavstart(iv+1));
            tempdFF2=fcordFF(behavend(iv)+1:behavstart(iv+1));
            minpk=fminpk((behavend(iv)+1):behavstart(iv+1)); %minimum height to for peak detection
            tracedif=ge(tempdFF,minpk);
            pks=tempdFF(tracedif);
            FP.outputs.boutbybout.analyses.anogenital.familiar.AUC.(FP.info.groups(ii))(iv,iii)=(trapz(tempdFF2)/(size(tempdFF2,1)*FP.info.samplingrate(1)));
            maxpkidx=find(cordFF==max(tempdFF),1,'first'); %find index of max peak. If multiple, it will select first
            FP.outputs.boutbybout.analyses.anogenital.familiar.MaxPeakHeight.(FP.info.groups(ii))(iv,iii)=max(tempdFF2); %max dFF
            FP.outputs.boutbybout.analyses.anogenital.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))(iv,iii)=time(maxpkidx); %time of max
            FP.outputs.boutbybout.analyses.anogenital.familiar.PeakFrequency.(FP.info.groups(ii))(iv,iii)=(size(pks,1)/(size(tempdFF,1)*FP.info.samplingrate(1)));
            FP.outputs.boutbybout.analyses.anogenital.familiar.MeanPeakHeight.(FP.info.groups(ii))(iv,iii)=mean(tempdFF2(tracedif));
        end
        
        %anogenital intruder
        %initialization
        temp1=find(diff(agistretch)==0);
        behavchange=setdiff(min(temp1):max(temp1),temp1); %find indices of changes in behavioural state
        behavlist=[1,behavchange,size(agfstretch,1)];
        behavstart=behavlist(1:2:end); %start of a behavioural sequence
        behavend=behavlist(2:2:end);  %end of a behavioural sequence
        if size(behavstart,2)>size(behavend,2) %safety check to ensure that arrays are of the same length
            behavstart=behavstart(1:size(behavend,2));
        else
            behavend=behavend(1:size(behavstart,2));
        end
        
        %batch through bouts 
        for iv=1:(size(behavstart,2)-1)
            tempdFF=idFF(behavend(iv)+1:behavstart(iv+1));
            tempdFF2=icordFF(behavend(iv)+1:behavstart(iv+1));
            FP.outputs.boutbybout.analyses.anogenital.intruder.AUC.(FP.info.groups(ii))(iv,iii)=(trapz(tempdFF2)/(size(tempdFF2,1)*FP.info.samplingrate(1)));
            minpk=iminpk((behavend(iv)+1):behavstart(iv+1)); %minimum height to for peak detection
            tracedif=ge(tempdFF,minpk);
            maxpkidx=find(cordFF==max(tempdFF),1,'first'); %find index of max peak. If multiple, it will select first
            FP.outputs.boutbybout.analyses.anogenital.intruder.MaxPeakHeight.(FP.info.groups(ii))(iv,iii)=max(tempdFF2); %max dFF
            FP.outputs.boutbybout.analyses.anogenital.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))(iv,iii)=time(maxpkidx); %time of max
            FP.outputs.boutbybout.analyses.anogenital.intruder.PeakFrequency.(FP.info.groups(ii))(iv,iii)=(size(pks,1)/(size(tempdFF,1)*FP.info.samplingrate(1)));
            FP.outputs.boutbybout.analyses.anogenital.intruder.MeanPeakHeight.(FP.info.groups(ii))(iv,iii)=mean(tempdFF2(tracedif));
        end
        
        %headtorso familiar
        %initialization
        temp1=find(diff(htfstretch)==0);
        behavchange=setdiff(min(temp1):max(temp1),temp1); %find indices of changes in behavioural state
        behavlist=[1,behavchange,size(agfstretch,1)];
        behavstart=behavlist(1:2:end); %start of a behavioural sequence
        behavend=behavlist(2:2:end);  %end of a behavioural sequence
        if size(behavstart,2)>size(behavend,2) %safety check to ensure that arrays are of the same length
            behavstart=behavstart(1:size(behavend,2));
        else
            behavend=behavend(1:size(behavstart,2));
        end
        
        %batch through bouts 
        for iv=1:(size(behavstart,2)-1)
            tempdFF=fdFF((behavend(iv)+1):behavstart(iv+1));
            tempdFF2=fcordFF(behavend(iv)+1:behavstart(iv+1));
            FP.outputs.boutbybout.analyses.headtorso.familiar.AUC.(FP.info.groups(ii))(iv,iii)=(trapz(tempdFF2)/(size(tempdFF2,1)*FP.info.samplingrate(1)));
            minpk=fminpk((behavend(iv)+1):behavstart(iv+1)); %minimum height to for peak detection
            tracedif=ge(tempdFF,minpk);
            maxpkidx=find(cordFF==max(tempdFF),1,'first'); %find index of max peak. If multiple, it will select first
            FP.outputs.boutbybout.analyses.headtorso.familiar.MaxPeakHeight.(FP.info.groups(ii))(iv,iii)=max(tempdFF2); %max dFF
            FP.outputs.boutbybout.analyses.headtorso.familiar.TimeMaxPeakHeight.(FP.info.groups(ii))(iv,iii)=time(maxpkidx); %time of max
            FP.outputs.boutbybout.analyses.headtorso.familiar.PeakFrequency.(FP.info.groups(ii))(iv,iii)=(size(pks,1)/(size(tempdFF,1)*FP.info.samplingrate(1)));
            FP.outputs.boutbybout.analyses.headtorso.familiar.MeanPeakHeight.(FP.info.groups(ii))(iv,iii)=mean(tempdFF2(tracedif));
        end
        
        %headtorso intruder
        %initialization
        temp1=find(diff(htistretch)==0);
        behavchange=setdiff(min(temp1):max(temp1),temp1); %find indices of changes in behavioural state
        behavlist=[1,behavchange,size(agfstretch,1)];
        behavstart=behavlist(1:2:end); %start of a behavioural sequence
        behavend=behavlist(2:2:end);  %end of a behavioural sequence
        if size(behavstart,2)>size(behavend,2) %safety check to ensure that arrays are of the same length
            behavstart=behavstart(1:size(behavend,2));
        else
            behavend=behavend(1:size(behavstart,2));
        end
        
        %batch through bouts 
        for iv=1:(size(behavstart,2)-1)
            tempdFF=idFF(behavend(iv)+1:behavstart(iv+1));
            tempdFF2=icordFF(behavend(iv)+1:behavstart(iv+1));
            FP.outputs.boutbybout.analyses.headtorso.intruder.AUC.(FP.info.groups(ii))(iv,iii)=(trapz(tempdFF2)/(size(tempdFF2,1)*FP.info.samplingrate(1)));
            minpk=iminpk((behavend(iv)+1):behavstart(iv+1)); %minimum height to for peak detection
            tracedif=ge(tempdFF,minpk);
            maxpkidx=find(cordFF==max(tempdFF),1,'first'); %find index of max peak. If multiple, it will select first
            FP.outputs.boutbybout.analyses.headtorso.intruder.MaxPeakHeight.(FP.info.groups(ii))(iv,iii)=max(tempdFF2); %max dFF
            FP.outputs.boutbybout.analyses.headtorso.intruder.TimeMaxPeakHeight.(FP.info.groups(ii))(iv,iii)=time(maxpkidx); %time of max
            FP.outputs.boutbybout.analyses.headtorso.intruder.PeakFrequency.(FP.info.groups(ii))(iv,iii)=(size(pks,1)/(size(tempdFF,1)*FP.info.samplingrate(1)));
            FP.outputs.boutbybout.analyses.headtorso.intruder.MeanPeakHeight.(FP.info.groups(ii))(iv,iii)=mean(tempdFF2(tracedif));
        end
    end
end

%clear remaining unnecessary variables
clear ii iii iv ix
clear time famtime inttime dFF cordFF fcordFF icordFF fdFF idFF
clear htf agf hti agi nchtf ncagf nchti ncagi htfstretch agfstretch htistretch agistretch
clear idxf1 idxf2 idxi1 idxi2 idx2
clear temp_time test minpk pks ftime itime temp1 tempdFF2 fminpk iminpk
clear behavchange behavend behavstart behavlist tempdFF locs maxpkidx tracedif

%process successful
disp('Process 5. Photometry Processing - Bout-by-Bout Information Complete');

%%  X. Save Structure to .mat File
%keep file for easy data access later

save('FP_inputs.mat','-struct','FP');

disp('Process X. Save Structure to a .mat File Complete');

%%  Y. Load .mat File
%   navigate to folder containing 'FP.mat'

load('FP_inputs.mat')
 w=whos;
 for ii=1:length(w)
     FP.(w(ii).name)=eval(w(ii).name); 
     clear w(ii).name
 end
 clear w info inputs outputs ii
 
 disp('Process Y. Load .mat File Complete');
 
%%  Z. Loading Without Behavioural Data
%  Run on a mouse-by-mouse basis
%  same processes as described above

%   Default information for each input is listed below
%
%   Analysis Information:
%       Z-Score Use:
%           0 if not necessary
%           1 if necessary
%           default: zscore = 0
%       %dF/F Use:
%           0 if not necessary
%           1 if necessary
%           default: percentdFF = 0
%       NOTE: ONLY SELECT ONE OF THESE OPTIONS
%
%       Limiting the time-range of the analysis:
%           0 if not necessary
%           1 if necessary
%           default: timelimit = 0
%       NOTE: Choosing this option does not limit the range of the period
%       being used for baseline correction or any kind of normalization.
%       Adjustments to the trace (Z-score, % change, and baseline
%       correction) will continue to be calculated for the entire trace.
%       Define limits to the time-range:
%           FP.info.timerange = [300 60];  [start duration]

%clear structure to start fresh
clear FP

%indicate whether z scores are wanted
FP.info.zscore = 0;
%indicate whether %dF/F is wanted
FP.info.percentdFF = 0;
%indicate whether analysis should be restricted to a specific time period
FP.info.timelimit = 0;
FP.info.timerange = [300 60];  %[start duration]

%peak detection rolling window time in seconds
FP.info.windowtime = 30;

%   Variables to be adjusted:
%       Number of frames for each section (see descriptions)
%       Individual Information:
%           Identify animal ID
%           Identify Doric output file
%           Number of frames in the dataset

%import raw data
prompt=strcat("Select Doric file");
disp(prompt);
[doricfile,pathname]=uigetfile(('*.csv'),prompt); %load doric file
doricdata=readmatrix(strcat(pathname,doricfile));
%trim if exceeds maximum trial duration (default = 180705, aka 25mins
%@120fps)
if size(doricdata,1)>180705 %#frames
   doricdata=doricdata(1:180705,:); %#frames
else
   %do nothing
end


%organize raw data
%adjust the number of frames (default = 180705, aka 25mins @120fps) as needed
FP.inputs.doric.time=doricdata(:,1); %timestamps from doric output
FP.inputs.doric.sig=doricdata(:,2); %calcium-dependent signal from doric output
FP.inputs.doric.iso=doricdata(:,4); %isosbestic sigal from doric output

%   The generation of delta F/F data from the photometry inputs
%   
%   User-inputted baseline period
%       default: baseline = [60 180]
%           first value is the number of seconds to trim from the start of
%           the trial, to minimize the influence of handling stress on the
%           baseline value
%           second value is the duration of the baseline period in seconds
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

%initialization
FP.info.baseline = [60 180];  %[start duration]

%full trace
time=FP.inputs.doric.time;
sig=FP.inputs.doric.sig;
iso=FP.inputs.doric.iso;
temp_x=1:length(iso);
temp_x=temp_x'; %transpose data frame
isofit=fit(temp_x,iso,'exp2'); %fit isosbestic data with a biexponential decay
sigfit=robustfit(isofit(temp_x),sig); %linearly scale biexponential decay to the raw 470 data using robust fit
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);
dFF=(sig-sigfit)./sigfit;   %delta F/F
        
%add framerate to info
FP.info.samplingrate=time(2)-time(1);

%rolling window frame number
FP.info.windowframes=round(FP.info.windowtime/FP.info.samplingrate);
        
%baseline adjustment
[~,idx1]=min(abs(time-(FP.info.baseline(1)))); %identify beginning of baseline
[~,idx2]=min(abs(time-((FP.info.baseline(2)+FP.info.baseline(1))))); %identify end of baseline
baseline=dFF(idx1:idx2); %baseline period
test=dFF(idx2:end); %test period
meanbaseline=mean(baseline); %calculate mean baseline dFF
mintest=min(test); %identify minimum dFF during the test segment
adjust=diff([mintest,meanbaseline]); %calculate the difference between baseline mean and test minimum
cordFF=dFF+abs(adjust); %add the absolute value of difference to every point in the trace
        
%z-score
%if selected, z-score will be applied to the analysis
%variables cordFF and dFF will no longer represent raw dF/F traces
%and will instead be zscores without baseline correction
   if FP.info.zscore(1) == 1
       cordFF=zscore(dFF);
       dFF=zscore(dFF);
   else
        %proceed with non-z-scored data
   end
        
% %dF/F
%if selected, analyses will be conducted using %dF/F
%variable dFF will represent the raw traces while cordFF will
%represent the %dF/F output
    if FP.info.percentdFF(1) == 1
       cordFF=((dFF-meanbaseline)./meanbaseline)*100; %percent change in dFF from baseline
    else
       %proceed with non-%dF/F data
    end

%limit time range if option has been selected
    if FP.info.timelimit(1) == 1
        [~,idx1]=min(abs(time-(FP.info.timerange(1)))); %identify beginning of baseline
        [~,idx2]=min(abs(time-((FP.info.timerange(2)+FP.info.timerange(1))))); %identify end of baseline
        cordFF=cordFF(idx1:idx2);
        dFF=dFF(idx1:idx2);
    else
        %proceed with full trace
    end
    
%area under the curve calculation
auc=trapz(cordFF); %area under corrected dFF curve
        
%peak detection and analysis
minpk=movmedian(dFF,FP.info.windowframes)+(2*movstd(dFF,FP.info.windowframes)); %minimum height to for peak 
%consideration is set to two standard deviations above the median 
%dFF value for the test segment of the recording
tracedif=ge(dFF,minpk);
pks=dFF(tracedif);
locs=cordFF(tracedif);
%since minimum peak height was calculated using the un-adjusted
%test portion, we must use the un-corrected dFF trace to find
%peaks. Their indicies can be found using position 2 of the
%findpeaks function
pkfreq=size(pks,1)/(FP.info.samplingrate(1)*size(time,1)); %total peak frequency
meanpk=mean(locs); %mean peak height
maxpkidx=find(dFF==max(dFF),1,'first'); %find index of max peak. If multiple, it will select first
maxpk=max(dFF); %max dFF
maxpktime=time(maxpkidx); %time of max

%input data into output structures
FP.outputs.traces.dFF_notcorrected=dFF;
FP.outputs.traces.correcteddFF=cordFF;
FP.outputs.analyses.AUC=auc;
FP.outputs.analyses.PeakFrequency=pkfreq;
FP.outputs.analyses.MeanPeakHeight=meanpk;
FP.outputs.analyses.MaxPeakHeight=maxpk;
FP.outputs.analyses.TimeOfMaxPeak=maxpktime;

%clear remaining variables
clear adjust auc baseline cordFF dFF doricdata doricfile idx1 idx2 iso isofit
clear meanbaseline meanpk minpk mintest pathname pkfreq pks prompt sig sigfit temp_x
clear test time locs maxpkidx maxpk maxpktime tracedif
