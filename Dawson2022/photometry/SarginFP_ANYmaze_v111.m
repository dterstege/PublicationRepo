%%                          SarginFP_ANYmaze.m
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

disp('Running SarginFP_ANYmaze');

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
%       Column of interest in ANYmaze data:
%           default: anymaze.columns = [4 5 6] (columns 4, 5, and 6)
%       Titles of column of interest in ANYmaze data:
%           default: anymaze.titles = ["left" "right" "middle"]
%       Total number of frames in the Doric files (framerate x duration)
%           default: fpframes = 108400
%       Total number of frames in the ANYmaze files (framerate x duration)
%           default: amframes = FP.info.amframes
%       Indicate the start and duration, if any, of a pickup window during
%       which data will be omitted
%           default: pickupwindow = [295 10];
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
FP.info.groups = ["control"]; %manually adjust

%input number of animals per group
FP.info.pergroup = [1]; %manually adjust

%input column of interest in the ANYmaze file
FP.info.anymaze.columns = [7 8]; %manually adjust

%input behaviours of interest in the ANYmaze file
FP.info.anymaze.titles = ["rightsniffing" "leftsniffing"]; %manually adjust

%input number of frames for full photometry trace (ie framerate x time)
FP.info.fpframes = 72280;

%input number of frames for full ANYmaze output (ie framerate x time)
FP.info.amframes = 4877;

%indicate start and duration of pickup window
FP.info.pickupwindow = [295 10];

%indicate whether z scores are wanted
FP.info.zscore = 0; %manually adjust

%indicate whether %dF/F is wanted
FP.info.percentdFF = 0; %manually adjust



%process successful
disp('Process 1. Initialization - Group Data Complete');

%%  02. Initialization - Subject Information
%   Somewhat time intensive
%   Requires manual selection of files
%   Default information for each input is listed below
%
%   Variables to be adjusted:
%       Individual Information:
%           Identify animal ID
%           Identify Doric output file
%           Identify SimBA output file

%prepare ANYmaze data import
FP.info.numberofbehav=size(FP.info.anymaze.columns,2);

%internal initialization
for ii = 1:FP.info.groupnum(1)
    %create structure arrays for raw doric timestamp inputs
    FP.inputs.doric.time.(FP.info.groups(ii))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
    %create structure arrays for raw doric calcium signal inputs
    FP.inputs.doric.sig.(FP.info.groups(ii))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
    %create structure arrays for raw doric isosbestic inputs
    FP.inputs.doric.iso.(FP.info.groups(ii))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
    %create structure array(s) for ANYmaze data
    for iii = 1:FP.info.numberofbehav(1)
        FP.inputs.anymaze.data.(FP.info.groups(ii)).(FP.info.anymaze.titles(iii))=zeros(FP.info.amframes,FP.info.pergroup(ii));
    end
end

%import and organize raw data for first group
for ii=1:FP.info.groupnum(1)
    for iii=1:FP.info.pergroup(ii)
        %import raw data
        prompt=strcat("Select Doric file for mouse #",num2str(iii)," in ",FP.info.groups(ii)," group:");
        disp(prompt);
        [doricfile,pathname]=uigetfile(('*.csv'),prompt); %load doric file
        doricdata=readmatrix(strcat(pathname,doricfile));
        %trim if exceeds maximum trial duration (default = FP.info.fpframes)
        if size(doricdata,1)>FP.info.fpframes %
            doricdata=doricdata(1:FP.info.fpframes,:); %
        else
            %do nothing
        end
        prompt=strcat("Select ANYmaze file for mouse #",num2str(iii)," in ",FP.info.groups(ii)," group:");
        disp(prompt);
        [amfile,pathname]=uigetfile(('*.csv'), prompt); %load simba familiar file
        amdata=readmatrix(strcat(pathname,amfile));
        
        amtrim=FP.info.amframes-size(amdata,1);
        
        %organize raw data
        FP.inputs.doric.time.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,1); %timestamps from doric output
        FP.inputs.doric.sig.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,2); %calcium-dependent signal from doric output
        FP.inputs.doric.iso.(FP.info.groups(ii))(1:size(doricdata,1),iii)=doricdata(:,4); %isosbestic sigal from doric output
        for iv = 1:FP.info.numberofbehav(1)
            FP.inputs.anymaze.data.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))((amtrim+1):FP.info.amframes,iii)=amdata(1:(FP.info.amframes-amtrim),FP.info.anymaze.columns(iv));
        end
    end
end


%clear remaining unnecessary variables
clear ii iii iv prompt pathname doricfile amfile
clear doricdata amdata amtrim

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
        FP.outputs.behaviourindependent.traces.dFF_notcorrected.(FP.info.groups(ii))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
        FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
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
        rawdFF=(sig-sigfit)./sigfit;   %delta F/F
        
        %add framerate to info
        FP.info.samplingrate=time(2)-time(1);
        
        %pickup window bounds
        [~,pick1]=min(abs(time-(FP.info.pickupwindow(1)))); %identify beginning of pickup
        [~,pick2]=min(abs(time-((FP.info.pickupwindow(2)+FP.info.pickupwindow(1))))); %identify end of pickup
        
        %baseline adjustment
        [~,idx1]=min(abs(time-(FP.info.baseline(1)))); %identify beginning of baseline
        [~,idx2]=min(abs(time-((FP.info.baseline(2)+FP.info.baseline(1))))); %identify end of baseline
        if pick1<idx2
            idx2=pick1;
        else
            %do nothing
        end
        baseline=rawdFF(idx1:idx2); %baseline period (start of baseline to end of baseline OR start of pickup - whichever comes first)
        test=rawdFF(pick2:end); %test period (end of pickup to end of recording)
        meanbaseline=mean(baseline); %calculate mean baseline dFF
        mintest=min(test); %identify minimum dFF during the test segment
        adjust=diff([mintest,meanbaseline]); %calculate the difference between baseline mean and test minimum
        rawcordFF=rawdFF+abs(adjust); %add the absolute value of difference to every point in the trace
        
        %z-score
        %if selected, z-score will be applied to the analysis
        %variables cordFF and dFF will no longer represent raw dF/F traces
        %and will instead be zscores without baseline correction
        if FP.info.zscore(1) == 1
            rawcordFF=zscore(rawdFF);
            rawdFF=zscore(rawdFF);
        else
            %proceed with non-z-scored data
        end
        
        % %dF/F
        %if selected, analyses will be conducted using %dF/F
        %variable dFF will represent the raw traces while cordFF will
        %represent the %dF/F output
        if FP.info.percentdFF(1) == 1
            rawcordFF=((dFF-meanbaseline)./meanbaseline)*100; %percent change in dFF from baseline
        else
            %proceed with non-%dF/F data
        end
        
        %duplicate raw traces for pickupwindow removal
        dFF=rawdFF;
        cordFF=rawcordFF;
        
        %exclude pickup window
        time(pick1:pick2)=[];
        dFF(pick1:pick2)=[];
        cordFF(pick1:pick2)=[];
        
        %area under the curve calculation
        auc=trapz(cordFF); %area under corrected dFF curve
        
        %peak detection and analysis
        minpk=median(dFF)+(2*std(dFF)); %minimum height to for peak 
        %consideration is set to two standard deviations above the median 
        %dFF value for the test segment of the recording
        [pks,locs]=findpeaks(dFF,'MinPeakHeight',(minpk),'MinPeakDistance',0);
        %since minimum peak height was calculated using the un-adjusted
        %test portion, we must use the un-corrected dFF trace to find
        %peaks. Their indicies can be found using position 2 of the
        %findpeaks function
        pkfreq=size(pks,1)/(FP.info.samplingrate(1)*size(time,1)); %total peak frequency
        meanpk=mean(cordFF(locs)); %mean peak height
        maxpkidx=find(dFF==max(dFF),1,'first'); %find index of max peak. If multiple, it will select first
        maxpk=max(dFF); %max dFF
        maxpktime=time(maxpkidx); %time of max

        %input data into output structures
        FP.outputs.behaviourindependent.traces.dFF_notcorrected.(FP.info.groups(ii))(:,iii)=rawdFF; %saved for visualization purposes only
        FP.outputs.behaviourindependent.traces.correcteddFF.(FP.info.groups(ii))(:,iii)=rawcordFF; %saved for visualization purposes only
        FP.outputs.behaviourindependent.traces.pickupdFF_notcorrected.(FP.info.groups(ii))(:,iii)=dFF;
        FP.outputs.behaviourindependent.traces.pickupcorrecteddFF.(FP.info.groups(ii))(:,iii)=cordFF;
        FP.outputs.behaviourindependent.analyses.AUC.(FP.info.groups(ii))=auc;
        FP.outputs.behaviourindependent.analyses.PeakFrequency.(FP.info.groups(ii))=pkfreq;
        FP.outputs.behaviourindependent.analyses.MeanPeakHeight.(FP.info.groups(ii))=meanpk;
        FP.outputs.behaviourindependent.analyses.MaxPeakHeight.(FP.info.groups(ii))=maxpk;
        FP.outputs.behaviourindependent.analyses.TimeMaxPeakHeight.(FP.info.groups(ii))=maxpktime;
        
    end
end


%clear remaining unnecessary variables
clear ii iii sig iso temp_x isofit sigfit dFF mintest adjust test cordFF
clear baseline idx1 idx2 meanbaseline time minpeak pks auc pkfreq meanpk
clear minpk locs maxpkidx maxpk maxpktime pick1 pick2 rawcordFF rawdFF

%process successful
disp('Process 3. Photometry Processing - Behaviour Independent Complete');

%%  4. Photometry Processing - Behaviour Epoch Dependent
%   
%   Integration of behavioural information into photometry analyses
%   Must run all previous processes before running this process
%
%   User-defined onset and duration of events
%       event.onset = [0 3600] (the onset time in seconds of the introduction
%       event)
%       event.duration = [600 600] (the duration time in seconds of the
%       introduction event)
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
FP.info.event.onset = 300; %start
FP.info.event.duration = 300; %duration

%batch processing
for ii=1:FP.info.groupnum(1)

        %create structure arrays for photometry outputs
        for iv = 1:FP.info.numberofbehav(1)
            FP.outputs.behaviourdependent.behaviour.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(FP.info.fpframes,FP.info.pergroup(ii));
            FP.outputs.behaviourdependent.dFF.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.duration/FP.info.samplingrate(1))),FP.info.pergroup(ii)); %placeholder array
            FP.outputs.behaviourdependent.analyses.AUC.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(1,FP.info.pergroup(ii));
            FP.outputs.behaviourdependent.analyses.PeakFrequency.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(1,FP.info.pergroup(ii));
            FP.outputs.behaviourdependent.analyses.MeanPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(1,FP.info.pergroup(ii));
            FP.outputs.behaviourdependent.analyses.MaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(1,FP.info.pergroup(ii));
            FP.outputs.behaviourdependent.analyses.TimeMaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=zeros(1,FP.info.pergroup(ii));
        end
        
        
    for iii=1:FP.info.pergroup(ii)
        %load data
        dFF=FP.outputs.behaviourindependent.traces.pickupdFF_notcorrected.(FP.info.groups(ii))(:,iii); %load not corrected dFF
        time=FP.inputs.doric.time.(FP.info.groups(ii))(:,iii);
        cordFF=FP.outputs.behaviourindependent.traces.pickupcorrecteddFF.(FP.info.groups(ii))(:,iii); %load corrected dFF
        
        %pickup window bounds
        [~,pick1]=min(abs(time-(FP.info.pickupwindow(1)))); %identify beginning of pickup
        [~,pick2]=min(abs(time-((FP.info.pickupwindow(2)+FP.info.pickupwindow(1))))); %identify end of pickup
        
        %exclude pickup window
        time(pick1:pick2)=[];
        
        %identify time segments of behavioural epochs
        [~,idx1]=min(abs(time-(FP.info.event.onset(1)))); %identify beginning of introduction event
        [~,idx2]=min(abs(time-((FP.info.event.duration(1)+FP.info.event.onset(1))))); %identify end of introduction event
        ecordFF=cordFF(idx1:idx2); %corrected dFF during the introduction event
        etime=time(idx1:idx2); %photometry time during introduction event
        
         
        
        for iv=1:FP.info.numberofbehav
        
        temp=FP.inputs.anymaze.data.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(:,iii); %load anymaze column

        %find indices of behaviours during epochs
        %stretch and align datasets
        %sets arrays to the same length and gap fills
        % *** NOTE: WILL NOT WORK PROPERLY IF THE FRAMERATE OF THE
        % BEHAVIOURAL RECORDING EXCEEDS THAT OF THE PHOTOMETRY RECORDING
        % ***
        eventtime=(FP.info.event.onset(1):(FP.info.event.duration(1)/size(temp,1)):(FP.info.event.onset(1)+FP.info.event.duration(1))); %dummy time series
        eventtime=eventtime(2:end)'; %trim starting zero
        tempstretch=zeros(size(etime,1),1); %generate sink to contain the new dataframes
        for v=1:size(etime,1) %batch one row at a time
            temp_time=etime(v);   
            [~,ix]=min(abs(eventtime-temp_time)); %find index to align closest timestamps
            tempstretch(v)=temp(ix,:); %apply index to head/torso
        end

        %find indices of bouts of behaviour of interest
        tempstretch=(tempstretch==1); %turns this into a logical array
        
        %save stretched binaries for visualization
        FP.outputs.behaviourdependent.behaviour.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))((idx1:idx2),iii)=tempstretch;
        
        %use incides of behaviours to isolate behaviour-specific photometry
        %signal
        %corrected
        temp=ecordFF(tempstretch);
        %not corrected
        nctemp=ecordFF(tempstretch);
        
        %write not corrected to structure for visualiaztion
        FP.outputs.behaviourdependent.analyses.dFF.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(1:size(nctemp,1),iii)=nctemp;
        
        %initialize peak analyses
        minpk=median(dFF)+(2*std(dFF)); %minimum height to for peak 
        
        %analyze and input data
        %familiar
        %head/torso
        FP.outputs.behaviourdependent.analyses.AUC.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(iii)=trapz(temp)/(size(temp,1)*(time(2)-time(1)));
        [pks,locs]=findpeaks(temp,'MinPeakHeight',(minpk),'MinPeakDistance',0);
        maxpkidx=find(cordFF==max(temp),1,'first'); %find index of max peak. If multiple, it will select first
        FP.outputs.behaviourdependent.analyses.MaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(iii)=max(temp); %max dFF
        FP.outputs.behaviourdependent.analyses.TimeMaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(iii)=time(maxpkidx); %time of max
        FP.outputs.behaviourdependent.analyses.PeakFrequency.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(iii)=size(pks,1)/(size(temp,1)*(time(2)-time(1)));
        FP.outputs.behaviourdependent.analyses.MeanPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(iii)=mean(cordFF(locs));
        end
    end
end


%clear remaining unnecessary variables
clear ii iii iv v ix
clear time famtime inttime dFF cordFF fcordFF icordFF
clear htf agf hti agi nchtf ncagf nchti ncagi htfstretch agfstretch htistretch agistretch
clear idxf1 idxf2 idxi1 idxi2 idx2 pick1 pick2
clear temp_time test minpk pks ftime itime locs maxpkidx
clear ag agstretch ecordFF etime eventtime ht htstretch idx1 ncag ncht nctemp temp tempstretch

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
            
        %create structure arrays for photometry outputs
        for iv = 1:FP.info.numberofbehav(1)
            FP.outputs.boutbybou.analyses.AUC.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.onset(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
            FP.outputs.boutbybou.analyses.PeakFrequency.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.onset(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
            FP.outputs.boutbybou.analyses.MeanPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.onset(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
            FP.outputs.boutbybou.analyses.MaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.onset(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
            FP.outputs.boutbybou.analyses.TimeMaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))=NaN((round(FP.info.event.onset(1)/FP.info.samplingrate(1))),FP.info.pergroup(ii));
        end
        
    for iii=1:FP.info.pergroup(ii)
        
        %load data
        dFF=FP.outputs.behaviourindependent.traces.pickupdFF_notcorrected.(FP.info.groups(ii))(:,iii); %load not corrected dFF
        time=FP.inputs.doric.time.(FP.info.groups(ii))(:,iii);
        cordFF=FP.outputs.behaviourindependent.traces.pickupcorrecteddFF.(FP.info.groups(ii))(:,iii); %load corrected dFF
        
        %pickup window bounds
        [~,pick1]=min(abs(time-(FP.info.pickupwindow(1)))); %identify beginning of pickup
        [~,pick2]=min(abs(time-((FP.info.pickupwindow(2)+FP.info.pickupwindow(1))))); %identify end of pickup
        
        %exclude pickup window
        time(pick1:pick2)=[];
        
        %identify time segments of introduction event
        [~,idx1]=min(abs(time-(FP.info.event.onset(1)))); %identify beginning of the event
        [~,idx2]=min(abs(time-((FP.info.event.duration(1)+FP.info.event.onset(1))))); %identify end of the event
        ecordFF=cordFF(idx1:idx2); %corrected dFF during the event
        edFF=cordFF(idx1:idx2); %dFF during the event
        etime=time(idx1:idx2); %photometry time during the event
        
        for iv=1:FP.info.numberofbehav
            
            temp=FP.inputs.anymaze.data.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(:,iii); %load anymaze column
        
            %find indices of behaviours during epochs
            %stretch and align datasets
            %sets arrays to the same length and gap fills
            % *** NOTE: WILL NOT WORK PROPERLY IF THE FRAMERATE OF THE
            % BEHAVIOURAL RECORDING EXCEEDS THAT OF THE PHOTOMETRY RECORDING
            % ***
            eventtime=(FP.info.event.onset(1):(FP.info.event.duration(1)/size(temp,1)):(FP.info.event.onset(1)+FP.info.event.duration(1))); %dummy time series
            eventtime=eventtime(2:end)'; %trim starting zero
            tempstretch=zeros(size(etime,1),1); %generate sink to contain the new dataframes
            for v=1:size(etime,1) %batch one row at a time
                temp_time=etime(v);   
                [~,ix]=min(abs(eventtime-temp_time)); %find index to align closest timestamps
                tempstretch(v)=temp(ix,:); %apply index to head/torso
            end

            %find indices of bouts of behaviour of interest
            tempstretch=(tempstretch==1); %turns this into a logical array
        
            %initialize peak analyses
            minpk=median(dFF)+(2*std(dFF)); %minimum height to for peak 
        
            %find indices of when behaviours change
            %assumes that the first frame will be a 0, meaning that there is no
            %interaction occuring

            %initialization
            temp1=find(diff(tempstretch)==0);
            behavchange=setdiff(min(temp1):max(temp1),temp1); %find indices of changes in behavioural state
            behavlist=[1,behavchange,size(tempstretch,1)];
            behavstart=behavlist(1:2:end); %start of a behavioural sequence
            behavend=behavlist(2:2:end);  %end of a behavioural sequence
            if size(behavstart,2)>size(behavend,2) %safety check to ensure that arrays are of the same length
                behavstart=behavstart(1:size(behavend,2));
            else
                behavend=behavend(1:size(behavstart,2));
            end
        
            %batch through bouts 
            for v=1:(size(behavstart,2)-1)
                tempdFF=edFF((behavend(v)+1):behavstart(v+1));
                tempdFF2=ecordFF(behavend(v)+1:behavstart(v+1));
                FP.outputs.boutbybou.analyses.AUC.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(v,iii)=(trapz(tempdFF)/(size(tempdFF,1)*FP.info.samplingrate(1)));
                [pks,locs]=findpeaks(tempdFF,'MinPeakHeight',(minpk),'MinPeakDistance',0);
                maxpkidx=find(cordFF==max(tempdFF),1,'first'); %find index of max peak. If multiple, it will select first
                FP.outputs.boutbybou.analyses.MaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(v,iii)=max(tempdFF); %max dFF
                FP.outputs.boutbybou.analyses.TimeMaxPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(v,iii)=time(maxpkidx); %time of max
                FP.outputs.boutbybou.analyses.PeakFrequency.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(v,iii)=(size(pks,1)/(size(tempdFF,1)*FP.info.samplingrate(1)));
                FP.outputs.boutbybou.analyses.MeanPeakHeight.(FP.info.groups(ii)).(FP.info.anymaze.titles(iv))(v,iii)=mean(tempdFF2(locs));
            end
        end
    end
end

%clear remaining unnecessary variables
clear ii iii iv v ix
clear time famtime inttime dFF cordFF fcordFF icordFF fdFF idFF
clear htf agf hti agi nchtf ncagf nchti ncagi htfstretch agfstretch htistretch agistretch
clear idxf1 idxf2 idxi1 idxi2 idx2 temp tempstretch
clear temp_time test minpk pks ftime itime temp1 tempdFF2
clear behavchange behavend behavstart behavlist tempdFF locs maxpkidx
clear ag agstretch ecordFF edFF etime eventtime ht htstretch idx1 pick1 pick2

%process successful
disp('Process 5. Photometry Processing - Bout-by-Bout Information Complete');

%%  X. Save Structure to .mat File
%keep file for easy data access later

save('FP.mat','-struct','FP');

disp('Process X. Save Structure to a .mat File Complete');

%%  Y. Load .mat File
%   navigate to folder containing 'FP.mat'

load('FP.mat')
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

%indicate start and duration of pickup window
FP.info.pickupwindow = [295 10];
%indicate whether z scores are wanted
FP.info.zscore = 0;
%indicate whether %dF/F is wanted
FP.info.percentdFF = 0;
%indicate whether analysis should be restricted to a specific time period
FP.info.timelimit = 0;
FP.info.timerange = [300 60];  %[start duration]

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
%trim if exceeds maximum trial duration (default = FP.info.fpframes, aka 15mins
%@120fps)
if size(doricdata,1)>FP.info.fpframes %
   doricdata=doricdata(1:FP.info.fpframes,:); %
else
   %do nothing
end


%organize raw data
%adjust the number of frames (default = FP.info.fpframes, aka 15mins @120fps) as needed
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


%pickup window bounds
[~,pick1]=min(abs(time-(FP.info.pickupwindow(1)))); %identify beginning of pickup
[~,pick2]=min(abs(time-((FP.info.pickupwindow(2)+FP.info.pickupwindow(1))))); %identify end of pickup
        
%baseline adjustment
[~,idx1]=min(abs(time-(FP.info.baseline(1)))); %identify beginning of baseline
[~,idx2]=min(abs(time-((FP.info.baseline(2)+FP.info.baseline(1))))); %identify end of baseline
if pick1<idx2
    idx2=pick1;
else
    %do nothing
end
baseline=dFF(idx1:idx2); %baseline period
test=rawdFF(pick2:end); %test period
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

%exclude pickup window
time(pick1:pick2)=[];
dFF(pick1:pick2)=[];
cordFF(pick1:pick2)=[];    
    
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
minpk=median(dFF)+(2*std(dFF)); %minimum height to for peak 
%consideration is set to two standard deviations above the median 
%dFF value for the test segment of the recording
[pks,locs]=findpeaks(dFF,'MinPeakHeight',(minpk),'MinPeakDistance',0);
locs=cordFF(locs);
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
clear test time locs maxpkidx maxpk maxpktime pick1 pick2
