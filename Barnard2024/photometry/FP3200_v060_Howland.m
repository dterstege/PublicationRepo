%%                          FP3002_v060_Howland.m
%
%   Fiber Photometry Analysis
%
%   Written for the lab of John Howland - University of Saskatchewan
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
%       Neurophotometrics Fiber Photometry Output (.csv)
%       ANYmaze Outputs (.csv)
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.1
%
%   V1.0.1  - Basic Processing
%
%   Created 08/28/2024 Dylan Terstege (https://www.github.com/dterstege)
%   Epp Lab, University of Calgary (https://www.epplab.com)
%   Contact: dylan.terstege@ucalgary.ca

clearvars

%%  1. Photometry Processing
%
%   DESCRIPTION
%       The purpose of this process is mainly for data
%       loading/initialization
%
%       If synchronizing with ANYmaze, this process does not consider the 
%       time when any TTL pulse occurred
%
%       This trace will include any baseline/post-testing periods
%
%       All data stored to "FP" structure element
%
%       Run on individual animals
%
%   INPUTS
%       Decide whether data should be analyzed as a z-score delta F/F (1)
%       or just as delta F/F (0)
        FP.inputs.zscore = 1;                   %set to 0 or 1

%       Specify LED States used in the recording
%       NOTE: Always have 415 (1) first with 470 (2) and/or 560 (4) after
        FP.inputs.LEDstates = [1 2];            %LED states in file
        FP.inputs.LEDID = ["iso" "GCaMP7f"];     %LED state IDs in file
%       NOTE: ALWAYS iso FIRST
%
%   OUTPUTS
%   All outputs from this process can be found in the structure tree under
%   FP.outputs.process1
%       fullrecording: data regardless of TTL state
%           Timestamp from start of the Bonsai recording
%           Includes baseline and post-testing periods
%       ttlon: data between TTL pulses (between when ANYmaze is turned on
%       and turned off)
%           Timestamp synched to TTL input
%       ttloff: data before (and after) the ANYmaze trial
%           Timestamps relative to the start of the Bonsai recording 
%
%                           load data
disp("Select parent folder");                   %identify parent folder
FP.ID=uigetdir();                               %set directory

path=fullfile(FP.ID,'/NPM.csv');                %input data
data=readmatrix(path);                          %extract data

for ii = 1:size(FP.inputs.LEDstates,2)
    temp=find(data(:,3)==FP.inputs.LEDstates(ii));
    isolated=data(temp,:);
    FP.inputs.(FP.inputs.LEDID(ii)).systemtime=isolated(:,2);
    FP.inputs.(FP.inputs.LEDID(ii)).computertime=isolated(:,4);
    FP.inputs.(FP.inputs.LEDID(ii)).rawfluor=isolated(:,5); %the index at this line should correspond to the green ROI of a given fibre (ex: 5 for G0)

end

fptime=FP.inputs.(FP.inputs.LEDID(1)).systemtime(:,1);
fptime=fptime-fptime(1);                        %normalize
FP.inputs.fpfr=round(1/fptime(2),0);            %frame rate of FP

iso=FP.inputs.(FP.inputs.LEDID(1)).rawfluor;
sig=FP.inputs.(FP.inputs.LEDID(2)).rawfluor;

%make sure inputs are the same lengths
n=min(numel(iso),numel(sig));
iso=iso(1:n);
sig=sig(1:n);

%filter noise
iso=filloutliers(iso,"center","movmedian",3);
sig=filloutliers(sig,"center","movmedian",3);

temp_x=1:length(iso);                       %note -- using a temporary x variable due to computational constraints in matlab
temp_x=temp_x';                             %transpose
isofit=fit(temp_x,iso,'exp2');              %fit isosbestic data with biexponential decay

sigfit=robustfit(isofit(temp_x),sig);       %scale biexponential decay to the raw 470 data using robust fit
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

dFF=(sig-sigfit)./sigfit;                   %delta F/F

%               z-score or omit
if FP.inputs.zscore==1
    %convert data to z-score dF/F
    dFF=zscore(dFF);
else
    %leave data as is
end

%lowpass - testing this out
fs=1e3;
dFF=lowpass(dFF,20,fs);

FP.outputs.process1.fullrecording.dFF=dFF;
FP.outputs.process1.fullrecording.time=fptime;

clearvars -except FP

disp('Process 1. Photometry Processing Complete');

%%  2. Photometry Segmentation
%
%   DESCRIPTION
%       Splits photometry trace based on segments of the test
%
%       User defines these segments. For most tasks, we will only have a
%       single segment. However, we can create as many segments as needed
%
%       Here, we will examine dF/F across the entire segment
%
%   INPUTS
%       Specify the number/duration (in seconds) of each segment of interest
%       For example: for 2 bouts that are 5 mins each, specify [300 300]
        FP.inputs.segment = [300 300];
%        
%       Specify the start time (in seconds) of each segment of interest
%       For example: a bout stating when the test starts and a second bout
%       starting at the 6 minute mark [0 360]
        FP.inputs.segmenttime = [0 360];

%       Specify the identifier for each bout
%       These should align with the names used in the BEHAVIOUR.xlsx file
        FP.inputs.segmentID = ["S" "T"];
%
%   OUTPUTS
%   All outputs from this process can be found in the structure tree under
%   FP.outputs.process2
%       dFF: dF/F by LED state
%           This dF/F has been processed using a low pass filter, similar
%           to what is used by the DORIC system. 
%       time: timstamp vector for photometry values during the user defined
%       segments
%

systemtime=FP.inputs.iso.systemtime;
systemtime=systemtime-systemtime(1);     %normalize to start at 0

for ii=1:size(FP.inputs.segment,2)
    %find bounds
    IOtime=FP.inputs.segmenttime(ii);
    [~,ix]=min(abs(systemtime-IOtime));  %find index to align closest timestamps
    dur=FP.inputs.segment(ii);
    dur=FP.inputs.fpfr*dur;
    
    iso=FP.inputs.(FP.inputs.LEDID(1)).rawfluor(ix:(ix+dur-1));
    sig=FP.inputs.(FP.inputs.LEDID(2)).rawfluor(ix:(ix+dur-1));
    
    %make sure inputs are the same lengths
    n=min(numel(iso),numel(sig));
    iso=iso(1:n);
    sig=sig(1:n);
    
    %filter noise
    iso=filloutliers(iso,"center","movmedian",3);
    sig=filloutliers(sig,"center","movmedian",3);
    
    temp_x=1:length(iso);                       %note -- using a temporary x variable due to computational constraints in matlab
    temp_x=temp_x';                             %transpose
    isofit=fit(temp_x,iso,'exp2');              %fit isosbestic data with biexponential decay
    
    sigfit=robustfit(isofit(temp_x),sig);       %scale biexponential decay to the raw 470 data using robust fit
    sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);
    
    dFF=(sig-sigfit)./sigfit;                   %delta F/F
    
    %               z-score or omit
    if FP.inputs.zscore==1
        %convert data to z-score dF/F
        dFF=zscore(dFF);
    else
        %leave data as is
    end
    
    %filter noise
    dFF=filloutliers(dFF,"center","movmedian",3);
    
    %lowpass - testing this out
    fs=1e3;
    dFF=lowpass(dFF,20,fs);
    
    FP.outputs.process2.(FP.inputs.segmentID(ii)).dFF=dFF;
    time=(0:(1/FP.inputs.fpfr):(FP.inputs.segment(ii)))';
    time=time(1:dur);
    FP.outputs.process2.(FP.inputs.segmentID(ii)).time=time;
    
end
clearvars -except FP

disp('Process 2. Photometry Segmentation Complete');

%%  3. Behaviour Integration
%
%   Integration of behavioural information into photometry analyses
%   Must run all previous processes before running this process
%
%   User-defined variables:
%       timebeforebehav = 2; time in seconds to consider before zone entry
%       behavtime = 5; time in seconds to consider after zone entry.  FP
%           data beyond this time will be omitted.  Bouts shorter than this
%           will not be considered.
%
%   Operational Definitions:
%       dFF:
%           Isosbestic data is fit to a biexponential decay
%           This biexponential decay is linearly scaled to the raw 470 data
%           (raw = sig, fitted = sigfit)
%           dFF = (sig-sigfit)./sigfit
%
%   NaN used as a placeholder for behaviour-specific dF/F traces. The trace
%   is limited to any values prior to the first NaN occurence

%   INPUTS
%       Specify the position of the novel odor
        FP.inputs.behaviour.novelposition = 5;
        
%       Specify minimum duration of interaction (in seconds)
        FP.inputs.behaviour.mindur = 0.5;
        
%       Specify duration prior to interation to be considered (in seconds)
%       NOTE: timeseries must have at least this much time before first
%       interaction (can omit first interaction in data file otherwise)
        FP.inputs.behaviour.timeprior = 2;
        
%       Specify duration after interation onset to be considered (in seconds)
%       NOTE: timeseries must have at least this much time after the onset
%       of the interaction (can omit last interaction in data file
%       otherwise)
        FP.inputs.behaviour.timeafter = 2;

%       Specify peak detection criteria for this process
        FP.inputs.peakdetection.process3 = ["stdev", 0.5]; %set first value to "setdev" or "mad"
        
%   OUTPUTS
%   All outputs from this process can be found in the structure tree under
%   FP.outputs.process3
%       dFF: dF/F by LED state and behaviour type
%           This dF/F is lifted directly from what was processed in process
%           2 and therefore has been processed by the same low pass filter
        
%load data
bpath=fullfile(FP.ID,'/BEHAVIOUR.xlsx');                           %Behaviour file
phase=readcell(bpath,'Range','E:E');
FP.inputs.behaviour.phase=phase(4:end);
start_time=readmatrix(bpath,'Range','F:F');
FP.inputs.behaviour.start_time=start_time;
end_time=readmatrix(bpath,'Range','G:G');
FP.inputs.behaviour.end_time=end_time;
odor=readcell(bpath,'Range','I:I');
FP.inputs.behaviour.odor=odor(4:end);
position=readmatrix(bpath,'Range','J:J');
FP.inputs.behaviour.position=position;
FP.inputs.behaviour.boutduration=FP.inputs.behaviour.end_time-FP.inputs.behaviour.start_time;

systemtime=FP.inputs.iso.systemtime;
systemtime=systemtime-systemtime(1);     %normalize to start at 0

%identify bouts meeting minimum duration
idx=find(FP.inputs.behaviour.boutduration>=FP.inputs.behaviour.mindur);
start_time=FP.inputs.behaviour.start_time(idx);
end_time=FP.inputs.behaviour.end_time(idx);
phase=FP.inputs.behaviour.phase(idx);
odor=FP.inputs.behaviour.odor(idx);
position=FP.inputs.behaviour.position(idx);

%storage vector
framesbefore=FP.inputs.behaviour.timeprior*FP.inputs.fpfr;
framesafter=FP.inputs.behaviour.timeafter*FP.inputs.fpfr;

outbin=nan((framesbefore+1+framesafter),size(idx,1));
peakfrequency=nan(size(idx,1),1);
meanpeakheight=nan(size(idx,1),1);
dynamicrange=nan(size(idx,1),1);
maximumvalue=nan(size(idx,1),1);
AUC=nan(size(idx,1),1);

%duration
bouttime=FP.inputs.behaviour.timeprior+FP.inputs.behaviour.timeafter+(1/FP.inputs.fpfr);

%sorting and analyses
for ii=1:size(idx,1)      
    if char(phase(ii)) == 'T' %make sure this aligns with the later of your two segments
        [~,startframe]=min(abs(systemtime-start_time(ii)));  %find index to align closest timestamps        
        [~,endframe]=min(abs(systemtime-end_time(ii)));  %find index to align closest timestamps
        startframe=startframe-(FP.inputs.segmenttime(2)*FP.inputs.fpfr);
    else
        [~,startframe]=min(abs(systemtime-start_time(ii)));  %find index to align closest timestamps
        [~,endframe]=min(abs(systemtime-end_time(ii)));  %find index to align closest timestamps
    end
    
    dFF=FP.outputs.process2.(char(phase(ii))).dFF;
    trace=dFF((startframe-framesbefore):(startframe+framesafter));
    outbin(:,ii)=trace;
    
    %peaks detection threshold
    if strcmp(FP.inputs.peakdetection.process3(1),"stdev")
        thr=str2double(FP.inputs.peakdetection.process3(2));
        medval=median(dFF);
        correction=std(dFF);
        thresh=medval+thr*correction;
    elseif strcmp(FP.inputs.peakdetection.process3(1),"mad")
        thr=str2double(FP.inputs.peakdetection.process3(2));
        medval=median(dFF);
        correction=mad(dFF);
        thresh=medval+thr*correction;
    else
        disp("Error: Please review peak detection input. Must be either stdev or mad.");
        return
    end
    
    %peak detection and analyses
    [pks,locs]=findpeaks(trace,'MinPeakHeight',thresh);
    peakfrequency(ii)=size(pks,1)/bouttime; %total peak frequency
    meanpeakheight(ii)=mean(pks); %mean peak height
    dynamicrange(ii)=abs(max(trace)-min(trace)); %dynamic range of the trace
    maximumvalue(ii)=max(trace); %maximum value in the trace
    AUC(ii)=trapz(trace); %area under the curve    
end

%storing data
FP.outputs.process3.traces=outbin;
FP.outputs.process3.analysistable=table(phase,odor,position,peakfrequency,meanpeakheight,dynamicrange,maximumvalue,AUC);

%clearvars -except FP

disp('Process 3. Behaviour Integration Complete');

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