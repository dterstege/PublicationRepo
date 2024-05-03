%%                          FP3002_v060.m
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
%       Neurophotometrics Fiber Photometry Output (.csv)
%       ANYmaze Outputs (.csv)
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.1
%
%   V1.0.1  - Basic Behavioural Alignment. Assumes consistent framerate
%           - Correlating FP signals using a rolling window approach
%
%   Created 01/07/2024 Dylan Terstege (https://www.github.com/dterstege)
%   Epp Lab, University of Calgary (https://www.epplab.com)
%   Contact: dylan.terstege@ucalgary.ca

%%  1. Photometry Processing
%
%   DESCRIPTION
%       Very simple analysis
%       Delta F/F
%       All data stored to "FP" structure element
%       Run on individual animals
%
%   INPUTS
%       Decide whether data should be analyzed as a z-score delta F/F (1)
%       or just as delta F/F (0)
        FP.inputs.zscore = 0;                   %set to 0 or 1

%       Specify LED States used in the recording
%       NOTE: Always have 415 (1) first with 470 (2) and/or 560 (4) after
        FP.inputs.LEDstates = [1 2];            %LED states in file
        FP.inputs.LEDID = ["iso" "GCaMP"];     %LED state IDs in file
%       NOTE: ALWAYS iso FIRST
%
%   OUTPUTS
%       fullrecording: data regardless of TTL state
%           Timestamp from start of the Bonsai recording
%       ttlon: data between TTL pulses (between when ANYmaze is turned on
%       and turned off)
%           Timestamp synched to TTL input
%       ttloff: data before (and after) the ANYmaze trial
%           Timestamps relative to the start of the Bonsai recording 
%
%                           load data
%
disp("Select parent folder");                   %identify parent folder
FP.ID=uigetdir();                               %set directory

path=fullfile(FP.ID,'/NPM.csv');                 %input data
data=readmatrix(path);                          %extract data

for ii = 1:size(FP.inputs.LEDstates,2)
    temp=find(data(:,3)==FP.inputs.LEDstates(ii));
    isolated=data(temp,:);
    FP.inputs.(FP.inputs.LEDID(ii)).systemtime=isolated(:,2);
    FP.inputs.(FP.inputs.LEDID(ii)).computertime=isolated(:,4);
    FP.inputs.(FP.inputs.LEDID(ii)).rawfluor=isolated(:,5);
end

fptime=FP.inputs.(FP.inputs.LEDID(ii)).systemtime(:,1);
fptime=fptime-fptime(1);                        %normalize
FP.inputs.fpfr=round(1/fptime(2),0);            %frame rate of FP

if size(FP.inputs.LEDstates,2) == 2
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
    
    FP.outputs.fullrecording.dFF=dFF;
    FP.outputs.fullrecording.time=fptime;
else
    iso=FP.inputs.(FP.inputs.LEDID(1)).rawfluor;
    first_sig=FP.inputs.(FP.inputs.LEDID(2)).rawfluor;
    second_sig=FP.inputs.(FP.inputs.LEDID(3)).rawfluor;
    
    %make sure inputs are the same lengths
    n=min(min(numel(iso),numel(first_sig)),numel(second_sig));
    iso=iso(1:n);
    first_sig=first_sig(1:n);
    second_sig=second_sig(1:n);
    
    %filter noise
    iso=filloutliers(iso,"center","movmedian",3);
    first_sig=filloutliers(first_sig,"center","movmedian",3);
    second_sig=filloutliers(second_sig,"center","movmedian",3);
    
    temp_x=1:length(iso);                       %note -- using a temporary x variable due to computational constraints in matlab
    temp_x=temp_x';                             %transpose
    isofit=fit(temp_x,iso,'exp2');              %fit isosbestic data with biexponential decay
    
    %first
    sigfit=robustfit(isofit(temp_x),first_sig); %scale biexponential decay to the raw 470 data using robust fit
    sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

    dFF=(first_sig-sigfit)./sigfit;             %delta F/F
    
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
    
    FP.outputs.fullrecording.(FP.inputs.LEDID(2)).dFF=dFF;
    FP.outputs.fullrecording.time=fptime;
    
    %second
    sigfit=robustfit(isofit(temp_x),second_sig); %scale biexponential decay to the raw 470 data using robust fit
    sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

    dFF=(second_sig-sigfit)./sigfit;             %delta F/F
    
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
    
    
    FP.outputs.fullrecording.(FP.inputs.LEDID(3)).dFF=dFF;
end

clearvars -except FP

disp('Process 1. Photometry Processing Complete');

%%  2. Photometry Segmentation
%
%   DESCRIPTION
%       Splits photometry trace based on TTL inputs
%
%   INPUTS
%       Specify the number/duration (in seconds) of each segment of interest
%       For example: for 2 bouts that are 5 mins each, specify [300 300]
        FP.inputs.segment = [300];

%       Specify the identifier for each bout
        FP.inputs.segmentID = ["YMaze"];
%
%   OUTPUTS
%

%load data
path=fullfile(FP.ID,'/IO.csv');                 %input data
data=readmatrix(path);                          %extract data

systemtime=FP.inputs.iso.systemtime;

for ii=1:size(FP.inputs.segment,2)
    %find bounds
    IOtime=data(ii,4);
    [~,ix]=min(abs(systemtime-IOtime));  %find index to align closest timestamps
    dur=FP.inputs.segment(ii);
    dur=FP.inputs.fpfr*dur;
    
    if size(FP.inputs.LEDstates,2) == 2
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
        
        FP.outputs.(FP.inputs.segmentID(ii)).dFF=dFF;
        time=(0:(1/FP.inputs.fpfr):(FP.inputs.segment(ii)))';
        time=time(1:dur);
        FP.outputs.(FP.inputs.segmentID(ii)).time=time;
        
    elseif size(FP.inputs.LEDstates,2) == 3
        iso=FP.inputs.(FP.inputs.LEDID(1)).rawfluor(ix:(ix+dur-1));
        first_sig=FP.inputs.(FP.inputs.LEDID(2)).rawfluor(ix:(ix+dur-1));
        second_sig=FP.inputs.(FP.inputs.LEDID(3)).rawfluor(ix:(ix+dur-1));

        %make sure inputs are the same lengths
        n=min(min(numel(iso),numel(first_sig)),numel(second_sig));
        iso=iso(1:n);
        first_sig=first_sig(1:n);
        second_sig=second_sig(1:n);

        %filter noise
        iso=filloutliers(iso,"center","movmedian",3);
        first_sig=filloutliers(first_sig,"center","movmedian",3);
        second_sig=filloutliers(second_sig,"center","movmedian",3);
        
        
        temp_x=1:length(iso);                       %note -- using a temporary x variable due to computational constraints in matlab
        temp_x=temp_x';                             %transpose
        isofit=fit(temp_x,iso,'exp2');              %fit isosbestic data with biexponential decay

        %first
        sigfit=robustfit(isofit(temp_x),first_sig); %scale biexponential decay to the raw 470 data using robust fit
        sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

        dFF=(first_sig-sigfit)./sigfit;             %delta F/F

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
        
        FP.outputs.(FP.inputs.segmentID(ii)).(FP.inputs.LEDID(2)).dFF=dFF;
        time=(0:(1/FP.inputs.fpfr):(FP.inputs.segment(ii)))';
        time=time(1:dur);
        FP.outputs.(FP.inputs.segmentID(ii)).time=time;
        
        sigfit=robustfit(isofit(temp_x),second_sig); %scale biexponential decay to the raw 470 data using robust fit
        sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

        dFF=(second_sig-sigfit)./sigfit;             %delta F/F

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

        FP.outputs.(FP.inputs.segmentID(ii)).(FP.inputs.LEDID(3)).dFF=dFF;
        time=(0:(1/FP.inputs.fpfr):(FP.inputs.segment(ii)))';
        time=time(1:dur);
        FP.outputs.(FP.inputs.segmentID(ii)).time=time;
        
    end
end
clearvars -except FP

disp('Process 2. Photometry Segmentation Complete');

%%  3. Behaviour Integration - written specific for ymaze
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
%       specify the number of columns containing ANYmaze data which you
%       would like to align
        FP.inputs.ANYmaze.UserSpecs.numbehavs = 3;
%       specify which columns to consider in the ANYmaze data file (time
%       excluded)
        FP.inputs.ANYmaze.UserSpecs.columns = [2 3 4];
%       specify what behaviour each column represents
        FP.inputs.ANYmaze.UserSpecs.behavs = ["xpos" "ypos" "middle"];
%       specify the alignment method for each column
%           linear: linear interpolation between known datapoints
%           nn: nearest neighbours. Works well with binary data
        FP.inputs.ANYmaze.UserSpecs.alignments = ["linear" "linear" "nn"];
        
%load data
ampath=fullfile(FP.ID,'/ANYmaze.csv');                              %ANYmaze file
am=readmatrix(ampath);                                              %extract ANYmaze
am(isnan(am))=0;                                                    %replace empty columns with zero

%align timestamps of behaviour
for ii=2%:size(FP.inputs.segment,2) %ignore complex maze for now
    fptime=FP.outputs.(FP.inputs.segmentID(ii)).time;
    amtime=am(1:end,1);
    
    for iii = 1:FP.inputs.ANYmaze.UserSpecs.numbehavs
        temp=am(:,FP.inputs.ANYmaze.UserSpecs.columns(iii));
        FP.inputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(iii)).raw=temp;
        tempstretch=nan(length(fptime),1);                              %storage array
        if strcmp(FP.inputs.ANYmaze.UserSpecs.alignments(iii),"nn")
            for iv=1:size(fptime,1)                                    %batch one row at a time
                temp_time=fptime(iv);   
                [~,ix]=min(abs(amtime-temp_time));                      %find index to align closest timestamps
                tempstretch(iv)=temp(ix,:);                            %apply index to behaviour
            end
            FP.outputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(iii)).stretched=tempstretch;
        elseif strcmp(FP.inputs.ANYmaze.UserSpecs.alignments(iii),"linear")
            tempstretch=interp1(amtime,temp,fptime);                    %linear interpolation
            FP.outputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(iii)).stretched=tempstretch;
        else
            disp('Please ensure that alignment type is correctly specified');
        end
    end
    
    if size(FP.inputs.LEDstates,2) == 2
        dFF=FP.outputs.(FP.inputs.segmentID(ii)).dFF;
        mid=FP.outputs.ANYmaze.middle.stretched;
        temp=logical(mid);
        dFF_on=dFF(temp);
        dFF_off=dFF(~temp);
        time_on=size(fptime(temp),1)*(1/FP.inputs.fpfr);
        time_off=size(fptime(~temp),1)*(1/FP.inputs.fpfr);
        FP.outputs.ANYmaze.middle.AUC.on=(trapz(dFF_on)/time_on);
        FP.outputs.ANYmaze.middle.AUC.off=(trapz(dFF_on)/time_off);
        FP.outputs.ANYmaze.middle.RANGE.on=abs(max(dFF_on)-min(dFF_on));
        FP.outputs.ANYmaze.middle.RANGE.off=abs(max(dFF_off)-min(dFF_off));
        
    elseif size(FP.inputs.LEDstates,2) == 3
        dFF=FP.outputs.(FP.inputs.segmentID(ii)).(FP.inputs.LEDID(2)).dFF;
        mid=FP.outputs.ANYmaze.middle.stretched;
        temp=logical(mid);
        dFF_on=dFF(temp);
        dFF_off=dFF(~temp);
        time_on=size(fptime(temp),1)*(1/FP.inputs.fpfr);
        time_off=size(fptime(~temp),1)*(1/FP.inputs.fpfr);
        FP.outputs.ANYmaze.middle.AUC.on=(trapz(dFF_on)/time_on);
        FP.outputs.ANYmaze.middle.AUC.off=(trapz(dFF_on)/time_off);
        FP.outputs.ANYmaze.middle.RANGE.on=abs(max(dFF_on)-min(dFF_on));
        FP.outputs.ANYmaze.middle.RANGE.off=abs(max(dFF_off)-min(dFF_off));
        
    end
    %do something
end

clearvars -except FP

disp('Process 3. Behaviour Integration Complete');

%%  X. Save Structure to .mat File
%keep file for easy data access later

save('NPM.mat','-struct','FP');

disp('Process X. Save Structure to a .mat File Complete');

%%  Y. Load .mat File
%   navigate to folder containing 'FP.mat'

load('NPM.mat')
 w=whos;
 for ii=1:length(w)
     FP.(w(ii).name)=eval(w(ii).name); 
     clear w(ii).name
 end
 clear w info inputs outputs ii
 
 disp('Process Y. Load .mat File Complete');