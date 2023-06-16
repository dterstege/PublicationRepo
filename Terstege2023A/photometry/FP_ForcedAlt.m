%%                              FP_ForcedAlt.m
%
%   Fiber Photometry Analysis - Forced Alternation
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
%       415 - Neurophotometrics Fiber Photometry Output (.csv)
%       470 - Neurophotometrics Fiber Photometry Output (.csv)
%       ANYmaze Outputs (.csv)
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.1
%
%   V1.0.1  - Basic Behavioural Alignment. Assumes consistent framerate
%
%   Created 05/27/2022 Dylan Terstege (https://www.github.com/dterstege)
%   Epp Lab, University of Calgary (https://www.epplab.com)
%   Contact: dylan.terstege@ucalgary.ca

clearvars

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
        FP.inputs.zscore = 0;                                       %set to 0 or 1
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
disp("Select parent folder");                                       %identify parent folder
FP.ID=uigetdir();                                                   %set directory

isopath=fullfile(FP.ID,'/415.csv');                                 %isosbestic file
greenpath=fullfile(FP.ID,'/470.csv');                               %470 signal file

iso=readmatrix(isopath);                                            %extract isosbestic
green=readmatrix(greenpath);                                        %extract 470 signal

fptime=iso(:,2);                                                    %time series
fptime=fptime-fptime(1);                                            %normalize
FP.inputs.fpfr=round(1/fptime(2),0);                                %frame rate of FP
flag=iso(:,7);                                                      %TTL flag
g_iso=iso(:,9);                                                     %iso fluorescence
green=green(:,9);                                                   %470 fluorescence

%               make all vectors the same length
%
%recordings often result in an incongruent numbers of frames across
%channels depending on which channel was being recorded from when the trial
%ended
if length(g_iso)>length(green)                                      %if iso is longer than green
    g_iso=g_iso(1:length(green));                                   %make iso same size as green
    fptime=fptime(1:length(green));                                 %make time same size as green
    flag=flag(1:legnth(green));                                     %make flag same size as green
else
    green=green(1:length(g_iso));                                   %make green same size as iso
end

%              storing variables to a data structure
%
FP.inputs.time=fptime;                                              %full time
FP.inputs.flags=flag;                                               %full flags
FP.inputs.isosbestic415=iso;                                        %full iso
FP.inputs.signal470=green;                                          %full 470

%                       overall dF/F calculations
%
%470
temp_x=1:length(g_iso);                                             %note -- using a temporary x variable due to computational constraints in matlab
temp_x=temp_x';                                                     %transpose
isofit=fit(temp_x,g_iso,'exp2');                                    %fit isosbestic data with biexponential decay

sigfit=robustfit(isofit(temp_x),green);                             %scale biexponential decay to the raw 470 data using robust fit
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

green_dFF=(green-sigfit)./sigfit;                                   %delta F/F


%               z-score or omit
if FP.inputs.zscore==1
    %convert data to z-score dF/F
    green_dFF=zscore(green_dFF);                                    %optional conversion of dF/F to z-score dF/F
else
    %leave data as is
end


%               storing variables to a data structure
%
FP.outputs.fullrecording.time=fptime;                               %full time
FP.outputs.fullrecording.dFF470=green_dFF;                          %full 470 dF/F

%               data during flag states
%
uniqueflags=unique(flag);                                           %identify all flag states
ttl_off=find(flag==uniqueflags(1));                                 %index before TTL (before trial)
ttl_on=find(flag==uniqueflags(2));                                  %index after TTL (during trial)

%non trial data
fptime_pretrial=fptime(ttl_off);                                    %timestamps of recordings during ttl off cycle
iso_pretrial=iso(ttl_off,:);                                        %isosbestic during ttl off cycle
green_pretrial=green(ttl_off,:);                                    %470 channel during ttl off cycle

%470
temp_x=1:length(iso_pretrial);                                      %note -- using a temporary x variable due to computational constraints in matlab
temp_x=temp_x';                                                     %transpose
isofit=fit(temp_x,iso_pretrial(:,1),'exp2');                        %fit isosbestic data with biexponential decay

sigfit=robustfit(isofit(temp_x),green_pretrial(:,1));               %scale biexponential decay to the raw 470 data using robust fit
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);                          

green_dFF=(green_pretrial(:,1)-sigfit)./sigfit;                     %delta F/F


%               z-score or omit
if FP.inputs.zscore==1
    %convert data to z-score dF/F
    green_dFF=zscore(green_dFF);                                    %optional conversion of dF/F to z-score dF/F
else
    %leave data as is - dF/F
end

%               storing variables to a data structure
%
FP.outputs.ttloff.time=fptime_pretrial;                             %full time
FP.outputs.ttloff.dFF470=green_dFF;                                 %full 470 dF/F

%trial data
fptime_trial=fptime(ttl_on);                                        %timestamps of recordings during ttl on cycle
iso_trial=iso(ttl_on,:);                                            %isosbestic during ttl on cycle
green_trial=green(ttl_on,:);                                        %470 channel during ttl on cycle

%470
temp_x=1:length(iso_trial);                                         %note -- using a temporary x variable due to computational constraints in matlab
temp_x=temp_x';                                                     %transpose
isofit=fit(temp_x,iso_trial(:,1),'exp2');                           %fit isosbestic data with biexponential decay

sigfit=robustfit(isofit(temp_x),green_trial(:,1));                  %scale biexponential decay to the raw 470 data using robust fit
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

green_dFF=(green_trial(:,1)-sigfit)./sigfit;                        %delta F/F

%               z-score or omit
if FP.inputs.zscore==1
    %convert data to z-score dF/F
    green_dFF=zscore(green_dFF);
else
    %leave data as is - dF/F
end

%               storing variables to a data structure
%
FP.outputs.ttlon.time=fptime_trial;                                 %full time
FP.outputs.ttlon.dFF470=green_dFF;                                  %full 470 dF/F

clearvars -except FP                                                %clean variable space

disp('Process 1. Photometry Processing Complete');

%%  2. ANYmaze Alignment
%
%   *NOTE* MUST BE RAN AFTER PROCESS 1
%
%   DESCRIPTION
%       Very basic alignment of ANYmaze data to the photometry data
%       Requires user to input specific vairables in the below INPUT
%       section
%       Requires that Process 1. Photometry Processing has previously been
%       ran (FP structure element must be in MATLAB Workspace)
%
%   INPUTS
%       specify the number of columns containing ANYmaze data which you
%       would like to align
        FP.inputs.ANYmaze.UserSpecs.numbehavs = 16;
%       specify which columns to consider in the ANYmaze data file (time
%       excluded)
        FP.inputs.ANYmaze.UserSpecs.columns = [4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19];
%       specify what behaviour each column represents
        FP.inputs.ANYmaze.UserSpecs.behavs = ["Novel" "FamiliarA" "FamiliarB" "Middle"...
            "NoveltoFamiliarAstart" "NoveltoFamiliarAstop"...
            "NoveltoFamiliarBstart" "NoveltoFamiliarBstop"...
            "FamiliarAtoNovelstart" "FamiliarAtoNovelstop"...
            "FamiliarAtoFamiliarBstart" "FamiliarAtoFamiliarBstop"...
            "FamiliarBtoNovelstart" "FamiliarBtoNovelstop"...
            "FamiliarBtoFamiliarAstart" "FamiliarBtoFamiliarAstop"];
%       specify the alignment method for each column
%           linear: linear interpolation between known datapoints
%           nn: nearest neighbours. Works well with binary data
        FP.inputs.ANYmaze.UserSpecs.alignments = ["nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn" "nn"];
%
%   OUTPUTS
%       All stretched data to be saved to the output branch of the FP
%       structure element
        
%               load data
ampath=fullfile(FP.ID,'/ANYmaze.csv');                              %ANYmaze file
am=readmatrix(ampath);                                              %extract ANYmaze
am(isnan(am))=0;                                                    %replace empty columns with zero
fptime=FP.outputs.ttlon.time;                                       %load FP time
endtime=fptime(end);                                                %final timestamp of FP
amfr=endtime/(length(am));                                          %framerate of ANYmaze
amtime=(0:amfr:endtime)';                                           %generate time series
amtime=amtime(1:length(am));                                        %trim starting zero


%               working with behavioural data itself
for ii = 1:FP.inputs.ANYmaze.UserSpecs.numbehavs
    temp=am(:,FP.inputs.ANYmaze.UserSpecs.columns(ii));
    FP.inputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(ii)).raw=temp;
    tempstretch=nan(length(fptime),1);                              %storage array
    if strcmp(FP.inputs.ANYmaze.UserSpecs.alignments(ii),"nn")
        for iii=1:size(fptime,1)                                    %batch one row at a time
            temp_time=fptime(iii);   
            [~,ix]=min(abs(amtime-temp_time));                      %find index to align closest timestamps
            tempstretch(iii)=temp(ix,:);                            %apply index to behaviour
        end
        FP.outputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(ii)).stretched=tempstretch;
    elseif strcmp(FP.inputs.ANYmaze.UserSpecs.alignments(ii),"linear")
        tempstretch=interp1(amtime,temp,fptime);                    %linear interpolation
        FP.outputs.ANYmaze.(FP.inputs.ANYmaze.UserSpecs.behavs(ii)).stretched=tempstretch;
    else
        disp('Please ensure that alignment type is correctly specified');
    end
end

clearvars -except FP

disp('Process 2. ANYmaze Alignment Complete');

%%  Process 3. Forced Alternation Specific Analyses
%

%init
green=FP.outputs.ttlon.dFF470;                                      %load fp signal

%Photometry Analyses in Each Arm
%NOVEL
%   All Novel
am=FP.outputs.ANYmaze.Novel.stretched;                              %load all novel entries (binary)
am=logical(am);                                                     %make novel column able to index into other columns
green_am=green(am,:);                                               %index novel column into fp to get list of fp values while in novel
auc_green=(trapz(green_am))/(size(green_am,1)/(FP.inputs.fpfr));    %normalized by seconds in zone
FP.outputs.Analyses.NovelArm.AllTime.AUC.Green=auc_green;           %area under the curve while in the novel zone
%   First Novel
am_on=strfind(am',[0 1]);                                           %find index of first entry into the novel zone -- NOTE: WILL THROW ERROR IF RECORDING STARTS WITH ANIMAL IN THIS ZONE
am_off=strfind(am',[1 0]);                                          %find index of first novel zone exit
if am_off(1)<am_on(1)
    am_off=am_off(2:end);                                           %bring together the indices of these events
else
    %do nothing
end
green_am=green(am_on(1):am_off(1),:);                               %fp data during first novel zone entry
auc_green=(trapz(green_am))/(size(green_am,1)/(FP.inputs.fpfr));    %normalized by seconds in zone
FP.outputs.Analyses.NovelArm.FirstVisit.AUC.Green=auc_green;        %area under the curve during the first entry into the novel zone
%FAMILIAR
%   All Familiar
amA=FP.outputs.ANYmaze.FamiliarA.stretched;                         %indices of Familiar A
amB=FP.outputs.ANYmaze.FamiliarB.stretched;                         %indices of Familiar B
amM=FP.outputs.ANYmaze.Middle.stretched;                            %indices of middle
am=amA+amB+amM;                                                     %add columns
am(am>1)=1;                                                         %just to be safe, in the rare instance that the nearest neighbours framerate stretching results in an overlap, set the summed 2 to a 1
am=logical(am);                                                     %make this column able to index
green_am=green(am,:);                                               %identify f data during all familiar entries
auc_green=(trapz(green_am))/(size(green_am,1)/(FP.inputs.fpfr));    %normalized by seconds in zone
FP.outputs.Analyses.Familiar.AllTime.AUC.Green=auc_green;           %area under the curve in the familiar arms

%ZONE ENTRY ANALYSES
%init
amA=FP.outputs.ANYmaze.FamiliarA.stretched;                         %Familiar A
amB=FP.outputs.ANYmaze.FamiliarB.stretched;                         %Familiar B
amN=FP.outputs.ANYmaze.Novel.stretched;                             %Novel
amA_in=strfind(amA',[0 1]);                                         %indices of entering Familiar A
amA_out=strfind(amA',[1 0]);                                        %indices of leaving Familiar A
amB_in=strfind(amB',[0 1]);                                         %indices of entering Familiar B
amB_out=strfind(amB',[1 0]);                                        %indices of leaving Familiar B
amN_in=strfind(amN',[0 1]);                                         %indices of entering Novel
amN_out=strfind(amN',[1 0]);                                        %indices of leaving Novel

%First Familiar - Novel Entry
%set windows - seconds
FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore = 2;            %modify as needed - time to consider prior to zone entry
FP.inputs.Analyses.UserSpecs.FirstFamNov.timeafter = 2;             %modify as needed - time to consider after zone entry
%analysis
green_FirstFamNov=green((amN_in(1)-((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)):(amN_in(1)+((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)));
FP.outputs.Analyses.ZoneEntry.FirstFamNov.traces.green=green_FirstFamNov;       %trace
FP.outputs.Analyses.ZoneEntry.FirstFamNov.AUC.green=trapz(green_FirstFamNov);   %AUC

%Mean Familiar - Familiar Entry
%set windows - seconds
FP.inputs.Analyses.UserSpecs.FamFam.timebefore = 2;                 %modify as needed - time to consider prior to zone entry
FP.inputs.Analyses.UserSpecs.FamFam.timeafter = 2;                  %modify as needed - time to consider after zone entry
%analysis
vectorlength=(FP.inputs.Analyses.UserSpecs.FamFam.timebefore*FP.inputs.fpfr)+(FP.inputs.Analyses.UserSpecs.FamFam.timeafter*FP.inputs.fpfr)+1;
famBfamA_green=nan(vectorlength,size(amA_in,2));
for ii=1:size(amA_in)
    [valfam,locfam]=min(abs(amA_in(ii)-amB_out));
    [valnov,locnov]=min(abs(amA_in(ii)-amN_out));
    if valfam>valnov
        %entry is novel to familiar
        %do nothing
    else
        %entry is familiarB to familiarA
        famBfamA_green(:,ii)=green((amA_in(ii)-((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)):(amA_in(ii)+((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)));
    end
end
famAfamB_green=nan(vectorlength,size(amB_in,2));
for ii=1:size(amB_in)
    [valfam,locfam]=min(abs(amB_in(ii)-amA_out));
    [valnov,locnov]=min(abs(amB_in(ii)-amN_out));
    if valfam>valnov
        %entry is novel to familiar
        %do nothing
    else
        %entry is familiarA to familiarB
        famAfamB_green(:,ii)=green((amB_in(ii)-((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)):(amB_in(ii)+((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)));
    end
end
famfam_green=horzcat(famAfamB_green,famBfamA_green);
FP.outputs.Analyses.ZoneEntry.FamFam.traces.green=famfam_green;     %trace - familiar to familiar

%mean novel - familiar entries
%set windows - seconds
FP.inputs.Analyses.UserSpecs.NovFam.timebefore = 2;                 %modify as needed - time to consider prior to entry
FP.inputs.Analyses.UserSpecs.NovFam.timeafter = 2;                  %modify as needed - time to consider after entry
%analysis
vectorlength=(FP.inputs.Analyses.UserSpecs.NovFam.timebefore*FP.inputs.fpfr)+(FP.inputs.Analyses.UserSpecs.NovFam.timeafter*FP.inputs.fpfr)+1;
novfamA_green=nan(vectorlength,size(amA_in,2));
for ii=1:size(amA_in)
    [valfam,locfam]=min(abs(amA_in(ii)-amB_out));
    [valnov,locnov]=min(abs(amA_in(ii)-amN_out));
    if valfam>valnov
        %entry is novel to familiar A
        novfamA_green(:,ii)=green((amA_in(ii)-((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)):(amA_in(ii)+((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)));
    else
        %entry is familiarB to familiarA
        %do nothing
    end
end
novfamB_green=nan(vectorlength,size(amB_in,2));
for ii=1:size(amB_in)
    [valfam,locfam]=min(abs(amB_in(ii)-amA_out));
    [valnov,locnov]=min(abs(amB_in(ii)-amN_out));
    if valfam>valnov
        %entry is novel to familiar
        novfamB_green(:,ii)=green((amB_in(ii)-((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)):(amB_in(ii)+((FP.inputs.Analyses.UserSpecs.FirstFamNov.timebefore)*FP.inputs.fpfr)));
        
    else
        %entry is familiarA to familiarB
        %do nothing
    end
end
novfam_green=horzcat(novfamB_green,novfamA_green);
FP.outputs.Analyses.ZoneEntry.NovFam.traces.green=novfam_green;


clearvars -except FP

disp('Process 3. Forced Alternation Specific Analyses Complete');

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
 