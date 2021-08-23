%%                              FP_Part2.m
%
%   Fiber Photometry Analysis - Part 2
%   
%   It is recommended that users read through the documentation in full
%   prior to using this analysis. This applies to each section in the code
%
%   All outputs are accessed through the structure element "FP" in the
%   MATLAB workspace
%
%                       *GENERAL INFORMATION*
%
%   Behavioural Paradigm: Contextual Fear Conditioning.  Neurophotometrics
%       photometry recording paired with ANYmaze behavioural tracking
%
%                        *ANALYSIS INFORMATION*
%
%   Inputs:
%       Neurophotometrics Photometry Output (.csv)
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
%       "behaviourindependent"
%
%   Operational Definitions:
%       dFF:
%           Isosbestic data is fit to a biexponential decay
%           This biexponential decay is linearly scaled to the raw 470 data
%           (raw = sig, fitted = sigfit)
%           dFF = (sig-sigfit)./sigfit
%
%
%                           *ASSUMPTIONS*
%                      *IMPORTANT - PLEASE READ*
%
%   Like most analyses, this code was developed with a number of
%   assumptions about the organization of data files. The script assumes
%   that the data is stored in files of very specific names. Under a parent
%   folder, a file named '415_.csv' should contain the isosbestic data from
%   the neurophotometrics system. A file named '470_.csv' should contain
%   the photometry signal. Runs each mouse individually. Assumes specific
%   TTL flags in the neurophotometrics outputs - adjust these as necessary.
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.1


%run block-wise using "Run Section" - jump to last section to save at any time

%Created 05/19/2021 Dylan Terstege
%Epp Lab, University of Calgary
%Contact: dylan.terstege@ucalgary.ca



%%  1. Initialization
%mandatory step

%load data
disp("Select parent folder");
FP.parent=uigetdir();

isopath=strcat(FP.parent,'/415_.csv'); %isosbestic file
sigpath=strcat(FP.parent,'/470_.csv'); %signal file

iso=readmatrix(isopath); %extract isosbestic
sig=readmatrix(sigpath); %extract signal

temp=iso(:,3)<20;           %pre TTL flag is 17 (415) or 18 (470)
hciso=iso(temp,4);          %homecage
iso(temp,:)=[];             %post TTL flag jumps to 273 (415) or 274 (470) - indicates when mouse is in CFC chamber
fptime=iso(:,2);            %time column
fptime=fptime-fptime(1);    %adjusts time so that first block is 0
iso=iso(:,4);               %only the  isosbestic column
temp=sig(:,3)<20;           %pre TTL flag is 17 (415) or 18 (470)
hcsig=sig(temp,4);          %homecage
sig(temp,:)=[];             %post TTL flag jumps to 273 (415) or 274 (470)
sig=sig(:,4);               %only the  signal column

if length(iso)>length(sig)  %make all vectors the same length
    iso=iso(1:length(sig));
    fptime=fptime(1:length(sig));
else
    sig=sig(1:length(iso));
end

if length(hciso)>length(hcsig)  %make all vectors the same length
    hciso=hciso(1:length(hcsig));
else
    hcsig=hcsig(1:length(hciso));
end

%trial signal
temp_x=1:length(iso);
temp_x=temp_x';
isofit=fit(temp_x,iso,'exp2');

sigfit=robustfit(isofit(temp_x),sig);
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

dFF=(sig-sigfit)./sigfit;   %delta F/F

%homecage signal
temp_x=1:length(hciso);
temp_x=temp_x';
isofit=fit(temp_x,hciso,'exp2');

sigfit=robustfit(isofit(temp_x),hcsig);
sigfit=isofit(temp_x)*sigfit(2)+sigfit(1);

hcdFF=(hcsig-sigfit)./sigfit;   %delta F/F

FP.time=fptime; %time
FP.signal=sig; %signal
FP.isosbestic=iso; %isosbestic
FP.dFF=dFF; %delta F over F
FP.HC=hcdFF; %delta F over F during home cage recording

clear dFF fptime iso isopath isofit parent sig sigfit sigpath temp temp_x hcsig hciso hcdFF

disp('Process 1. Initialization Complete');

%%  2.  ANYmaze integration

%align behavioural timestamps with FP
%will duplicate behavioural values to gapfill through Nearest Neighbours
%       no interpolation
%this keeps binary for freezing data

disp(FP.parent);
disp("Select ANYmaze output");
[behavfile,behavpath]=uigetfile('.csv');
behavfile=strcat(behavpath,behavfile);
behavtime=readmatrix(behavfile,'OutputType','duration');
behavtime=behavtime(:,1);
behavtime=seconds(behavtime);
behav=readmatrix(behavfile);
behav(:,1)=behavtime;
behavstretch=zeros(length(FP.dFF),size(behav,2));

for ii=1:length(FP.dFF)
    temp_time=FP.time(ii);   %align closest timestamps - Nearest Neighbour
    [~,ix]=min(abs(behavtime-temp_time));
    behavstretch(ii,:)=behav(ix,:);
end

FP.freeze=behavstretch(:,5); %only considering binary freezing column - adjust as necessary

clear ii ix behav behavfile behavpath behavstretch behavtime temp_time

disp('Process 2. ANYmaze Integration Complete');