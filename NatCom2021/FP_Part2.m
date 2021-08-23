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
%   Run block-wise
%   Jump to 2nd last section (Section X) to save data structure at any time
%   Load previously saved data structures using last section (Section Y)
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
%           Baseline correction:
%               Trial is split into test and baseline segments
%               Mean dFF value during baseline is calculated (X frames)
%               Minimum dFF value from the test segment is calculated
%               The difference between these two values is then added to
%               every dFF value during the test segment
%               This effectively raises the minimum dFF value from the test
%               segment to the mean baseline segment value
%       peaks:
%           Peaks defined as any point at which the dFF trace makes a
%           downward deflection after rising
%           Peak detection restricted to peaks exceeding a height of two
%           standard deviations above the median dFF value
%           Mean peak height is relative to 0
%
%                           *ASSUMPTIONS*
%                      *IMPORTANT - PLEASE READ*
%
%   Like most analyses, this code was developed with a number of
%   assumptions about the organization of data files. The script assumes
%   that the data is stored in an intermediate file 'test_summary.xlsx'.
%   See attached 'test_summary.xlsx' file for example of data organization.
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
%adjust groups to suit the data set

%load data
%MUST NAVIGATE TO FOLDER WITH 'test_summary.xlsx' FIRST
FP.is.full=readmatrix('test_summary.xlsx','Sheet','rawiso');
FP.is.full(:,17)=[];      %remove spacing column
FP.r.full=readmatrix('test_summary.xlsx','Sheet','rawsig');
FP.r.full(:,17)=[];      %remove spacing column
FP.dFF.full=readmatrix('test_summary.xlsx','Sheet','dFF');
FP.dFF.full(:,17)=[];     %remove spacing column
FP.HC.full=readmatrix('test_summary.xlsx','Sheet','HCdFF');
FP.HC.full(:,17)=[];       %remove spacing column
FP.F.full=readmatrix('test_summary.xlsx','Sheet','F');
FP.F.full(:,17)=[];       %remove spacing column

%split into groups - Adjust based on group-specific columns
FP.is.sed.full=FP.is.full(:,1:5); %group 1 range in xlsx sheet (5 mice)
FP.is.run.full=FP.is.full(:,7:11); %group 2 range in xlsx sheet (5 mice)
FP.r.sed.full=FP.r.full(:,1:5);
FP.r.run.full=FP.r.full(:,7:11);
FP.dFF.sed.full=FP.dFF.full(:,1:5);
FP.dFF.run.full=FP.dFF.full(:,7:11);
FP.F.sed.full=FP.F.full(:,1:5);
FP.F.run.full=FP.F.full(:,7:11);
FP.HC.sed.full=FP.HC.full(:,1:5);
FP.HC.run.full=FP.HC.full(:,7:11);

disp('Process 1. Initialization Complete');

%%  2. Minimum dF/F Homecage to Zero
%adjusts trace so that the minimum dFF from the homecage recording are set
%to zero. All other values in the trace are scaled accordingly
%calculates AUC and mean peak height
%total number of peaks and mean peak height
%calculated on a mouse-by-mouse basis

%initialize
%generate output arrays
FP.outputs.AUC.full.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.min1.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.min2.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.min3.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.min4.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.min5.sed=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.peaks.sed.number=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.peaks.sed.height=zeros(size(FP.dFF.sed.full,2),1);
FP.outputs.AUC.full.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.AUC.min1.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.AUC.min2.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.AUC.min3.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.AUC.min4.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.AUC.min5.run=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.peaks.run.number=zeros(size(FP.dFF.run.full,2),1);
FP.outputs.peaks.run.height=zeros(size(FP.dFF.run.full,2),1);

%metrics from control group while freezing
FP.outputs.freeze.sed.AUCnorm=zeros(size(FP.rs.sed.full,2),1);
FP.outputs.freeze.sed.peakrate=zeros(size(FP.rs.sed.full,2),1);
FP.outputs.freeze.sed.peakheight=zeros(size(FP.rs.sed.full,2),1);

%metrics from control group while moving
FP.outputs.move.sed.AUCnorm=zeros(size(FP.rs.sed.full,2),1);
FP.outputs.move.sed.peakrate=zeros(size(FP.rs.sed.full,2),1);
FP.outputs.move.sed.peakheight=zeros(size(FP.rs.sed.full,2),1);

%metrics from running group while freezing
FP.outputs.freeze.run.AUCnorm=zeros(size(FP.rs.run.full,2),1);
FP.outputs.freeze.run.peakrate=zeros(size(FP.rs.run.full,2),1);
FP.outputs.freeze.run.peakheight=zeros(size(FP.rs.run.full,2),1);

%metrics from running group while moving
FP.outputs.move.run.AUCnorm=zeros(size(FP.rs.run.full,2),1);
FP.outputs.move.run.peakrate=zeros(size(FP.rs.run.full,2),1);
FP.outputs.move.run.peakheight=zeros(size(FP.rs.run.full,2),1);

%CTRL
for ii=1:size(FP.F.sed.full,2)
    idxf=(FP.F.sed.full(:,ii)==1); %freezing
    idxf=idxf(1:6000);
    idxm=(FP.F.sed.full(:,ii)==0); %movement
    idxm=idxm(1:6000);
    hc=FP.HC.sed.full(1:3600,ii); %identify homecage period to be used as baseline (1:3600 frames)
    hc=mean(hc);
    full=FP.dFF.sed.full(1:6000,ii); %identify testing period (1:6000 frames)
    mindFF=min(full);
    adjust=abs(abs(mindFF)-abs(hc));
    full=full+adjust;
    freeze=full(idxf);    
    timef=size(freeze,1)*0.05; %framerate of the photometry recording (ex, 0.05s/frame)
    move=full(idxm);    
    timem=size(move,1)*0.05;m%framerate of the photometry recording (ex, 0.05s/frame)
    FP.outputs.freeze.sed.AUCnorm(ii)=trapz(freeze)/timef;
    FP.outputs.move.sed.AUCnorm(ii)=trapz(move)/timem;
    maxdFF=median(full)+(2*std(full));
    FP.outputs.AUC.full.sed(ii)=trapz(full(1:6000));
    FP.outputs.AUC.min1.sed(ii)=trapz(full(1:1200));
    FP.outputs.AUC.min2.sed(ii)=trapz(full(1201:2400));
    FP.outputs.AUC.min3.sed(ii)=trapz(full(2401:3600));
    FP.outputs.AUC.min4.sed(ii)=trapz(full(3601:4800));
    FP.outputs.AUC.min5.sed(ii)=trapz(full(4801:6000));
    [pks,~,~,~]=findpeaks(full(1:6000),'MinPeakHeight',maxdFF);
    FP.outputs.peaks.sed.number(ii)=size(pks,1);
    FP.outputs.peaks.sed.height(ii)=mean(pks);
    [pks,~,~,~]=findpeaks(freeze,'MinPeakHeight',maxdFF);
    FP.outputs.freese.sed.peakrate(ii)=size(pks,1)/timef;
    FP.outputs.freeze.sed.peakheight(ii)=mean(pks);
    [pks,~,~,~]=findpeaks(move,'MinPeakHeight',maxdFF);
    FP.outputs.move.sed.peakrate(ii)=size(pks,1)/timem;
    FP.outputs.move.sed.peakheight(ii)=mean(pks);
end

%RUN
for ii=1:size(FP.F.run.full,2)
    idxf=(FP.F.run.full(:,ii)==1); %freezing
    idxf=idxf(1:6000);
    idxm=(FP.F.sed.full(:,ii)==0); %movement
    idxm=idxm(1:6000);
    hc=FP.HC.run.full(1:3600,ii); %identify homecage period to be used as baseline (1:3600 frames)
    hc=mean(hc);
    full=FP.dFF.run.full(1:6000,ii); %identify testing period (1:6000 frames)
    mindFF=min(full);
    adjust=abs(abs(mindFF)-abs(hc));
    full=full+adjust;
    freeze=full(idxf);    
    timef=size(freeze,1)*0.05; %framerate of the photometry recording (ex, 0.05s/frame)
    move=full(idxm);    
    timem=size(move,1)*0.05; %framerate of the photometry recording (ex, 0.05s/frame)
    FP.outputs.freeze.run.AUCnorm(ii)=trapz(freeze)/timef;
    FP.outputs.move.run.AUCnorm(ii)=trapz(move)/timem;
    maxdFF=median(full)+(2*std(full));
    FP.outputs.AUC.full.run(ii)=trapz(full(1:6000,:));
    FP.outputs.AUC.min1.run(ii)=trapz(full(1:1200));
    FP.outputs.AUC.min2.run(ii)=trapz(full(1201:2400));
    FP.outputs.AUC.min3.run(ii)=trapz(full(2401:3600));
    FP.outputs.AUC.min4.run(ii)=trapz(full(3601:4800));
    FP.outputs.AUC.min5.run(ii)=trapz(full(4801:6000));
    [pks,~,~,~]=findpeaks(full(1:6000),'MinPeakHeight',maxdFF);
    FP.outputs.peaks.run.number(ii)=size(pks,1);
    FP.outputs.peaks.run.height(ii)=mean(pks);
    [pks,~,~,~]=findpeaks(freeze,'MinPeakHeight',maxdFF);
    FP.outputs.freese.run.peakrate(ii)=size(pks,1)/timef;
    FP.outputs.freeze.run.peakheight(ii)=mean(pks);
    [pks,~,~,~]=findpeaks(move,'MinPeakHeight',maxdFF);
    FP.outputs.move.run.peakrate(ii)=size(pks,1)/timem;
    FP.outputs.move.run.peakheight(ii)=mean(pks);
end
clear ii hc full mindFF adjust maxdFF pks
clear freeze idxf idxm move timef timem

disp('Process 2. dF/F Minimum to Mean Homecage dF/F');