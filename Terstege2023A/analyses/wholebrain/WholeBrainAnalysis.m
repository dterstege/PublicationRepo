%%                          WholeBrainAnalysis.m
%   Whole Brain Analysis
%   
%   It is recommended that users read through the documentation in full
%   prior to using this analysis. This applies to each section in the code.
%   Sections are to be ran sequentially.
%
%   All outputs can be accessed through the structure element "WB" in the
%   MATLAB workspace
%
%                       *GENERAL INFORMATION*
%
%   WholeBrainAnalysis.m was developed for the compilation and analysis of
%   datasets collected using the Whole Brain software suite developed by
%   Daniel Furth (http://www.wholebrainsoftware.org) and incorporates a 
%   modification using scripts and plugins 
%   (https://github.com/dterstege/CavalieriPointMask) developed by 
%   Dylan Terstege.
%
%   Whole Brain is a neuroanatomical information system.  Neuroanatomical
%   data obtained from microscope images is encoded and stored in
%   stereotactic coordinate form within the Allen Mouse Brain Atlas.
%   The software suite can be accessed at http://www.wholebrainsoftware.org
%
%
%   Additional Whole Brain scripts can be obtained from Dylan Terstege upon
%   request.
%
%   Can be saved and re-loaded at any time using Processes X and Y, found
%   at the end of this script.
%
%   NOTE: Avoid using spaces in any file, folder, or variable names. MATLAB
%   doesn't handle this well and in many cases the program will throw an
%   error in response.
%
%   INPUTS:
%       "cells.csv" files: naming structure is not of particular
%       importance, but these files should contain the number of segmented
%       LABELS per region. This is a raw Whole Brain output. In the example
%       of c-Fos mapping, this is the number of c-fos+ cells per region.
%       "grids.csv" files: naming structure is not of particular
%       importance, but these files hsoudl contain the number of segmented
%       POINTS per region. This is a raw Whole Brain output. These points
%       are counted from the Cavalier Point Masks
%
%   Requires "WB_Atlas.mat" to be in the MATLAB path
%
%                       *DEVELOPER INFORMATION*
%
%   Version 1.0.4
%
%   Check for outliers before running correlations.  Only correlates
%   non-outliers
%
%   Option to threshold by p value percentile rank
%
%
%   Created 11/06/2021 Dylan Terstege (https://github.com/dterstege)
%   Epp Lab, University of Calgary (https://epplab.com)
%   Contact: dylan.terstege@ucalgary.ca

welcomemessage = ['         *** Running Whole Brain Analysis *** ', newline, '                    Version 1.0.1',newline, newline, 'This program was written and developed by Dylan Terstege'];
disp(welcomemessage);
clear welcomemessage

%%  01. Initialization
%   During this process, data will be loaded into the MATLAB structure
%   element 'WB' for later use.  This step in the analysis is reasonable
%   hands-on and requires the user to point the script to several
%   directories and input animal IDs.  The results of this data loading
%   can be found under the 'inputs' element within 'WB'
%
%   Requires "WB_Atlas.mat" to be in the MATLAB path
%
%   There are several variables which the user should adjust prior to
%   starting the analysis. These are listed below
%
%   Group Information:
%       Input number of groups
%           default: groupnum = 2
%       Input group identifiers
%           default: groups = ["CTRL" "MWT"]
%       Input number of animals per group
%           default: pergroup = [10 10]
%       Input individual IDs
%           defaults:
%           IDs.CTRL = ["DT42br" "DT44br" "DT45w" "DT63" "DT64" "DT65" "DT71" "DT76" "DT92" "DT93"];
%           IDs.MWT = ["DT11" "DT12" "DT12bl" "DT12br" "DT13" "DT14" "DT15" "DT16w" "DT17br" "DT23"];
%
%   Image Specifications:
%       Input area factor
%           default: areafactor = 0.56286606



%clear structure to start fresh
clear WB

%input number of groups
WB.info.groupnum = 2;
%input group identifiers
WB.info.groups = ["CONTROL" "GiDREADD"]; %manually adjust, no spaces
%input number of animals per group
WB.info.pergroup = [10 10]; %manually adjust
%input individual IDs
%make sure group names appear in same order as inputted in the group
%identifiers variable
WB.info.IDs.CONTROL = ["DT33" "DT34" "DT35" "DT51" "DT52" "DT53" "DT54" "DT55" "DT75" "DT95"];
WB.info.IDs.GiDREADD = ["DT11" "DT12" "DT13" "DT14" "DT71" "DT72" "DT91" "DT92" "DT93" "DT94"];

%area factor
WB.info.areafactor = 0.56286606; % ((dist/pix)*(pix/point))^2
%load full whole brain atlas
WB.info.fullatlas = load('WB_Atlas.mat');

%load whole brain output files
for ii=1:WB.info.groupnum(1)
    for iii=1:WB.info.pergroup(ii)
        %import raw data
        %import cells
        prompt=strcat("Select parent folder for mouse ",(WB.info.IDs.(WB.info.groups(ii))(iii))," in ",WB.info.groups(ii)," group:");
        disp(prompt);
        parentpath=uigetdir; %identify parent directory
        cellpath=fullfile(parentpath,'cell');
        celldir=dir(cellpath); %set directory
        cellcsv={celldir(:).name}; %scan directory
        cellcsv(ismember(cellcsv,{'.','..'}))=[]; %filter
        cellcsv=cellcsv(~startsWith(cellcsv,'._')); %removes dot underscore files
        csvindex=contains(cellcsv,'.csv'); %identify csv
        cellcsv=cellcsv(csvindex); %list csv files
        cellcount=numel(cellcsv); %count csv
        WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).cells=NaN(numel(WB.info.fullatlas.WB_Atlas),cellcount);
        for iiii=1:cellcount
            data=readcell(char(fullfile(cellpath,cellcsv(iiii)))); %load cell file
            regions=data(2:end,2); %identify regions in file
            counts=data(2:end,3);
            [~,k]=intersect(WB.info.fullatlas.WB_Atlas(2:end),regions); %find index of regions in file with full atlas
            WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).cells(k,iiii)=cell2mat(counts);
        end
        
        %import grids
        gridpath=fullfile(parentpath,'grid');
        griddir=dir(gridpath); %set directory
        gridcsv={griddir(:).name}; %scan directory
        gridcsv(ismember(gridcsv,{'.','..'}))=[]; %filter
        gridcsv=gridcsv(~startsWith(gridcsv,'._')); %removes dot underscore files
        csvindex=contains(gridcsv,'.csv'); %identify csv
        gridcsv=gridcsv(csvindex); %list csv files
        gridcount=numel(gridcsv); %count csv
        WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).areas=NaN(numel(WB.info.fullatlas.WB_Atlas),gridcount);
        for iiii=1:gridcount
            data=readcell(char(fullfile(gridpath,gridcsv(iiii)))); %load grid file
            regions=data(2:end,2); %identify regions in file
            area=cell2mat(data(2:end,3));
            area=area*(WB.info.areafactor);
            [~,k]=intersect(WB.info.fullatlas.WB_Atlas,regions); %find index of regions in file with full atlas
            WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).areas(k,iiii)=area;
        end
    end
end
%process successful

%clean up variable space
clearvars -except WB

%process complete
disp('Process 1. Initialization Complete');

%%  02.  Regional Label Density
%   During this process, the density of the segmented label will be
%   calculated using a user-defined atlas organization.  This process is
%   completely hands-off once the atlas has been created and added to the
%   MATLAB path.  Its name should be inputed to the script in place of the
%   default atlas.
%
%   Requires "User_Atlas_98reg.mat", or whatever the user names it, to be in the
%   MATLAB path
%
%   User must have the WB structure in their workspace prior to starting
%   this process

%load user atlas
WB.info.useratlas = load('User_Atlas_90reg.mat'); %adjust as needed
regnum=((numel(WB.info.useratlas.User_Atlas_90reg))/3)-1; %number of regions %adjust as needed
regiondims=cell2mat(WB.info.useratlas.User_Atlas_90reg(2:end,2:3)); %identifies where regions start and end %adjust as needed
regnames=WB.info.useratlas.User_Atlas_90reg(2:end,1); %adjust as needed
WB.info.useratlas.regnames=regnames;

for ii=1:WB.info.groupnum(1)
    %create output arrays
    WB.outputs.(WB.info.groups(ii)).totalcells=NaN(regnum,WB.info.pergroup(ii)); %NaN array
    WB.outputs.(WB.info.groups(ii)).totalarea=NaN(regnum,WB.info.pergroup(ii)); %NaN array
    WB.outputs.(WB.info.groups(ii)).density=NaN(regnum,WB.info.pergroup(ii)); %NaN array
    WB.outputs.(WB.info.groups(ii)).overalldens=NaN(WB.info.pergroup(ii),1);
    WB.outputs.(WB.info.groups(ii)).overallarea=NaN(WB.info.pergroup(ii),1);
    for iii=1:WB.info.pergroup(ii)
        cells=nansum(WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).cells,2); %total number of segmented labels in each region
        areas=nansum(WB.inputs.(WB.info.groups(ii)).(char(WB.info.IDs.(WB.info.groups(ii))(iii))).areas,2); %total area of each region
        for iiii=1:regnum %sequentially run through regions in custom atlas organization
            temp=regiondims(iiii,:); %isolate dimensions of region
            tempcell=sum(cells(temp(1):temp(2),1)); %isolate and sum row(s) within region bounds
            temparea=sum(areas(temp(1):temp(2),1)); %isolate and sum row(s) within region bounds
            tempdensity=tempcell/temparea; %calculate density
            WB.outputs.(WB.info.groups(ii)).totalcells(iiii,iii)=tempcell; %save cells
            WB.outputs.(WB.info.groups(ii)).totalarea(iiii,iii)=temparea; %save area
            WB.outputs.(WB.info.groups(ii)).density(iiii,iii)=tempdensity; %save density
        end
        %corrections
        idx=((WB.outputs.(WB.info.groups(ii)).totalarea)==0);
        WB.outputs.(WB.info.groups(ii)).totalcells(idx)=0;
        WB.outputs.(WB.info.groups(ii)).density(idx)=nan;
        
        %brain-wide totals
        WB.outputs.(WB.info.groups(ii)).overallarea(iii)=sum(WB.outputs.(WB.info.groups(ii)).totalarea(isfinite(WB.outputs.(WB.info.groups(ii)).totalarea)),1);
        temp=sum(WB.outputs.(WB.info.groups(ii)).totalcells(isfinite(WB.outputs.(WB.info.groups(ii)).totalcells)),1);
        WB.outputs.(WB.info.groups(ii)).overalldens(iii)=temp/(WB.outputs.(WB.info.groups(ii)).overallarea(iii));
    end
end

%clean up variable space
clearvars -except WB

%process complete
disp('Process 2. Regional Label Density Complete');

%% 03.   Basic Network Analyses
%   During this process, pairwise correlation matrices will be generated
%   from the segmented label densities from 
%
%   Required packages in MATLAB path:
%       fdr_bh.m (REF)
%       degrees_und.m (REF)
%
%   Outputs:
%       - Correlation Matrices
%               Pearson's R
%           p-value matrix
%       - Binarized Adjacency Matrix
%       - Survival curve of functional connections with increased Pearson's
%       threshold
%
%   User must have the WB structure in their workspace prior to starting
%   this process
%
%   NOTE: 'NaN' values in density matrix will have major implications on the
%   correlation matrices.  Be sure to consider these values very carefully
%   during interpretation of results. NaN occurs when a region was not
%   registered for an animal

%user input variables
%thresholding
%   p value percentile rank or raw, (1=percentile; 0=raw)
    WB.info.networks.threshold.type = 0;
%   percentile to consider
    WB.info.networks.threshold.percentile = 75;
%   significance threshold for pariwise correlations
    WB.info.networks.threshold.pthresh = 0.05;
%   false discovery rate threshold for pariwise correlations
    WB.info.networks.threshold.fdrthresh = 0.05;
%   cutoff for binarizing matrix
    WB.info.networks.threshold.binary = 0.8;
%correlation strength survival
%   minimum value to look at
    WB.info.networks.corrstrength.startR=0.00;
%   maximum value to look at    
    WB.info.networks.corrstrength.lastR=1.00;
%   resolution of analysis
    WB.info.networks.corrstrength.stepSize=0.01;



%analysis
for ii=1:WB.info.groupnum(1)
    %load files
    dens=WB.outputs.(WB.info.groups(ii)).density;
    dens(isnan(dens))=0;
    
    %check for outliers in density data
    outliers=isoutlier(dens','grubbs');
    WB.outputs.(WB.info.groups(ii)).outliers=outliers;
    %dens(outliers')=NaN;
    
    %correlate data
    [cormat,pVal]=corr(dens','Type','Pearson','rows','pairwise'); %pearson pairwise correlation coefficient
    cormat(isnan(cormat))=0;
    
    %save un-thresholded data
    WB.outputs.(WB.info.groups(ii)).networks.correlationmatrix.unthresholded=cormat;
    WB.outputs.(WB.info.groups(ii)).networks.correlationmatrix.pvalues=pVal;
    
    %optional: save correlation matrix heatmap
    %correlation matrix image
    figure('visible','off');
    clims=[-1 1];
    imagesc(cormat,clims);
    colorbar;
    colormap(redbluecmap);
    set(gcf,'Position',get(0,'Screensize'));
    axis('square');
    set(gca,'YTick',[]);
    set(gca,'XTick',[]);
    figname=strcat(WB.info.groups(ii),'_correlation.png');
    saveas(gcf,figname); %makes graph image
    close gcf
    
    %mean R
    if ~WB.info.networks.threshold.type
        cormat(pVal>=(WB.info.networks.threshold.pthresh))=NaN;
        fdr=fdr_bh(pVal,(WB.info.networks.threshold.fdrthresh));
        cormat(fdr==0)=NaN;
        WB.outputs.(WB.info.groups(ii)).networks.meanR=(mean(cormat,'omitnan'))';
    elseif WB.info.networks.threshold.type
        temp_pVal=triu(pVal);
        temp_pVal=temp_pVal(:);
        temp_pVal(temp_pVal==0)=100;
        percentile_pVal=prctile(temp_pVal,(100-WB.info.networks.threshold.percentile),"all");
        cormat(pVal>=(percentile_pVal))=NaN;
        WB.outputs.(WB.info.groups(ii)).networks.meanR=(mean(cormat,'omitnan'))';
    end
    
    %threshold
    cormat(isnan(cormat))=0;
    WB.outputs.(WB.info.groups(ii)).networks.correlationmatrix.thresholded.weighted=cormat;
    cormat(cormat~=0)=1;
    WB.outputs.(WB.info.groups(ii)).networks.correlationmatrix.thresholded.binary=cormat;
    
    %optional: save adjacency matrix
    %adjacency matrix image
    %figure('visible','off');
    %clims=[-1 1];
    %imagesc(cormat,clims);
    %colorbar;
    %colormap(redbluecmap);
    %set(gcf,'Position',get(0,'Screensize'));
    %axis('square');
    %set(gca,'YTick',[]);
    %set(gca,'XTick',[]);
    %figname=strcat(WB.info.groups(ii),'_adjacency.png');
    %saveas(gcf,figname); %makes graph image
    %close gcf
    
    %correlation strength survival
    numsteps=((WB.info.networks.corrstrength.lastR)-(WB.info.networks.corrstrength.startR))/(WB.info.networks.corrstrength.stepSize);
    corsurvivalcurve=zeros(numsteps,1);
    tempR=WB.info.networks.corrstrength.startR;
    for iii=1:numsteps
        temp=cormat;
        temp(temp<tempR & temp>(-(tempR)))=0; %threshold
        temp(temp>tempR | temp<(-(tempR)))=1; %threshold
        deg=degrees_und(temp);
        totalconnections=sum(deg)/2; %number of connections
        corsurvivalcurve(iii)=totalconnections;
        tempR=tempR+(WB.info.networks.corrstrength.stepSize); %adjust threshold to next step
    end
    WB.outputs.(WB.info.groups(ii)).networks.corrstrengthsurvival=corsurvivalcurve;
    
    %binarize to undirected adjacency matrix
    cormat(cormat<(WB.info.networks.threshold.binary) & cormat>-(WB.info.networks.threshold.binary))=0;
    cormat(cormat>(WB.info.networks.threshold.binary) | cormat<-(WB.info.networks.threshold.binary))=1;
    
    %for only negatives
    %cormat(cormat>-(WB.info.networks.threshold.binary))=0;
    %cormat(cormat<-(WB.info.networks.threshold.binary))=1;
    
    
    %optional: save adjacency matrix circle plot
    %G=graph(cormat,'upper','omitselfloops');
    %figure('visible','off');h=plot(G,'layout','circle'); h.MarkerSize=2;
    %regnames=WB.info.useratlas.(char(fieldnames(WB.info.useratlas)))(2:end,1);
    %labelnode(h,(1:size(cormat,1)),regnames);
    %set(gcf,'Position',get(0,'Screensize'));
    %axis('square');
    %set(gca,'YTick',[]);
    %set(gca,'XTick',[]);
    %figname=strcat(WB.info.groups(ii),'_circleplot.png');
    %saveas(gcf,figname); %makes graph image
    %close gcf
    
    %save graph
    G=graph(cormat,'upper','omitselfloops');
    G.Nodes.Name=WB.info.useratlas.regnames;
    WB.outputs.(WB.info.groups(ii)).networks.G=G;
    
    %general network parameters
    %local efficiency
    le=efficiency_bin(cormat,1);
    %clustering coefficient
    cc=clustering_coef_bu(cormat);
    %betweenness centrality
    bc=betweenness_bin(cormat);

    %save outputs
    WB.outputs.(WB.info.groups(ii)).networks.clusteringcoeff=cc;
    WB.outputs.(WB.info.groups(ii)).networks.localefficiency=le;
    WB.outputs.(WB.info.groups(ii)).networks.betweennesscentrality=bc;
    
end

%clean up variable space
clearvars -except WB

%process complete
disp('Process 3. Basic Network Analyses Complete');

%%  X. Save Structure to .mat File
%keep file for easy data access later

save('WB','-struct','WB');

disp('Process X. Save Structure to a .mat File Complete');

%%  Y. Load .mat File
%   navigate to folder containing 'WB.mat'

load('WB_90reg_GiDREADD.mat')
 w=whos;
 for ii=1:length(w)
     WB.(w(ii).name)=eval(w(ii).name); 
     clear w(ii).name
 end
 
%clean up variable space
clearvars -except WB
 
disp('Process Y. Load .mat File Complete');