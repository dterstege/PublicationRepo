% Generate some example data
x = randn(100,1);  % x position data (by time series)
x = abs(x*100)+1;
y = randn(100,1);  % y position data (by time series)
y = abs(y*100)+1;
z = randn(100,1);  % signal by position (by time series)

%% load real data
[filename,foldername]=uigetfile('*.mat', 'Select a mat file');
fullname=fullfile(foldername, filename);
load(fullname);

x=outputs.ANYmaze.xpos.stretched;
y=outputs.ANYmaze.ypos.stretched;
z=outputs.YMaze.GCaMP.dFF;

%clean up variable space
clearvars -except x y z

%% PLOT
%XY to whole numbers
x = round(x);
y = round(y);

%can't be at pixel 0, so we shift everything by 1
x = x+1;
y = y+1;

% 3D matrix for storage
c = nan(max(x),max(y),size(z,1));

% Sorting
for ii = 1:size(z,1)
    c(x(ii),y(ii),ii) =z (ii);
end

% Average signal across repeated spatial occurences
c = mean(c,3,'omitnan');
c(isnan(c)) = 0;

% Spatial downsampling
c = imresize(c, [170 140], 'bilinear'); %approx 1/3

% Rotate
c = rot90(c);

% Gaussian filter
c = imgaussfilt(c,3);

%cleaning
c(c == 0) = NaN;
c(47:139,1:93) = NaN;
c(50:139,80:100) = NaN;

% Visualize
figure
imagesc(c);
colormap(jet);
axis xy
colorbar
axis([-10 200 -60 180])

%clean up variable space
clearvars -except c

%% batch to group
IDs=["DT023" "DT024" "DT031" "DT073" "DT083"];
%["DT023" "DT024" "DT025" "DT031" "DT041" "DT043" "DT044" "DT051" "DT072" "DT073" "DT083" "DT093" "DT094" "DT095" "DT101" "DT112"]'
session="W1";

out=nan(140,170,size(IDs,1));

for ii = 1:size(IDs,1)
    filename=strcat(IDs(ii),"_",session,".mat");
    load(filename);
    
    x=outputs.ANYmaze.xpos.stretched;
    y=outputs.ANYmaze.ypos.stretched;
    z=outputs.ComplexMaze.GCaMP.dFF; %change depending on maze
    
    x = round(x);
    y = round(y);
    
    x = x+1;
    y = y+1;
    
    x = fillmissing(x,'nearest');
    y = fillmissing(y,'nearest');
    
    % 3D matrix for storage
    c = nan(max(x),max(y),size(z,1));

    % Sorting
    for iii = 1:size(z,1)
        c(x(iii),y(iii),iii) =z (iii);
    end

    % Average signal across repeated spatial occurences
    c = mean(c,3,'omitnan');
    c(isnan(c)) = 0;

    % Spatial downsampling
    c = imresize(c, [170 140], 'bilinear');

    % Rotate
    c = rot90(c);

    % Gaussian filter
    c = imgaussfilt(c,3);
    
    out(:,:,ii)=c;
end


%mask 
load('ymask.mat')
out = out + Ymask;

group=out;
out=mean(out,3,'omitnan');

figure
imagesc(out,'AlphaData', 1-isnan(out));
colormap(jet);
axis xy
colorbar
axis([-10 200 -60 180])
ax=gca;
ax.CLim = [-0.002 0.005];

%clean up variable space
clearvars -except out group

%% Compare groups
load('week1_group.mat')
w1=group;

load('week2_group.mat')
w2=group;

load('week3_group.mat')
w3=group;

[h2,p2,~,stats2] = ttest2(w1,w2, 'dim',3,'tail', 'both');
[h3,p3,~,stats3] = ttest2(w1,w3, 'dim',3,'tail', 'both');

filtered_p2 = p2;
filtered_p2(filtered_p2>=0.05)=NaN;
figure
imagesc(filtered_p2,'AlphaData', 1-isnan(filtered_p2));
colormap(jet);
axis xy
colorbar
axis([-10 200 -60 180])


filtered_p3 = p3;
filtered_p3(filtered_p3>=0.05)=NaN;
figure
imagesc(filtered_p3,'AlphaData', 1-isnan(filtered_p3));
colormap(jet);
axis xy
colorbar
axis([-10 200 -60 180])

