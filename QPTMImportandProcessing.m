TrackMateQPTable = uigetfile;
%Built off the output of Trackmate/QP table combiner. Assumes the same
%variables will be present in the same order
% Add table to be analyzed. Attempts to remove non numeric variables
paras = removevars(TrackMateQPTable,{'MANUAL_SPOT_COLOR','Bin','MANUAL_EDGE_COLOR','LABEL_edgeTable'}); %[output:04613fa7]
% paras = removevars(TrackmateQPTable,{'MANUAL_SPOT_COLOR','Bin','MANUAL_EDGE_COLOR'});
% paras = removevars(TrackmateQPTable,{'MANUAL_SPOT_COLOR','MANUAL_EDGE_COLOR','LABEL_edgeTable'});
l1 = paras.Properties.VariableNames; % Saves what variables are which for later
paras = table2array(paras);  %Changes format and sorts based on track number and frame (Sorting probably not necessary)
% paras = str2double(paras);
sortedpcp = sortrows(paras,[42,13]);
varnum=size(sortedpcp,2); %Setting up dimensions of storage space
tracks = sortedpcp(:,3);
tracknum = length(unique(tracks));
frames = max(sortedpcp(:,13));
sortedtrack=zeros(tracknum,varnum,frames); %Sorts into :Rows= Track ID, Columns = Frame, Depth = Variable 
for i =1:length(sortedpcp)
    sortedtrack(1+sortedpcp(i,3),:,(sortedpcp(i,13))) = sortedpcp(i,:);
end
sortedtrack = permute(sortedtrack,[1,3,2]);
sortedtrack = sortedtrack(any(sortedtrack ~= 0,[2 3]),:,:); %Attempts to remove any tracks that are exclusively zero
sortedtrack = changem(sortedtrack,NaN);
%%
%Uses the forwardmost 50 cells to estabish a front tissue boundary, saves
%distance to tissue front as as a new variable at depth zval
zval = size(sortedtrack,3)+1;
sortedtrack(:,:,zval)=0;
% sortedtrack(:,:,9) = sortedtrack(:,:,9)*-1;
for i=1:size(sortedtrack,2)
    tempsort =sortrows(sortedtrack(:,i,9));
    meanext(i) = mean(tempsort(1:50));
    sortedtrack(:,i,zval)=meanext(i)-sortedtrack(:,i,9);
end
%%
% % Alternative import for preformatted data
sortedtrack =  wtcombined; 
zval=size(sortedtrack,3)
%%
%%Histograms spread of the x positions of the cells
delx = sortedtrack(:,2:end,9)-sortedtrack(:,1:end-1,9); %Calculates change in x per cell
disp(['Number of tracks = ' num2str(length(sortedtrack))])
%%

%Sets boundaries for binning
% Currently divides the entire cell field into 25 "rows", each containing
% 1/25th of the cells 
% 
Numberofrows =25;
bounds=zeros(size(sortedtrack,2),Numberofrows);
for i=1:size(sortedtrack,2)
    for j=0:Numberofrows
        bounds(i,j+1)=quantile(sortedtrack(:,i,zval),j/Numberofrows);
    end
end
%Sanity Check Plots
plot(bounds)
xlabel('Frame')
ylabel('X-position')
hold off
clf
BoundsrelativetoBack = bounds(:,:)-bounds(:,1);
plot([.5:.5:size(BoundsrelativetoBack,1)/2],BoundsrelativetoBack(:,1:25)/1.53)
ylabel("X-position relative to back of frame, um")
xlabel("Hour")
clf
BoundWidths = bounds(:,2:end)-bounds(:,1:end-1);
plot([.5:.5:size(BoundWidths,1)/2],BoundWidths/1.53);
axis([0,25,0,150])
ylabel("Width of bin, um")
clf

%%

% Assigns each cell a bin based on previous cut offs. 
BinInd= zeros(size(sortedtrack,1),size(sortedtrack,2));
for j=1:size(sortedtrack,2)
    for i=1:25  
    BinInd(:,j)= BinInd(:,j)+(sortedtrack(:,j,zval) > bounds(j,i) & sortedtrack(:,j,zval)< bounds(j,i+1))*i;
    end
end
histogram(BinInd(:,10))
xlabel('Bin assigned (If N/A value, assigned to 0)')
ylabel('# of traces')
%%
%Sanity check plot: The x positions included in each bin- possible failure
%case that the tissue is too long or short and the binning over or
%undershoots accordingly
scatter(sortedtrack(:,23,zval),BinInd(:,23))
ylabel('Bin ID')
xlabel('Cell X Position')
clf
%%
Numberofrows =25;
delx = sortedtrack(:,2:end,9)-sortedtrack(:,1:end-1,9); %Calculates change in x per cell

%Calculating the means of major stats for each bin
Speedmeans=[]; Polmeans=[]; Xmeans=[]; Tissuemeans=[]; Delmeans=[];EllipAspRatioMeans=[];
for i=1:size(sortedtrack,2)
    for j=1:Numberofrows
        Polmeans(j,i) = mean(sortedtrack(BinInd(:,i)==j,i,4),'omitmissing');
        Xmeans(j,i) = mean(sortedtrack(BinInd(:,i)==j,i,37),'omitmissing');
        Speedmeans(j,i) = mean(sortedtrack(BinInd(:,i)==j,i,73),'omitmissing');
        EllipAspRatioMeans(j,i) = mean(sortedtrack(BinInd(:,i)==j,i,29),'omitmissing');
        AreaMeans(j,i) = mean(sortedtrack(BinInd(:,i)==j,i,30),'omitmissing');
        Tissuemeans(j,i)=tissuealign(sortedtrack(BinInd(:,i)==j,i,4),sortedtrack(BinInd(:,i)==j,i,35));
    end
end
%Same thing, but for the change in X position 
for i=1:size(sortedtrack,2)-1
    for j=1:25
Delmeans(j,i) = mean(-1*delx(BinInd(:,i)==j,i),'omitmissing');
end
end
%%
%Speed heatmap
h=heatmap(Speedmeans(1:25,[2:36])'.*2/1.53);
%Conversion factor pixels/halfhour -> um/ hour is 2/1.53
h.GridVisible = 'off';
Ax = gca;
title('Speed, um/hr')
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.XLabel= 'Position along tissue';
h.YLabel= 'Time (hr)';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%%
%Area Heatmap 
h=heatmap(AreaMeans(1:25,[2:36])'*(1/1.53)^2);
%Conversion factor pixels^2 -> um^2 is (1/1.53)^2
h.GridVisible = 'off';
Ax = gca;
title('Area, um^2')
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.XLabel= 'Position along tissue';
h.YLabel= 'Time (hr)';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
h.ColorLimits = [400 900];
%%
h=heatmap(Delmeans(1:25,[2:36])');
h.GridVisible = 'off';
Ax = gca;
title('Change in X-position')
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

%%
h=heatmap(Tissuemeans(1:25,[2:36])');
h.ColorLimits=[0.05 0.3];  %Color Limits
h.GridVisible = 'off';
Ax = gca;
title('Normalized Vector Average Polarity')

Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


%%

recell = num2cell(2:48);
recell2 = num2cell(1:0.5:24);%setting up the heatmap axes
rows = num2cell(([1:25]));
h=heatmap(Xmeans');
h.ColorLimits=[0.06 0.11];  %Color Limits
h.GridVisible = 'off';
Ax = gca;
title('X-Magnitude Means')
h.YDisplayData = recell;
h.YDisplayLabels = recell2;
h.XDisplayData = rows;
h.YDisplayData = string(2:36);
h.XLabel= 'Position along tissue';
h.YLabel= 'Time (hr)';
h.XDisplayLabels = string(25:-1:1);
%%
%Now setting up for Markov models
%Set stringency of thresholding
%Based on testing: Bulk: Prim = 2, Sec = 1.5
InitialThreshold = 2;
SecondaryThreshold = 1.5;
%%
%Histogram the mean x-displacment and attempts to fit to a double guassian,
%assuming smaller curve is nonmigratory cells
%Sets border between migratory and nonmigratory cells at 2 SDs above mean
%of smaller curve
h= histogram(Delmeans,'BinWidth',0.3,'Normalization','pdf','DisplayStyle','stairs');
counts = h.Values;
binEdges = h.BinEdges;
clf
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
[fite,gof] = fit(binCenters',counts','gauss2')%,'StartPoint',[0.5,3,1,0.1,15,3],'Lower',[0,0.1,0,0.05,1,0])
hold on
plot(fite,binCenters,counts)
xline(fite.b1+InitialThreshold*fite.c1)
hold off
clf
delborder = fite.b1+InitialThreshold*fite.c1 %pixels/half hour

subplot(1,2,1)
truths = Delmeans(1:25,1:36)<delborder;
h=heatmap(double(~truths)');
title('Migratory Cells')
colorbar off
h.YDisplayData = string(1:36);
h.XDisplayData = string([1:25]);
h.XDisplayLabels = string(25:-1:1);
subplot(1,2,2)
h=heatmap((Delmeans(1:25,2:36))'*2/1.53);
title('X-displacement')
h.ColorLimits = [0,30];
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.GridVisible = 'off';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clf
%%
hold on
h= histogram(Delmeans(truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.15);
histogram(Delmeans(~truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.5);
legend("Speed<Thresh",'Speed>Thresh')
title('X-Displacment')
hold off 
clf
%%
%Setting threshold for Xcomponent magnitude
hold on
h= histogram(Xmeans(truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.001);
histogram(Xmeans(~truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.001);
legend("Speed<Thresh",'Speed>Thresh')
title('X-component Magnitude')
hold off 

counts = h.Values;
binEdges = h.BinEdges;
clf
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
[fite,gof] = fit(binCenters',counts','gauss1','StartPoint',[20,.1,0.1])
hold on
plot(fite,binCenters,counts)
xline(fite.b1+SecondaryThreshold*fite.c1)
hold off
clf
xborder = fite.b1+SecondaryThreshold*fite.c1

subplot(1,2,1)
xtruths = Xmeans(1:25,1:36)<xborder;
h=heatmap(double(~xtruths)');
title('Plarized Cells');
colorbar off
h.YDisplayData = string(1:36);
h.XDisplayData = string([1:25]);
h.XDisplayLabels = string(25:-1:1);
subplot(1,2,2)
h=heatmap((Xmeans(1:25,2:36))');
title('X Mag Means')
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.GridVisible = 'off';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clf

%%
%Setting threshold for Normalized vector average polarity
hold on
h= histogram(Tissuemeans(truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.01);
histogram(Tissuemeans(~truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.01);
legend("Speed<Thresh",'Speed>Thresh');
title("Tissue Alignment");
hold off 

counts = h.Values;
binEdges = h.BinEdges;
clf
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
[fite,gof] = fit(binCenters',counts','gauss1','StartPoint',[20,.1,0.1])
hold on
plot(fite,binCenters,counts)
xline(fite.b1+SecondaryThreshold*fite.c1)
hold off
clf
tissueborder = fite.b1+SecondaryThreshold*fite.c1

subplot(1,2,1)
ttruths = Tissuemeans(1:25,1:36)<tissueborder;
h=heatmap(double(~ttruths)');
title('Polarized Cells')
colorbar off
h.YDisplayData = string(1:36);
h.XDisplayData = string([1:25]);
h.XDisplayLabels = string(25:-1:1);
subplot(1,2,2)
h=heatmap((Tissuemeans(1:25,2:36))');
title('TA Means')
% h.ColorLimits = [0,30];
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.GridVisible = 'off';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clf

%%
%Setting threshold for Polarity magnitude
hold on
h= histogram(Polmeans(truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.002);
histogram(Polmeans(~truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.002);
legend("Speed<Thresh",'Speed>Thresh')
title('Polarity Magnitude')
hold off 

counts = h.Values;
binEdges = h.BinEdges;
clf
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
[fite,gof] = fit(binCenters',counts','gauss1','StartPoint',[20,.1,0.1])
hold on
plot(fite,binCenters,counts)
xline(fite.b1+SecondaryThreshold*fite.c1)
hold off
clf
polborder = fite.b1+SecondaryThreshold*fite.c1

subplot(1,2,1)
ptruths = Polmeans(1:25,1:36)<polborder;
h=heatmap(double(~ptruths)');
title('Polarized Cells');
colorbar off
h.YDisplayData = string(1:36);
h.XDisplayData = string([1:25]);
h.XDisplayLabels = string(25:-1:1);
subplot(1,2,2)
h=heatmap((Polmeans(1:25,2:36))');
title('Pol Means')
% h.ColorLimits = [0,30];
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
h.GridVisible = 'off';
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clf

%%
%Setting threshold for Speed
hold on
h= histogram(Speedmeans(truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.2);
histogram(Speedmeans(~truths),'DisplayStyle','stairs','Normalization','percentage','BinWidth',0.2);
legend("Speed<Thresh",'Speed>Thresh')
title('Speed')
hold off 

counts = h.Values;
binEdges = h.BinEdges;
clf
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
[fite,gof] = fit(binCenters',counts','gauss1','StartPoint',[20,5,1])
hold on
plot(fite,binCenters,counts)
xline(fite.b1+SecondaryThreshold*fite.c1)
hold off
clf
Speedborder = fite.b1+InitialThreshold*fite.c1
%%

function TA = tissuealign(Polmag,Polang)

%inputs
%Currently built to take a matrix w/ rows as individual cells and columns
%as frame 
%Polarity magnitude 
% Polarity angle (in radians)
%Actual Computation 
totmag = sum(Polmag,1,'omitmissing'); %Compute Maximum Length
%Set up to transform from 180 domain to 360 domain -> this is the most
%rigorous way to do this I can think of. It will systematically
%underestimate alignment, but everything is "correct" from an axial
%perspective
%Compute x and y components in new space
ycomp2 = Polmag.*sin(2*Polang);
xcomp2 = Polmag.*cos(2*Polang);
%Compute vector sum and length of resulting vector
xtot = sum(xcomp2,1,'omitmissing');
ytot = sum(ycomp2,1,'omitmissing'); 
leggy = sqrt((xtot).^2+(ytot).^2);
%Normalize the length of the vector sum vs maximum length
TA = leggy./totmag;
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:04613fa7]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Incorrect number or types of inputs or outputs for function removevars."}}
%---
