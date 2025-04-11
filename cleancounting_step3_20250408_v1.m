% time for step 3: take in the thresholds and get the site occupancies
% (binarize latcounts)

%load('step2_20250331.mat');
%load('step2_20250409.mat');
%Stat.LatCount is 7651x1x20


%% need the ring lattice occupation stuff. What goes into this: latcountsout(sites by 1 by number of shots by number of files), matrix A (rings by files)
%latcountsin = reshape(latcountsout,2); % now it is sites by shots by files
%%
RingLatOccup=zeros(NumSite,1,numShots,numfiles);
RingLatUnoccup=zeros(NumSite,1,numShots,numfiles);
LatCountin= zeros(NumSite,1,numShots,numfiles);

% make the threshold matrix go from 1 to 7651 rather than from 1 to 50
% Inputs:
% thresholds: 50 x 3 (rings x data files)
% sitesPerRing: a vector of length 50 containing number of sites per ring

sitesPerRing = ringPop;  % Fill this with all 50 values!
% Make sure length(sitesPerRing) == 50

%numSites = sum(sitesPerRing);          % Should be 7651
%numFiles = size(thresholds, 2);        % Should be 3

thresholdMatrix = zeros(numSites, numfiles);  % Preallocate output
siteIndex = 1;

for ring = 1:length(sitesPerRing)
    nSites = sitesPerRing(ring);
    idx = siteIndex:(siteIndex + nSites - 1);     % Index range for this ring
    thresholdMatrix(idx, :) = repmat(A(ring, :), nSites, 1);  % Copy threshold
    siteIndex = siteIndex + nSites;
end

% Assume:
% data is 7651 x 1 x 20 x 3
% thresholdMatrix is 7651 x 3

% Step 1: Expand thresholdMatrix to match size of data
% We want it to become 7651 x 1 x 20 x 3
data = latcountsout;
thresholdExpanded = reshape(thresholdMatrix, [numSites, 1, 1, numfiles]);
thresholdExpanded = repmat(thresholdExpanded, [1, 1, numShots, 1]);  % Now same size as data

% Step 2: Compare data to thresholds
RingLatOccup   = data > thresholdExpanded;   % Logical matrix: 1 if count > threshold
RingLatUnoccup = data < thresholdExpanded;   % Logical matrix: 1 if count < threshold

%% Filling fraction vs. ring for each file

%it seems like this worked up until the last 2 data sets in the "goodforlooping" folder.

image = 1;


hexloccups2 = squeeze(sum(RingLatOccup,3)); % this is the sum of the lattice occupancies across each data set
totatomnum = sum(hexloccups2,1)/numShots;
%comparetotnum = numSites(1)*fill_mean;
 % total atom number vs time
figure;
%c=1:numfiles;
%h=scatter(labels,totatomnum,[],c,"filled")
scatter(labels,totatomnum,"filled")
colormap(gca,"parula")
%set(h,{"MarkerFaceColor"},get(h,'Color'));
xlabel('Time (ms)');
ylabel('Atom number (r = 50 ROI, 7651 sites)')
%%
RingOccups2=zeros(hexrad,numfiles);
RingFill2 = zeros(hexrad,numfiles);
% this array is 1951 by 13 (number of sites by number of data sets)
 for ii = 1:numfiles
     for j = 1:hexrad
         RingOccups2(j,ii) = sum(hexloccups2(ringidx{:,j},ii))/numShots;
         RingFill2(j,ii) = RingOccups2(j,ii)/ringPop(j);
     end
 end

 %% % and std filling
% Assume: RingLatOccup = 7651 × 20 × 14 (after squeeze)
RingLatOccup2 = squeeze(RingLatOccup);
sitesPerRing = ringPop;  % length 50

numRings = hexrad;
%numShots = size(RingLatOccup, 2);  % 20
numFiles = numfiles;  % 14

% Preallocate output: std deviation per ring per file
ringStd = zeros(numRings, numFiles);
fillingperShot = zeros(numShots,numfiles);

siteIdx = 1;

for r = 1:numRings
    nSites = sitesPerRing(r);
    idx = siteIdx:(siteIdx + nSites - 1);

    % Extract data for this ring: size nSites × 20 × 14
    ringData = RingLatOccup2(idx, :, :);

    % Compute filling per shot: mean across sites
    fillingPerShot = squeeze(mean(ringData, 1));  % size: 20 × 14

    % Compute std dev across shots, for each file
    ringStd(r, :) = std(fillingPerShot, 0, 1);  % 1 × 14

    siteIdx = siteIdx + nSites;
end
%%
ringsem = ringStd/sqrt(numShots);
 %%

 % with that, I have filling vs. ring for all 13 times. make a plot with
 % hold time and filling vs ring plots...
ringnum = 1:hexrad;
 figure;
 hold on;
 colors = parula(numfiles);
 chosenfiles = [1;8;13];
% colors = lines(13);
 for p=1:length(chosenfiles)
     i = chosenfiles(p);
     plot(ringnum,RingFill2(:,i),'Color',colors(i,:),'LineWidth',2,'DisplayName',num2str(labels(i)));
 end
 xlabel('Ring')
 ylabel('Filling fraction')
 title('Filling fraction vs ring for different times, blue curve 20240725 data');
 legend('show')
 grid on;

 ringfill5ms = RingFill2(:,1);
 %%
 ringnum = 1:hexrad;

figure;
hold on;

% Choose PRL-like color palette
colors = parula(numfiles);
chosenfiles = [1; 8; 13];

% Plot each dataset
for p = 1:length(chosenfiles)
    i = chosenfiles(p);
    plot(ringnum, RingFill2(:,i), ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s', num2str(labels(i))));
end

% Axes labels (with units if relevant)
xlabel('Ring number', 'FontSize', 12, 'FontName', 'Helvetica');
ylabel('Filling fraction', 'FontSize', 12, 'FontName', 'Helvetica');

% Axis settings
set(gca, ...
    'FontSize', 12, ...
    'FontName', 'Helvetica', ...
    'LineWidth', 1, ...
    'Box', 'on', ...
    'TickDir', 'in', ...
    'TickLength', [0.015 0.015]);

% Optional: fine-tune axis limits or ticks
% xlim([0 hexrad+1]);
% ylim([0 1]);

legend('show', ...
    'Box', 'off', ...
    'FontSize', 12, ...
    'Location', 'northeast');

% grid off; % PRL figures often don't have grids unless essential

% Optional: set figure size for export (in inches)
% set(gcf, 'Units', 'Inches', 'Position', [1, 1, 3.4, 3]);  % ~8.6cm width
% set(gcf, 'PaperPositionMode', 'auto');  % Ensures correct export size

% Set paper size to match figure size for clean export
set(gcf, 'Units', 'Inches');
figPos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperSize', [figPos(3), figPos(4)]);
set(gcf, 'PaperPosition', [0, 0, figPos(3), figPos(4)]);
set(gcf, 'PaperPositionMode', 'manual');

% Export as PDF
print(gcf, 'filling_vs_ring_v2', '-dpdf', '-r300');

%% With error bars
figure;
hold on;

colors = parula(numfiles);   % Or your favorite color map
chosenfiles = [1; 8; 13];

for p = 1:length(chosenfiles)
    i = chosenfiles(p);

    errorbar(ringnum, RingFill2(:,i), ringsem(:,i), ...
        'Color', colors(i,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s', num2str(labels(i))));
end

xlabel('Ring number', 'FontSize', 12, 'FontName', 'Helvetica');
ylabel('Filling fraction', 'FontSize', 12, 'FontName', 'Helvetica');

set(gca, ...
    'FontSize', 12, ...
    'FontName', 'Helvetica', ...
    'LineWidth', 1.2, ...
    'TickDir', 'out', ...
    'TickLength', [0.015 0.015]);

box on;
legend('show', 'FontSize', 12, 'Location', 'northeast');

% Optional: export to PDF
% print(gcf, 'filling_vs_ring_with_errorbars', '-dpdf', '-r300');
% Set paper size to match figure size for clean export
set(gcf, 'Units', 'Inches');
figPos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperSize', [figPos(3), figPos(4)]);
set(gcf, 'PaperPosition', [0, 0, figPos(3), figPos(4)]);
set(gcf, 'PaperPositionMode', 'manual');

% Export as PDF
print(gcf, 'filling_vs_ring_v3', '-dpdf', '-r300');

 %% actually try doing this in a nested loop: total filling in each time as a function of the ROI
figure;
for i = 1:15;
    numsiteloop = siteSum(i);
    fillingloop = sum(RingOccups2(1:i,:),1)/numsiteloop;
    plot(labels,fillingloop);
    hold on
end
%%
%%
figure;
hold on;
legendLabels = cell(1, 15); % Preallocate a cell array for legend labels

for i = 1:30
    numsiteloop = siteSum(i);
    fillingloop = sum(RingOccups2(1:i, :), 1) / numsiteloop;
    plot(labels, fillingloop);
    
    legendLabels{i} = sprintf('Number of sites: %d', numsiteloop); % Store legend label
end

legend(legendLabels); % Apply the legend after the loop
title('Filling fraction vs. time, blue detuned OP for different ROI sizes')
xlabel("Time (ms)")
ylabel("Filling Fraction")
hold off;

%%
%% make another one which is total atom number as a function of ROI Size
figure;
hold on;
legendLabels = cell(1, 15); % Preallocate a cell array for legend labels

for i = 1:50
    numsiteloop = siteSum(i);
    atomsloop = sum(RingOccups2(1:i, :), 1);
    plot(labels, atomsloop);
    legendLabels{i} = sprintf('Number of sites: %d', numsiteloop); % Store legend label
end

legend(legendLabels); % Apply the legend after the loop
title('Total atom number vs. time, blue detuned OP for different ROI sizes')
xlabel("Time (ms)")
ylabel("Total Atom Number")
hold off;

%% Now take the total atom number vs time and plot it alongside the side CCD data (which I will load in here)
load('redist_lifetime_processed.mat');

numvstime = sum(RingOccups2(:,:),1);
t = 1:3000;
scalingit = 1.8;
sidenumvstimefit = (236*exp(-t/48.6) + 849*exp(-t/5766))*scalingit;

startnumtot = validnewY(1)*1.8;

figure;
plot(labels,numvstime,'r--o','MarkerSize',10);
hold on;
plot(validnewX,validnewY*scalingit,'b--o','MarkerSize',10);
hold on;
plot(t,sidenumvstimefit,'LineWidth',2);
legend('Lower CCD atom number in radius 50 ROI (7651 sites)', ...
    'Side CCD atom number multiplied by 1.8','Fit on side CCD num')
xlabel('Time (ms)')
ylabel('Atoms')

%% compare histograms
thresholdsout=A;
% unpack latcountsout and thresholdsout and plot the histograms vs. ring
% for all of the 13 data sets.
% binwidth=10;
% mus = zeros(1,hexrad);
% 
% for ii = 1:numfiles
%     figure;
%     sgtitle(num2str(labels(ii)));
%     hold on;
%     for j=1:hexrad
%         LatCount=reshape(latcountsout(ringidx{:,j}, 1, :,ii), [], 1);
%         subplot(5,(hexrad/5),j)
%         MaxCount = max(latcountsout(:,:,2,ii));
%         MinCount = min(latcountsout(:,:,2,ii));
%         histogram(LatCount,'BinEdges', MinCount:binwidth:MaxCount, 'EdgeColor','none');
%         % histfit(LatCount, 50);
%         xlabel('Counts')
%         ylabel('Frequency')
%         xlim([-300, MaxCount])
%         ylim([0 30]);
%         % fitting
%         %pd = fitdist(LatCount,'Normal');
%         %mus(j) = pd.mu;
%         xline(thresholdsout(j,ii),'LineWidth',1)
%         hold on
%         title('ring rad = '+ string(j));
% end
% 
% end
%% compare the thresholds
thresholdsout=A;
figure;
 hold on;
 colors = parula(numfiles);
% colors = lines(13);
 for i = 1:numfiles
     plot(RingNum,thresholdsout(:,i),'Color',colors(i,:),'LineWidth',2,'DisplayName',num2str(labels(i)));
 end
 xlabel('Ring')
 ylabel('Threshold')
 title('Threshold vs ring for different times, blue curve 20240725 data');
 legend('show')
 grid on;

%% make a plot of the filling vs. time in different parts of the trap
% center filling = sum(RingLatOccup, first 217 sites)
% slope filling = sum(RingLatOccup, sites corresponding to ring 25 to ring
% 29


