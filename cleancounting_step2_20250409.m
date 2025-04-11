% Step 2 of the cleancounting script. Take in the latcountsout from the
% step 1 workspace. Plot the histograms, perform first round of fitting.
%load('step1_20250328.mat'); % this is just on the first 3 data sets
load('step1_20250409.mat'); % this is all of the data sets
%%
% unpack latcountsout and thresholdsout and plot the histograms vs. ring
% for all of the data sets.
binwidth=10;
mus = zeros(1,hexrad);

for ii = 1:numfiles
    %figure;
    sgtitle(num2str(labels(ii)));
    hold on;
    for j=1:hexrad
        LatCount=reshape(latcountsout(ringidx{:,j}, 1, :,ii), [], 1);
        subplot(5,(hexrad/5),j)
        MaxCount = max(latcountsout(:,:,2,ii));
        MinCount = min(latcountsout(:,:,2,ii));
        histogram(LatCount,'BinEdges', MinCount:binwidth:MaxCount, 'EdgeColor','none');
        % histfit(LatCount, 50);
        xlabel('Counts')
        ylabel('Frequency')
        xlim([-300, MaxCount])
        ylim([0 30]);
        hold on
        title('ring rad = '+ string(j));
end
    
end

%% Fit the thresholds for the first time...
%% Fitting test count histogram 20250218 trying again- summed distributions to fit better at larger rings
% here is the code that creates figures that explain what's going on here.
% adapt the bulk of this code into getThreshold4

% 20250306 I am trying ring-dependent guesses for the Gauss and Maxwell
% dists
% i am using linear fits of the fit parameters at different rings to get
% linear equations to feed the guess
files1 = 1:numfiles;
lookatfiles = length(files1);
thresholdsgaussmax = zeros(hexrad,lookatfiles); % this is fine if we're doing one file at a time.
thresholdsgauss =zeros(hexrad,lookatfiles);
thresholdschosen = zeros(hexrad,lookatfiles);
%lookatfiles = 1;
for ii=1:lookatfiles; % this specifies the file we are looking at. 
    % Need to get the actual full 20 shot histogram out
% testing with 2 Gaussians
rows = 5;
cols = 10;
files2 = files1(ii);
%figure;
%t = tiledlayout(rows, cols, 'TileSpacing', 'compact');


%ring=hexrad;
for j = 1:hexrad
    ring = j;
    image=1;

% ring = 45;
% image = 1;

% Selecting count data indexed by ring
MaxCount = max(reshape(latcountsout(:, image, j), [], 1));
MinCount = min(reshape(latcountsout(:, image, j), [], 1));
%LatCount=reshape(Stat.LatCount(ringidx{:,ring}, image, j), [], 1);
LatCount = reshape(latcountsout(ringidx{:,j},1,:,ii),[],1);

% Plotting to see bare histogram
% figure
% histogram(LatCount,'BinEdges', 0:binwidth:MaxCount, 'EdgeColor','none');

% Fitting histograms
[y,edges] = histcounts(LatCount,100);
binwidth = edges(2)-edges(1);
x = edges(1:end-1)+binwidth/2;
[maxy,peakIdx]= max(y);
peakX = x(peakIdx);
stdX = std(LatCount);
com = 1/sum(y)*sum(x.*y);
thresh=com;

ring1 = 5;

foc = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1,x(1),binwidth,0, 0, 10],...
    'Upper',[2*max(y),thresh,thresh-min(x),10, 100000, 10000],...
    'StartPoint',[3.86*ring+3.7,-30.2*ring+1145,-5.5*ring+275,-5*10^-7*ring+0.0000573, -50*ring+1629, -8.55*ring+824]);

gaussEqn = 'a0*exp(-(x-a1)^2/(2*a2^2))';

sumgaussmax = 'a0*exp(-(x-a1)^2/(2*a2^2))+b0*(max((x-b1),0))^2*exp(-(x-b1)^2/(2*b2^2))';

sumgauss = 'a0*exp(-(x-a1)^2/(2*a2^2))+b0*exp(-(x-b1)^2/(2*b2^2))';

ftc = fittype(sumgaussmax,...
    'independent',{'x'},...
    'dependent',{'y'},...
    'coefficients',{'a0','a1','a2','b0','b1','b2'},...
    'options',foc);


lengthc = length(x);
[cfit,cgof]=fit(x',y',ftc);


% Compute fitted parameters
xc = 1:6000;
paramsC = coeffvalues(cfit); % [a0, a1, a2,b0,b1,b2]
gauss1c = paramsC(1).*exp(-(xc-paramsC(2)).^2/(2.*paramsC(3)^2));
%gauss2c = paramsC(4).*exp(-(xc-paramsC(5)).^2/(2*paramsC(6)^2));
maxc = paramsC(4).*(max(xc-paramsC(5),0)).^2.*exp(-(xc-paramsC(5)).^2/(2*paramsC(6)^2));
gauss2c = maxc;

% set threshold by integrating
cumgauss1 = cumsum(gauss1c,"reverse");
cumgauss2 = cumsum(gauss2c,"forward");
% figure
% plot(xc, cumgauss1)
% hold on
% plot(xc,cumgauss2)

% test threshold: where these cumulative fns cross
diff = cumgauss1-cumgauss2;
% figure
% plot(diff)

target_value = 1.75E-3; % set to 2.5 sigma for Gaussian
crossidx = find(diff(1:end-1) > target_value & diff(2:end) < target_value);

Threshold1 = crossidx;


% contx = xc; afit = gauss1c; bfit = gauss2c;
% % overlap region
% idxoa = (contx >=Threshold1 & contx<= contx(end));
% idxob = (contx <= Threshold1 & contx >= contx(1));
% % Calculate imaging fidelity
% intatot = trapz(xc,afit);
% inta = trapz(xc(idxoa),afit(idxoa));
% intbtot =trapz(xc,bfit);
% intb = trapz(xc(idxob),bfit(idxob));
% Fidelity = 1-(inta+intb)/(intatot+intbtot-inta-intb);


% fod with guesses for a solo Gaussian
fod = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[1,x(1),binwidth],...
    'Upper',[2*max(y),thresh,thresh-min(x)],...
    'StartPoint',[max(y),thresh/2,(thresh-min(x))/2]);

gaussEqn = 'a0*exp(-(x-a1)^2/(2*a2^2))';

ftd = fittype(gaussEqn,...
    'independent',{'x'},...
    'dependent',{'y'},...
    'coefficients',{'a0','a1','a2'},...
    'options',fod);


lengthd = length(x);
[dfit,dgof]=fit(x',y',ftd);

% Compute fitted parameters
xd = 1:6000;
paramsD = coeffvalues(dfit); % [a0, a1, a2]
gauss1d = paramsD(1).*exp(-(xc-paramsD(2)).^2/(2.*paramsD(3)^2));

target_value = 1.75E-3; % set to 2.5 sigma for Gaussian
%target_value =0.1; % set to 2.5 sigma for Gaussian
%target_value = 0 ; % now that i'm using 2 gaussians
% crossidxd = find(gauss1d(1:end-1) > target_value & gauss1d(2:end) < target_value);
% if crossidxd ==[]
%     crossidxd = 3*paramsD(3);
% end

crossidxd = find(gauss1d(1:end-1) > target_value & gauss1d(2:end) < target_value);

if isempty(crossidxd)
    crossidxd = 3 * paramsD(3); % Ensure paramsD(3) exists
end

if isempty(crossidx)
    crossidx = 0;
end


thresholdsgauss(j,ii) = crossidxd;
thresholdsgaussmax(j,ii) = crossidx;

if dgof.rmse <= cgof.rmse
    Threshold1 = crossidxd;
else 
    Threshold1 = crossidx;
end

% if Threshold1<=275
%     Threshold1 = 275;
% end

thresholdschosen(j,ii)=Threshold1;
% Plotting histogram fitted to a single Gaussian

% nexttile;
% 
% histogram(LatCount, 'BinEdges', MinCount:binwidth:MaxCount, 'EdgeColor', 'none');
% hold on
% sgtitle(num2str(labels(ii)));
% plot(dfit,'k');% 'HandleVisibility', 'off'); % No legend
% plot(cfit,'r'); %, 'HandleVisibility', 'off'); % No legend
% plot(xc, gauss1c, 'g','HandleVisibility', 'off');
% plot(xc, gauss2c, 'g', 'HandleVisibility', 'off');
% %plot(xd, gauss1d, 'b','HandleVisibility', 'off');
% 
% xline(crossidx, 'HandleVisibility', 'off');
% xline(Threshold1, 'r', 'HandleVisibility', 'off');
% xlim([-1000 4000])
% title(num2str(j));
% 
% legend off; % Explicitly remove legends
% 
% drawnow;

end

end

%%
% plot the thresholds out: the one from gaussmax, the one from single
% gauss, and the one from thresholdschosen
thresholdsgaussmax1 = thresholdsgaussmax(:,1);
thresholdsgauss1 = thresholdsgauss(:,1);
thresholdschosen1 =thresholdschosen(:,1);
rings = 1:hexrad;
figure;
plot(rings,thresholdschosen1,'o');
hold on
plot(rings,thresholdsgauss1,'k');
hold on
plot(rings,thresholdsgaussmax1,'b');


%%
% add some code to bring the threshold up to the tail value if it dips

A = thresholdschosen;

% Step 1: Set any element lower than 275 to 275
A(A < 350) = 350;
% 
% for ii = 1:numfiles
%     for j = 25:hexrad
%         if A(j,ii)<-31*j+1265
%             A(j,ii) = -31*j+1265
%         end
% 
%     end
% end

% Step 2: Check for outliers and adjust them
for col = 1:size(A, 2)
    for row = 2:size(A, 1) % Start from the second row
        if abs(A(row, col) - A(row-1, col)) > 375
            % Set the value to the previous value minus 100
            A(row, col) = A(row-1, col) - 100;
        end
    end
end

% Display the adjusted array
disp('Adjusted array:');
disp(A);
A_avg = mean(A,2);
%%

% Step 2: Flag elements that differ from the element above by more than 100
% Create a logical array of the same size as A
flagged = false(size(A));

% Loop through each column
for col = 1:size(A, 2)
    for row = 2:size(A, 1) % Start from the second row
        if abs(A(row, col) - A(row-1, col)) > 375
            flagged(row, col) = true;
        end
    end
end

% Display the flagged elements (optional)
[rowIdx, colIdx] = find(flagged);
for i = 1:length(rowIdx)
    fprintf('Outlier found at (%d, %d): %d\n', rowIdx(i), colIdx(i), A(rowIdx(i), colIdx(i)));
end

%%
figure;
%scatter(rings, thresholdschosen, 50, 'r', 'filled'); % 50 is the marker size
hold on
%plot(rings, thresholdsgauss, 'k');  % Black line
scatter(rings, thresholdsgaussmax, 'b','DisplayName','Threshold from Gaussian and Maxwell');  % Blue line
plot(rings, A,'g','DisplayName','Smoothed Threshold'); % add one with the adjusted array
plot(rings,A_avg,'m','LineWidth',2,'DisplayName','Averaged across times smoothed threshold');
xlabel('Ring Number (lattice sites)');
ylabel('Threshold (counts per site)');
%legend('Threshold from Gaussian and Maxwell','Smoothed Threshold','Averaged across times smoothed threshold')
legend;

%% And plot some exemplary histograms with the thresholds on them... chosen, smoothed, etc, etc

figure;
rows = 5;
cols = 3;
files2 = files1(ii);
%figure;
%t = tiledlayout(rows, cols, 'TileSpacing', 'compact');
ii = 1;
for j = 25:39
LatCount = reshape(latcountsout(ringidx{:,j},1,:,ii),[],1); % j is ring, ii is file.
nexttile;
histogram(LatCount, 'BinEdges', MinCount:binwidth:MaxCount, 'EdgeColor', 'none');
hold on
sgtitle(num2str(labels(ii)));
% plot(dfit,'k');% 'HandleVisibility', 'off'); % No legend
% plot(cfit,'r'); %, 'HandleVisibility', 'off'); % No legend
% plot(xc, gauss1c, 'g','HandleVisibility', 'off');
% plot(xc, gauss2c, 'g', 'HandleVisibility', 'off');
%plot(xd, gauss1d, 'b','HandleVisibility', 'off');

xline(A(j,ii), 'k','HandleVisibility', 'off');
xline(thresholdschosen(j,ii), 'r', 'HandleVisibility', 'off');
xlim([-1000 4000])
title(num2str(j));

legend off; % Explicitly remove legends

end

% What leaves this script is A- the smoothed thresholds. It is size (# of
% rings x # of files)