%% Cleaned up version of the counting code. 
% Step 1: Loop through, get latcountsout arranged by ring for each file
% run this and save the workspace. Load the workspace into the next script.
%% Loop through folder of data, label them with the redist time
%Bg = load('C:\Users\qmspc\Documents\MATLAB\Test_app_analysis\full analysis test\data\CleanBackground_20240327_HSSpeed=2_VSSpeed=1.mat');
Bg = load('CleanBackground_20240327_HSSpeed=2_VSSpeed=1.mat')% datafolder="C:\Users\qmspc\Desktop\NewLabData\2024\06 June\20240627\amgcp_0.5_freqscan\"; % June 27 am gcp 0.5 data
datafolder = "/Users/Lauren/Documents/MATLAB/Counting Atoms 2025/test2/" % new file path for blue redist data
% previously ran for test1
imgfiles=dir(datafolder+"*.mat");
numfiles=length(imgfiles);
labels=[];
subImgNum=1; % this is FK1 data so that works.
fill_mean=zeros(numfiles,1);
fill_std = zeros(numfiles,1);

hexrad = 50

ROI_shape='Hexagon'; % Current options: 'Rectangular', 'Hexagon'
RingNum=1:hexrad;


if ROI_shape == 'Hexagon'
     % Site is now arranged in rings
    ringPop = zeros(1, hexrad);
                for j = 1:hexrad
                    siteSum(j)=1;
                    for i=1:j
                        siteSum(j) = siteSum(j)+6*i;
                    end
                    ringPop(j)= 6*j;
                end
    ringPop(1)=7;
    ringidx{:,1}=1:7;
    for i = 2:hexrad
        ringidx{:,i}= siteSum(i-1)+1:siteSum(i);
    end
end
numSites = siteSum(hexrad)
%%;
% 1951 for 25
numShots = 20;
latcountsout = zeros(numSites,1,numShots,numfiles);

RingFill = zeros(hexrad,numfiles);
% loop through the folder and get all of the filling
% fractions as a fn of ring

for ii=1:numfiles
    load(datafolder+imgfiles(ii).name)
    labels(ii)=str2num(cell2mat(regexp(imgfiles(ii).name,'\d?\.?\d+','match')))-13754.55;
    numFK = Data.Andor19330.Config.NumFrames; %Fast kinetic number (e.g. 2, 4, 8)
    numShots=size(Data.Andor19330.Image,3); %number of shots in Data
    SubImg = zeros(numFK, 2);
    for i = 1:numFK
        SubImg(i, 1) = 1+(i-1)*1024/numFK;
        SubImg(i, 2) = i*1024/numFK;
    end
    Data.Andor19330.SubImg=SubImg;
    % Fitting entire sample to Gaussian
[xcenter, xwidth, ycenter, ywidth,xlimits,ylimits] = funFitGaussXY(Data, Bg);

% Extract lattice vectors, generate lattice sites, group counts by site
[Stat,Site,NumSite,Lat,MeanSum,MeanBg,BgOffset] = getStatFull_edit(Data,Bg,xcenter,ycenter,ROI_shape,hexrad,ringidx);
latcountsout(:,:,:,ii) = Stat.LatCount;

end
%%
 %%
 shot = 1;
 PT3= Site*Lat.V+Stat.LatOffset(shot,:)+[(numFK-1)*(1024/numFK),0];
%PT4 = Site(20:37,:)*Lat.V+Stat.LatOffset(shot,:)+[(numFK-1)*(1024/numFK),0];
figure
imagesc(double(Data.Andor19330.Image(:,:,shot))-mean(Bg.Data.Andor19330.Image,3));daspect([1 1 1]); caxis([0, 200]);colorbar; hold on; 
%scatter(PT(:,2),PT(:,1));
hold on;
% scatter(PT2(:,2),PT2(:,1), "filled");
scatter(PT3(:,2),PT3(:,1), "filled", 'or');
%scatter(PT4(:,2),PT4(:,1), "filled", 'or');

%% save the workspace to use in step 2
filename = ['step1_a_' datestr(now, 'yyyymmdd') '.mat'];
save(filename);


%% funFitGaussXY

function [xc, wx, yc, wy,xlimits,ylimits] = funFitGaussXY(Data, Bg, varargin)
    numFK=Data.Andor19330.Config.NumFrames;
    fkShift=floor(1024/numFK);
    d=mean(Data.Andor19330.Image((1024-fkShift+1):1024,1:1024,:),3)-mean(Bg.Data.Andor19330.Image((1024-fkShift+1):1024,1:1024,:),3);
    % d=mean(Data.Andor19330.Image((1024-fkShift+1):1024,1:1024,:),3); % no background subtraction

    [szY, szX] = size(d);

    xvec = linspace(1,szX,szX); % x coords if one is not provided
    yvec = linspace(1,szY,szY); % y coords if one is not provided
    
    xsum = sum(d,1); % integrated X profile (summed over rows)
    ysum = sum(d,2)'; % integrated Y profile (summed over cols)
    
    [xfit,xgof]=fitGaussX(xvec, xsum);
    [yfit,ygof]=fitGaussX(yvec, ysum);

    wx=xfit.w; % horizontal 1/e^2 radius
    wy=yfit.w; % vertical 1/e^2 radius
    xc=xfit.xc; % x0,y0 ellipse centre coordinates
    yc=yfit.xc;
    sx=xfit.w/2; % horizontal 1sigma
    sy=yfit.w/2; % horizontal 1sigma
    t=-pi:0.01:pi;
    xw=xc+wx*cos(t);
    yw=yc+wy*sin(t);
    xs=xc+sx*cos(t);
    ys=yc+sy*sin(t);

% hold off
    function [gXYfit, gXYgof] = fitGaussX(xvec, xsum)
    gaussXYEqn = 'a0*exp(-2*((x-xc)^2/w^2))+c'; % This fitted width is the 1/e^2 waist, NOT 1sigma.
    
    foLow = [0.001*min(xsum),xvec(1),0.01*(xvec(end)-xvec(1)),min(xsum)];
    foUp = [max(xsum),xvec(end),10*(xvec(end)-xvec(1)),max(xsum)];
    foStart = [mean(xsum),xvec(floor(length(xvec)/2)),0.25*(xvec(end)-xvec(1)),min(xsum)];
    
    fo = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',foLow,...
        'Upper',foUp,...
        'StartPoint',foStart);
    
    ft = fittype(gaussXYEqn,...
        'independent',{'x'},...
        'dependent',{'y'},...
        'coefficients',{'a0','xc','w','c'},...
        'options',fo);
    
    [gXYfit, gXYgof] = fit(xvec', xsum',ft);
    end

    xlimits=[xc-3*sy xc+3*sy];
    ylimits=[yc-3*sx yc+3*sx];


end

%% getStatFull_edit

function [Stat,Site,NumSite,Lat,MeanSum,MeanBg,BgOffset] = getStatFull_edit(Data,Bg,xc,yc,ROI_shape,hexrad,ringidx,Options)
    arguments
        Data (1,1) struct
        Bg (1,1) struct
        xc (1,1) {mustBeNumeric}
        yc (1,1) {mustBeNumeric}
        ROI_shape
        hexrad
        ringidx

        % Options for loading calibration
        Options.LatMode = "Lower"

        % Options for PSF width for deconvolution kernel generation
        Options.Sigma = 1.86;
        
        % Options for setting thresholds automatically
        Options.ThresholdMode (1,:) double = 0
        Options.MaxIgnore (1,1) int16 = 1
        Options.AbsMax (1,1) double = 10000
        Options.AbsMin (1,1) double = -300
        Options.ThresholdMin (1,1) double = 200
        Options.ThresholdMax (1,1) double = 4200

        % Options for showing graphics
        Options.ShowMeanSig (1,1) logical = false
        Options.ShowFFTPeak (1,1) logical = false
        Options.ShowUpdate (1,1) logical = false
    end

    tic

    BgCenter = [50,20];
    BgR = 2;
    
    % Gaussian + Lorentzian parameters  a*exp(-(x-b)^2/(2*c^2))+d/((x-b)^2+e^2)
    funca=112;
    funcb=0;
    funcc=1.608;
    funcd=197.4;
    funce=4.015;
    funcPSF = @(x,y) (funca*exp(-(x.^2+y.^2)/(2*funcc^2))+funcd./((x.^2+y.^2)+funce^2))/(funca*2*pi*funcc+funcd*pi^2/funce);

    %funcPSF = @(x,y) exp(-0.5*(x.^2+y.^2)/Options.Sigma^2)/(2*pi*Options.Sigma^2);
    %funcPSF = @(x,y) 0.03;
    %funcPSF = @(x,y) exp(-0.5*(x.^2+y.^2)/Options.Sigma^2)/(sqrt(2*pi*Options.Sigma^2));


    % Initialize
    numFK = Data.Andor19330.Config.NumFrames;
    switch numFK
        case 1
            Lat = loadCalibration(Options.LatMode,[yc,xc]);
            %Stat.LatFormat = {'Rectangular',-4:4,-4:4};
            Stat.LatFormat= {'Hexagon',hexrad};
        case {2,4}
            Lat = loadCalibration(Options.LatMode,[yc,xc]);
            %Stat.LatFormat = {'Rectangular',-10:10,-10:10};
            Stat.LatFormat= {'Hexagon',hexrad};
            %Stat.LatFormat={'Ring',ringrad}
        case 8
            Lat = loadCalibration(Options.LatMode,[60,540]);
            Stat.LatFormat = {'Rectangular',-4:4,-8:8};
    end
    
    [Site,NumSite] = prepareSite(Stat.LatFormat);
    [~,YPixels,NumSave] = size(Data.Andor19330.Image);
    XSize=YPixels/numFK;
    % XSize = Data.SubImg(1,2)-Data.SubImg(1,1)+1; 
    
    % Get pre-calibration with averaged data
    [Lat,MeanBg,BgOffset,RFFT,FFTPeakFit,MeanSum] = getPreCalibration(Data,Bg,Lat);
    RFFT=[100 100];
    % if Options.ShowMeanSig
    %     figure(Units="normalized",OuterPosition=[0.1,0.1,0.6,0.8],Name="Mean signal")
    %     imagesc(MeanSum)
    %     checkCalibration(Lat,Site,acquireMarkerSize,'w')
    %     daspect([1 1 1])
    %     colorbar
    % end
    if Options.ShowFFTPeak
        showFFTPeakFit(FFTPeakFit)
    end

    % Get deconvolution pattern for background counts
    BgCSite = round((BgCenter-Lat.R)/Lat.V);
    [BgSite,NumBg] = prepareSite({'Rectangular',BgCSite(1)+(-BgR:BgR),BgCSite(2)+(-BgR:BgR)});
    BgDeconv = getDeconv(Lat,funcPSF,BgSite,XSize,YPixels);
%     assignin("base","BgDeconv",BgDeconv)

    % Perform analysis on each image
    % if Options.ShowUpdate
    %     FigMain = figure(Units="normalized", ...
    %         Name="Signal images",OuterPosition=[0,0,1,1]);
    % end
    LatNew = Lat;
    Stat.LatCount = zeros(NumSite,numFK,NumSave);
    Stat.LatOffset = nan(NumSave,2);
    Stat.BgCount = zeros(NumBg,numFK,NumSave);

    for i = 1:NumSave
        % Get background-subtracted image
        Signal = double(Data.Andor19330.Image(:,:,i))-MeanBg-BgOffset;
        % figure; imagesc(Signal); % this is to check
        % Get new offset calibration
        LatNew = getCalibration(Signal,RFFT,LatNew,'Full',numFK);
        Stat.LatOffset(i,:) = LatNew.R;
        
        % Get deconvolution pattern for the new offset
        Deconv = getDeconv(LatNew,funcPSF,Site,XSize,YPixels);
       
        % Get lattice site counts for each site
        Stat.LatCount(:,:,i) = getCount(Signal,Data.Andor19330.SubImg,Deconv);
        Stat.BgCount(:,:,i) = getCount(Signal,Data.Andor19330.SubImg,BgDeconv);
    end

end
function showUpdate(Signal,Count,Site,Lat)
    % Parameters
    XR = 100;
    YR = 100;

    XPixels = size(Signal,1);
    NumSubImg = size(Count,2);
    MarkerSize = acquireMarkerSize;

    XRange = round(Lat.R(1))+(-XR:XR);
    YRange = round(Lat.R(2))+(-YR:YR);

    switch NumSubImg
        case {1,2}
            DimX = NumSubImg;
            DimY = 1;
        case {4,8}
            DimX = NumSubImg/2;
            DimY = 2;
    end

    for i = 1:NumSubImg
        XRangeBox = XRange+(i-1)*(floor(XPixels/NumSubImg));
        LatNew = Lat;
        LatNew.R = Lat.R+[(i-1)*(floor(XPixels/NumSubImg)),0];

        % subplot(DimX,DimY,i)
        % imagesc(YRange,XRangeBox,Signal(XRangeBox,YRange))
        % checkCalibration(LatNew,Site,MarkerSize,'w')
        % title(sprintf('Image %d',i))
        % daspect([1 1 1])
        % colorbar
    end
end

%%  loadCalibration

function Lat = loadCalibration(Mode,varargin)
    switch Mode
        case 'Upper'
%             load('LatUpperCCD_20210614.mat','Lat')
            load("Lat_CCD19330_20230216_Atoms","Lat")
        case 'Lower'
            load('precalibration_lowerccd_20240522.mat','Lat')
    end
    
    if nargin>1
        Lat.R = varargin{1};
    end
end
%% prepareSite

function [Site,NumSite] = prepareSite(LatFormat)
    switch LatFormat{1}
        case 'Rectangular'
            % Format: {[X1,X2,...,XEnd]},{[Y1,Y2,...,YEnd]}
            [Y,X] = meshgrid(LatFormat{3},LatFormat{2});
            Site = [X(:),Y(:)];
            NumSite = length(LatFormat{2})*length(LatFormat{3});
        case 'MultipleRectangular'
            % Format: {[X11,...,X1End;X21,...,X2End;,,,]},{Y}
        case 'MaskedRectangular'
            % Format: {[X1,X2,...,XEnd]},{[Y1,Y2,...,YEnd]},{[1,0,1,...,1]}
        case 'Square'
            % WIP
        case 'Hexagon'
            rbig=LatFormat{2};
        
            for r = 1:rbig
                siteSum(r)=1;
                for i=1:r
                    siteSum(r) = siteSum(r)+6*i;
                end
                
                x=-r:r;
                [a, b]=meshgrid(x,x);
                sitetemp=[a(:),b(:)];
                
                % Generating hexagon coordinates in lattice space
                for i=0:(r-1)
                    for j=1:(r-i)
                        % disp([r-i, -(r-j+1)])
                        % disp([-(r-i), r-j+1])
                        if ismember([r-i, -(r-j+1)],sitetemp,"rows")==true || ismember([-(r-i), r-j+1],sitetemp,"rows")==true
                            [a,b]=ismember([r-i, -(r-j+1)],sitetemp,"rows");
                            sitetemp(b,:)=[];
                            [c,d]=ismember([-(r-i), r-j+1],sitetemp,"rows");
                            sitetemp(d,:)=[];
                            
                        end
                    end
                end
                Site{:,:,r}=sitetemp;
            end

            %Filling in rings
            ring{:,:,1} =Site{:,:,1};
            for i = 2:rbig
                % ring{:,:,i} = zeros(1,siteSum(i));
                outer = Site{:,:,i};
                inner = Site{:,:,i-1};
                ringtemp = outer(~ismember(outer, inner, 'rows'), :);
                ring{:,:,i}= ringtemp;
            end
            rings = ring;

            % Initialize an empty array to store the concatenated result
            concatenated_rings = [];
            
            % Loop through each page (3rd dimension) of the cell array
            for i = 1:numel(rings)
                concatenated_rings = [concatenated_rings; rings{i}];
            end
            Site = concatenated_rings;
            NumSite = numel(Site)/2;
    
    end
end

%% getPreCalibration

function [Lat,MeanBg,BgOffset,RFFT,FFTPeakFit,MeanSum] = getPreCalibration(Data,Bg,Lat)
    
    numFK = Data.Andor19330.Config.NumFrames;
    MeanBg = mean(Bg.Data.Andor19330.Image,3);
    MeanSig = mean(Data.Andor19330.Image,3)-MeanBg;
    BgOffset = cancelOffset(MeanSig,numFK);
    
    % Calculate RFFT based on the ROI center position
    YPixels = size(Bg.Data.Andor19330.Image,2);
    XSize=YPixels/numFK;
    Corner=[1,1;
        1,YPixels;
        XSize,1;
        XSize,YPixels];

    % RFFT = min(round(min(abs(Corner-Lat.R)))-12,[200 200]);
    RFFT = [100 100];
    if any(RFFT<0)
        error('Lattice center is too close to the image edge!')
    elseif any(RFFT<5)
        warning('Lattice center is too close to the image edge! RFFT = %d',min(RFFT))
    end

    % Print out old calibration results
    V1 = Lat.V(1,:);
    V2 = Lat.V(2,:);
    fprintf('\nOld lattice calibration:\n')
    fprintf('\tV1=(%4.2f, %4.2f),\t|V1|=%4.2fpx\n',V1(1),V1(2),norm(V1))
    fprintf('\tV2=(%4.2f, %4.2f),\t|V2|=%4.2fpx\n',V2(1),V2(2),norm(V2))
    fprintf('\tAngle=%4.2f deg\n\n',acosd(V1*V2'/(norm(V1)*norm(V2))))

    % Get calibration
    [Lat,FFTPeakFit,MeanSum] = getCalibration(Data.Andor19330.Image(:,:,1),RFFT,Lat,'Full',numFK);

    % Print out new calibration results
    printCalibration(Lat)

end
%% getCalibration

function [Lat,FFTPeakFit,SignalSum] = getCalibration(Signal,RFFT,Lat,Mode,varargin)
    % Parameters
    LatChangeThreshold = 0.002;
    CalBkgMin = 20;
    CalBkgMax = 1000;
    RFit = 7;
    % NumSubImg = 1;
    
    if nargin>4
        numFK = varargin{1};
        CalBkgMin = CalBkgMin*sqrt(numFK);
        CalBkgMax = CalBkgMax*numFK;
    end

    XSize = floor(size(Signal,1)/numFK);
    SignalSum = zeros(XSize,size(Signal,2));
    for i = 1:numFK
        SignalSum = SignalSum+double(Signal((i-1)*XSize+(1:XSize),:));
    end
    [SignalBox,FFTX,FFTY] = prepareBox(SignalSum,round(Lat.R),RFFT);

    if strcmp(Mode,'Full')
        PeakPosInit = (2*RFFT+1).*Lat.K+RFFT+1;
        [PeakPos,FFTPeakFit] = signalFFT(SignalBox,PeakPosInit,RFit);
        
        Lat.K = (PeakPos-RFFT-1)./(2*RFFT+1);
        Lat.V = (inv(Lat.K(1:2,:)))';

        VDis = vecnorm(Lat.V'-Lat.V')./vecnorm(Lat.V');
        if any(VDis>LatChangeThreshold)
            warning('off','backtrace')
            warning('Lattice vector length changed significantly by %.2f%%.',...
                100*(max(VDis)))
            warning('on','backtrace')
            A = [confint(FFTPeakFit{1}{1},0.95);...
                confint(FFTPeakFit{2}{1},0.95);...
                confint(FFTPeakFit{3}{1},0.95)];
            A = (A([2,4,6],[5,6])-A([1,3,5],[5,6]));
            A = A./(vecnorm(Lat.K')'*(2*RFFT+1));
            if all(A<LatChangeThreshold)
                save(sprintf('NewLatCal_%s.mat',datestr(now,'yyyymmdd')),'Lat')
                fprintf('\nNew lattice calibration saved\n')
            end
        end
    else
        FFTPeakFit = cell(0);
    end
    
    % Extract lattice center coordinates from phase at FFT peak
    [Y,X] = meshgrid(FFTY,FFTX);
    Phase = zeros(1,2);
    SignalModified = SignalBox;
    SignalModified(SignalBox<CalBkgMin | SignalBox>CalBkgMax) = 0;
    for i = 1:2
        PhaseMask = exp(-1i*2*pi*(Lat.K(i,1)*X+Lat.K(i,2)*Y));
        Phase(i) = angle(sum(PhaseMask.*SignalModified,'all'));
    end
    Lat.R = (round(Lat.R*Lat.K(1:2,:)'+Phase/(2*pi))-1/(2*pi)*Phase)*Lat.V;
end

function [PeakPos,FFTPeakFit] = signalFFT(Data,PeakPosInit,RFit)
       
    PeakPos = PeakPosInit;
    FFTPeakFit = cell(1,3);
    
    FFT = abs(fftshift(fft2(Data)));    
    for i = 1:3       
        XC = round(PeakPosInit(i,1));
        YC = round(PeakPosInit(i,2));
        PeakX = XC+(-RFit(1):RFit(1));
        PeakY = YC+(-RFit(end):RFit(end));
        PeakData = FFT(PeakX,PeakY);
        
        % Fitting FFT peaks
        [PeakFit,GOF,X,Y,Z] = fit2DGaussian(PeakData,-RFit(1):RFit(1),-RFit(end):RFit(end)); 
        FFTPeakFit{i} = {PeakFit,[X,Y],Z,GOF};

        if GOF.rsquare<0.5
            PeakPos = PeakPosInit;
            warning('off','backtrace')
            warning('FFT peak fit might be off (rsquare=%.2f), not updating.',...
                GOF.rsquare)
            warning('on','backtrace')
            return
        else
            PeakPos(i,:) = [PeakFit.u0,PeakFit.v0]+[XC,YC];
        end
    end
end

%% printCalibration

function printCalibration(Lat)
    V1 = Lat.V(1,:);
    V2 = Lat.V(2,:);
    V3 = V1+V2;
    fprintf('\nNew lattice calibration:\n')
    fprintf('\tV1=(%4.2f, %4.2f),\t|V1|=%4.2fpx\n',V1(1),V1(2),norm(V1))
    fprintf('\tV2=(%4.2f, %4.2f),\t|V2|=%4.2fpx\n',V2(1),V2(2),norm(V2))
    fprintf('\tV3=(%4.2f, %4.2f),\t|V3|=%4.2fpx\n',V3(1),V3(2),norm(V3))
    fprintf('\tAngle<V1,V2>=%4.2f deg\n',acosd(V1*V2'/(norm(V1)*norm(V2))))
    fprintf('\tAngle<V1,V3>=%4.2f deg\n\n',acosd(V1*V3'/(norm(V1)*norm(V3))))
end
%% prepareBox

function [DataBox,DataX,DataY] = prepareBox(Data,RC,R)
    DataX = RC(1)+(-R(1):R(1));
    DataY = RC(2)+(-R(end):R(end));
%     disp(DataX);
%     disp(RC(1));
%     disp(R);
    DataBox = Data(DataX,DataY);

end

%% fit2DGaussian

function [Fit,GOF,X,Y,Z] = fit2DGaussian(Signal,varargin)
    [XSize,YSize] = size(Signal);

    switch nargin
        case 3
            XRange = varargin{1};
            YRange = varargin{2};
            XSize = XRange(end)-XRange(1);
            YSize = YRange(end)-YRange(1);
        case 2
            XRange = varargin{1};
            YRange = varargin{1};
            XSize = XRange(end)-XRange(1);
            YSize = XSize;
        case 1
            XRange = 1:XSize;
            YRange = 1:YSize;
        otherwise
            error("Wrong number of inputs")
    end

    [Y,X,Z] = prepareSurfaceData(YRange,XRange,Signal);
    Max = max(Signal(:))+1;
    Min = min(Signal(:))-1;
    Diff = Max-Min;

    % Define 2D Gaussian fit type
    PeakFT = fittype('a*exp(-0.5*((u-u0)^2/b^2+(v-v0)^2/c^2))+d',...
                    'independent',{'u','v'},...
                    'coefficient',{'a','b','c','d','u0','v0'});
    PeakFO = fitoptions(PeakFT);

    PeakFO.Upper = [5*Diff,XSize,YSize,Max,XRange(end),YRange(end)];
    PeakFO.Lower = [0,0,0,Min,XRange(1),YRange(1)];
    PeakFO.StartPoint = [Diff,XSize/10,YSize/10,Min, ...
        (XRange(1)+XRange(end))/2,(YRange(1)+YRange(end))/2];
    PeakFO.Display = "off";

    [Fit,GOF] = fit([X,Y],Z,PeakFT,PeakFO);
end
%% getDeconv

function  [Deconv,DecPat] = getDeconv(Lat,funcPSF,Site,XLim,YLim)
    
    % Default deconvolution parameters
    PSFR = 10;
    RPattern = 30;
    Factor = 5;
    LatRLim = 2;
    RDeconv = 15;
    Threshold = 0.01;
    
    % Initialize
    NumSite = size(Site,1);
    Deconv = cell(NumSite,3);
    
    % Get the deconvolution pattern
    DecPat = matDeconv(Lat,funcPSF,PSFR,RPattern,Factor,LatRLim);
    %DecPat = ones(size(DecPat));
    %DecPat = DecPat/sum(DecPat,"all");
    
    % For each lattice site, find corresponding pixels and weights
    for i = 1:NumSite

        % Convert lattice X-Y index to CCD space coordinates
        Center = Site(i,:)*Lat.V+Lat.R;
        
        % Find X and Y range of pixels
        CX = round(Center(1));
        CY = round(Center(2));

        XMin = max(CX-RDeconv,1);
        XMax = min(CX+RDeconv,XLim);
        YMin = max(CY-RDeconv,1);
        YMax = min(CY+RDeconv,YLim);
        
        % Generate a list of pixel coordinates
        [Y,X] = meshgrid(YMin:YMax,XMin:XMax);
        XList = X(:);
        YList = Y(:);
        
        % Find the distance to site center
        XDis = XList-Center(1);
        YDis = YList-Center(2);
       
        % Find the cooresponding index in the deconvolution pattern
        XDeconv = round(Factor*(XDis+RPattern))+1;
        YDeconv = round(Factor*(YDis+RPattern))+1;
        XYDeconv = XDeconv+(YDeconv-1)*(2*Factor*RPattern+1);
        
        % Assign weights
        WeightList = DecPat(XYDeconv);
        Index = abs(WeightList)>Threshold;     
        Deconv{i,1} = XList(Index);
        Deconv{i,2} = YList(Index);
        Deconv{i,3} = WeightList(Index);
    end
    
end

function Pattern = matDeconv(Lat,funcPSF,PSFR,RPattern,Factor,LatRLim)
% Calculate deconvolution pattern by inverting (-LatRLim:LatRLim) PSF

    % Number of sites
    NumSite = (2*LatRLim+1)^2;

    % Number of pixels
    NumPx = (2*Factor*RPattern+1)^2;

    M = zeros(NumSite,NumPx);
    if funcPSF(PSFR,PSFR)>0.001
        warning(['Probability density at edge is significant = %.4f\nCheck' ...
            ' PSFR (radius for calculating PSF spread)'],funcPSF(PSFR,PSFR))
    end
    
    % For each lattice site, find its spread into nearby pixels
    for i = -LatRLim:LatRLim
        for j = -LatRLim:LatRLim
            
            % Site index
            Site = (i+LatRLim+1)+(j+LatRLim)*(2*LatRLim+1);

            % Lattice site coordinate
            Center = [i,j]*Lat.V;

            % Convert coordinate to magnified pixel index
            CXIndex = round(Factor*(Center(1)+RPattern))+1;
            CYIndex = round(Factor*(Center(2)+RPattern))+1;

            % Range of pixel index to run through
            xMin = CXIndex-PSFR*Factor;
            xMax = CXIndex+PSFR*Factor;
            yMin = CYIndex-PSFR*Factor;
            yMax = CYIndex+PSFR*Factor;

            % Go through all pixels and assign the spread
            x = xMin:xMax;
            y = yMin:yMax;
            Pixel = x'+(y-1)*(2*Factor*RPattern+1);
            [YP,XP] = meshgrid((y-1)/Factor-RPattern,(x-1)/Factor-RPattern);
            M(Site,Pixel) = funcPSF(XP(:)-Center(1),YP(:)-Center(2))/Factor^2;
            
        end
    end

    % Convert transfer matrix to deconvolution pattern
    MInv = (M*M')\M;
    Pattern = reshape(MInv(round(NumSite/2),:),sqrt(NumPx),[]);

    % Re-normalize deconvolution pattern
    Area = abs(det(Lat.V));
    %disp(Area);
    Pattern = Area/(sum(Pattern,"all")/Factor^2)*Pattern;
    %Pattern= Pattern/sum(Pattern,"all");
end

%% cancelOffset

function [BgOffset,STD] = cancelOffset(Signal,numFK,Name)
    arguments
        Signal
        numFK
        Name.YBgSize (1,1) double = 100
    end

    % NumSubImg = size(SubImg,1);
    [XPixels,YPixels] = size(Signal);
    XSize=YPixels/numFK;

    BgOffset = zeros(XPixels,YPixels);
    STD = zeros(numFK,2);

    if YPixels<2*Name.YBgSize+200
        warning('Not enough space to cancel background offset!')
        return
    end

    YRange1 = 1:Name.YBgSize;
    YRange2 = YPixels+((1-Name.YBgSize):0);
    for i = 1:numFK
        % XRange = SubImg(i,1):SubImg(i,2);
        XRange = 1+(i-1)*XSize:XSize+(i-1)*XSize;
        BgBox1 = Signal(XRange,YRange1);
        BgBox2 = Signal(XRange,YRange2);
        [XOut1,YOut1,ZOut1] = prepareSurfaceData(XRange,YRange1',BgBox1');
        [XOut2,YOut2,ZOut2] = prepareSurfaceData(XRange,YRange2',BgBox2');
        XOut = [XOut1;XOut2];
        YOut = [YOut1;YOut2];
        ZOut = [ZOut1;ZOut2];
        try
        XYFit = fit([XOut,YOut],ZOut,'poly11');
        catch
            assignin("caller","XOut",XOut)
            assignin("caller","YOut",YOut)
            assignin("caller","ZOut",ZOut)
            assignin("caller","XRange",XRange)
            assignin("caller","YRange1",YRange1)
            assignin("caller","YRange2",YRange2)
            assignin("caller","BgBox1",BgBox1)
            assignin("caller","BgBox2",BgBox2)
        end
        
        % Background offset canceling with fitted plane
        BgOffset(XRange,:) = XYFit.p00+XYFit.p10*XRange'+XYFit.p01*(1:YPixels);

        BgBoxNew1 = BgBox1-BgOffset(XRange,YRange1);
        BgBoxNew2 = BgBox2-BgOffset(XRange,YRange2);
        STD(i,:) = [std(BgBoxNew1(:)),std(BgBoxNew2(:))];
    end
    
    warning('off','backtrace')
    if any(BgOffset>2)
        warning('Noticable background offset: %4.2f',max(BgOffset(:)))
    end
    if any(STD>6)
        %warning('Noticable background noise: %4.2f',max(STD))
    end
    warning('on','backtrace')
end

%% getCount
function Count = getCount(Signal,XStart,Deconv)
    XPixels = size(Signal,1);
    NumSite = size(Deconv,1);
    NumSubImg = size(XStart,1);
    Count = zeros(NumSite,NumSubImg);

    for i = 1:NumSite
        List = Deconv{i,1}+XPixels*(Deconv{i,2}-1);
        if size(List,1)>0
            for j = 1:NumSubImg
                Count(i,j) = Deconv{i,3}'*Signal(List+XStart(j)-1); %only
                %line in original code
%                 Count(i,j) = ones(size(Deconv{i,3}'))*Signal(List+XStart(j)-1);
            end
        end
    end
end

