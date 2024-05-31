function [index_vals,velocity,meanvel,stdvel,badvals] = LSPIV_Parallel_Multi_Box( LineScanTable,numCores,pixel_Size,windowSize,Vessel_Type)
% LSPIV with parallel processing enabled.
%
% For additional information, please see corresponding manuscript:
%
% 'Line-Scanning Particle Image Velocimetry: an Optical Approach for
% Quantifying a Wide Range of Blood Flow Speeds in Live Animals'
% by Tyson N. Kim, Patrick W. Goodwill, Yeni Chen, Steven M. Conolly, Chris
% B. Schaffer, Dorian Liepmann, Rong A. Wang
% 
% PWG 3/28/2012


%numWorkers    = 4;  % number of workers on this machine.  
                     % depends on number of processors in your machine
                     % A safe starting point is typically 4,
                     % MATLAB supports up to 12 local workers in 
                     % R2011b.
                     %
                     % If you have trouble, you can access matlabpool
                     % directly: e.g. try typing "matlabpool 12" for 12
                     % workers.

% Parameters to improve fits

warning('off','all')
warning

maxGaussWidth = 100;  % maximum width of peak during peak fitting

% Judge correctness of fit
numstd        = 3;  %num of stdard deviation from the mean before flagging

pixel_size = pixel_Size;
numWorkers = numCores;
windowsize = windowSize; %in # scans, this will be converted to velocity points
                      %if one scan is 1/2600 s, then windowsize=2600 means
                      %a 1 second moving window.  Choose the window size
                      %according to experiment.

%%  settings
% Ask user for setting 
speedSetting = Vessel_Type;

if speedSetting == 1   % CAPILLARY SETTING
    numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 5;
elseif speedSetting == 2   % ARTERY SETTING
    numavgs       = 100;  %up to 100 (or more) for noisy or slow data
    skipamt       = 25;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 1;
elseif speedSetting == 3   % USER SETTING
    %disp('settings are hard coded in the script, see script!');
    numavgs       = 100;  %up to 200 (or more) for troublesome data. However
                          %you will lose some of the info in the peaks and
                          %troughs
    skipamt       = 3;   %if it is 2, it skips every other point.  3 = skips 2/3rds of points, etc.
    shiftamt      = 1;
end


%% Import the data from a multi-frame tif and make into a single array
%  The goal is a file format that is one single array,
%  so modify this section to accomodate your raw data format.
%
%  This particular file format assumes

imageLines = LineScanTable;
%%

%change 1 to 2
%% choose where in the image to process
%figure(2)
%imagesc(imageLines(1:size(imageLines,1),:))
%colormap('gray')
%set(gcf,'position',[10,10,size(imageLines,2),size(imageLines,1)])

%title('Select the boundaries of the region of interest 1/2');
%[X1,Y1] = ginput(1);
%line([X1 X1],[1 size(imageLines,2)]);

%title('Select the boundaries of the region of interest 2/2');
%[X2,Y2] = ginput(1);
%line([X2 X2],[1 size(imageLines,2)]);
%refresh
%pause(.01);

%startColumn   = round(min(X1, X2));      % Defines what part of the image we perform LSPIV on.
%endColumn     = round(max(X1, X2));

startColumn = 1;
endColumn = floor(size(imageLines,2));

%% startup parallel processing
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local',numWorkers)
else
    disp('Matlabpool Already Detected');
end

tic

%% minus out background signal (PWG 6/4/2009)
disp('DC correction')
DCoffset = sum(imageLines,1) / size(imageLines,1);

test = repmat(DCoffset,size(imageLines,1),1);

imageLinesDC = double(imageLines) - test;

%% do LSPIV correlation
disp('LSPIV begin');

scene_fft  = fft(imageLinesDC(1:end-shiftamt,:),[],2);
test_img   = zeros(size(scene_fft));
test_img(:,startColumn:endColumn)   = imageLinesDC(shiftamt+1:end, startColumn:endColumn);
test_fft   = fft(test_img,[],2);
W = 1./sqrt(abs(scene_fft)) ./ sqrt(abs(test_fft)); % phase only

LSPIVresultFFT = scene_fft .* conj(test_fft) .* W; 
LSPIVresult = ifft(LSPIVresultFFT,[],2);
disp('LSPIV complete');

toc

%% find shift amounts
disp('Find the peaks');
velocity = [];
maxpxlshift = round(size(imageLines,2)/2)-1;


index_vals = skipamt:skipamt:(size(LSPIVresult,1) - numavgs);
numpixels = size(LSPIVresult,2);
velocity  = nan(size(index_vals));
amps      = nan(size(index_vals));
sigmas    = nan(size(index_vals));
goodness  = nan(size(index_vals));


%% iterate through
sprintf('Length: %.15g',length(index_vals));

parfor index = 1:length(index_vals)
    %quick
    if mod(index_vals(index),100) == 0
        fprintf('line: %d\n',index_vals(index))
    end
    
    LSPIVresult_AVG   = fftshift(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1)) ...
                                      / max(sum(LSPIVresult(index_vals(index):index_vals(index)+numavgs,:),1));
    
    % find a good guess for the center
    c = zeros(1, numpixels);
    c(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift) = ...
        LSPIVresult_AVG(numpixels/2-maxpxlshift:numpixels/2+maxpxlshift);
    [maxval, maxindex] = max(c);
    
    % fit a guassian to the xcorrelation to get a subpixel shift
    options = fitoptions('gauss1');
    options.Lower      = [0    numpixels/2-maxpxlshift   0            0];
    options.Upper      = [1e9  numpixels/2+maxpxlshift  maxGaussWidth 1];
    options.StartPoint = [1 maxindex 10 .1];
    [q,good] = fit((1:length(LSPIVresult_AVG))',LSPIVresult_AVG','a1*exp(-((x-b1)/c1)^2) + d1',options);
    
    %save the data
    velocity(index)  = (q.b1 - size(LSPIVresult,2)/2 - 1)/shiftamt;
    amps(index)      = q.a1;
    sigmas(index)    = q.c1;
    goodness(index)  = good.rsquare;
end
%% find possible bad fits
toc

% Find bad velocity points using a moving window 
pixel_windowsize = round(windowsize / skipamt);

badpixels = zeros(size(velocity));
for index = 1:1:length(velocity)-pixel_windowsize
    pmean = mean(velocity(index:index+pixel_windowsize-1)); %partial window mean
    pstd  = std(velocity(index:index+pixel_windowsize-1));  %partial std 
    
    pbadpts = find((velocity(index:index+pixel_windowsize-1) > pmean + pstd*numstd) | ...
                   (velocity(index:index+pixel_windowsize-1) < pmean - pstd*numstd));

    badpixels(index+pbadpts-1) = badpixels(index+pbadpts-1) + 1; %running sum of bad pts
end
badvals  = find(badpixels > 0); % turn pixels into indicies
goodvals = find(badpixels == 0);

meanvel  = abs(mean(velocity(goodvals))*pixel_size*windowsize); %overall mean
stdvel   = abs(std(velocity(goodvals))*pixel_size*windowsize);  %overall std

Switch_to_mm = false;
if abs(meanvel) > 1000
   meanvel = meanvel/1000;
   stdvel = stdvel/1000;
   Switch_to_mm = true;
end

end

