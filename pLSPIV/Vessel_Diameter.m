%% Calculate Blood Vessel Width Based on Full-Width at Half-Maximum method
%Input: Single 8-bit tif image
%TO USE: Draw straight line perpendicular to blood vessel
%User drawn lines saved to image
%Calculates vessel width based off user defined pixel_size
%Output: Saves vessel width to excel file

%Made by Christian Shizuo Burns 2019
%Get video
clear all
clc
close all

%turn off warnings
warning('off','all')
warning

%Set pixel_size as um/pixel
pixel_size = 0.8;
windowsize = 30;

%Select image to analyze
disp('import raw data');
[fname,pathname]=uigetfile('*.TIF','pick a linescan file');%loads file
if fname == 0; beep; disp('Cancelled'); end
fullFileName = fullfile(pathname,fname);

%Create output excel file
newExcelFile = strcat(fullFileName(1:end-4),'_PlotData','.xlsx');

%% tif image stackuc 
InfoImage=imfinfo(fullFileName);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
imageBitDepth=InfoImage(1).BitDepth;
if imageBitDepth == 8; FinalImage=zeros(nImage,mImage,NumberImages,'uint8'); end
if imageBitDepth == 16; FinalImage=zeros(nImage,mImage,NumberImages,'uint16'); end

TifLink = Tiff(fullFileName, 'r');
FinalImage=TifLink.read();
TifLink.close();

snapShot = FinalImage;

%Display input image
imshow(snapShot);

drawnLine = 'YES'; 

%% Get straight line for vessel width
linefinished = 'NO';
linei = 1;
                        
if strcmpi(drawnLine,'Yes')
   disp('Select Vessel Width');
   imshow(snapShot);
   while strcmpi(linefinished,'NO')
        linehFH(linei) = drawline;
        while linehFH(linei).Selected == 0
             pause;
        end
        linefinished = questdlg('Finished?', ...
            'confirmation', ...
            'YES', 'NO', 'UNDO', 'NO');
        if strcmpi(linefinished, 'UNDO')
            %closes window
            delete(linehFH(linei))
            linefinished = 'NO';
        else
            linexy{linei} = linehFH(linei).Position;
            lineInitialPosition{linei,1} = [floor(linexy{linei}(1,1)) floor(linexy{linei}(1,2))];  
            linelabellingNumber{linei} = int2str(linei);
            linei = linei+1;
        end
   end
end

%% Display
TempOutPutImage = snapShot;

%label vessel width
if strcmpi(drawnLine,'YES')
    TempOutPutImage = insertText(TempOutPutImage,cell2mat(lineInitialPosition),linelabellingNumber,'FontSize',10,'BoxColor','blue');
end

%shuffle line
if strcmpi(drawnLine,'YES')
   for q =1:linei-1
       lineFinalPosition{q}(1) = linexy{q}(1,1);
       lineFinalPosition{q}(2) = linexy{q}(1,2);
       lineFinalPosition{q}(3) = linexy{q}(2,1);
       lineFinalPosition{q}(4) = linexy{q}(2,2);
   end
end

%add line width
if strcmpi(drawnLine,'YES')
   for q=1:linei-1
      TempOutPutImage = insertShape(TempOutPutImage,'Line',lineFinalPosition{q},'Color','blue'); 
   end
end

%% Show all positions
imshow(TempOutPutImage);

fname1 = erase(fullFileName, '.tif');
snapshotName = strcat(fname1, '_lines.png');

imwrite(TempOutPutImage,snapshotName,'png');

%Width Line
VesWidth_X = {};
VesWidth_Y = {};

if strcmpi(drawnLine,'YES')
   linenumberSlices = size(linexy); 
   for j=1:linenumberSlices
       VesWidth_X{j} = 0;
       VesWidth_Y{j} = 0; 
   end
end

%% Get improfile line
%improfile of vessel width
if strcmpi(drawnLine,'YES')
    for j=1:linenumberSlices(2)
        textprogressbar('VesselScan: ');
        for x = 1:NumberImages
            VesselLineScan{j}(x,:) = improfile(FinalImage(:,:,x),linexy{j}(:,1),linexy{j}(:,2));
        
            textprogressbar((x/NumberImages)*100);
        end
        
        VesWidth_X{j} = linexy{j}(:,1);
        VesWidth_Y{j} = linexy{j}(:,2);
        textprogressbar('end');
    end
end

%% Plot values
%Calculate vessel width
if strcmpi(drawnLine,'YES')
    for VesselNumber = 1:linenumberSlices(2)
        
        for lineNumber = 1:NumberImages
            linelength = size(VesselLineScan{VesselNumber}(1,:));
            linelength2 = linspace(1,linelength(2),linelength(2));
            x = linelength2(1,:);
            y = VesselLineScan{VesselNumber}(lineNumber,:);
            VesselLineWidth = fwhm(x,y); 
        end
        
        %save excel data
        excel_array = zeros(2,1);

        excel_array(1,1:end) = VesselNumber;
        excel_array(2,1:end) = VesselLineWidth;
        
        writematrix(excel_array,newExcelFile,'Sheet',VesselNumber);
        
        fprintf('\nMean Vessel Width %0.2f\n',nanmean(VesselLineWidth)*pixel_size);
        
    end
end
