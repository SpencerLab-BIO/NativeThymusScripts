%--------------------------------------------------------------------------
%Pseudo line-scanning particle image velocimetry (p3LSPIV)
%--------------------------------------------------------------------------
%Used to calculate blood flow in individual blood vessels where the blood
%lumen has been labelled with a fluorescent dye which provides negative
%contrast against the RBC moving through the vessel.
%Instead of taking a line scan as an input, a pseudo line scan is generated
%by getting the pixel intensities of a line drawn across a video. A pseudo
%line scan is a 2D image where x is the width of the scan and y is the
%frames of the video.
%--------------------------------------------------------------------------
%User can select from one of 4 scans: Straight line scan, Handdrawn line
% %scan, box scan and vessel width scan.

%Straight line scan - A user specifies a straight line roi which is used to
%generate a pseudo line scan.

%Handdrawn line scan - user draws a line roi of interest. Used to follow
%blood flow that turns.

%Box scan - performs a series of parallel line scans across a box roi
%interest specified by the user. Allows for the multiple parallel line
%scans to be performed in a single vessel.

%Vessel Width line scan - straight line scan drawn perpendicular to the
%blood vessel. Used to calculate vessel width via FWHM methdod.

%Made by Christian Shizuo Burns 2019
%Edited by Iqbal Bhandal

%Developed on matlab 19a
%All matlab additions included

%Modified from script found here:
%Kim, T. N., Goodwill, P. W., Chen, Y., Conolly, S. M., Schaffer, C. B., Liepmann, D., & Wang, R. A. (2012). 
%Line-scanning particle image velocimetry: an optical approach for quantifying a wide range of blood flow speeds in live animals. 
%PloS one, 7(6), e38590. https://doi.org/10.1371/journal.pone.0038590
%Enter at your own risk

%Clear workspace
clear all
clc
close all

%turn off warnings
%CSB: I don't know what's causing them
warning('off','all')
warning

%% Initialize Variables

%IDK what this was for
%I think it waas suppose to prevent redudant user input about these
%variables, but it doesn't seem to be used
fpscheck = false;

%um per pixel
%Old
%pixel_size = 0.3672;
%As of January 2020
pixel_size = 0.31;

%Number of cores which will be used for parallel processing
numWorkers = 8;

%Frame rate
windowsize = 30;

%number of points scanned in rectangle
num_points = 10;

%% Active Imagej if live video is selected

%Download MIJ to run Imagej from matlab
%https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab
%Add MIJJ for java path
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar';
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar';

%Ask user if video should be used
%ImgJ = questdlg('Video?',...
%               'confirmation',...
%                'YES','NO','NO');
ImgJ = 'No';
%Start imagej
if strcmpi(ImgJ,'YES'); MIJ.start; end

%Ask user to select video tif file
disp('import raw data');
[fname,pathname]=uigetfile('*.TIF','pick a linescan file');%loads file
if fname == 0; beep; disp('Cancelled'); end

%File name and location
fullFileName = fullfile(pathname,fname);

%Display user selected video
if strcmpi(ImgJ,'YES'); ij.IJ.open(fullFileName); end

%Create empty excel file to be used to store data related to calculated
%blood flow and vessel width
newExcelFile = strcat(fullFileName(1:end-4),'_PlotData','.xlsx');

%% Load Tif video stack into memory

%Get information about tif video
%mImage = width of image
%nImage = height of image
%NumberImages = number of frames
%imageBitDepth = bit depth 
%Bit depth = 8bit or 16bit
InfoImage=imfinfo(fullFileName);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
imageBitDepth=InfoImage(1).BitDepth;

%Create empty arrays equal to tif video 
if imageBitDepth == 8; FinalImage=zeros(nImage,mImage,NumberImages,'uint8'); end
if imageBitDepth == 16; FinalImage=zeros(nImage,mImage,NumberImages,'uint16'); end

%Load tif into memory
%FinalImage = loaded tif video
TifLink = Tiff(fullFileName, 'r');
textprogressbar('Image Download: ');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
   
   textprogressbar((i/NumberImages)*100);
end

TifLink.close();
textprogressbar('end');

%Display first frame from video
snapShot = FinalImage(:,:,1);
imshow(snapShot);

%Set variables
%finished -> used to determine when handdrawn line aquisition is complete
finished = 'NO';
%Initialize index
i = 1;
%Create list of color options
col = ['g' 'b' 'r' 'y' 'm' 'c'];

%% Hemodynamics calculation option
%Calculate shear rate
%If RBC density is known, can calculate additional parameters
%Check this paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5318670/

%Ask if hemodynamics will be calculated
HemoDyna = questdlg('Calculate Hemodynamics?',...
                            'confirmation',...
                            'YES','NO','NO');
if strcmpi(HemoDyna,'YES')                     
    HemoMsg = msgbox("To calculate the hemodynamics the psuedo line and Vessel Width line need to chronologically be the same");
end

%% Get Line Options
%   -Straight Line - User specifies a straight line
%   -Handdrawn Line - User draws line correspond to pseudo line scan
%   -Box Scan - User specifies a box area - parallel straight line drawn
%   from the top of the box to the bottom will be used as straight lines
%   -Blood Vessel Width - User specifies straight line perpendicular to
%   blood vessel. Calculates width via FWHM method
%   -Rolling Lymphocyte - smooths intensity plot to prevent interference
%   from background RBC

%Get users scan type
%Use Ctrl to select multiple options
list = {'Straight','Handdrawn','Box','Width'};
[indx,tf] = listdlg('ListString',list);
LineOptList = [1 2 3 4];

for t=1:length(indx)
    LineOptSelect = ismember(LineOptList,indx);
end

%Set variables for use of straight, handdrawn, box and width scans 
SLine = false;
HLine = false;
drawnRectangle = 'NO';
drawnLine = 'NO';


if LineOptSelect(1)
   SLine = true; 
end
if LineOptSelect(2)
   HLine = true; 
end
if LineOptSelect(3)
    drawnRectangle = 'YES'; 
end
if LineOptSelect(4)
   drawnLine = 'YES'; 
end

%% Get Straight Line Positions
%Get handdrawn positions from user
slinefinished = 'NO';
slinei = 1;

%Ask user to draw straight line on snapshot(first frame of video)
%Line can be adjusted after initially drawnn.
%Click on line and press enter once done.
%Prompt will ask if user is finished
%Yes - precede to next scan type or to plotting data
%No - ask user to repeat line scan selection
%Undo - redo previous line

%Output:
%slinexy - positional daata of line
%slineInitialPosition - first point of line - used to place label
%slinelabellingNumber - corresponds to order of line
%slinei = current index in array 

if SLine
   disp('Select Straight Lines');
   imshow(snapShot);
   while strcmpi(slinefinished,'NO')
        slinehFH(slinei) = drawline;
        while slinehFH(slinei).Selected == 0
             pause;
        end
        slinefinished = questdlg('Finished?', ...
            'confirmation', ...
            'YES', 'NO', 'UNDO', 'NO');
        if strcmpi(slinefinished, 'UNDO')
            %closes window
            delete(slinehFH(slinei))
            slinefinished = 'NO';
        else
            slinexy{slinei} = slinehFH(slinei).Position;
            slineInitialPosition{slinei,1} = [floor(slinexy{slinei}(1,1)) floor(slinexy{slinei}(1,2))];  
            slinelabellingNumber{slinei} = int2str(slinei);
            slinei = slinei+1;
        end
   end
end

%% Get freehand handdrawn positions
%Get freehald position data
%Same as Straight line but line can be curved
%Matlab records freehand data as a series of lines of specific length

%Output:
%xy = positional data of handdrawn line
%FinalPostion = position at the end of the line used for labelling
%labellingNumber = number used to label position
%i = current indext

while strcmpi(finished,'NO') && HLine
  disp('Select Handdrawn Lines');
  imshow(snapShot);
  hFH(i) = drawfreehand('Closed',false);
  
  FinalColor{i} = col(mod(i,6)+1);
  %setColor(hFH(i),FinalColor{i});
  
  finished = questdlg('Finished?', ...
      'confirmation', ...
      'YES', 'NO', 'UNDO', 'NO');
  if strcmpi(finished, 'UNDO')
      %closes window
      delete(hFH(i))
      finished = 'NO';
  else
      xy{i}  = hFH(i).Position;
      
      FinalPosition{i,1} = [floor(xy{i}(end,1)) xy{i}(end,2)];
      labellingNumber{i} = int2str(i); 
      
      i = i + 1;
  end
end

%% Get drawn rectangle
%Get user input for box scan.
%User specifies a box area to be examinged. Box size and orientation can be
%changed after box is drawn.
%A series of imaginary parallel lines are drawn from the top of the box to
%the bottom and used to generate line scans.

%Variables related to user input in box scan
rectanglefinished = 'NO';
reci = 1;
             
%Ask user to draw a box. Once drawn the box can be moved around, resized
%and rotated.

%Output:
%rectanglexy = position of the box saved as coordinate positions. Check matlab
%website
%recVerticesTest = Position of verticies of box. Used to calculate edge
%recInitialPosition = position used to place label
%reclabellingNumber = number label that will be used
%reci = array position

if strcmpi(drawnRectangle,'YES')
    disp('Select Box Scan')
    imshow(snapShot);
    while strcmpi(rectanglefinished,'NO')
        rectanglehFH(reci) = drawrectangle('Rotatable',true);
  
        while rectanglehFH(reci).Selected == 0
             pause;
        end
        rectanglefinished = questdlg('Finished?', ...
            'confirmation', ...
            'YES', 'NO', 'UNDO', 'NO');
        if strcmpi(rectanglefinished, 'UNDO')
            %closes window
            delete(rectanglehFH(reci))
            rectanglefinished = 'NO';
        else
            rectanglexy{reci}  = rectanglehFH(reci).Position;
      
            recVerticesTest{reci} = rectanglehFH(reci).Vertices;
            
            recInitialPosition{reci,1} = [floor(recVerticesTest{reci}(1,1)) floor(recVerticesTest{reci}(1,2))];
            reclabellingNumber{reci} = int2str(reci); 
      
            reci = reci + 1;
        end
    end
end

%% Get straight line for vessel width
%Get position data from straight line used to calculate vessel width
%Gets data in similiar manner to previous straight line

linefinished = 'NO';
linei = 1;

%drawnLine = questdlg('Draw Vessel Width?',...
%                            'confirmation',...
%                            'YES','NO','NO');
                        
%Output:
%linexy = position of line
%lineInitialPosition = position of label for line
%lineLabellingNumber = label for line
%linei = line index value

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

%% Get Straight Line Position and label number on snapshot (first frame)
TempOutPutImage = snapShot;

if SLine
    TempOutPutImage = insertText(TempOutPutImage,cell2mat(slineInitialPosition),slinelabellingNumber,'FontSize',10,'BoxColor','blue');
end

%% Get position of free hand line and label number
%Overlays number label for each line on the image
%Overlays line position for each line onto the image


%Add freehand line
if HLine
    
    %Insert text over iamge
    TempOutPutImage = insertText(TempOutPutImage,cell2mat(FinalPosition),labellingNumber,'FontSize',12,'BoxColor','magenta');
    
    %Insert Line over image
    %shuffle together
    %Convert 2 column x,y data to single column xy
    for g = 1:i-1
        for u = 1:size(xy{g}(:,1))
            TempLinePosition{g}(2*u-1) = xy{g}(u,1);
            TempLinePosition{g}(2*u) = xy{g}(u,2);
        end
    end
    
    %add free hand lines
    for q = 1:i-1
        TempOutPutImage = insertShape(TempOutPutImage,'Line',TempLinePosition{q}(1,:),'LineWidth',5,'Color','magenta');
    end
end

%% get rectangle position and label number
if strcmpi(drawnRectangle,'YES')
    TempOutPutImage = insertText(TempOutPutImage,cell2mat(recInitialPosition),reclabellingNumber,'FontSize',10,'BoxColor','red');
end

%% label vessel width
if strcmpi(drawnLine,'YES')
    TempOutPutImage = insertText(TempOutPutImage,cell2mat(lineInitialPosition),linelabellingNumber,'FontSize',10,'BoxColor','blue');
end

%% Get position of line and add to the snapshot(first frame)
%Add straight line outline, box outline, and vessel width outline 
%Positions of straight line, rectangle and vessel line need to be formated
%to be added to image

%shuffle Straight line
if SLine
   for q =1:slinei-1
       slineFinalPosition{q}(1) = slinexy{q}(1,1);
       slineFinalPosition{q}(2) = slinexy{q}(1,2);
       slineFinalPosition{q}(3) = slinexy{q}(2,1);
       slineFinalPosition{q}(4) = slinexy{q}(2,2);
   end
end

%shuffle rectangle together
if strcmpi(drawnRectangle,'YES')
    for q = 1:reci-1
        recFinalPosition{q}(1) = recVerticesTest{q}(1,1); 
        recFinalPosition{q}(2) = recVerticesTest{q}(1,2); 
        recFinalPosition{q}(3) = recVerticesTest{q}(2,1); 
        recFinalPosition{q}(4) = recVerticesTest{q}(2,2); 
        recFinalPosition{q}(5) = recVerticesTest{q}(3,1); 
        recFinalPosition{q}(6) = recVerticesTest{q}(3,2); 
        recFinalPosition{q}(7) = recVerticesTest{q}(4,1); 
        recFinalPosition{q}(8) = recVerticesTest{q}(4,2); 
        recFinalPosition{q}(9) = recVerticesTest{q}(1,1); 
        recFinalPosition{q}(10) = recVerticesTest{q}(1,2); 
        
        %top
        recx{q}(1,:) = linspace(recVerticesTest{q}(1,1),recVerticesTest{q}(4,1),num_points);
        recy{q}(1,:) = linspace(recVerticesTest{q}(1,2),recVerticesTest{q}(4,2),num_points);
        %bottom
        recx{q}(2,:) = linspace(recVerticesTest{q}(2,1),recVerticesTest{q}(3,1),num_points);
        recy{q}(2,:) =  linspace(recVerticesTest{q}(2,2),recVerticesTest{q}(3,2),num_points);
        
    end
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

%add straight line outline
if SLine
   for q=1:slinei-1
      TempOutPutImage = insertShape(TempOutPutImage,'Line',slineFinalPosition{q},'Color','blue'); 
   end
end

%add box outline
if strcmpi(drawnRectangle,'YES')
    for q = 1:reci-1
        TempOutPutImage = insertShape(TempOutPutImage,'Polygon',recFinalPosition{q},'Color','red');
    end
end

%add line width outline
if strcmpi(drawnLine,'YES')
   for q=1:linei-1
      TempOutPutImage = insertShape(TempOutPutImage,'Line',lineFinalPosition{q},'Color','blue'); 
   end
end

%% Show all positions
%Display image with labels and lines overlayed
%Saves image to folder with video + _lines

imshow(TempOutPutImage);

fname1 = erase(fullFileName, '.tif');
snapshotName = strcat(fname1, '_lines.png');

imwrite(TempOutPutImage,snapshotName,'png');

%% Initialize parameters related to velocity calculation
%Create frame works which will be filled with data related to intensity of
%pixels across the various line scans

snumberSlices = [0 0];
numberSlices = [0 0];
recnumberSlices = [0 0];

simagelines = 0;

%get number of positions 
%Straight Line
sVesWidth_X = {};
sVesWidth_Y = {};

%Handdrawn Line
if SLine
   slinenumberSlices = size(slinexy); 
   for j=1:slinenumberSlices
       sVesWidth_X{j} = 0;
       sVesWidth_Y{j} = 0; 
   end
end

%Handdrawn Line
if HLine
    numberSlices = size(xy);
end

%Box scan
if strcmpi(drawnRectangle,'YES')
    recnumberSlices = size(rectanglexy);
end

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
%Uses improfile function from matlab to get the line scan from the input
%video. Get the pixel values corresponding to the line scans specified by
%the user. Saves line scans to a 2D array with x being the width of the
%line and y being frames.

%Improfile records the pixel values of the data that is in a line directly 
%in between two points. Improfile starts from the opposite corners of the 
%positions specified by as the start and end of the line.

%Output:
%sVesselLineScan = 2D line scan corresponding the pixel values specified by
%the user
%sVesWidth_X = doesn't seemed to be used. Potentially to check if position
%is filled
%sVesWidth_Y = same as above

%LineScan1 = 2D line scan corresponding the pixel values specified by
%the user by the handdrawn line
%LineScan = modified version of handdrawn line scan which has removed 
%duplicate columns which can be generated by the handdrawn line scan

%recLineScan = 2D X N line scan corresponding the pixel values specified by
%the user by the box line. If multiplied by the number of parallel lines 
%scanned across box.

%VesselLineScan = 2D line scan corresponding to pixel values specified by
%the user by the vessel width line.


%improfile of Straight Vessel Scan
if SLine
    for j=1:slinenumberSlices(2)
        textprogressbar('StraightVesselScan: ');
        for x = 1:NumberImages
            stempVesselLineScan{j}(x,:) = improfile(FinalImage(:,:,x),slinexy{j}(:,1),slinexy{j}(:,2));
        
            textprogressbar((x/NumberImages)*100);
        end
        
        sVesselLineScan{j} = stempVesselLineScan{j};
        
        sVesWidth_X{j} = slinexy{j}(:,1);
        sVesWidth_Y{j} = slinexy{j}(:,2);
        textprogressbar('end');
    end
end

%improfile of VesselScan
if HLine
    for j=1:numberSlices(2)
        %get intensity
        textprogressbar('LineScan: ');
        for x = 1:NumberImages
            LineScan1{j}(x,:) = improfile(FinalImage(:,:,x),xy{j}(:,1),xy{j}(:,2));

            textprogressbar((x/NumberImages)*100);
        end
        LineX{j} = xy{j}(:,1);
        LineY{j} = xy{j}(:,2);

        textprogressbar('end');
    end

    %Compress Linescan Data
    for j=1:numberSlices(2)
        textprogressbar('LineScan Compression: ');
        realPos = 1; 

        for x = 1:(length(LineScan1{j}(1,:))-1)
            if LineScan1{j}(:,x) == LineScan1{j}(:,x+1)
            else
                LineScan{j}(:,realPos) = LineScan1{j}(:,x);
                realPos = realPos+1;
                textprogressbar(x/(length(LineScan1{j}(1,:))-1)*100);
            end
        end
        textprogressbar('end');
    end
end

%get improfile of rectangles
if strcmpi(drawnRectangle,'YES')
    %though number of boxes
    for j=1:recnumberSlices(2)
        %get positions of lines
        for k=1:num_points
            %get intensity
            textprogressbar('BoxScan: ');
            for x = 1:NumberImages
                recLineScan{j}{k}(x,:) = improfile(FinalImage(:,:,x),recx{j}(:,k),recy{j}(:,k));
                textprogressbar((x/NumberImages)*100);
            end
            textprogressbar('end');
        end
    end
end

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
%Data is handed to one of 3 functions to calculate velocity or vessel width
%LSPIV_Parallel_Multi = modified program used to calculate velocity from 2D
%line scan. Used with Straight and Handdrawn line scans. Assumes X = width
%and y = frames. Allows user to specify end boundaries of line scan which
%will be analyzed. Cannot scan last numavgs frames of the image due to 
%averaging required. Full analysis will scan to the end of the image by 
%reducing numavgs to not extend past the end of the line scan. Reduced 
%accuracy region marked in blue. Final data values is saved to excel sheet.

%LSPIV_Parallel_Multi_Box = Same as previous function but boundaries are 
%automatically set to the ends of the line scan. Used with the box scan.
%Does not have full analysis.

%FWHM = program used to calculate vessel with via The Full-Width at
%Half-Maximum method. Download location can be found here: https://www.mathworks.com/matlabcentral/fileexchange/10590-fwhm

%Absolute value of data is plotted.
%Region of the line scan used to calculate velocity drawn above graph of
%velocity over time.

%Output:
%imageLines = 2D line scan corresponding to the user selected region of
%analysis.
%startcolumn = position of one user specified boundary of line scan.
%endcolumn = position of end user specified boundary.
%index_vals = time positions corresponding to points were velocity was
%analyzed.
%windowsize = frame rate
%velocity = array of velocitys corresponding to index_vals time points
%pixel_size = size of pixels in um. um/pixel
%badvals = velocity measurements that are significantly outside expected values
%meanvel = mean of velocities
%stdev = standard deviation of velocities
%Switch_to_mm = converts the y axis from um/sec to mm/sec when velocity
%excedes 1000 um/sec
%numworkers = number of cores used for parallel processing. Redundant
%fpscheck = IDK. possible related to not inputing data multiple times, but
%not used otherwise.
%full_analysis = Enable full analysis of line scan data. Normally, the
%LSPIV scripts will not analyze the end of the data. Numavgs is used by the
%LSPIV scripts to average numavgs number of frames to reduce background
%noise when calculating velocity. To overcome this when full analysis is
%engaged, numavgs is reduced towards the end of the line scan to not exceed
%the end of the line scan. This reduces accuracy as numavgs is reduced.
%Uncertain areas are plotted in blue.
%index_length = number of time points to be graphed
%numavgs = number of frames which are averaged to reduced background when
%calculating velocity.
%skipmat = point spacing. Potentially divide windowsize to get better time frame - NOT CURRENTLY USED
%rolling = whether or not rolling lymphocyte was imaged. Used to turn on skipmat - NOT CURRENTLY
%USED

rolling = 0;

%Straight Line Graph
if SLine
     for simagelines = 1:slinenumberSlices(2)
        
         %get line scan
        [simageLines,sstartColumn,sendColumn,sindex_vals,swindowsize,svelocity,spixel_size,sbadvals, smeanvel,sstdvel,sSwitch_to_mm,snumWorkers, fpscheck,sfull_analysis,sindex_length,snumavgs,skipamt,rolling] = LSPIV_Parallel_Multi(sVesselLineScan{simagelines}, nImage, mImage, fpscheck, pixel_size, numWorkers, windowsize);
        
        %Detemine figure number which will not override snapshot
        sFinalFigureTitle = simagelines + 2;
        figure(sFinalFigureTitle);

        %IDK. Doesn't seem to be used
        srealVel(simagelines, 1:length(svelocity)) = svelocity(1,1:length(svelocity));

        %save excel data
        %Create empty data sheet and time and velocity
        ssize_index_vals = size(sindex_vals);

        excel_array = zeros(2,ssize_index_vals(2));

        excel_array(1,1:end) = (sindex_vals/swindowsize);
        excel_array(2,1:end ) = (svelocity*spixel_size*(swindowsize));
        
        excel_array(3,1) = sindex_length;

        %Save data to excel sheet
        writematrix(excel_array,newExcelFile,'Sheet',sFinalFigureTitle-2);

        %imagelines
        subplot(8,1,1:2);

        %%%% linescan of selected data
        %Rotates selected line scan 90 degrees to match the velocity
        %positions
        ssize1 = sendColumn - sstartColumn;
        if sfull_analysis
            ssize2 = NumberImages;
        else
            ssize2 = NumberImages-snumavgs;
        end
        simgtmp = zeros([ssize1+1 ssize2]); % to enable BW and color simultaneously
        simgtmp(:,:,1) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn));
        simgtmp(:,:,2) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn)); 
        simgtmp(:,:,3) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn));
        imagesc(simgtmp/max(max(max(simgtmp))))
    
        %Save Image
        fname1 = erase(fullFileName, '.tif');
        snapshotName = strcat(fname1, "_",int2str(simagelines),'_Straight_Line_scan.jpeg');
        
        soloScan = uint8(zeros([ssize2 ssize1+1]));
        
        soloScan(:,:,1) = simageLines(1:ssize2,sstartColumn:sendColumn);
        soloScan(:,:,2) = simageLines(1:ssize2,sstartColumn:sendColumn);
        soloScan(:,:,3) = simageLines(1:ssize2,sstartColumn:sendColumn);
        
        imwrite(soloScan,snapshotName,'jpeg');
        %End Save Image
        
        t = title(simagelines);

        ylabel('[pixels]');

        %velocity plot
        %Switches between displaying at mm and um scale depending on values
        %Full analysis enabled will display uncertain data in blud and
        %accurate data in green.
        %Solid line = mean of data points
        %Dotted lines = 4 * standard deveiation
        if sSwitch_to_mm == true
            %subplot(numberSlices(2)*2,1,(imagelines*2))
            subplot(8,1,4:5)
            
            if sfull_analysis
                plot(sindex_vals(1:sindex_length)/swindowsize, abs((svelocity(1:sindex_length)*spixel_size*(swindowsize))/1000),'g.');
                hold all
                plot(sindex_vals(sindex_length:end)/swindowsize, abs((svelocity(sindex_length:end)*spixel_size*(swindowsize))/1000),'b.');
                plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*(swindowsize)/1000), 'ro');
                disp('Full analysis');
            else
                plot(sindex_vals/swindowsize, abs((svelocity*spixel_size*(swindowsize))/1000),'.');
                hold all
                plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*(swindowsize)/1000), 'ro');
            end
            
            hold off
            xlim([sindex_vals(1)/swindowsize sindex_vals(end)/swindowsize]);
            supperbound = smeanvel+sstdvel*4;
            slowerbound = smeanvel-sstdvel*4;
            
            if slowerbound < 0
                slowerbound = 0;
            end
            
            ylim([slowerbound supperbound]);
            title('Fitted Pixel Displacement');
            ylabel({'Velocity'; '[mm/sec]'});
        else
            subplot(8,1,4:5)
            
           if sfull_analysis
                plot(sindex_vals(1:sindex_length)/swindowsize, abs((svelocity(1:sindex_length)*spixel_size*(swindowsize))),'g.');
                hold all
                plot(sindex_vals(sindex_length:end)/swindowsize, abs((svelocity(sindex_length:end)*spixel_size*(swindowsize))),'b.');
                plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*swindowsize), 'ro');
                disp('Full analysis');
           else
                plot(sindex_vals/swindowsize, abs((svelocity*spixel_size*(swindowsize))),'.');
                hold all
                plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*(swindowsize)), 'ro');
           end
            
            hold off
            xlim([sindex_vals(1)/swindowsize sindex_vals(end)/swindowsize]);
            supperbound = smeanvel+sstdvel*4;
            slowerbound = smeanvel-sstdvel*4;
            
            if slowerbound < 0
                slowerbound = 0;
            end
            
            ylim([slowerbound supperbound]);
            title('Fitted Pixel Displacement');
            ylabel({'Velocity'; '[um/sec]'});
        end
        xlabel('Time [sec]');

        sh = line([sindex_vals(1)/swindowsize sindex_vals(end)/swindowsize], [smeanvel smeanvel]);
        set(sh, 'LineStyle','--','Color','k');
        sh = line([sindex_vals(1)/swindowsize sindex_vals(end)/swindowsize], [(smeanvel+sstdvel) (smeanvel+sstdvel)]);
        set(sh, 'LineStyle','--','Color',[.5 .5 .5]);
        sh = line([sindex_vals(1)/swindowsize sindex_vals(end)/swindowsize], [(smeanvel-sstdvel) (smeanvel-sstdvel)]);
        set(sh, 'LineStyle','--','Color',[.5 .5 .5]);
        
        %Plots cutoff velocity, mean velocity and stdev velocity. 
        %Cutoff velocity is the max value which could probably be
        %calculated with sloped lines. Exceeding reduces accuracy of velocity.
        %This occurs when RBC travel so fast they exit the line scan area between frame.
        %Relies on the length of the line scan being analyzed, the frame 
        %rate and pixel size.
        subplot(8,1,7:8);
        sZer = (((sendColumn - sstartColumn)* (swindowsize))* spixel_size);
        sxt = [0 0.5];
        syt = [1 1];
        sxz = [0 0.5];
        syz = [0.6 0.6];
        sxx = [0 0.5];
        syx = [0.2 0.2];
        sstr = {'Cutoff Velocity [um/sec] =', sZer};
        sstr1 ={'Mean  Velocity [um/sec] =', smeanvel};
        sstr2 ={'Stdev  Velocity [um/sec] =', sstdvel};
        text(sxt,syt,sstr)
        hold all
        text(sxz,syz,sstr1)
        text(sxx,syx,sstr2)
        set(gca,'visible','off')
        hold off

        sboolWdthX = isempty(sVesWidth_X);
        sboolWdthY = isempty(sVesWidth_Y);
       
        if sSwitch_to_mm == true
            fprintf('\nMean  Velocity %0.2f [mm/sec]\n', smeanvel);
            fprintf('Stdev Velocity %0.2f [mm/sec]\n', sstdvel);
        else
            fprintf('\nMean  Velocity %0.2f [um/sec]\n', smeanvel);
            fprintf('Stdev Velocity %0.2f [um/sec]\n', sstdvel);
        end

        fprintf('Cutoff Velocity %0.2f [um/sec]\n', sZer);
    end
end

%% Handdrawn Line Graph
if HLine
    for imagelines = 1:numberSlices(2)
        %get line scan
        [imageLines,startColumn,endColumn,index_vals,windowsize,velocity,pixel_size,badvals, meanvel,stdvel,Switch_to_mm,numWorkers,fpscheck,full_analysis,index_length,numavgs,skipamt,rolling] = LSPIV_Parallel_Multi(LineScan{imagelines}, nImage, mImage, fpscheck, pixel_size, numWorkers, windowsize);
        FinalFigureTitle = simagelines + imagelines + 2;
        figure(FinalFigureTitle);

        realVel(imagelines, 1:length(velocity)) = velocity(1,1:length(velocity));

        %save excel data
        size_index_vals = size(index_vals);

        excel_array = zeros(2,size_index_vals(2));

        excel_array(1,1:end) = (index_vals/windowsize);
        excel_array(2,1:end ) = (velocity*pixel_size*windowsize);

        excel_array(3,1) = index_length;
        
        writematrix(excel_array,newExcelFile,'Sheet',FinalFigureTitle-2);

        %imagelines
        %subplot(numberSlices(2)*2,1,(imagelines*2)-1)
        subplot(8,1,1:2);

        %%%% linescan of selected data
        size1 = endColumn - startColumn;
        if full_analysis
            size2 = NumberImages;
        else
            size2 = NumberImages-numavgs;
        end
        imgtmp = zeros([size1+1 size2]); % to enable BW and color simultaneously
        imgtmp(:,:,1) = rot90(imageLines(1:size2,startColumn:endColumn));
        imgtmp(:,:,2) = rot90(imageLines(1:size2,startColumn:endColumn)); 
        imgtmp(:,:,3) = rot90(imageLines(1:size2,startColumn:endColumn));
        imagesc(imgtmp/max(max(max(imgtmp))))

        %Save Image
        fname1 = erase(fullFileName, '.tif');
        snapshotName = strcat(fname1, "_",int2str(imagelines),'_Hand_Line_scan.jpeg');
        
        soloScan = uint8(zeros([size2 size1+1]));
        
        soloScan(:,:,1) = imageLines(1:size2,startColumn:endColumn);
        soloScan(:,:,2) = imageLines(1:size2,startColumn:endColumn);
        soloScan(:,:,3) = imageLines(1:size2,startColumn:endColumn);
        
        imwrite(soloScan,snapshotName,'jpeg');
        %End Save Image
        
        t = title(imagelines);

        ylabel('[pixels]');

        %velocity
        if Switch_to_mm == true
            %subplot(numberSlices(2)*2,1,(imagelines*2))
            subplot(8,1,4:5)
           
            
            if full_analysis
                plot(index_vals(1:index_length)/windowsize, abs((velocity(1:index_length)*pixel_size*windowsize)/1000),'g.');
                hold all
                plot(index_vals(index_length:end)/windowsize, abs((velocity(index_length:end)*pixel_size*windowsize)/1000),'b.');
                plot(index_vals(badvals)/windowsize, abs(velocity(badvals)*pixel_size*windowsize/1000), 'ro');
                disp('Full analysis');
            else
                plot(index_vals/windowsize, abs((velocity*pixel_size*windowsize)/1000),'.');
                hold all
                plot(index_vals(badvals)/windowsize, abs(velocity(badvals)*pixel_size*windowsize/1000), 'ro');
            end
        
            
            hold off
            xlim([index_vals(1)/windowsize index_vals(end)/windowsize]);
            upperbound = meanvel+stdvel*4;
            lowerbound = meanvel-stdvel*4;
            
            if lowerbound < 0
                lowerbound = 0;
            end
            
            ylim([lowerbound upperbound]);
            title('Fitted Pixel Displacement');
            ylabel({'Velocity'; '[mm/sec]'});
        else
            subplot(8,1,4:5)
            
            
            if full_analysis
                plot(index_vals(1:index_length)/windowsize, abs((velocity(1:index_length)*pixel_size*windowsize)),'g.');
                hold all
                plot(index_vals(index_length:end)/windowsize, abs((velocity(index_length:end)*pixel_size*windowsize)),'b.');
                plot(index_vals(badvals)/windowsize, abs(velocity(badvals)*pixel_size*windowsize), 'ro');
                disp('Full analysis');
            else
                plot(index_vals/windowsize, abs((velocity*pixel_size*windowsize)),'.');
                hold all
                plot(index_vals(badvals)/windowsize, abs(velocity(badvals)*pixel_size*windowsize), 'ro');
            end
            
            
            hold off
            xlim([index_vals(1)/windowsize index_vals(end)/windowsize]);
            upperbound = meanvel+stdvel*4;
            lowerbound = meanvel-stdvel*4;
            
            if lowerbound < 0
                lowerbound = 0;
            end
            
            ylim([lowerbound upperbound]);
            title('Fitted Pixel Displacement');
            ylabel({'Velocity'; '[um/sec]'});
        end
        xlabel('Time [sec]');

        h = line([index_vals(1)/windowsize index_vals(end)/windowsize], [meanvel meanvel]);
        set(h, 'LineStyle','--','Color','k');
        h = line([index_vals(1)/windowsize index_vals(end)/windowsize], [(meanvel+stdvel) (meanvel+stdvel)]);
        set(h, 'LineStyle','--','Color',[.5 .5 .5]);
        h = line([index_vals(1)/windowsize index_vals(end)/windowsize], [(meanvel-stdvel) (meanvel-stdvel)]);
        set(h, 'LineStyle','--','Color',[.5 .5 .5]);

        subplot(8,1,7:8);
        Zer = (((endColumn - startColumn)* windowsize)* pixel_size);
        xt = [0 0.5];
        yt = [1 1];
        xz = [0 0.5];
        yz = [0.6 0.6];
        xx = [0 0.5];
        yx = [0.2 0.2];
        str = {'Cutoff Velocity [um/sec] =', Zer};
        str1 ={'Mean  Velocity [um/sec] =', meanvel};
        str2 ={'Stdev  Velocity [um/sec] =', stdvel};
        text(xt,yt,str)
        hold all
        text(xz,yz,str1)
        text(xx,yx,str2)
        set(gca,'visible','off')
        hold off

        boolWdthX = isempty(VesWidth_X);
        boolWdthY = isempty(VesWidth_Y);
       
        if Switch_to_mm == true
            fprintf('\nMean  Velocity %0.2f [mm/sec]\n', meanvel);
            fprintf('Stdev Velocity %0.2f [mm/sec]\n', stdvel);
        else
            fprintf('\nMean  Velocity %0.2f [um/sec]\n', meanvel);
            fprintf('Stdev Velocity %0.2f [um/sec]\n', stdvel);
        end

        fprintf('Cutoff Velocity %0.2f [um/sec]\n', Zer);
    end
end

%% Plot boxscan data in 3d picture 
if strcmpi(drawnRectangle,'YES')
    %number of boxes being graphed
    for recNumber = 1:recnumberSlices(2)
        %set values for: windowsize, num_workers, pixel_size, bloodflow type, automate window,  
        str = {'capillary', 'artery', 'user'};
        [speedSetting,v] = listdlg('PromptString','Select a configuration',...
                'SelectionMode','single',...
                'ListString',str);
        if v == 0; beep; disp('Cancelled'); end
        %create new figure per box simagelines + imagelines
        figure(simagelines + recNumber+numberSlices(2)+2);
        
        %number of slices per box
        for recSlices = 1:num_points
            %LSPIV
            [recindex_vals,recvelocity,boxMean{recNumber}(recSlices),boxSTDev{recNumber}(recSlices),recbadvals] = LSPIV_Parallel_Multi_Box(recLineScan{recNumber}{recSlices},numWorkers,pixel_size,windowsize,speedSetting);
            %add 3D plot
            %vectors must be the same length
            
            %save to excel file
            size_index_vals = size(recindex_vals);
            
            if recSlices == 1
                excel_array = zeros(recSlices*2,size_index_vals(2));
            end
            
            excel_array((recSlices*2)-1,1:end) = (recindex_vals/windowsize);
            excel_array(recSlices*2,1:end) = (recvelocity*pixel_size*windowsize);
            
            %limit plots to stdev being 20% of mean
            if (boxSTDev{recNumber}(recSlices) / boxMean{recNumber}(recSlices)) < 1
                recSlices2 = repmat(recSlices,[size(recindex_vals),1]); 
                
                plot3(recSlices2,recindex_vals/windowsize,abs(recvelocity*pixel_size*windowsize),'Color','b')
                hold all
                
                xlabel('Box Slice')
                ylabel('Time [sec]')
                zlabel({'Velocity'; '[um/sec]'})
                
                recSlices2 = repmat(recSlices,[size(recbadvals),1]);
                plot3(recSlices2,recindex_vals(recbadvals)/windowsize,abs(recvelocity(recbadvals)*pixel_size*windowsize),'ro')
                grid on
            end
        end
        hold off
        writematrix(excel_array,newExcelFile,'Sheet',simagelines + recNumber+numberSlices(2));
    end
end

%% Calculate vessel width
if strcmpi(drawnLine,'YES')
    for VesselNumber = 1:linenumberSlices(2)
        
        for lineNumber = 1:NumberImages
            linelength = size(VesselLineScan{VesselNumber}(1,:));
            linelength2 = linspace(1,linelength(2),linelength(2));
            x = linelength2(1,:);
            y = VesselLineScan{VesselNumber}(lineNumber,:);
            VesselLineWidth{VesselNumber}(lineNumber) = fwhm(x,y); 
        end
        
        %plot mean line width per second
        NumberSeconds = NumberImages/windowsize;
        
        for h = 1:floor(NumberSeconds)
            pointValue(h) = mean(VesselLineWidth{VesselNumber}((windowsize*(h-1))+1:(windowsize*h)))*pixel_size;
        end
        
        figure(simagelines + recnumberSlices(2)+numberSlices(2)+2+VesselNumber)
        
        widthlength = linspace(1,floor(NumberSeconds),floor(NumberSeconds));
        plot(widthlength,pointValue)
        
        xlabel('Time [sec]')
        ylabel('Vessel Diameter [um]')
        
        %save excel data
        excel_array = zeros(2,floor(NumberSeconds));

        excel_array(1,1:end) = widthlength;
        excel_array(2,1:end) = pointValue;
        
         writematrix(excel_array,newExcelFile,'Sheet',simagelines + recnumberSlices(2)+numberSlices(2)+VesselNumber);
        
        fprintf('\nMean Vessel Width %0.2f\n',nanmean(VesselLineWidth{VesselNumber})*pixel_size);
        fprintf('\nStDev Vessel Width %0.2f\n',std(VesselLineWidth{VesselNumber})*pixel_size);
        
        if strcmpi(HemoDyna,'YES')                     
            avgWdth(VesselNumber) = mean(VesselLineWidth{VesselNumber}*pixel_size);
            avgVel(VesselNumber) = abs(mean(realVel(VesselNumber,:))*pixel_size*windowsize);
            shearRate(VesselNumber) = (8*avgVel()/avgWdth());
            
            fprintf('Average Velocity %0.2f [um/sec]\n', avgVel(VesselNumber));
            fprintf('Average Width %0.2f [um/sec]\n', avgWdth(VesselNumber));
            fprintf('Shear Rate %0.2f [1/sec]\n', shearRate(VesselNumber));
        end
    end
end

%Close imagej and video
if strcmpi(ImgJ,'YES')
    MIJ.closeAllWindows;
    MIJ.exit;
end