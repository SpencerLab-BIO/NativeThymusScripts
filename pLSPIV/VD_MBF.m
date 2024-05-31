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
windowsize = 120;
fpscheck = false;
numWorkers = 4;

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

SLine = true;

if SLine
   slinenumberSlices = 50; 
   for j=1:slinenumberSlices
       sVesWidth_X{j} = 0;
       sVesWidth_Y{j} = 0; 
   end
end

%Straight Line Graph
if SLine
     for simagelines = 1:slinenumberSlices
        
        %Pass info
        [simageLines,sstartColumn,sendColumn,sindex_vals,swindowsize,svelocity,spixel_size,sbadvals, smeanvel,sstdvel,sSwitch_to_mm,snumWorkers, sfpscheck,sfull_analysis,sindex_length,snumavgs,skipamt] =  VD_MBF_LSPIV(snapShot, nImage, mImage, fpscheck, pixel_size, numWorkers, windowsize);
        sFinalFigureTitle = simagelines + 2;
        figure(sFinalFigureTitle);

        srealVel(simagelines, 1:length(svelocity)) = svelocity(1,1:length(svelocity));

        %save excel data
        ssize_index_vals = size(sindex_vals);

        excel_array = zeros(2,ssize_index_vals(2));

        excel_array(1,1:end) = (sindex_vals/swindowsize);
        excel_array(2,1:end ) = (svelocity*spixel_size*(swindowsize));

        %Save data to excel sheet
        writematrix(excel_array,newExcelFile,'Sheet',sFinalFigureTitle-2);

        %imagelines
        %subplot(numberSlices(2)*2,1,(imagelines*2)-1)
        subplot(8,1,1:2);

        %%%% linescan of selected data
        ssize1 = sendColumn - sstartColumn;
        ssize2 = nImage-skipamt;
        ssize2 = nImage-100;
        simgtmp = zeros([ssize1+1 ssize2]); % to enable BW and color simultaneously
        simgtmp(:,:,1) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn));
        simgtmp(:,:,2) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn)); 
        simgtmp(:,:,3) = rot90(simageLines(1:ssize2,sstartColumn:sendColumn));
        imagesc(simgtmp./max(max(max(simgtmp))))

        t = title(simagelines);

        ylabel('[pixels]');

        %velocity
        if sSwitch_to_mm == true
            %subplot(numberSlices(2)*2,1,(imagelines*2))
            subplot(8,1,4:5)
            plot(sindex_vals/swindowsize, abs((svelocity(1:end-100)*spixel_size*(swindowsize))/1000),'.');
            hold all
            plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*(swindowsize)/1000), 'ro');
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
            plot(sindex_vals(1:end-33)/swindowsize, abs(svelocity(1:end-33)*spixel_size*(swindowsize)),'.');
            hold all
            plot(sindex_vals(sbadvals)/swindowsize, abs(svelocity(sbadvals)*spixel_size*(swindowsize)), 'ro');
            hold off
            xlim([sindex_vals(1)/swindowsize sindex_vals(end-33)/swindowsize]);
            
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

        subplot(8,1,7:8);
        sZer = (((sendColumn - sstartColumn)* (swindowsize/skipamt))* spixel_size);
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
        break;
     end 
end


