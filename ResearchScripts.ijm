//Set image scale to 1x1 pixels
//Add point to image
//Set back to 0.31 x 0.31um
function addScaleBar(){
	height = getHeight();
	width = getWidth();
	
	cornerH = 8;
	cornerW = 8;
	
	BarW = 161;
	BarH = 15;
	
	makePoint(width-(cornerW+BarW),height-(cornerH+BarH), "small yellow hybrid");
}

//Get % Area covered by blood vessels
function volumeAnalysis(){
	dir = getDir("Input Folder");
	
	filelist = getFileList(dir) 
	setBatchMode(true);
	for (i = 0; i < lengthOf(filelist); i++) {
	    if (endsWith(filelist[i], ".tif")) { 
	        open(dir + File.separator + filelist[i]);
	        title = getTitle();
	        run("Auto Threshold", "method=Otsu white");
	        setOption("BlackBackground", false);
			run("Convert to Mask");
			//run("Invert LUT");
			run("Analyze Particles...", "pixel summarize");
			close(title);
	    } 
	}
	setBatchMode("exit and display");
}

//Gets random positions within 3D object(pixel = 255)
//Will not directly overlap position
//Saves positions to folder
//Files needed:
//Downscaled image -> Image of the downscaled blood vessel network. FOV from Image will be saved to folder.
//Mask -> Mask corresponding to the area of original image that will have a position randomly selected from.
function getAnalysisROI(){
	dir = getDirectory("Choose a Input folder with Stitched Slices");
	File.makeDirectory(dir+File.separator+"VesselROI");
	
	n = 60;
	getVoxelSize(width, height, depth, unit);
	Dimension = 300;
	dimensions = Dimension/width;
	
	l = getWidth();
	w = getHeight();
	
	z = getSliceNumber();
	h = nSlices;
	count = 0;
	
	Imagelist = getList("image.titles");
	for(i=0;i<Imagelist.length;i++){
		if(Imagelist[i].contains("Mask")){
			MaskTitle = Imagelist[i];
		}
		else{
			OriginalTitle = Imagelist[i];
		}
	}
	
	while (true){
		run("ROI Manager...");
		setBatchMode(true);
		selectWindow(MaskTitle);
		x1 = random()*l;
		y1 = random()*w;
		z1 = floor(random()*h)+1;
		setSlice(z1);
		pixelt = getValue(x1, y1-25);
		pixelr = getValue(x1+25, y1);
		pixell = getValue(x1-25, y1);
		pixelb = getValue(x1, y1+25);
		//Check if position overlaps with previus positions
		//if(pixel == 255 && checkOverlap(x1,y1,z1,dimensions)){
		if(pixelt == 255 && pixelr == 255 && pixell == 255 && pixelb == 255 && checkOverlap(x1,y1,z1,dimensions)){
			selectWindow(OriginalTitle);
			setSlice(z1);
			makeRectangle(x1-(Dimension/2), y1-(Dimension/2), dimensions, dimensions);
			roiManager("Add");
			run("Duplicate...", " ");
			title = getTitle();
			saveAs(dir+File.separator+"VesselROI"+File.separator+"Image"+count+".tif");
			close(title);
			count++;
		}
		if(count==n){
			run("Select All");
			roiManager("Save", dir+File.separator+"VesselROI"+File.separator+"RoiSet.zip");
			break;
		}
	}
	setBatchMode("exit and display");
}

//Checks if ROI overlaps with any of the previous positions in the ROI manager
function checkOverlap(x,y,z,dimension){
	n = roiManager('count');
	returner = false;
	if(n == 0){
		return true;
	}
	for (i = 0; i < n; i++) {
		roiManager('select', i);
		Roi.getPosition(channel, slice, frame);
		Roi.getContainedPoints(xpoints, ypoints);
		
		makeRectangle(x, y, dimension, dimension);
		for(j=0;j<xpoints.length;j++){
			if(Roi.contains(xpoints[j], ypoints[j]) && abs(z-slice)<3 ){
				return false;
			}
		}
	}
	return true;
}

//Assumes Distance Map generated from 3D Suite 3D Distance Map. 
function getDistanceVolumeFromDistanceMap(){
	radius = getNumber("Radius?", 100);
	
	//Stack.getStatistics(voxelCount, mean, min, max, stdDev);
	getStatistics(area, mean, min, max, std, histogram);
	getVoxelSize(width, height, depth, unit);
	
	print(min);
	print(max);
	//print(voxelCount);
	
	setThreshold(1, radius);
	setOption("BlackBackground", false);
	run("Convert to Mask");
}

//Get channels from stitched slices
//Open each channel
//Adjust pixel parameters
//downscale image by 2
//save stack to scaled folder
function getChannelandDownscale(){
	dir = getDirectory("Choose a Input folder with Stitched Slices");
	File.makeDirectory(dir+File.separator+"scaled");
	//import ch1
	run("Image Sequence...", "dir=["+dir+"] filter=c1 sort");
	//set pixel values
	Stack.setXUnit("um");
	run("Properties...", "pixel_width=0.31 pixel_height=0.31 voxel_depth=3");
	//downscale image
	run("TransformJ Scale", "x-factor=0.5 y-factor=0.5 z-factor=0.5 interpolation=Linear preserve");
	//save image
	saveAs(dir+File.separator+"scaled"+File.separator+"ch1.tif");
	close("*");
	run("Collect Garbage");
	
	run("Image Sequence...", "dir=["+dir+"] filter=c2 sort");
	//set pixel values
	Stack.setXUnit("um");
	run("Properties...", "pixel_width=0.31 pixel_height=0.31 voxel_depth=3");
	//downscale image
	run("TransformJ Scale", "x-factor=0.5 y-factor=0.5 z-factor=0.5 interpolation=Linear preserve");
	//save image
	saveAs(dir+File.separator+"scaled"+File.separator+"ch2.tif");
	close("*");
	run("Collect Garbage");
	
	run("Image Sequence...", "dir=["+dir+"] filter=c3 sort");
	//set pixel values
	Stack.setXUnit("um");
	run("Properties...", "pixel_width=0.31 pixel_height=0.31 voxel_depth=3");
	//downscale image
	run("TransformJ Scale", "x-factor=0.5 y-factor=0.5 z-factor=0.5 interpolation=Linear preserve");
	//save image
	saveAs(dir+File.separator+"scaled"+File.separator+"ch3.tif");
	close("*");
}

//Segment Vessel Data with Labkit
//Assumes folder contains stack of individual tif files which correspond to single 3D zplane
//Saves files to new folder
function segmentThymus(){
	dir = getDirectory("Choose a Input folder");
	filelist = getFileList(dir);
	
	ending = "tif";
	
	File.makeDirectory(dir+File.separator+"Processed Folder");
	 
	//setBatchMode(true);
	for (i = 0; i < lengthOf(filelist); i++) {
		if (endsWith(filelist[i], ending)) { 
			open(dir + File.separator + filelist[i]);
			title = getTitle();

			
			run("Segment Image With Labkit", "input=["+title+"] segmenter_file=["+dir+"v1.classifier] use_gpu=true");
			selectWindow("output");
			run("Duplicate...", "duplicate");
			run("Convert to Mask");
			run("Fill Holes");	
			run("Invert");
			
			saveAs(".tiff", dir+"Processed Folder"+File.separator+title);
			close("*");
			run("Collect Garbage");
		}
	}
	//setBatchMode("exit and display");
}

