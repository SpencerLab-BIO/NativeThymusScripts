//pixelRatio();
velocity();

//Get Velocity Measurements
function velocity(){
	FrameRate1 = FirstCheck();
	getRoiDistance(FrameRate1);
}

function FirstCheck(){
	//Get Velocity Measurements
	checkPixelSetting();
	
	Width = getWidth();
	Height = getHeight();
	
	FrameRate = 30;
	
	if(Width == 512){
		FrameRate = 60;
	}
	if(Width == 256){
		FrameRate = 120;
	}
	return FrameRate;
}


function getRoiDistance(FrameRate1){	
	//Get # of ROI regions
	n = roiManager('count');
	for (i = 0; i < n; i++) {
	    roiManager('select', i);
	    roiManager("Measure");
	}
	//Get Distance Measurements
	XArray = newArray(n);
	YArray = newArray(n);
	
	for (i = 0; i < n; i++) {
		XArray[i] = Table.get("X", i);
		YArray[i] = Table.get("Y", i);
	}
	calcVelocity(XArray,YArray,FrameRate1);
}

function calcVelocity(XArray,YArray,FrameRate){
	DistanceArray = newArray(XArray.length-1);
	Sum = 0;
	for(j = 0;j < XArray.length-1; j++){
		DistanceArray[j] = Math.sqrt(Math.sqr(XArray[j+1]-XArray[j])+Math.sqr(YArray[j+1]-YArray[j]));
		Sum = DistanceArray[j] + Sum;
	}
	Average = Sum/(XArray.length-1);
	print(Average*FrameRate);
}

function checkPixelSetting(){
	run("Clear Results");
	print("\\Clear");
	
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(pixelWidth == 1 && pixelHeight == 1){
		Stack.setXUnit("um");
		run("Properties...", "pixel_width=0.3100000 pixel_height=0.3100000 voxel_depth=2.0000000");
	}
}

function pixelRatio(){
	channel1 = 3;
	channel2 = 2;
	 
	Stack.setChannel(channel1);
	run("Measure");
	Stack.setChannel(channel2);
	run("Measure");

	results1 = Table.get("Mean", 0);
	results2 = Table.get("Mean", 1);
	print(results1);
	print(results2);
	print(results2/results1);

	close("Results");
}
