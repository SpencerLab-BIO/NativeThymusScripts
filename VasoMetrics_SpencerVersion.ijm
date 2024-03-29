// Developed by Konnor McDowell at the Seattle Children's Research Institute, Summer 2019
// Updated to include movie support, general bug fixes, and general improvements by Dr. Pablo Blinder's lab, June 2021 (thanks!!!!)

//Image handling edited by Christian Burns and Ruth Verriender at the Spencer Lab, UC Merced, Winter 2023

test1();

function test1(){
	count = 0;
	originalFileName = getTitle();

	while(true){
		//waitForUser("Go to Z Plane");
		
		Stack.getPosition(a1, b2, c3); //channels, slices, frames
		run("Make Substack...", "  slices="+b2+"");
		FileName = getTitle();
		print("Vessel:"+(count+1));
	
		Vaso(count);

		close(FileName);
		selectWindow(originalFileName);
	
		count=count+1;
		//if(count > 3){
			//checkContinue();
			//count = 0;
		//}
		roiManager("show all without labels");
		checkContinue();
		
	}
	print("FINISHED");
}

function test2(){
	count = 0;
	while (true) {
		print("Calling Dialog");
		callDialog();
		print("Dialog Called");
	
		print("Calling Vaso Loop");
		Vaso();
		print("Exit Vaso");
		
		count=count+1;
		if(count > 10){
			checkContinue();
			count = 0;
		}
	}
}

function callDialog(){
	waitForUser("Measure");
	
}

function checkContinue(){
	Dialog.create("Continue?");
	UserResp = Dialog.addString("Continue?", "Yes");
	Dialog.show();
}

function Vaso(count){
	//clearAll();
	
	setOption("ExpandableArrays", true);

	// Prepare for operation by getting information about the image
	//originalFileName = getInfo("image.filename");
	originalFileName = getTitle();
	getDimensions(width, height, channels, slices, frames); // x and y are width and height, slices is z?
	getPixelSize(pixelUnit, pixelWidth, pixelHeight);
	//print("Got info");

	//print("originalFileName:",originalFileName);
	//print("slices:",slices);
	//print("frames:",frames);

	//Point where Max intensity projection is taken
	 //"projection=[Max Intensity] all"); //run MIP over intervals (selectable?) instead of all
	 //Move Outside Vaso Loop=====================
	//Stack.getPosition(a1, b2, c3); //channels, slices, frames
	//run("Make Substack...", "  slices="+b2+"");
	//originalFileName = getTitle();
	//print("Got Substack");
	//=============================================
	
	//Check if ROI is populated
	useOld = false;
	//if (roiManager("count") > 0) {
	//	useOld = !getBoolean("ROI already exists. Use the existing ROI for measurements?", "Reset ROI", "Use Existing ROI");
	//	if(!useOld) roiManager("reset"); 
	//}

	//Clear old result //Modify to store results from each MIP
	//if (nResults > 0) showMessageWithCancel("Action Required","Results table is populated. Proceeding will clear the results.");
	//run("Clear Results");
	
	//Ask user for line
	if (!useOld) {
		//print("test1");
		setTool("polyline");
		waitForUser("Please draw through line.\nRight click when selecting the final point.\nClick OK when finished.");
		if (selectionType() == -1) exit();
		if (selectionType() != 6 && selectionType() != 5) exit("Through line must be a line or polyline.");
		getSelectionCoordinates(x, y);
	}
	//print("Got Lines");

	
	startingSlice = getSliceNumber();
	meanFwhms = newArray;
	meanFwhmSTDEVs = newArray;
	maxFWHM = 0;
	
	if (!useOld) {
		//print("test2");
		// Determine the settings for the cross-lines
		cLength = getCLineLength(x, y, height);
		cSpace = getNumber("Enter the desired distance between cross-lines (" + pixelUnit + "):", 5 * pixelWidth) / pixelWidth;
		
		// Calculate any excess distance on the through-line so that the cross-lines can be centered
		totalLength = 0;
		for (i = 0; i < x.length - 1; i++) totalLength += sqrt(pow(x[i+1]-x[i],2) + pow(y[i+1]-y[i],2));
		excessLength = totalLength % cSpace;
	}
 	else {
		print("test")
		roiManager("select", 0);
		getSelectionCoordinates(cXs, cYs);
		cLength = sqrt(pow((cXs[1] - cXs[0]), 2) + pow((cYs[1] - cYs[0]), 2)) / 2;
	}
	selectWindow (originalFileName); 

	//Get measurements from original file
	//slices = frames;
	print("Getting Measurements");
	for (slice = 1; slice <= slices; slice++) {
		//print("SLice: "+slices);
		showProgress(slice, slices);
		fwhms = newArray;

		if (!useOld) {
			// Parse through the line segments and create the cross-lines
			if (slice == 1) {
				lastX = 0;
				lastY = 0;
				for (i = 0; i < x.length - 1; i++) {
					//print("A loop started");
					slope = (y[i+1] - y[i]) / (x[i+1] - x[i]);
					inverseSlope = -1 / slope;
					yInt = y[i] + (-1 * slope * x[i]);
		
					// Move the starting point based on the excess so that the cross-lines are centered
					if (i == 0) {
						val = getIntersection(x[i], y[i], excessLength / 2, slope, yInt);
						if (x[i+1] > x[i]) {
							movingX = maxOf(val[0], val[1]);
						} else {
							movingX = minOf(val[0], val[1]);
						}
					}
		
					// Move the starting points on line-segments following the first so that equal cross-line distance is maintained
					if (i > 0) {
						val = getIntersection(lastX, lastY, cSpace, slope, yInt);
						if (x[i+1] > x[i]) {
							movingX = maxOf(val[0], val[1]);
						} else {
							movingX = minOf(val[0], val[1]);
						}
					}
			
					while (movingX <= maxOf(x[i], x[i+1]) && movingX >= minOf(x[i], x[i+1])) {
						//print("In another loop");
						movingY = slope * movingX + yInt;
						invYInt = movingY + (-1 * inverseSlope * movingX);
						val = getIntersection(movingX, movingY, cLength, inverseSlope, invYInt);
				
						// The Cross Line:
						makeLine(maxOf(val[0], val[1]), inverseSlope * maxOf(val[0], val[1]) + invYInt, minOf(val[0], val[1]), inverseSlope * minOf(val[0], val[1]) + invYInt);
						Roi.setStrokeColor("red");
						roiManager("add");
						roiManager("select", roiManager("count") - 1);
						roiManager("rename", (count+1)+"-"+roiManager("count"));
		
						// Move the current x value
						lastX = movingX;
						lastY = movingY;
						val = getIntersection(movingX, movingY, cSpace, slope, yInt);
						if (x[i+1] > x[i]) {
							movingX = maxOf(val[0], val[1]);
						} else {
							movingX = minOf(val[0], val[1]);
						}
					}
					roiManager("show all without labels");
				}
			}
	
		}
	
		// Parse through the cross-lines and obtain the FWHM values
		profiles = getProfilesForSlice(slice);
		fwhms = getFWHMFromProfiles(profiles, (cLength * 2) + 1, pixelWidth);
	
		nonNanFwhms = newArray;
		for (i = 0; i < fwhms.length; i++) if (!isNaN(fwhms[i])) nonNanFwhms[nonNanFwhms.length] = fwhms[i];
		Array.getStatistics(nonNanFwhms, min, max, mean, stdDev);
		//Get# of Results measurements
		//Set Results to next one
		numData = nResults();

		//setResult("Frame", slice-1, slice);
		//Stack.getUnits(A, B, C, Time, Value);
		//setResult("Mean (" + pixelUnit + ")", slice-1, mean);
		//setResult("SD", slice-1, stdDev);
		setResult("Frame", numData, slice);
		Stack.getUnits(A, B, C, Time, Value);
		setResult("Mean (" + pixelUnit + ")", numData, mean);
		setResult("SD", numData, stdDev);
		
		
		for (i = 0; i < fwhms.length; i++) {
			if ((i+1) <= roiManager("count")) {
				//if (fwhms[i] == 0) setResult("Line " + (i+1), slice-1, NaN);
				//else setResult("Line " + (i+1), slice-1, fwhms[i]);
				if (fwhms[i] == 0) setResult("Line " + (i+1), numData, NaN);
				else setResult("Line " + (i+1), numData, fwhms[i]);
			}
		}
		//print("Some loop ended");
	
		// Set the maximum FWHM measurement to be used in the heatmap
		if (max > maxFWHM) maxFWHM = max;
	
		if (mean == NaN || mean == (1/0) || mean == -(1/0)) {
			mean = 0;
			stdDev = 0;
		}
		meanFwhms[meanFwhms.length] = mean;
		meanFwhmSTDEVs[meanFwhmSTDEVs.length] = stdDev;
	}
	setSlice(startingSlice);
	updateResults();
	
	if (slices > 1) {
		xValues = Array.slice(Array.getSequence(slices + 1), 1);
		xLabel = "Frame Number";	
		Array.getStatistics(meanFwhms, min, maxY, mean, stdDev);
		Array.getStatistics(xValues, minX, maxX, mean, stdDev);
		Plot.create("Plot of Results", xLabel, "Mean FWHM (" + pixelUnit + ")");
		Plot.setLimits(minX, maxX, 0, maxY);

		Plot.add("connected circle", xValues, meanFwhms);
		Plot.add("error bars", xValues, meanFwhmSTDEVs);
		Plot.show();
	}
}

function getIntersection(centerX, centerY, radius, slope, yInt) {
	a = pow(slope, 2) + 1;
	b = 2 * ((slope * yInt) - (slope * centerY) - centerX);
	c = pow(centerY, 2) - pow(radius, 2) + pow(centerX, 2) - (2 * yInt * centerY) + pow(yInt, 2);
	x1 = (-1 * b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
	x2 = (-1 * b - sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);
	return newArray(x1, x2);
}

function index(a, value) { 
	for (i=0; i<a.length; i++) if (a[i]==value) return i; 
	return -1; 
} 

function getProfilesForSlice(sliceNum) {
	profiles = newArray;
	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select", i);
		// Obtain the intensity profile and normalize it
		getDimensions(width, height, channels, slices, frames);
		
		if (slices > 1) {
			setSlice(sliceNum);
		} else if (frames > 1) {
			Stack.setFrame(sliceNum);
		}

		profile = getProfile();
		Array.getStatistics(profile, min, max, mean, stdDev);
		for (c=0; c<profile.length; c++) profile[c] = (profile[c] - min) / (max - min);
		
		profiles = Array.concat(profiles, profile);
	}
	return profiles;
}

function getFWHMFromProfiles(profiles, profileLength, pixelScale) {
	fwhms = newArray;
	for (i = 0; i < (profiles.length / profileLength); i++) {
		profile = Array.slice(profiles, i * profileLength, (i+1) * profileLength - 1);
		// Obtain the FWHM value for this profile
		Array.getStatistics(profile, min, max, mean, stdDev);
		halfMax = max / 2;
		

	intersects = getYIntersects(halfMax, profile);
		Array.getStatistics(intersects, min, max, mean, stdDev);
	
		fwhm = (max - min) * pixelScale;

		// Determine the derivative of this profile and use it to expand the bounds of FWHM height
		derivative = newArray;
		for (x = 0; x < (profile.length - 1); x++) derivative[derivative.length] = profile[x+1] - profile[x];
		copiedProfile = Array.copy(profile);
		Array.sort(copiedProfile);
		median = (copiedProfile[floor((copiedProfile.length - 1) / 2)] + copiedProfile[round((copiedProfile.length - 1) / 2)]) / 2;
		intersects = getYIntersects(0, derivative);
		leftIntChange = 0;
		rightIntChange = profile.length - 1;
		for (x = 0; x < intersects.length; x++) {
			if (intersects[x] > leftIntChange && min - intersects[x] > 0) leftIntChange = intersects[x];
			if (intersects[x] < rightIntChange && intersects[x] - max > 0) rightIntChange = intersects[x];
		}

		// Using the adjusted intersects, recalculate the half max and find the x distance
		Array.getStatistics(profile, min, max, mean, stdDev);
		halfMax = (max - minOf(profile[leftIntChange], profile[rightIntChange])) / 2;
		//halfMax = max / 2;
		fwhm = fwhmFromProfile(profile, halfMax, ((rightIntChange + leftIntChange) / 2), false) * pixelScale;		
		fwhms[fwhms.length] = fwhm;
	}
	return fwhms;
}

function getCLineLength(x, y, imgHeight) {
	if (getBoolean("Select cross-line length method", "Automatically Calculate Length", "Manually Enter Length")) {
		getDimensions(width, height, channels, slices, frames);
		if (slices > 1 || frames > 1) {
			run("Z Project...", "projection=[Max Intensity] all");
		} else {
			run("Duplicate...", " ");
		}
		
		run("Median...", "radius=10 stack");
		maxPeakDist = 0;
		for (i = 0; i < x.length - 1; i++) {
			invM = -1 / ((y[i+1] - y[i]) / (x[i+1] - x[i]));
			invB = ((y[i+1] + y[i])/2) + (-1 * invM * ((x[i+1] + x[i]) / 2));

			if (invM > 0) {
				makeLine(0, invB, (imgHeight - invB) / invM, imgHeight);
			} else {
				makeLine(0, invB, (0 - invB) / invM, 0);	
			}
			
			profile = getProfile();
			dist = sqrt(pow(((x[i+1] + x[i]) / 2) - 0, 2) + pow(((y[i+1] + y[i])/2) - invB, 2));
			sample = Array.slice(profile,dist-10,dist+10);
			Array.getStatistics(sample, min, max, mean, stdDev);
			halfMax = max / 2;
			fwhm = fwhmFromProfile(profile, halfMax, dist, true);

			peakDist = fwhm + (0.65 * fwhm);
			if (peakDist > maxPeakDist) maxPeakDist = peakDist;
		}
		close();

		if (maxPeakDist == 0 || maxPeakDist > 80) {
			return getNumber("Automatic Line Length Calculation Failed. Please input length for crosslines (in pixels).", 20) / 2;
		} else {
			return maxPeakDist / 2;	
		}
	} else {
		return getNumber("Please input length for crosslines (in pixels).", 20) / 2
	}
}


function fwhmFromProfile(profile, targetY, vesselCenterX, minDist) {
	intersects = getYIntersects(targetY, profile);
	if (targetY < 0.2) return NaN;
	if (intersects.length < 2) return NaN;
	
	leftX = intersects[0];
	rightX = intersects[intersects.length - 1];
	for (x = 0; x < intersects.length; x++) {
		if (minDist && intersects[x] < vesselCenterX && vesselCenterX - intersects[x] < vesselCenterX - leftX) leftX = intersects[x];
		if (minDist && intersects[x] > vesselCenterX && intersects[x] - vesselCenterX < rightX - vesselCenterX) rightX = intersects[x];
		if (!minDist && intersects[x] < vesselCenterX && vesselCenterX - intersects[x] > vesselCenterX - leftX) leftX = intersects[x];
		if (!minDist && intersects[x] > vesselCenterX && intersects[x] - vesselCenterX > rightX - vesselCenterX) rightX = intersects[x];
	}

	return rightX - leftX;
}

function getYIntersects(targetY, fx) {
	intersects = newArray;
	for (c=0; c < (fx.length - 1); c++) {
		profileSlope = (fx[c+1] - fx[c]);
		profileYInt = fx[c] + (-1 * profileSlope * c);
		xInt = (targetY - profileYInt) / profileSlope;
		if (xInt >= c && xInt <= (c+1)) {
			intersects[intersects.length] = xInt;
		}
	}
	return intersects;
}

function clearAll(){ //At what point is this function used? (When are results cleared)
	roiManager("reset"); 
	//run("Clear Results");
}


