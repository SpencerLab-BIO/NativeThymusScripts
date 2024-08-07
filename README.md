# NativeThymusScripts
Collection of scripts used for analysis of data in "Intravital Two-Photon Microscopy of the Native Thymus"

# Scripts and Usage:
`Intravital two-photon microscopy of the native mouse thymus.xlsx` - Excel file containing raw data used to generate figures in paper.

`VasoMetric_SpencerVersion.ijm` - Modified version of the [VasoMetric script](https://github.com/mcdowellkonnor/ResearchMacros). Modified Script no longer applies a max intensity projection to an image stack. Intended to be used while scrolling through video. No modifications made to vessel diameter calculations. 

`Manual_BloodFlow_Velocity.ijm` - Script used to manually calculate blood flow. Open video and place ROI in the approximate centroid of a negatively contrast labelled RBC for each frame. The script calculates velocity based on the displacement of the centroid across different frames. 

`ResearchScripts.ijm` - Collection of functions used for analysis:

- `addScaleBar()` - Used to properly position scale bar in same spot relative to bottom right corner.
- `volumeAnalysis()` - Get Percentage Area covered by blood vessels. Apply to ROIs generated by `getAnalysisROI()` function.
- `getAnalysisROI()` - Gets random positions within 3D object(pixel = 255). Will not directly overlap position. Saves positions to folder. Files needed: `Downscaled image` -> Image of the downscaled blood vessel network. FOV from Image will be saved to folder. `Mask Image` -> Mask corresponding to the area of original image that will have a position randomly selected from.
- `checkOverlap(x,y,z,dimension)` - called in `getAnalysisROI()` function. Returns `true` if ROI does not overlaps with other ROI in ROIManager. Returns `false` if ROI overlaps in XY within 3 slices of previous ROI.
- `getDistanceVolumeFromDistanceMap()` - Returns mask corresponding to some distance from edge of Distance Map. Assumes Distance Map generated from [3D Suite 3D Distance Map](https://mcib3d.frama.io/3d-suite-imagej/plugins/Binary/3D-Distance-Map-EVF/).
- `getChannelandDownscale()` - Get channels from stitched slices. Open each channel. Adjust pixel parameters. Downscale image by 2. Save stack to scaled folder.
- `segmentThymus()` - Segment Vessel Data with Labkit. Assumes folder contains stack of individual tif files which correspond to single 3D zplane. Saves files to new folder.


`pLSPIV` - Main Function: Pseudo-Line Scanning Particle Image Velocimetry  

Fork of original research paper which allows for the analysis of pseudo line scans: [LSPIV](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0038590)  

Installation: Remember to install Imagej plugin. Installation instructions can be found [here](https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab)  

`Multi_Blood_Flow` - Main function: Takes a blood flow video file in multiff format and calculated blood vessel flow velocity and diameter.  

Parameters:  

`Vessel_Diameter` - Get Vessel Diameter based on [FWHM](https://www.mathworks.com/matlabcentral/fileexchange/10590-fwhm) of a single image   


`VD_MBF` - Plot velocity from input line scan image

`VC_MBF_LSPIV` - Calculate velocity passed from `VD_MBF`  


`LSPIV_Parallel_Multi` - Gets velocity specified by hand drawn line 
- Called after hand drawn line scan is generated at the start of `Multi_Blood_Flow` script
  
`LSPIV_Parallel_Multi_Box` - Get velocity of line specified by box scan
- Called after box line scan is generated by `Multi_Blood_Flow` script
  
`fwhm` - Full-Width Half-Measure: Get width of blood vessel by looking at area with significant intensity
- Called at the end of the `Multi_Blood_Flow` script
  
`textprogressbar` - Print out a progress bar as program iterates through image files
- Called in `Multi_Blood_Flow` script as files are loaded and analyzed 


# 3D Model:
3D Model of Adhesion Thymus Holder
- `Adhesion Holder Design.pdf` - Diagram of Adhesion Thymus Holder.
- `Roundedv5.scad` - Openscad model of the Adhesion Thymus Holder.
- `Roundedv6.scad` - Openscad model of the Adhesion Thymus Holder. THe size of the ring is automatically included on the side of the holder. The ring of the holder will also adjust its position when the size is changed
