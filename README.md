# NativeThymusScripts
Collection of scripts used for analysis of data in "Intravital Two-Photon Microscopy of the Native Thymus"

# Scripts and Usage:
_VasoMetric_SpencerVersion.ijm_ - Modified version of the [VasoMetric script](https://github.com/mcdowellkonnor/ResearchMacros). Modified Script no longer applies a max intensity projection to an image stack. Intended to be used while scrolling through video. No modifications made to vessel diameter calculations. 

_Manual_BloodFlow_Velocity.ijm_ - Script used to manually calculate blood flow. Open video and place ROI in the approximate centroid of a negatively contrast labelled RBC for each frame. The script calculates velocity based on the displacement of the centroid across different frames. 
