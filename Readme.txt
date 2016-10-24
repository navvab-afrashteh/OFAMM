Optical-Flow Analysis Toolbox 
in Matlab (Registered) for investigating 
the spatiotemporal dynamics of 
Mesoscale brain activity (OFAMM) Version 1.0

Copyright 2016, Navvab Afrashteh, Samsoon Inayat, Mostafa Mohsenvand, Majid H. Mohajerani

Provided under the GNU General Public License, version 3 (GPL-3.0) 


Main Folder: OFAMM
contains toolbox m files and subfolders

Sample Data Folder: SampleData
contains simulated and real wide-field imaging data

What does the sample data contain?
In each folder, "ImgSeq.mat" file contains the image sequence and "Mask.mat" contains the mask to be used for the respective image sequence.
After executing optical flow algorithms on the ImgSeq.mat, the toolbox creates "uvResults.mat" file in which velocity vector fields are stored
Applying source/sink analysis onto the velocity vector fields gives "SourceSinkResults.mat" file

RoiSet.zip files were generated with ImageJ containing the pixels of interest for which trajectories and temporal velocities are determined

How to use the OFAMM toolbox?
Before using the toolbox, “OFAMM” folder and its subfolders need to be added to Matlab path.
The toolbox contains a graphical user interface (GUI) which will show up by running “OFAMM_v1.m”.
The first step includes loading of the preprocessed image sequence as well as the binary mask image (if applicable) into the GUI memory space. The loaded image sequence can be visualized 
with or without the mask (controlled by a checkbox) in the main axes shown in the GUI and using a sliding bar or “Play” pushbutton. The user now can select a region of interest or pixels 
to view the signal value (e.g. average value over pixels for the region of interest) over time in a separate figure window. Next, optical flow analysis can be performed to estimate velocity 
vector fields. The parameters for running optical flow algorithms can be set prior to running them using edit boxes in “Optical Flow Parameters” panel. The “Run Optical Flow Analysis” 
button will run the checked optical flow methods on the image sequence while showing progress bars for individual methods. The resulting velocity vector fields will be stored in the 
“uvResults.mat” file (saved in the same directory as the image sequence) with variables inside it named as “uv” followed by the respective method e.g. uvCLG variable stores the velocity 
vector fields estimated by the CLG method. The user can then visualize the vector fields by checking the “Vector Field” checkbox. Additionally, for a region of interest or an individual 
pixel, the user can view the distribution (histogram) of instantaneous speeds and directions in a separate figure window by pushing the “Instantaneous Speeds/Directions” button. By using 
the “Temporal Speeds/Trajectory Lengths” button for the selected region of interest or pixels, the user can view the trajectories of pixels in the main axes plotted on top of the current 
image sequence while in a separate figure window, the distribution of temporal speeds and length of trajectories can be seen. For finding trajectories of selected pixels and their temporal 
speeds, the user can enter starting and ending frames whose default values are first and last frames of the image sequence. The results of finding pixel trajectories and temporal speeds 
can be saved by pressing the “Save” button in the “Speed Histogram and Angle Rose Plot” panel and will be stored in a “.mat” file using the name specified in the box next to “Save” button. 
The variable stored in this file is of “Structure” type and each element of the structure includes the self-explanatory named variables containing pixel location, trajectory (streamline), 
and speed. After the estimation of velocity vector fields, “Run Source/Sink Analysis” button, will identify sources and sinks and store results in “SourceSinkResults.mat” file saved in the 
same directory as the image sequence. In this file, the variable names are also self-explanatory and includes the locations of sources and sinks in space time, their sizes and strengths, 
and the contour information that can be visualized to see shape of a source or sink. “Sources/Sinks” checkbox will allow the visualization of these results in the main axes. 

Reference:
“Optical-flow analysis toolbox for characterization of spatiotemporal dynamics in mesoscale brain activity”, Navvab Afrashteh, Samsoon Inayat, Mostafa Mohsenvand, Majid H. Mohajerani.
