# pendelluft
<i>An algorithm for pendelluft detection in EIT imaging</i>

PhD Prof. John A. Sims; BME Victor Fidelis; Undergrad Mateus Ramirez;

This algorithm's aim is to automatically detect Pendelluft instances in an EIT video recording. I picked on this project somewhen around October/2022 and it was initially coded entirely in MATLAB. My main goal from then to March/2023 was to understand the code then port and adapt it to Python.

### What has been done? (March 2023)
1. All codes ported succesfully from MATLAB to Python
    - with the exception of "HSoptical.mat", since we will use Farneback Optical Flow Calculation
2. The codes loaded both example videos ("Normal" and "Pendelluft") and processed them simultaneously. It was changed.
    - the codes now can open an user specified file and save its results in ".mat" files for further and latter analysis
    - Instead of simultaneously processing pigs N and P, this code will process only the chosen file (from whatever pig or being it may come from)
    - It seems that Principal.m didn't save its highly precious calculated values in a specific file. That was changed in here.
3. Most of the codes are commented, so it is easier for everyone to understand what is happening
4. Added some script tools for debugging, isolated image plotting and video conversion (explained below)
5. Brief comparisons between Farneback and Horn & Schunck methods for Optical flow calculation (images located in Images/Comparisons/2023 - 3. March)

### So... what does the code do now?
(March 2023)
When the codes are run to completion, a folder is created for storing images of the Motion Vector fields U and V plotted against its respectives frames. Next to the end two main images are generated, one containing 4 plots each for Divergence-Free and Curl-Free maximum and minimum potentials, and the second containg two plots of the Φ values for both kinds of potentials. In the end the Φ values are displayed with regards to the lung side.

## How to run it?
0. Get ready the txt video files (normal and pendelluft pigs samples are provided in "txt_video_files/") and install the libraries used in the codes
1. Go to Scripts/ and run "Code_1.py"
    - Choose the txt file
    - It may take a while to render all the images depending on the number of frame the video has (can be disabled by editing the code)
2. Run "Code_2.py"
    - Choose the "Images_U_V.mat" file generated from step 1. It is in the folder with the same name as the previous txt file
    - This really will take some time
3. Run "Code_3.py"
    - Choose the file "potenciaisE_W.mat" located in the same folder as "Images_U_V.mat"
    - This will be faster than step 2
    - Save the plots if you wish
### Do I need to run everything again to see the last images?
NO! Throughout the codes' executions, various '.mat' files are saved so we can access them afterwards. If we want to generate the last images again, we can use the "Loader.py" script in "scripts\tools\\". Check out "Tools" section below.
### Tools
<i> Some tools are included with the scripts in scripts\tools\ that allow us to check some specific things</i>
- Loader.py: loads "calculos_freqmap.mat" and "calculos_phi.mat" to recreate all the outputs generated in "Code_3.py"
  - Divergence-free and Curl free potentials maximum and minimum plots (4 plots)
  - Φ plotting
  - Φ values by potential and by lung side
- ftsm.py: compares HS and Farneback methods by plotting U and V vectors from two '.mat' files onto the same frame image
- UV_trial.py: loads an user chosen '.mat' file and plots (saves too) the selected frames with its Motion Vector Fields
- video_maker_from_array.py: generates a video file in greyscale from images contained in a '.mat' file

### TODOs:
- todo

Comparison between Python Matplotlib pyplot's quivers and MATLAB's:<br>
![image](https://user-images.githubusercontent.com/126175949/222615900-4bb1d8f1-ec9c-4962-b73b-47ae11a98b54.png)
![image](https://user-images.githubusercontent.com/126175949/222615874-d1802376-6c63-47bb-986f-97ef7d824b39.png)
