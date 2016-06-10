README for the CoBaL walking code

I. General Workflow
Data is stored in files, either .csv for raw data or .mat for processed data. Each script 
loads data from previous processing steps and saves output in new files. There is a 
certain order in which scripts have to be called. All scripts will work on data in the
current directory. 

II. File naming convention
File names contain relevant information separated by underscores. The order is 
1. date, 
2. subject ID, 
3. trial type, 
4. trial number, 
5. data type, 
e.g. 20160324_JXG_calibration_002_markerTrajectories.
The date convention YYYYMMDD, but that is not assumed or enforced in the code. 
The subject ID and trial types can be any textual strings.
The trial number is a number string with at least three digits.
The data type is a description of what kind of data the file holds, lower case with new 
words starting with upper case.
There are also summary files, containing data from many trials. These drop the trial type 
and number, e.g. 20160518_HR_resultsMarker.mat.

III. Processing steps
1. General processing
importCsv 
    scans for any .csv and attempts to import the contained data. Will store different 
    data types in different files, e.g. marker, forcplate, emg, labview.
saveSubjectInfoToFile(<height>, <weight>, <gender>)
    saves the provided info and some other parameters to subjectInfo.mat
processForceplateData
    takes forceplate data from either the labview data or the analog forceplate data if 
    available. Data is filtered and transformed into the A_cw frame.
filterEmgData
    applies a sequence of filter and rectification to the EMG data
findStepEvents
stepEventGui(<trial number>) (to validate and correct step events)
findRelevantDataStretches
    
2. For stimulus response paradigms
analyzeStimulusResponse
plotStimulusResponse

3. For unperturbed walking paradigms
analyzeGaitParameters
plotGaitParameters

4. For inverse kinematics and dynamics
createModel
markerToAngleTrajectories
calculateKinematicVariables
visualizeKinematics




