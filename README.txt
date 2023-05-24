Date: 5/1/23

Data Flow:
Raw EEG during HFO and nonHFO events are preprocessed through 80-500 Hz filtering and smoothed (root mean square).  For each patient, the HFO rate along with the functional connectivity networks are first characterized with their centralities computed.  For the HFO rates and centralities, the critical resection percentage (CReP) and centrality ranks grouped into different percentiles for each electrode type are extracted.

Software requirements: 
Matlab R2021a Version 9.10 or above (no other toolboxes necessary)

Instructions:
1. Get data from: https://doi.org/10.7302/n9rt-sc45

2. Unzip "Data.zip" file.
	/BKGEEG		- Preprocessed background (nonHFO) EEG data (.mat)
	/CombinedData	- Final centrality and hfo features for all patients (.mat)
	/HFOEEG		- Preprocessed HFO specific EEG data (.mat)
	/PatientData	- Patient level data (.mat)
	
2. Download and move all "UMHS-xxxx-elecxx.mat" files into the /HFOEEG folder.

3. Download up-to-date code from GitHub repository (above) and move to same directory as the /Data folder to run.
For a given patient, change the patient identifier number under %% CHANGE ME %%.

4. Run "driver.m" file to run entire dataprocessing pipeline.  For more details on pipeline, see "driver.m."


Details on data file contents:
/Data/BKGEEG/UMHS-xxxx.mat
- Background EEG data that is filtered (80-500 Hz) and smoothed (root mean square RMS)

/Data/HFOEEG/UMHS-xxxx-elecXX.mat
- Up to 500 samples of EEG data during which a HFO was detected on the specified electrode. Like background data this is also filtred and smoothed.

/Data/CombinedData/combinedPatientData.mat
- All centrality and HFO rate features for all patients. Includes critical resection percentage (CReP) and centrality ranks grouped into different percentiles for each electrode type.

/Data/PatientData/UMHS-xxxx.mat
- Patient level data including HFO information, electrode information, patient surgical outcome, and centrality.  The outputs of 'getCReP.m' and 'getPercentileRanks.m' for CReP information and centrality ranks grouped into different percentiles are respectively saved under the variables 'crepTable' and 'prtileTable.'


Variable Abbreviations and Disambiguations:
SOZ: Seizure onset zone electrodes
RV: Resected volume electrodes
BOTH: SOZ & RV
REST: Neither SOZ nor RV
perSOZinRV (SOZRVpart): Percentage of SOZ that is in the RV
[-High] subclass: SOZRVpart >= .8
[-Low] subclass: SOZRVpart < .8