## Preprocessing scripts

Run scripts in order (01 --> 06) for each session. By default, scripts expect all raw data files in `../data_raw`:

- behavioral data from MonkeyLogic (*.bhv or *.bhv2); could be multiple files if needed to restart task during recording
	- e.g.: amntprob4x4_George_2021-02-11.bhv2
- probe coordinates (*.csv)
	- e.g.: George00_rec18_02112021_elec_info_clean.csv
- spike times and waveforms (sampled 40 kHz, *.pl2); could be multiple files
	- e.g.: George00_rec18_02112021-units.pl2
- LFP (sampled 1 kHz, *.pl2)
	- e.g.: George00_rec18_0211_2021-LFP.pl2


#### Convert raw data into useable .mat formats

1. `process01_bhv.m`: extract useful task variables and timestamps
2. `process02_neural.m`: smooth spike times into firing rates; bandpass filter LFP magnitudes; cut up both to align to task events (pictures onset, choice) and smooth with useful boxcars (200ms, 100ms, 20ms)
3. `process03_unitlocations.m`: calculate (ML, AP, DV) coordinates for each contact

#### Do linear regressions for all neurons

4. `process04_unitregressions.m`: model firing rate ~ trial type + lever direction + max value (+ trial number); use 100 ms windows

#### Decode value from OFC or CN sessions

5. `process05_decode_value.m`: use all neurons + LFP magnitudes (either OFC or CN channels) to train value decoder (forced trials only) and apply weights to free trials; threshold into states

#### Decoder direction from CN sessions

6. `process06_decode_direction.m`: use all neurons + LFP magnitudes (CN channels only) to train direction decoder (correct free trials) and apply weights to held-out free trials

#### Rearrange data for to generate figures

Raw data, regresion outputs, and decoder outputs are organized by subject/experiment/session on the lab synology server. All figure generating scripts assume this organization. Remember to copy the outputs from these preprocessing steps into the same directory structure, or update the figure scripts.

