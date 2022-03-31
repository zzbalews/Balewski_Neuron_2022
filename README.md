# Balewski_Neuron_2022
Code for Balewski et al. (2022) Neuron


## Data availability
Any data, code, and additional information supporting this work are available from the lead contact, Joni Wallis, on request.

#### Raw data for each session:
- anatomical coordinates for recording sites
- behavioral data (MonkeyLogic: *.bhv2)
- spike times and waveforms (40 kHz, Plexon: *.pl2)
- LFP (1 kHz, Plexon: *.pl2)

#### Data after preprocessing (all *.mat):
- reformatted behavioral data
- spike times, aligned to task events (200ms, 100ms, and 20ms boxcar smoothing)
- bandpass filtered LFP magnitudes, aligned to task events (200ms, 100ms, and 20ms boxcar smoothing)
- single unit regression outputs
- value decoding (free trials only): decoding strength (% change from baseline posterior probability) and states for chosen, unchosen, and unavailable values
- direction decoding (free trials only): decoding strength (posterior probability) and staets for chosen and unchosen direction

## Code organization
 - `info`: meta data for each recording session
 - `functions`: general matlab functions, used for all analyses
 - `preprocessing`: scripts and functions to preprocess raw data into useful files to - generate paper figures
 - `results_anatomy`: scripts and functions to make anatomy figures (Fig. 2a, 5a, S1)
 - `results_behavior`: scripts and functions to make behavioral figures (Fig. 1)
 - `results_units`: scripts and functions to make single unit regression figures (Fig. 2b-e, 5b-e)
 - `results_decoding`: scripts and functions to make direction and value decoding figures (Fig. 3, 4, 6, 7, S2)
