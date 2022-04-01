# Balewski_Neuron_2022
Code for Balewski et al. (2022) Neuron


## Data availability
Any data, code, and additional information supporting this work are available from the lead contact, Joni Wallis, on request.

#### Raw data for each session: (default: expect in `data_raw`)
- probe coordinates
- behavioral data (MonkeyLogic: *.bhv2)
- spike times and waveforms (40 kHz, Plexon: *.pl2)
- LFP (1 kHz, Plexon: *.pl2)

#### Data after preprocessing (all *.mat): (default: expect in `data_processed`)
- reformatted behavioral data
- spike times, aligned to task events (200ms, 100ms, and 20ms boxcar smoothing)
- bandpass filtered LFP magnitudes, aligned to task events (200ms, 100ms, and 20ms boxcar smoothing)

#### Single unit regressions and population decoding output (default: expect in `output_regressions` and `output_decoding`)
- single unit regression outputs
- value decoding (free trials only): decoding strength (% change from baseline posterior probability) and states for chosen, unchosen, and unavailable values
- direction decoding (free trials only): decoding strength (posterior probability) and staets for chosen and unchosen direction

#### Expected data directory structure from server
More info here.

## Code organization

 - `info`: meta data for each recording session
 - `functions`: general matlab functions, used for all analyses
 - `data_*`: placeholders for expected data locations (contact Joni for data)
 - `preprocessing`: scripts and functions to convert raw behavior and neural data into useful formats, then do single neuron regressions and decoding; each script runs on a single session (need to update vars at the top for different sessions); produces all outputs necessary to generate paper figures in `figures_*` directories
 - `output_*`: placeholder for outputs of regressions and decoding completed by preprocessing steps
 - `figures_anatomy`: scripts and functions to make anatomy figures (Fig. 2a, 5a, S1)
 - `figures_behavior`: scripts and functions to make behavioral figures (Fig. 1)
 - `figures_units`: scripts and functions to make single unit regression figures (Fig. 2b-e, 5b-e)
 - `figures_decoding`: scripts and functions to make direction and value decoding figures (Fig. 3, 4, 6, 7, S2)

## To generate paper figures:
1. Download data (contact Joni Wallis)
2. To start from raw data: move all raw data for each session to `data_raw` and complete all `preprocessing` steps; move outputs to expected data directory structure
3. To start from regression and decoding outputs: run scripts in all `figures_*`
