# Balewski_Neuron_2022
Code for Balewski et al. (2022) Neuron


## Data
Any data, code, and additional information supporting this work are available from the lead contact, Joni Wallis, on request.

#### Expected data directory structure from server
```
.
├── [subject]
|   ├── OFC+CN
|   |   ├── [session ID]
|   |   |   ├── raw
|   |   |   |   ├── *.bhv or *.bhv2: raw behavioral data
|   |   |   |   ├── *elec_info_clean.csv: probe coordinate info
|   |   |   |   ├── *sorting_notes.xlsx: spike sorting notes
|   |   |   |   ├── *units.pl2: spike times and waveforms (40 kHz)
|   |   |   |   ├── *LFP.pl2: LFP (1 kHz)
|   |   |   |   ├── [session ID].pl2: raw pl2 file, including continuous spk channels (only some sessions)
|   |   |   ├── processed
|   |   |   |   ├── [bhv file]_clean.mat: processed trial info and timestamps
|   |   |   |   ├── *unit_MLAPDVloc.mat: contact coordinates (within subject space)
|   |   |   |   ├── *boxcar_[smoothing]_[task event]_[region].mat: boxcar smoothed firing rates and bandpassed LFP
|   |   |   |   ├── *slice_[smoothing]_[task event]_[region].mat: average firing rates and bandpassed LPF in time slice
|   |   |   ├── output_unit_regressions
|   |   |   |   ├── *expval_bin_alltrials_[task event]_[region].mat: single neuron regression output
|   |   |   ├── output_decoding
|   |   |   |   ├── decoding_direction_[region]_[task event]
|   |   |   |   |   ├── *free*.mat: decoder output
|   |   |   |   |   ├── *free*.png: quick visualization of decoder output (only some sessions)

```

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
1. Download data: contact Joni Wallis
2. To process raw data (one session at a time): 
	- move these files into `data_raw`: (from `raw` on server)
		1. *.bhv or *.bhv2: raw behavioral data
		2. *elec_info_clean.csv: probe coordinate info
		3. *units.pl2: spike times and waveforms (40 kHz)
		4. *LFP.pl2: LFP (1 kHz)
	- complete all `preprocessing` steps
	- move all files from `data_processed`, `output_unit_regressions`, and `output_decoding` to appropriate server location (or same spot in local copy of sever directory structure)
3. To generate figures (starting from regression and decoding outputs for all sessions):
	- complete all `figures_behavior` steps
	- complete all `figures_anatomy` steps
	- complete all `figures_units` steps
	- complete all `figures_decoding` steps





