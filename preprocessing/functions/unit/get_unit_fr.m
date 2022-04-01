function units = get_unit_fr(pl2intermediate,spkraw)

load(strjoin({pl2intermediate,spkraw},'/'),'raw_firing_rate','unit_names');

units = struct();
units.id = unit_names;
units.fr = raw_firing_rate;

end

