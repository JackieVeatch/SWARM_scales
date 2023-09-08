# SWARM_scales
code for identifying the physical and biological timescales of Palmer Deep Canyon during the SWARM 2020 field season using glider and HFR data

Physical timescale of Palmer Deep Canyon using project SWARM HFR data. Defined as the time is takes for the autocorrelation to pass the efold.

Biological timescales of Plamer Deep Canyon using project SWARM glider data. Glider was piloted to hold station at head of canyon. Defined as the average time of consecuatively measured phytoplankton patches (high integrated phytoplankton).
Note: change variable in line 100 to change definition of phytoplankton patch. I used integrated particle backscatter to negate effects of non photochemical quenching. Currently, chlorophyll_a ia coded into the script. I saw little difference in results switching between particle backscatter and chlorophyll a.
