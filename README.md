# SBSGEM_standalone
Standalone ROOT macros for analysis of GEM data.

Requirements: Working ROOT installation; decoded GEM data in standard ROOT Tree format. 

GEM_reconstruct.C: main clustering/hit reconstruction/tracking code; developed for analysis of GEM data from INFN cosmic ray test stand. Takes decoded GEM data as input (strips fired by module and six ADC samples from APV25), reconstructs 2D hits, finds straight-line tracks, produces diagnostic histograms and ROOT Tree. 

GEM_reconstruct_HallAtest.C: same as GEM_reconstruct.C, but specialized for analysis of Hall A GEM test data from 2016. Uses additional information from the calorimeter that was part of the test; computes simple sum of calorimeter ADC values.

GEM_align.C: Alignment code; takes the output of GEM_reconstruct.C as input, determines best set of translational (x0, y0, z0) and rotational (alpha_x, alpha_y, alpha_z) "yaw, pitch, roll" offests to minimize chi-squared of straight line tracks.

Usage example:

root [0] .L GEM_reconstruct.C+
root [1] GEM_reconstruct("../HitData/GEMfixNov18/Hit_run3805_*.root","configINFN.txt","temp.root");


General usage: GEM_reconstruct( infilename, configfilename, outfilename );

Here infilename gives the absolute or relative path name from the working directory to the decoded GEM data file in ROOT Tree format.

configfilename gives the name of the text configuration file, which defines the geometry and all relevant parameters for the clustering, hit reconstruction, and tracking.

outfilename is the desired name of the output ROOT file.

Similarly, for the alignment code:

root [0] .L GEM_align.C+
root [1] GEM_align(inputfilename, configfilename, outputfilename);

where inputfilename is the name of a root file that is the output of GEM_reconstruct containing the relevant reconstructed hit and track information, configfilename is the text configuration file, and outputfilename is the name of a text output file containing the new alignment parameters for each module. The output of GEM_align can be directly copied and pasted into the configuration file for GEM_reconstruct to repeat the tracking with the new alignment parameters.

Example configuration files are:

configINFN.txt: Example configuration file for INFN GEM four layer cosmic test data from Nov. 2018
configHallAtest.txt: Example configuration file for five layer UVA GEM Hall A beam test data from 2016.
configalign.txt: Example alignment configuration file for INFN four-layer cosmic test data.
configalignHallAtest.txt: Example alignment configuration file for five-layer UVA GEM Hall A beam test data 2016.

Detailed documentation of configuration parameters is forthcoming. In the mean time, read the source code.




