Notes on CAM-Chem cloud chemistry code
source code: https://github.com/ESCOMP/CAM/blob/cam_development/src/chemistry/aerosol/mo_setsox.F90
Notes on Cloud Chemistry Science
Regarding the  HO2 (g) --> H2O2 (aq) reaction in setsox.F code 
This is a combination of HO2 (g) dissolving into the drop and its aqueous-phase self reaction to form H2O2 in the aqueous phase. This can be written as:
HO2 (g) ↔ HO2 (aq)
HO2 (aq) ↔ O2- (aq) + H+
HO2 (aq) + O2- (aq) → H2O2 (aq) 

I will write more later…. – Mary

Working through a potential bug with HO2(g) -> H2O2(aq)
Relevant parts of the code
https://github.com/NCAR/CAM-ACOM-dev/blob/c58c42e10e1e6f6520d10c30da8cac1ba8ee76f8/src/chemistry/aerosol/mo_setsox.F90#L236-L241
    real(r8), parameter :: kh0 = 9.e3_r8            ! HO2(g)          -> Ho2(a)
    real(r8), parameter :: kh1 = 2.05e-5_r8         ! HO2(a)          -> H+ + O2-
    real(r8), parameter :: kh2 = 8.6e5_r8           ! HO2(a) + ho2(a) -> h2o2(a) + o2
    real(r8), parameter :: kh3 = 1.e8_r8            ! HO2(a) + o2-    -> h2o2(a) + o2
    real(r8), parameter :: Ra = 8314._r8/101325._r8 ! universal constant   (atm)/(M-K)
    real(r8), parameter :: xkw = 1.e-14_r8          ! water acidity

https://github.com/NCAR/CAM-ACOM-dev/blob/c58c42e10e1e6f6520d10c30da8cac1ba8ee76f8/src/chemistry/aerosol/mo_setsox.F90#L683-L703
          kh4 = (kh2 + kh3*kh1/xph(i,k)) / ((1._r8 + kh1/xph(i,k))**2)
          ho2s = kh0*xho2(i,k)*patm*(1._r8 + kh1/xph(i,k))  ! ho2s = ho2(a)+o2-
          r1h2o2 = kh4*ho2s*ho2s                         ! prod(h2o2) in mole/L(w)/s

          if ( cloud_borne ) then
             r2h2o2 = r1h2o2*xl        &    ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  / const0*1.e+6_r8  &    ! correct a bug here ????
                  / xam
          else
             r2h2o2 = r1h2o2*xl  &          ! mole/L(w)/s   * L(w)/fm3(a) = mole/fm3(a)/s
                  * const0     &          ! mole/fm3(a)/s * 1.e-3       = mole/cm3(a)/s
                  / xam                   ! /cm3(a)/s    / air-den     = mix-ratio/s
          endif

          if ( .not. cloud_borne) then    ! this seems to be specific to aerosols that are not cloud borne
             xh2o2(i,k) = xh2o2(i,k) + r2h2o2*dtime         ! updated h2o2 by het production
          endif
Reference
Schwartz JGR 1984
Variables with units
kh0    ! Appears to be Henry’s Law constant for HO2 (HHO2 in Schwartz) (mol L-1 atm-1),
       ! although the value differs from the 1.2e-3 M atm-1 in Schwartz
kh1    ! Ka1 for HO2 (mol L-1) (K8 in Schwartz)
kh2    ! Second-order reaction rate constant (L mol-1 s-1) (k10a in Schwartz)
kh3    ! Second-order reaction rate constant (L mol-1 s-1) (k10b in Schwartz)
kh4    ! Overall rate constant for formation of H2O2 (L mol-1 s-1) (k* in Schwartz) see below
xph    ! [H+] (mol L-1)
ho2s   ! looks like the condensed-phase [HO2] (mol L-1) based on the effective Henry’s Law constant
       ! from Schwartz H* = HHO2(1 + K8 / [H+]) (mol L-1 atm-1) see below
r1h2o2 ! overall rate of production of [H2O2] (mol L-1 s-1) (d[H2O2]/dt in Schwartz) see below
r2h2o2 ! same as r2h2o2 but in (vmr s-1) where vmr = mol_h2o2 mol_air-1
xl     ! cloud liquid water concentration (L_h2o L_air-1)
xam    ! air density (molecule m-3)
const0 ! unknown conversion factor 1000 / AVOGADRO (units: ??? possibly: L mol molecule-1 m-3)

Useful term in Schwartz
[O2(-1)] = [HO2]+[O2-] (mol L-1) = [HO2](1 + K8 / [H+])

Overall [H2O2] production in Schwartz given by
d[H2O2]/dt = k*[O2(-1)]2
k* = ( k10a + k10b * ( K8 / [H+] ) ) / ( 1 + K8 / [H+] )
Factoring out units
kh4 = ( kh2 + kh3 * kh1 / xph(i,k) ) / ( ( 1._r8 + kh1 / xph(i,k) )**2 )
L mol-1 s-1 = (L mol-1 s-1 + L mol-1 s-1 mol L-1 L mol-1) / (1 + (mol L-1 L mol-1)**2) ! Looks right

ho2s = kh0 * xho2(i,k) * patm * ( 1._r8 + kh1 / xph(i,k) )
mol L-1 = mol L-1 atm_ho2-1 atm_ho2 atm_air-1 atm_air ( 1 + mol L-1 L mol-1 ) ! Looks right

r1h2o2 = kh4 * ho2s * ho2s
mol L-1 s-1 = L mol-1 s-1 mol L-1 mol L-1 ! Looks right

r2h2o2 = r1h2o2 * xl * const0 / xam
r2h2o2 = r1h2o2 * xl * ( 1e3 ) / AVOGADRO / xam
mol_h2o2 mol_air-1 s-1 = mol_h2o2 L_h2o-1 s-1 * L_h2o L_air-1 * ( 1e3 L_air m_air-3 ) * molecule_air-1 mol_air m_air3 molecule_air-1 ! this doesn’t look right

Would have expected this:
r2h2o2 = r1h2o2 * xl * ( 1e3 ) * AVOGADRO / xam
mol_h2o2 mol_air-1 s-1 = mol_h2o2 L_h2o-1 s-1 * L_h2o L_air-1 * ( 1e3 L_air m_air-3 ) * m_air3 molecule_air-1 * molecule_air mol_air-1

What could have happened?
Being off by a factor of 1e46 would be too noticeable to slip by. I’m wondering if this code never actually is triggered. It’s only on for cloud_borne == FALSE, which seems to indicate BAM aerosols. If with BAM aerosols, gas-phase HO2 is not present, then [H2O2] would never be affected by HO2 under any scenario. It could be this is why it’s turned off and flagged for a potential bug when cloud_borne == TRUE. Strangely, the unused block for r2h202 calculation when cloud_borne == TRUE, looks correct to me. [Update] I think that even if H2O2 is present with BAM, this rate will always be a factor of 1e46 too small and be essentially zero.


Answer-Changing Modifications
Here is the running list of code changes that are expected to change results by more than round-off errors.
Truncated Conversion Factors
AVOGADRO should be 6.02214076×1023 mol−1 (exact SI value)
BOLTZMANN should be 1.380649×10−23 J K-1 (exact SI value)
Pressure conversion should be 101325.0 Pa atm-1
GAS_CONSTANT should be AVOGADRO*BOLTZMANN J K-1 mol-1
[RESOLVED] Updating constants
CO2
CO2 mixing ratios from the model state should be used when available (turn off cloud chemistry otherwise)
[RESOLVED] Using CO2 from model state
Henry’s Law Parameters
Several parameters for calculating Henry’s Law constants may have been roughly converted from the values in the literature. Exact values for (e.g., pressure) conversions should be used.
I found references for most, but not all, of the parameter values (see notes here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L346-L459 )
Mary recommended looking in the JPL database
Simone mentioned that there are similar parameters being read from a NetCDF file in another part of the model that we might be able to use.
[RESOLVED] Using HLC parameters from land model (except HO2, which has a more detailed treatment)
Condensed-Phase Reactant Concentrations
In some cases, the Henry’s Law Constant (H) is used to calculate condensed phase reactant concentrations (e.g., [SO2]) and in other cases, the Effective HLC (Heff) is used, resulting the reactant concentration being the sum of all condensed-phase partitioning products (e.g., [SO2], [HSO3-], [SO32-].). These should be consistent.
See notes here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L703-L737 and here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L780-L796 
[RESOLVED] See notes in Outstanding Issues
Fixed Pressure for Condensed-Phase Reactions
When cloud_borne is FALSE, in several cases condensed phase reaction rates are calculated assuming a pressure of 1 atm. It is not clear why this is done.
See notes here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L703-L737 and here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L780-L796 
[RESOLVED] Using actual pressure instead of 1 atm
Incorrect accounting for SO4 production from H2O2
When cloud_borne is FALSE, [H2O2] is not properly reacted away under some conditions
See notes here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L749-L759 
[RESOLVED] Fixing bug
Inconsistent use of Ideal Gas Assumptions
Air density is an input to the cloud chemistry function, but in many of the rate/equilibrium calculations, various forms of the gas constant are used to calculate air density, which could differ from the input value. I suggest we remove the air density input and add comments indicating the cloud chemistry code assumes an ideal gas. (We can also reduce the number of separate times this calculation is done.)
[RESOLVED] Removing air density input, calculating it instead, and adding a note that the cloud chemistry function assumes an Ideal Gas.
Confusing/Inconsistent Logic for Keeping Concentrations Positive
The flux through the reaction pathways for H2O2 + S(IV) -> S(VI) and O3 + S(IV) -> S(VI) are roughly calculated as d[S(VI)]/dt*t, where t is the time step. Then, logic is used to ensure the flux doesn’t exceed any of the initial reactant concentrations, but the small number used to keep values positive differs from the one used elsewhere in the module, and the code can be written much more succinctly if consistent logic can be applied.
See alternative and original solutions here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L739-L779 
[RESOLVED] Using proposed logic


Collecting Snapshot Files for Cloud Chemistry
Common Steps for All Configurations
Things in purple will likely be different for different configurations.
Clone CAM and checkout the tag you want to run, and initialize the submodules
cd /glade/work/mattdawson
git clone https://github.com/NCAR/CAM-ACOM-dev.git
cd CAM-ACOM-dev
git checkout develop-add-chemistry-tests
./bin/git-fleximod update
Modify the source code to exclude the new cloud chemistry modules
Define the CPP flag DISABLE_CLOUDS at the top of the cloud_aqueous_chemistry module.
Set up the code coverage collection
Modify the ccs_config/machines/derecho/config_machines.xml file to include the lcov module
<modules>
  <command name="load">cesmdev/1.0</command>
  <command name="load">ncarenv/23.09</command>
  <command name="purge"/>
  <command name="load">craype</command>
  <command name="load">lcov</command>
</modules>
Modify the ccs_config/machines/cmake_macros/gnu.cmake file to include the coverage compiler flags
if (DEBUG)
  string(APPEND FFLAGS " -g -Wall -Og -fbacktrace -fprofile-arcs -ftest-coverage -ffpe-trap=zero,overflow -fcheck=bounds")
  string(APPEND LDFLAGS " -lgcov")
endif()
Initialize the submodels and run create a new case
export CASE_DIR=/glade/work/mattdawson/cases/cloud-chemistry/QPC6-f10_f10_mg37-cloud-chemistry
cd cime/scripts
./create_newcase --compset QPC6 --res f10_f10_mg37 --case $CASE_DIR --walltime 00:29:00 --mach derecho --project P19010000 --compiler gnu --queue develop --run-unsupported
After creating the case, run
cd $CASE_DIR
./xmlchange STOP_OPTION=nsteps
./xmlchange STOP_N=9
./xmlchange ROF_NCPL=48
./xmlchange GLC_NCPL=48
./xmlchange DEBUG=TRUE
./xmlchange COMPILER=gnu
./xmlchange DOUT_S=FALSE
./case.setup
Open up user_nl_cam w/ an editor and add these entries (note specific output variables for each configuration correspond to actual indices used in certain input/output arrays)
For all configurations, to use the new NetCDF file for Henry’s Law Constant parameters, include:
dep_data_file="/glade/work/mattdawson/cloud-chem-dev/dep_data_c20250417.nc"
BAM Configuration
mfilt=1,1,1,1,1,1
ndens=1,1,1,1,1,1
nhtfrq=1,1,1,1,1,1
write_nstep0=.false.
fincl2 = "cloud_aqh2so4_1_out","cloud_aqso4_1_out","cloud_aqso4_h2o2_out","cloud_aqso4_o3_out","cloud_cldfrc_in","cloud_cldnum_in","cloud_lwc_in","cloud_mbar_in","cloud_pdel_in","cloud_press_in","cloud_qcw_1_in","cloud_qcw_1_out","cloud_qcw_8_in","cloud_qcw_8_out","cloud_qcw_13_in","cloud_qcw_13_out","cloud_qcw_14_in","cloud_qcw_14_out","cloud_qcw_83_in","cloud_qcw_83_out","cloud_qcw_85_in","cloud_qcw_85_out","cloud_qcw_86_in","cloud_qcw_86_out","cloud_qin_1_in","cloud_qin_1_out","cloud_qin_8_in","cloud_qin_8_out","cloud_qin_13_in","cloud_qin_13_out","cloud_qin_14_in","cloud_qin_14_out","cloud_qin_83_in","cloud_qin_83_out","cloud_qin_85_in","cloud_qin_85_out","cloud_qin_86_in","cloud_qin_86_out","cloud_tfld_in","cloud_xhnm_in","cloud_xphlwc_out"
CARMA Configuration
mfilt=1,1,1,1,1,1
ndens=1,1,1,1,1,1
nhtfrq=1,1,1,1,1,1
write_nstep0=.false.
fincl2 = "cloud_aqh2so4_1_out","cloud_aqh2so4_2_out","cloud_aqh2so4_3_out","cloud_aqh2so4_4_out","cloud_aqh2so4_5_out","cloud_aqh2so4_6_out","cloud_aqh2so4_7_out","cloud_aqh2so4_8_out","cloud_aqh2so4_9_out","cloud_aqh2so4_10_out","cloud_aqh2so4_11_out","cloud_aqh2so4_12_out","cloud_aqh2so4_13_out","cloud_aqh2so4_14_out","cloud_aqh2so4_15_out","cloud_aqh2so4_16_out","cloud_aqh2so4_17_out","cloud_aqh2so4_18_out","cloud_aqh2so4_19_out","cloud_aqh2so4_20_out","cloud_aqh2so4_21_out","cloud_aqh2so4_22_out","cloud_aqh2so4_23_out","cloud_aqh2so4_24_out","cloud_aqh2so4_25_out","cloud_aqh2so4_26_out","cloud_aqh2so4_27_out","cloud_aqh2so4_28_out","cloud_aqh2so4_29_out","cloud_aqh2so4_30_out","cloud_aqh2so4_31_out","cloud_aqh2so4_32_out","cloud_aqh2so4_33_out","cloud_aqh2so4_34_out","cloud_aqh2so4_35_out","cloud_aqh2so4_36_out","cloud_aqh2so4_37_out","cloud_aqh2so4_38_out","cloud_aqh2so4_39_out","cloud_aqh2so4_40_out","cloud_aqso4_1_out","cloud_aqso4_2_out","cloud_aqso4_3_out","cloud_aqso4_4_out","cloud_aqso4_5_out","cloud_aqso4_6_out","cloud_aqso4_7_out","cloud_aqso4_8_out","cloud_aqso4_9_out","cloud_aqso4_10_out","cloud_aqso4_11_out","cloud_aqso4_12_out","cloud_aqso4_13_out","cloud_aqso4_14_out","cloud_aqso4_15_out","cloud_aqso4_16_out","cloud_aqso4_17_out","cloud_aqso4_18_out","cloud_aqso4_19_out","cloud_aqso4_20_out","cloud_aqso4_21_out","cloud_aqso4_22_out","cloud_aqso4_23_out","cloud_aqso4_24_out","cloud_aqso4_25_out","cloud_aqso4_26_out","cloud_aqso4_27_out","cloud_aqso4_28_out","cloud_aqso4_29_out","cloud_aqso4_30_out","cloud_aqso4_31_out","cloud_aqso4_32_out","cloud_aqso4_33_out","cloud_aqso4_34_out","cloud_aqso4_35_out","cloud_aqso4_36_out","cloud_aqso4_37_out","cloud_aqso4_38_out","cloud_aqso4_39_out","cloud_aqso4_40_out","cloud_aqso4_h2o2_out","cloud_aqso4_o3_out","cloud_cldfrc_in","cloud_cldnum_in","cloud_lwc_in","cloud_mbar_in","cloud_pdel_in","cloud_press_in","cloud_qcw_74_in","cloud_qcw_75_in","cloud_qcw_84_in","cloud_qcw_112_in","cloud_qcw_122_in","cloud_qcw_137_in","cloud_qcw_176_in","cloud_qin_74_in","cloud_qin_75_in","cloud_qin_84_in","cloud_qin_112_in","cloud_qin_122_in","cloud_qin_137_in","cloud_qin_176_in","cloud_qcw_74_out","cloud_qcw_75_out","cloud_qcw_84_out","cloud_qcw_112_out","cloud_qcw_122_out","cloud_qcw_137_out","cloud_qcw_176_out","cloud_qcw_1_out","cloud_qcw_11_out","cloud_qcw_21_out","cloud_qcw_31_out","cloud_qcw_41_out","cloud_qcw_51_out","cloud_qcw_61_out","cloud_qcw_71_out","cloud_qcw_81_out","cloud_qcw_91_out","cloud_qcw_101_out","cloud_qcw_111_out","cloud_qcw_121_out","cloud_qcw_131_out","cloud_qcw_141_out","cloud_qcw_151_out","cloud_qcw_161_out","cloud_qcw_171_out","cloud_qcw_181_out","cloud_qcw_191_out","cloud_qin_74_out","cloud_qin_75_out","cloud_qin_84_out","cloud_qin_112_out","cloud_qin_122_out","cloud_qin_137_out","cloud_qin_176_out","cloud_tfld_in","cloud_xhnm_in","cloud_xphlwc_out"

Submit the first run
qcmd -A P19010000 -- ./case.build
./case.submit
Processing Coverage Output
Follow the instructions to download and process the coverage data here: Running CAM with lcov on Derecho
Folder Locations
All cloned CAM code is in /glade/work/mattdawson/cloud-chemistry/
All CASE folders are in /glade/work/mattdawson/cases/cloud-chemistry/
All results are in /glade/derecho/scratch/mattdawson/

Publish Coverage Results
Copy coverage results to moffatt.cgd.ucar.edu:/project/webshare/projects/CLUBB-MF/gcov/cloud-chemistry

MAM5 Aerosols
Notes


Configuration data
./create_newcase --compset QPC6 --res f10_f10_mg37 --case $CASE_DIR --walltime 00:29:00 --mach derecho --project P19010000 --compiler gnu --queue develop --run-unsupported
Bulk Aerosols
Notes
 *****************************
 Cloud Chemistry scalar inputs
 *****************************
 ncol:           13
 lchnk:           37
 loffset:            3
 dtime:    1800.0000000000000     
 press dims:           16          26
 invariants dims:           13          26           4
 qcw dims:           13          26         103
 qin dims:           13          26         103
 setsox: pcols =           16
 setsox: pver =           26
 setsox: gas_pcnst =          103
 setsox: nfs =            4
 *********************************
 Cloud Chemistry output dimensions
 *********************************
 ncol:           13
 lchnk:           37
 qcw dims:           13          26         103
 qin dims:           13          26         103
 xphlwc dims:           13          26
 aqso4 dims:           13           1
 aqh2so4 dims:           13           1
 aqso4_h2o2 dims:           13
 aqso4_o3 dims:           13


Configuration data
./create_newcase --compset QPMOZ --res f10_f10_mg37 --case $CASE_DIR --walltime 00:29:00 --mach derecho --project P19010000 --compiler gnu --queue develop --run-unsupported
CARMA Aerosols
Notes
 *****************************
 Cloud Chemistry scalar inputs
 *****************************
 ncol:           13
 lchnk:           37
 loffset:            9
 dtime:    1800.0000000000000     
 press dims:           16          32
 invariants dims:           13          32           3
 qcw dims:           13          32         220
 qin dims:           13          32         202
 setsox: pcols =           16
 setsox: pver =           32
 setsox: gas_pcnst =          202
 setsox: nfs =            3
 setsox: adv_mass =    133.14134000000001        104.14260000000000        28.010400000000001        204.34260000000000        78.110399999999998        160.12219999999999        126.10860000000000        98.098200000000006        84.072400000000002        98.098200000000006        98.098200000000006        112.12400000000000        72.143799999999999        56.103200000000001        79.903999999999996        115.35670000000000        95.903400000000005        141.90894000000000        99.716849999999994        106.12080000000000        124.13500000000001        26.036799999999999        28.051600000000001        46.065800000000003        62.065199999999997        30.066400000000002        42.077399999999997        76.090999999999994        44.092199999999998        110.10920000000000        153.82180000000000        165.36450600000001        148.91021000000001        137.36750300000000        187.37531000000001        170.92101299999999        154.46671599999999        120.91320600000000        173.83380000000000        30.025200000000002        94.937200000000004        133.40230000000000        44.051000000000002        50.485900000000001        41.050939999999997        58.076799999999999        72.061400000000006        60.050400000000003        76.049800000000005        32.039999999999999        48.039400000000001        16.040600000000001        252.73040000000000        35.452700000000000        70.905400000000000        102.90420000000000        51.452100000000002        97.457639999999998        100.91685000000000        28.010400000000001        44.009799999999998        66.007205999999996        82.461502999999993        108.13560000000000        62.132399999999997        28.010400000000001        78.064599999999999        18.998403000000000        60.050400000000003        58.035600000000002        1.0074000000000001        2.0148000000000001        259.82361300000002        34.013599999999997        98.078400000000002        80.911400000000000        116.94800300000000        100.49370600000000        86.467905999999999        36.460099999999997        27.025140000000000        46.024600000000000        20.005803000000000        63.012340000000002        79.011740000000003        96.910799999999995        52.459499999999998        135.11493999999999        116.11239999999999        74.076200000000000        100.11300000000000        118.12720000000000        68.114199999999997        147.12594000000001        147.12594000000001        162.11794000000000        163.12533999999999        118.12720000000000        184.35020000000000        70.087800000000001        120.10080000000001        72.102599999999995        104.10140000000000        147.08474000000001        136.22839999999999        70.087800000000001        14.006740000000001        44.012880000000003        108.01048000000000        147.12594000000001        145.11114000000001        17.028939999999999        18.036339999999999        28.010400000000001        28.010400000000001        30.006139999999998        46.005540000000003        62.004939999999998        119.07434000000001        231.23954000000001        15.999400000000000        47.998199999999997        47.998199999999997        67.451499999999996        60.076400000000000        133.10014000000001        121.04794000000000        183.11774000000000        93.102400000000003        94.109800000000007        176.12160000000000        92.090400000000002        90.075599999999994        32.066000000000003        146.05641900000001        48.065399999999997        64.064800000000005        80.064200000000000        250.44499999999999        250.44499999999999        250.44499999999999        250.44499999999999        250.44499999999999        28.010400000000001        310.58240000000001        140.13440000000000        200.22600000000000        215.24014000000000        186.24140000000000        168.22720000000001        154.20140000000001        174.14800000000000        92.136200000000002        150.12600000000000        106.16200000000001        188.17380000000000        122.16140000000000        204.17320000000001        14.006740000000001        14.006740000000001        137.11220000000000        103.13520000000000        253.34819999999999        159.11480000000000        159.11480000000000        123.12760000000000        61.057800000000000        75.083600000000004        109.10180000000000        75.042400000000001        47.031999999999996        129.08959999999999        105.10880000000000        61.057800000000000        77.057199999999995        33.006200000000000        63.031399999999998        117.11980000000000        117.11980000000000        117.11980000000000        233.35579999999999        119.09340000000000        115.06380000000000        101.07920000000000        117.07859999999999        103.09399999999999        185.23400000000001        230.23213999999999        15.999400000000000        17.006799999999998        175.11420000000001        91.082999999999998        89.068200000000004        199.21860000000001        185.23400000000001        173.14060000000001        173.14060000000001        149.11859999999999        187.16640000000001        187.16640000000001        203.16579999999999        18.014199999999999     
 setsox: mwdry =    28.966000000000001     
 setsox: gravit =    9.8061600000000002     
 setsox: pi =    3.1415926535897931


 *********************************
 Cloud Chemistry output dimensions
 *********************************
 ncol:           13
 lchnk:           37
 qcw dims:           13          32         220
 qin dims:           13          32         202
 xphlwc dims:           13          32
 aqso4 dims:           13         220
 aqh2so4 dims:           13         220
 aqso4_h2o2 dims:           13
 aqso4_o3 dims:           13
CARMA aerosol Initialization data:
id_msa =           -1
id_h2so4 =           75
id_so2 =          137
id_h2o2 =           74
id_nh3 =          112
nbins =           40
nspec_max =           10
ncnst_tot =          220
nspec =           10          10          10          10          10          10          10          10          10          10          10          10          10          10          10          10          10          10          10          10           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1           1
bin_idx(           1 ,:) =            1           2           3           4           5           6           7           8           9          10
bin_idx(           2 ,:) =           11          12          13          14          15          16          17          18          19          20
bin_idx(           3 ,:) =           21          22          23          24          25          26          27          28          29          30
bin_idx(           4 ,:) =           31          32          33          34          35          36          37          38          39          40
…
bin_idx(          36 ,:) =          216           0           0           0           0           0           0           0           0           0
bin_idx(          37 ,:) =          217           0           0           0           0           0           0           0           0           0
bin_idx(          38 ,:) =          218           0           0           0           0           0           0           0           0           0
bin_idx(          39 ,:) =          219           0           0           0           0           0           0           0           0           0
bin_idx(          40 ,:) =          220           0           0           0           0           0           0           0           0           0
history_aerosol =  F




Configuration data
./create_newcase --compset QPCARMATS --res f10_f10_mg37 --case $CASE_DIR --walltime 00:29:00 --mach derecho --project P19010000 --compiler gnu --queue develop --run-unsupported

Bux-Fix PR
Makes answer-changing modifications to src/chemistry/mo_setsox.F90, in preparation for code refactoring.

Development branch: CAM-ACOM-dev/develop-cloud-bug-fix
Outstanding Issues
The configurations I’ve been using for MAM (QPC6) and BAM (QPMOZ) aerosols don’t appear to have CO2 in the model state.
Are there better configurations to use?
Does this affect the decision to turn off cloud chemistry when CO2 is not in the state?
Or could the name of CO2 be something other than “CO2”?
[RESOLVED] Used rad_cnst_get_gas() function as Francis suggested.
There is still the outstanding question of what the correct values for the S(IV) oxidation reactant concentrations are. (total condensed species concentrations, undissociated condensed species only, or something else)
See notes here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L703-L737 and here: https://github.com/NCAR/CAM-ACOM-dev/blob/bb9cb1faf32de6f95a24340862161f9294fe397a/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L780-L796 
There is a NCL script in /glade/u/home/barthm/ForCAM-SIMA/ that has the S(IV) oxidation reactant concentrations. I tried to use the same variable names as those in CAM-chem, but I removed the arrays to just scalars. The script calculates effective Henry’s law (heh2o2, heo3, heso2), their phase ratio (px), their rates of reaction (rah2o2, rao3), and production of sulfate (pso4). The rah2o2 and rao3 have values that match what is in Seinfeld and Pandis textbook (chapter 6 in my version) which is based on Hoffmann and Calvert (1985) and …. The pso4 equations include the fraction of S(IV) (f_hso3, f_so3). 

The script sets the T, p, liquid water content, and pH for doing the calculations. The pso4 results can be compared to the figure in Seinfeld and Pandis, although I use a T, p more appropriate for clouds. [Updated]
The Henry’s Law parameters have been updated to use the data from the NetCDF file Francis provided. The only one I left, for now, is the one for HO2. The NetCDF file only included the partitioning parameters, whereas the cloud chemistry code seems to include some condensed-phase reaction rate constants. Should I use the partitioning parameters for HO2 from the NetCDF file along with the condensed-phase reaction rate constants from the existing cloud chemistry code to calculate the rate of H2O2 production?
See notes here for the references I found for what’s in the current cloud chemistry code: https://github.com/NCAR/CAM-ACOM-dev/blob/64d81f9ee67f95f4254fd60c3737a555260cb5ff/src/chemistry/aerosol/cloud_aqueous_chemistry.F90#L465-L483
The reactions are also described in some comments in the original code: https://github.com/NCAR/CAM-ACOM-dev/blob/64d81f9ee67f95f4254fd60c3737a555260cb5ff/src/chemistry/aerosol/mo_setsox.F90#L267-L270 and https://github.com/NCAR/CAM-ACOM-dev/blob/64d81f9ee67f95f4254fd60c3737a555260cb5ff/src/chemistry/aerosol/mo_setsox.F90#L739-L759 

For this set of reactions, please set the following:
    real(r8), parameter :: kh0 = 690._r8            ! HO2(g)          -> Ho2(a)            reference: JPL 19-5
    real(r8), parameter :: kh1 = 1.6e-5_r8         ! HO2(a)          -> H+ + O2-             reference: JPL 19-5 
    real(r8), parameter :: kh2 = 8.3e5_r8           ! HO2(a) + ho2(a) -> h2o2(a) + o2               reference: JPL; Bielski et al. 1985
    real(r8), parameter :: kh3 = 9.7e7_r8            ! HO2(a) + O2- -> h2o2(a) + o2       reference: JPL; Bielski et al. 1985

Answering the question: Use the NetCDF file but add 1.6e-5 for the K1 dissociation constant.
[Updated]
Testing updated code
A clone of our ACOM fork of CAM with the develop-cloud-bug-fix branch checked out, and an update NetCDF file with the HLC parameters are in: /glade/work/mattdawson/cloud-chem-dev/
I gave the ncar group rwx to this folder, so you should be able to change anything. (Just note that the updated NetCDF file is not saved anywhere other than here, so there is no back up if it’s changed or deleted.)

Note that many of the updates to parameters could still contain typos (I am very prone to typos). I don’t really have a way of testing this function with the snapshot tests because the results are expected to change. Maybe we can go over the git diff linked above in our next meeting very carefully. I’m not sure how else to double check this.

I ran six scenarios (QPC6, QPMOZ, QPCARMATS) with f10_f10_mg37 with code from cam_development and with develop-cloud-bug-fix.

The “cases” are here: /glade/work/mattdawson/cases/cloud-chemistry/
The results are here: /glade/derecho/scratch/mattdawson/

If you would like to run other scenarios, just be sure to add this line to the user_nl_cam file:
dep_data_file="/glade/work/mattdawson/cloud-chem-dev/dep_data_c20250417.nc"


