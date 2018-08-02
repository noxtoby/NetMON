#!/usr/bin/env python
# coding: utf-8
__doc__ = """
ppmi_merger.py

A script to merge tables (CSV files) of interest from the 
Parkinson's Progression Markers Initiative dataset, available from 
the Laboratory Of NeuroImaging (LONI) Imaging Data Archive (IDA) 
at the University of Southern California:
http://ida.loni.usc.edu

Script developed by Neil P Oxtoby as part of his NetMON research project.
Funded by a grant from the Biomarkers Across Neurodegenerative Disease 
(BAND) program sponsored by
The Michael J Fox Foundation for Parkinson's Research, 
Alzheimer's Association, 
Alzheimer's Research UK, and
Weston Brain Institute.

Instructions:
  0. Download all of the PPMI CSV files into a single folder: 
     `PPMI_csv_path` (define yourself below)
  
  Set some variable names
"""
__author__ = "Neil P Oxtoby"
__copyright__ = "Copyright 2015-2018, Neil P Oxtoby"
__credits__ = ["Neil Oxtoby", "Coffee", "Biomarkers Across Neurodegenerative Disease"]
__license__ = "TBC"
__version__ = "1.0"
__maintainer__ = "Neil Oxtoby"
__email__ = "neil@neiloxtoby.com"
__status__ = "Pretty much done"
__python_version__ = "3.5.2"

# ## Notes on PPMI
#   - Page name (PAG_NAME) is table's key in the Data Dictionary (where PAG_NAME is present)
#   - Be aware of variable coding
#     - e.g., Gender can be
#       - text (Female/Male), or 
#       - numeric (Female = 0,1; Male = 2). Note that 0 is female of child-bearing potential.
#   - Diagnosis in lab data tables is sometimes APPRDX in clinical data tables
# 
# Useful files:
# - **Data Dictionary**
#   - `Data_Dictionary.csv`
# - CRF (Case Report Forms) - don't list the variable names
#   - `PPMI-CRF-All-in-One-AM13.pdf`
# - **Code List**: DECODE is explanation of score (CODE)
#   - `Code_List.csv`
# - PPMI derived variable definitions and score calculations
#   - `Derived_Variable_Definitions_and_Score_Calculations.csv`
# 
# Enrolled subjects: in both `SCREEN` and `RANDOM`, plus nonempty `ENROLLDT` in `RANDOM`
# 
# APPRDX (numbers as of early 2018):
#   - 1 PD (n=423) - recent diagnosis of PD (<2 yrs) who are not taking PD medications
#   - 2 HC (n=196) - without PD, >30 yrs old, no first-degree blood relative with PD
#   - 3 SWEDD (n=64) - PD subjects whose DaTscan does not show evidence of a dopaminergic deficit
#   - 4 Prodromal (n=65) - at-risk cohort (without PD), diagnosed of hyposmia or REM sleep behavior disorder (RBD)
#   - 5(6) Genetic Cohort PD (unaffected) - with (without) PD, having a genetic mutation in LRRK2, GBA, or SNCA
#   - 7(8) Genetic Registry PD (unaffected) - with (without) PD, having either:
#     - a genetic mutation in LRRK2, GBA, or SNCA or 
#     - a first-degree relative with a LRRK2, GBA, or SNCA mutation
# 
# Note that `CLINDX`\[PRIMDIAG\] = (1) Idiopathic PD; (2) Alzheimer's disease; (3) Chromosome-17 frontotemporal dementia; (4) Corticobasal degeneration; (5) Dementia with Lewy bodies; (6) Dopa-responsive dystonia; (7) Essential tremor; (8) Hemiparkinson/hemiatrophy syndrome; (9) Juv. autosomal recessive parkinsonism; (10) Motor neuron disease with parkinsonism; (11) Multiple system atrophy; (12) Neuroleptic-induced parkinsonism; (13) Normal pressure hydrocephalus; (14) Progressive supranuclear palsy; (15) Psychogenic illness; (16) Vascular parkinsonism; (17) No PD nor other neurological disorder; (18) Spinocerebellar Ataxia (SCA); (19) Other neurological disorder(s)
# 
# Visit dates are scattered. Here's some info on visits:
#   - ST is a Symptomatic Therapy visit where therapy started. Sometimes this unscheduled visit replaces a planned visit (see notes in `SIGNATURE` Signature_Form.csv)
#   - From `CODE_LIST`:
#     - Uxx are unscheduled visits (UT1 include T1/T2 + DaTscan followup);
#     - AVx are unscheduled telephone AV-133; 
#     - PW is premature withdrawal; 
#     - Txx telephone contact;
#     - V01-04: Months 3/6/9/12; V05-V12: six-monthly 18 to 60; V13-V15: months 72,84,96
#     - SCBL: screening/baseline combined;
#     - RS1: rescreen;
#     - Xxx: Transfer event
# 
# See http://www.ppmi-info.org/wp-content/uploads/2015/12/PPMI-data-access-final.mp4 for presentation with example of extracting alpha-synuclein data. Also from this presentation: beware ST visits (where symptomatic therapy started), not present in DATSCAN SBR because it's "been accounted for" (whatever that means).

# ### === Dependencies, Functions ===

import os,sys,glob
from datetime import datetime
import pandas as pd
import numpy as np
from scipy import stats

#* Only if you want to plot anything
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
plt.style.use('ggplot')
from matplotlib import rcParams # For changing default figure size

#*** Custom Functions
## Remove duplicates by taking the most recent RUNDATE
def remove_duplicates(df, ID='PATNO', keyColumns=['PATNO', 'EVENT_ID', 'TYPE', 'TESTNAME'], dateColumn='RUNDATE'):
    """
    Idea: Sort database by entries expected to be identical for repeated tests, 
          keeping only the latest result
    """
    # Sort the data frame
    df_sorted = df.sort_values(keyColumns + [dateColumn], ascending = False)
    df_len = len(df_sorted)
    # Key columns for identifying retests:
    df_sorted_keyCols = df_sorted[keyColumns]
    # Step through the sorted list, line by line.
    # If the columns match up, only keep the newest result.
    # rowsRetained stores the rows to be kept (in the sorted data frame).
    rows_retained = []
    rows_removed = []
    last_entry = []
    for row in range(0, df_len):
        entry = df_sorted_keyCols.iloc[row].tolist()
        if (entry != last_entry):
            rows_retained.append(row)
        else:
            rows_removed.append(row)
        last_entry = entry
    # Eliminate the duplicate/obsolete rows
    df_sorted = df_sorted.take(rows_retained, axis = 0)
    # Sort the database again, in a more user-friendly ascending format, and return:
    return df_sorted.sort_values(keyColumns + [dateColumn], ascending = True), rows_removed, rows_retained
#***

#* Convert super-tall (e.g., biomarkers stacked in same column) 
#* to semi-wide (biomarkers split, but still multiple rows per subject) format
def unstack_df(df_tall,idCols_Subject_Visit,tallCol_Name_Value):
    id_ = idCols_Subject_Visit[0]
    vis_ = idCols_Subject_Visit[1]
    tall_name = tallCol_Name_Value[0]
    tall_value = tallCol_Name_Value[1]
    v = df_tall[id_].values
    if type(v[0])==str:
        id_is_str = True
    else:
        id_is_str = False
    df_tall_pivot = df_tall.copy()
    df_tall_pivot['ID'] = df_tall_pivot[id_].map(str) + '_' + df_tall_pivot[vis_].map(str)
    df_tall_pivot = df_tall_pivot.pivot(index='ID',columns=tall_name,values=tall_value)
    df_tall_pivot = pd.DataFrame(df_tall_pivot.to_records())
    if id_is_str:
        df_tall_pivot[id_] = df_tall_pivot['ID'].map(lambda x: x.split('_')[0])
    else:
        df_tall_pivot[id_] = df_tall_pivot['ID'].map(lambda x: int(x.split('_')[0]))
    df_tall_pivot[vis_] = df_tall_pivot['ID'].map(lambda x: x.split('_')[1])
    df_tall_pivot.drop('ID',axis=1,inplace=True)
    # Move id and visit to be the first columns
    cols = df_tall_pivot.columns.isin([id_,vis_])
    df_tall_pivot = df_tall_pivot[df_tall_pivot.columns[cols].append(df_tall_pivot.columns[~cols])]
    return df_tall_pivot
#***

#*** Courtesy of Markus Konrad:
#    https://mkonrad.net/2016/04/16/cross-join--cartesian-product-between-pandas-dataframes.html
def df_crossjoin(df1, df2, **kwargs):
    """
    Make a cross join (cartesian product) between two dataframes by using a constant temporary key.
    Also sets a MultiIndex which is the cartesian product of the indices of the input dataframes.
    See: https://github.com/pydata/pandas/issues/5401
    :param df1 dataframe 1
    :param df1 dataframe 2
    :param kwargs keyword arguments that will be passed to pd.merge()
    :return cross join of df1 and df2
    """
    df1['_tmpkey'] = 1
    df2['_tmpkey'] = 1

    res = pd.merge(df1, df2, on='_tmpkey', **kwargs).drop('_tmpkey', axis=1)
    res.index = pd.MultiIndex.from_product((df1.index, df2.index))

    df1.drop('_tmpkey', axis=1, inplace=True)
    df2.drop('_tmpkey', axis=1, inplace=True)

    return res
#***

# ## === Manually load PPMI tables (CSV files) into pandas DataFrames ===
# Read tables using pandas.read_csv( )
#   - Drop uninteresting columns (`update_stamp`, `REC_ID`, etc.)
#   - Ensure `PATNO` is a string
#******** FIX ME **********
PPMI_csv_path = os.path.join('path','to','PPMI','Data')
PPMI_csv_path = os.path.join('/Users/noxtoby/Documents/Research/UCLPOND/Projects/201803-PPMI-EBM/Data','CSV')
#****** END FIX ME ********
PPMIMERGE_csv = os.path.normpath(os.path.join(PPMI_csv_path,'PPMIMERGE.csv')) # output file

#* This one needs special treatment
BIOLAB_csv = os.path.normpath(os.path.join(data_path,'Biospecimen_Analysis_Results_with_Elapsed_Times.csv'))

daysInAYear = 365.25
runDate = datetime.today().isoformat().replace('-','')[0:8]

#* Identify all the CSV files
PPMI_tables = glob.glob(os.path.join(PPMI_csv_path,'*.csv'))

#* Load data tables
PPMI_df = []
PPMI_df_names = []
# Naming scheme: CSV filename with certain characters removed
for k in range(0,len(PPMI_tables)):
    PPMI_df.append(pd.read_csv(PPMI_tables[k],low_memory=False))
    nom = os.path.basename(PPMI_tables[k])
    nom = nom.replace('.csv','').replace('-','_').replace('+','_').replace(',','') # remove hyphen, plus, comma
    if any('PATNO' == PPMI_df[k].columns):
        PPMI_df[k].PATNO = PPMI_df[k].PATNO.map(str)
    if any('SUBJECT_NUMBER' == PPMI_df[k].columns):
        PPMI_df[k].PATNO = PPMI_df[k].SUBJECT_NUMBER.map(str)
    PPMI_df_names.append(nom)
    #print('{0}\n - {1}'.format(PPMI_df_names[-1],PPMI_df[k].columns.values))

#* For decoding visit codes
Keys = ['PATNO','EVENT_ID']
Code_List_n = np.where([1 if n=="Code_List" else 0 for n in PPMI_df_names])[0][0]
df_Code_List = PPMI_df[Code_List_n]
EVENT_ID_CODE = df_Code_List['CODE'][df_Code_List.CDL_NAME=='EVENT_ID']
EVENT_ID_DECODE = df_Code_List['DECODE'][df_Code_List.CDL_NAME=='EVENT_ID']
EVENT_ID_CODEs, I1 = np.unique(EVENT_ID_CODE, return_index=True)
EVENT_ID_DECODEs = EVENT_ID_DECODE.iloc[I1].values
df_EVENT_ID = pd.DataFrame(data={'CODE':EVENT_ID_CODEs, 'DECODE':EVENT_ID_DECODEs})

# ***
# ##  Data preparation: tidying, renaming columns for consistency, adding dates, etc.
# 
# Keys for joining most tables: `['PATNO','EVENT_ID']`
# 
# **Beware**: joins will fail if keys are of different types (str and int, for example). 
# ### Tasks
# 1. Identify enrolled subjects - add demographics and some covariates:
#     - `RANDOM` + `SCREEN`: gender, ethnicity (`HISPLAT`), race (Amer/Asian/Afro/Pacific/White/NotSpecified), enrolment diagnosis
#     - `FAMHXPD`: family history of PD (converted to numbers of 1st-/2nd-degree relatives with PD)
#     - `MUTRSLT`: genetic category (GENECAT = 1,2,3: LRRK2,GBA,SNCA)
#     - `MUTTEST`: duration of PD at enrollment (PDDURYRS)
#     - `SOCIOECO`: dominant hand (`HANDED` = 1,2,3: I guess R,L,ambidextrous)
# 2. Tables that need dates (`INFODT`):
#   1. `BIOANALYS` - get from `LAB` (`Laboratory_Procedures_with_Elapsed_Times.csv`)
#     1. Also needs some columns renaming: `CLINICAL_EVENT` to `EVENT_ID`
#     2. Also needs to be pivoted: `TESTNAME` entries become columns, populated by `TESTVALUE`
# 3. All tables:
#   - Calculate years since baseline: convert INFODT fields to DATENUM
#   - Rename date column: add suffix (`PAG_NAME` if it exists)
#   - Harmonise variable encoding and naming, e.g., APPRDX in clinical tables is DIAGNOSIS in lab tables
# 4. Create PPMIMERGE blank table from Subjects`x`Visits
# 5. Merge DataFrames of interest using pandas.DataFrame.merge( )
#   - CSF: Alpha-synuclein, Hemoglobin, pTau, tTau, ABeta 1-42, Abeta 42; ApoE genotype 
#     - `BIOANALYS` (my pivot-table version)
#   - Demographics: Age at clinical visit, Years_bl, Medication (`PDMEDUSE`)
#   - Diagnostics - `CLINDX`\[INFODT, PRIMDIAG, DCRTREM, DCRIGID, DCBRADY\] or `PRIMDXPD`\[PRIMDIAG\]; potentially add other DX and features (DXFxxxxx and/or `PDFEAT`)
#   - Vital signs - `VITAL`\[WGTKG, HTCM, TEMPC, SYSSUP, DIASUP, HRSUP, SYSSTND, DIASTND, HRSTND\]: weight, height, temp, blood pressure, heart rate
#   - Neurological 
#     - `PECN` (Cranial nerves) \[CN1RSP,CN2RSP,CN346RSP,CN5RSP,CN7RSP,CN8RSP,CN910RSP,CN11RSP,CN12RSP\]
#     - `PENEURO` (General neurological: muscle strength, etc.) \[MSRARSP, MSRACM, MSLARSP, MSLACM, MSRLRSP, MSRLCM, MSLLRSP, MSLLCM, COFNRRSP, COFNRCM, COFNLRSP, COFNLCM, COHSRRSP, COHSRCM, COHSLRSP, COHSLCM, SENRARSP, SENRACM, SENLARSP, SENLACM, SENRLRSP, SENRLCM, SENLLRSP, SENLLCM, RFLRARSP, RFLRACM, RFLLARSP, RFLLACM, RFLRLRSP, RFLRLCM, RFLLLRSP, RFLLLCM, PLRRRSP, PLRRCM, PLRLRSP, PLRLCM\]
#   - Neuropsych
#     - `MOCA`\[MCATOT\]
#     - `SFT`\[DVS_SFTANIM, DVT_SFTANIM, VLTANIM, VLTVEG, VLTFRUIT\] (total 3: sum 1-3)
#     - `LNSPD`\[LNS_TOTRAW, DVS_LNS\] (total 21: sum 1-7.a+b+c)
#     - `HVLT`\[DVT_TOTAL_RECALL, DVT_DELAYED_RECALL, DVT_RETENTION, DVT_RECOG_DISC_INDEX\] (can split: immediate recall; delayed recall; delayed recognition)
#     - `LINEORNT`\[JLO_TOTRAW, JLO_TOTCALC, AGE_ASSESS_JLO, DVS_JLO_MSSA, DVS_JLO_MSSAE\] (odd numbered questions at first visit; even at second)
#     - `SDM`\[SDMTOTAL, DVSD_SDM, DVT_SDM\]
#   - Neurobehaviour:
#     - `GDSSHORT`\[GDSSATIS, GDSDROPD, GDSEMPTY, GDSBORED, GDSGSPIR, GDSAFRAD, GDSHAPPY, GDSHLPLS, GDSHOME, GDSMEMRY, GDSALIVE, GDSWRTLS, GDSENRGY, GDSHOPLS, GDSBETER\] (depression: 0-4 normal, 5-7 mild, 8-11 moderate, 12-15 severe)
#       - Convert to single depression score: sum(GDSSATIS, GDSGSPIR, GDSHAPPY, GDSALIVE, GDSENRGY == 0) + sum(GDSDROPD, GDSEMPTY, GDSBORED, GDSAFRAD, GDSHLPLS, GDSHOME, GDSMEMRY, GDSWRTLS, GDSHOPLS, GDSBETER == 1)
#       - Test is only considered valid if 12 or more answers are provided
#     - `STAI`\[STAIAD[1-40]\] (anxiety)
#       - Convert to a single score somehow (not sure how, as it's not in `Data_Acquisition_and_Usage_20140609.pdf`)
#     - `QUIPCS`\[TMGAMBLE, CNTRLGMB, TMSEX, CNTRLSEX, TMBUY, CNTRLBUY, TMEAT, CNTRLEAT, TMTORACT, TMTMTACT, TMTRWD, TMDISMED, CNTRLDSM\] (impulsive-compulsive disorders)
#       - Total score: sum of all
#     - `EPWORTH`\[ESS[1-8]\] Total is sum of 8 Qs, each worth 0-3 (sleepiness scale: biomarker cutpoint is normal <=9, sleepy >=10)
#     - (Not sure how to total the score, so I've omitted this. Perhaps convert to binary DX of RBD.) `REMSLEEP`\[\] (sleep behaviour disorder questionnaire)
#       - Total: sum(1-9) + ?
#       - Maybe also add `REMBHVDS`\[REMONSDT, REMONEST, RBDDXDT, RBDDXEST\] (date of RBD onset, date estimated?, date of RBD DX, date estimated?)
#   - Motor-symptoms: new UPDRS
#     - `NUPDRS1`[NP1COG, NP1HALL, NP1DPRS, NP1ANXS, NP1APAT, NP1DDS] - sum total
#     - `NUPDRS1p`[NP1SLPN, NP1SLPD, NP1PAIN, NP1URIN, NP1CNST, NP1LTHD, NP1FATG] - sum total
#     - `NUPDRS2p`[NP2SPCH, NP2SALV, NP2SWAL, NP2EAT, NP2DRES, NP2HYGN, NP2HWRT, NP2HOBB, NP2TURN, NP2TRMR, NP2RISE, NP2WALK, NP2FREZ] - sum total
#     - (prior to medication) `NUPDRS3`[NP3SPCH, NP3FACXP, NP3RIGN, NP3RIGRU, NP3RIGLU, PN3RIGRL, NP3RIGLL, NP3FTAPR, NP3FTAPL, NP3HMOVR, NP3HMOVL, NP3PRSPR, NP3PRSPL, NP3TTAPR, NP3TTAPL, NP3LGAGR, NP3LGAGL, NP3RISNG, NP3GAIT, NP3FRZGT, NP3PSTBL, NP3POSTR, NP3BRADY, NP3PTRMR, NP3PTRML, NP3KTRMR, NP3KTRML, NP3RTARU, NP3RTALU, NP3RTARL, NP3RTALL, NP3RTALJ, NP3RTCON] - sum total; Also `NUPDRS3`[NHY] is Hoehn and Yahr stage
#     - (post-medication) `PAG_NAME` = `NUPDRS3A` => additionally ANNUAL_TIME_BTW_DOSE_NUPDRS,PD_MED_USE
#     - `NUPDRS4`[NP4WDYSK, NP4DYSKI, NP4OFF, NP4FLCTI, NP4FLCTX, NP4DYSTN] - sum total
#   - Non-motor:
#     - Olfactory `UPSIT`[UPSITBK[1-4]] - sum total
#     - Autonomic `SCOPAAUT`[SCAUx]
#       - Six subscales (23A and 26 do not contribute to the score)
#         - Gastrointestinal: sum 1-7
#         - Urinary: sum 8-13 (beware coding 9 for "uses a catheter", which scores as 3 "often")
#         - Cardiovascular: sum 14-16
#         - Thermoregulatory: sum 17,18,20-21
#         - Pupillomotor: 19
#         - Sexual: 22+23 (male) or 24+25 (female) (beware 9 coding for not-applicable - scores 0)
#       - Total out of 23
#   - Imaging:
#     - `DATSCAN`[CAUDATE_R, CAUDATE_L, PUTAMEN_R, PUTAMEN_L] (occipital:striatal binding ratio, SBR in L&R caudate and putamen)
#     - AV133 (`AVIMAG`, `ind_av133_sbr`[RCAUD-S, RPUTANT-S, RPUTPOST-S, LCAUD-S, LPUTANT-S, LPUTPOST-S], `ind_av133_metadata`[scan_quality_rating_pet, 12_scan_quality_rating_mr]) 
#     - `MRI` volumes/atrophy
#     - `DTIROI` - FA and Eigenvalues (I think): see [PPMI_DTI_ROI_Methods_20160915.pdf](/Users/noxtoby/Documents/Research/UCLPOND/Projects/201803-PPMI-EBM/Data/CSV/PPMI_DTI_ROI_Methods_20160915.pdf), the six ROIs are manually-drawn ones in the L&R SN, and two reference regions in the L&R cerebral peduncle

# ### Data Preparation step 1: Identify enrolled subjects and demographics
#*** Enrolled subjects: in both SCREEN and RANDOM, plus non-null ENROLLDT in RANDOM table ***
df_RANDOM = PPMI_df[np.where([1 if n=="Randomization_table" else 0 for n in PPMI_df_names])[0][0]].copy()
df_RANDOM_enrolled = df_RANDOM.loc[df_RANDOM.ENROLLDT.notnull(),:]
df_SCREEN = PPMI_df[np.where([1 if n=="Screening___Demographics" else 0 for n in PPMI_df_names])[0][0]].copy()
# Inner join: RANDOM and SCREEN. Get APPRDX and race/ethnicity from SCREEN
df_RANDOM_enrolled = pd.merge(df_RANDOM_enrolled[['PATNO','BIRTHDT','GENDER','ENROLLDT','CONSNTDT']],
                              df_SCREEN[['PATNO','APPRDX',
                                         'HISPLAT','RAINDALS','RAASIAN','RABLACK','RAHAWOPI','RAWHITE','RANOS']].sort_values('PATNO',axis=0),
                              on='PATNO',suffixes=['_Random','_Screen'])
df_RANDOM_enrolled.rename(index=str, columns={'APPRDX':'APPRDX_SCREEN'})
#* HANDED,EDUCYRS  from SOCIOECO
df_SOCIOECO = PPMI_df[np.where([1 if n=="Socio_Economics" else 0 for n in PPMI_df_names])[0][0]].copy()
df_RANDOM_enrolled = pd.merge(df_RANDOM_enrolled,
                              df_SOCIOECO[['PATNO','HANDED','EDUCYRS']].sort_values('PATNO',axis=0),
                              on='PATNO',suffixes=['_Random','_SocioEco'])
PATNO = df_RANDOM_enrolled.PATNO
df_RANDOM_enrolled

# ### Data Preparation step 2: Tables that need dates
# #### 2.1 `BIOANALYS`
# 1. Change `CLINICAL_EVENT` to `EVENT_ID`
# 2. Add `INFODT` from `Laboratory_Procedures_with_Elapsed_Times`
# 3. Cleanup: prefer reprocessed Bioanalyses (remove old ones)
#  - (PPMI_CSF_baseline_analysis_Assay_variability_April_2014.pdf says to keep 2013 rather than 2011, but there are others that get reprocessed)
# 4. Convert `TESTVALUE` to numeric using a map: replaces text entries such as 'below detection limit' with the detection limit
# 5. Convert long format to semi-wide: one column per `TESTNAME`, populated by `TESTVALUE`
Bio_n = np.where([1 if n=="Biospecimen_Analysis_Results" else 0 for n in PPMI_df_names])[0][0]
df_BIO = PPMI_df[Bio_n].copy()
df_BIO.PI_INSTITUTION = df_BIO.PI_INSTITUTION.str.replace("’","") #* Annoying weirdly-encoded apostrophe
df_BIO.PATNO = df_BIO.PATNO.map(str)

#*** 1. Translate CLINICAL_EVENT to EVENT_ID
event_list = ['SC', 'BL', 'SCBL', 
              'V01', 'V02', 'V03', 'V04', 'V05', 'V06', 'V07', 'V08', 'V09',
              'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'PW', 'ST', 
              'U01', 'U02', 'U03', 'U04', 'U05', 'U06', 'UT1']
#** EVENT_ID descriptions
event_desc = []
for k in range(len(event_list)):
    event_desc.append(df_EVENT_ID.DECODE[df_EVENT_ID.CODE == event_list[k]].values[0])
event_dict = dict(zip(event_list, event_desc))
#** Discard all rows with non-recognized keys
df_BIO_ = df_BIO[(df_BIO['CLINICAL_EVENT'].isin(event_list))]
#** Rename the column
df_BIO_ = df_BIO_.rename(columns = {'CLINICAL_EVENT' : 'EVENT_ID'})
#** Translate (if they have a translation)
df_BIO_['EVENT_DESC'] = df_BIO_['EVENT_ID'].apply(lambda x : event_dict[x])

#*** 2. Add INFODT from Laboratory_Procedures_with_Elapsed_Times
#       (seems not required for EVENT_ID =='SC' because these were the
#        DNA test results; these dates do exist - in the VITAL table).
Lab_n = np.where([1 if n=="Laboratory_Procedures_with_Elapsed_Times" else 0 for n in PPMI_df_names])[0][0]
df_LAB = PPMI_df[Lab_n].copy()
df_LAB = df_LAB.drop(['REC_ID','F_STATUS','PAG_NAME','ORIG_ENTRY','LAST_UPDATE','QUERY','SITE_APRV'],axis=1)
#* Merge Bio to Lab
df_BIOLAB = df_BIO_.merge(df_LAB[Keys+['INFODT','PDMEDYN']],how='left',on=Keys,suffixes=('_BIO','_LAB'))

#*** 3. Remove earlier Bioanalyses that were reprocessed, and also remove rows: DNA, RNA
df_BIOLAB = df_BIOLAB.drop(['PI_INSTITUTION','PI_NAME','PROJECTID','update_stamp'],axis=1) # remove some columns for convenience
df_BIOLAB.TYPE = df_BIOLAB.TYPE.str.replace('Cerebrospinal fluid','CSF') # Recode to CSF
df_BIOLAB_ = df_BIOLAB[(~df_BIOLAB.TYPE.isin(['DNA','RNA'])) & 
                       (df_BIOLAB.TESTNAME != 'ApoE Genotype') ].copy() # remove DNA, RNA rows (except APOE genotype)
keyColumns = Keys + ['TYPE','TESTNAME']
dateColumn = 'RUNDATE'
df_BIOLAB_, rowsRemoved, rowsRetained = remove_duplicates( df_BIOLAB_, 'PATNO', keyColumns, dateColumn)

#*** 4. Convert TESTVALUE to numeric (removes 'below detection limit', etc.)
TESTVALUE_dict = {
    # CSF Hemoglobin 
    'below detection limit':20, 'below':20, '>20':20, '<20':20, 'above':12500,
    '>12500 ng/ml':12500,'>12500ng/ml':12500,
    # pTau
    '<8':8, 
    # Abeta 1-42
    '<200':200, '>1700':1700,
    # tTau
    '<80':80,
    # ApoE Genotype
    'e2/e2':'22', 'e2/e4':'24', 'e3/e2':'32', 'e3/e3':'33', 'e4/e3':'43', 'e4/e4':'44'}
df_BIOLAB_["TESTVALUE"].replace(TESTVALUE_dict, inplace=True)
s = '*'
df_BIOLAB_['key'] = df_BIOLAB_['EVENT_ID'].map(str) +s+ df_BIOLAB_['INFODT'].map(str) +s+ df_BIOLAB_['PDMEDYN'].map(str)

#*** 5. Convert long to semi-wide format: pivot on TESTNAME,TESTVALUE
df_BIOLAB_SemiWide = unstack_df(df_BIOLAB_,['PATNO','key'],['TESTNAME','TESTVALUE'])
df_BIOLAB_SemiWide['EVENT_ID'] = df_BIOLAB_SemiWide['key'].map(lambda x: x.split(s)[0])
df_BIOLAB_SemiWide['INFODT'] = df_BIOLAB_SemiWide['key'].map(lambda x: x.split(s)[1])
df_BIOLAB_SemiWide['PDMEDYN'] = df_BIOLAB_SemiWide['key'].map(lambda x: x.split(s)[2])
df_BIOLAB_SemiWide.drop('key',axis=1,inplace=True)
cols = df_BIOLAB_SemiWide.columns.isin(['PATNO','EVENT_ID','INFODT','PDMEDYN'])
df_BIOLAB_SemiWide = df_BIOLAB_SemiWide[df_BIOLAB_SemiWide.columns[cols].append(df_BIOLAB_SemiWide.columns[~cols])]
df_BIOLAB_SemiWide.to_csv(BIOLAB_csv)


# ### Data Preparation step 3: Prepare tables of interest
# - 3.1 Calculate total scores and subscores for various assessments
# 
# When merging later, rename `INFODT` columns: add suffix (e.g., `PAG_NAME` where it exists)
def calculate_scores_MOCA(df_MOCA):
    """
    Calculates MOCA subscores by summing the relevant questions:
    
    MOCA_visuospatial = MCAALTTM + MCACUBE + MCACLCKC + MCACLCKN + MCACLCKH
    MCOA_naming = MCALION + MCARHINO + MCACAMEL
    MCOA_attention = MCAFDS + MCABDS + MCAVIGIL + MCASER7
    MOCA_language = MCASNTNC + MCAVF
    MOCA_delayed_recall = MCAREC1 + MCAREC2 + MCAREC3 + MCAREC4 + MCAREC5
    MOCA_orientation = MCADATE + MCAMONTH + MCAYR + MCADAY + MCAPLACE + MCACITY
    
    Returns a list of added column names
    """
    cols_visuospatial = ['MCAALTTM','MCACUBE','MCACLCKC','MCACLCKN','MCACLCKH']
    cols_naming = ['MCALION','MCARHINO','MCACAMEL']
    cols_attention = ['MCAFDS','MCABDS','MCAVIGIL','MCASER7']
    cols_language = ['MCASNTNC','MCAVF']
    cols_recall = ['MCAREC1','MCAREC2','MCAREC3','MCAREC4','MCAREC5']
    cols_orientation = ['MCADATE','MCAMONTH','MCAYR','MCADAY','MCAPLACE','MCACITY']
    
    df_MOCA['MOCA_Visuospatial'] = df_MOCA[cols_visuospatial].sum(axis=1)
    df_MOCA['MOCA_Naming'] = df_MOCA[cols_naming].sum(axis=1)
    df_MOCA['MOCA_Attention'] = df_MOCA[cols_attention].sum(axis=1)
    df_MOCA['MOCA_Language'] = df_MOCA[cols_language].sum(axis=1)
    df_MOCA['MOCA_DelayedRecall'] = df_MOCA[cols_recall].sum(axis=1)
    df_MOCA['MOCA_Orientation'] = df_MOCA[cols_orientation].sum(axis=1)
    
    return ['MOCA_Visuospatial','MOCA_Naming','MOCA_Attention','MOCA_Language','MOCA_DelayedRecall','MOCA_Orientation']

def calculate_scores_GDSSHORT(df_GDSSHORT):
    """
    Calculates a single depression score:
    
      GDS_TOT = sum(GDSSATIS, GDSGSPIR, GDSHAPPY, GDSALIVE, GDSENRGY == 0)
              + sum(GDSDROPD, GDSEMPTY, GDSBORED, GDSAFRAD, GDSHLPLS, GDSHOME, GDSMEMRY, GDSWRTLS, GDSHOPLS, GDSBETER == 1)
    
    The GDS is valid only if 12 or more answers are provided. 
    For invalid tests, this function returns GDS_TOT = NA
    
    Returns a list of added column names
    """
    # Negative depression - depression indicated by 0 for these questions
    negative_questions = [ 'GDSSATIS', 'GDSGSPIR', 'GDSHAPPY', 'GDSALIVE', 'GDSENRGY' ]
    # Positive depression - depression indicated by 1 for these questions
    positive_questions = [ 'GDSDROPD', 'GDSEMPTY', 'GDSBORED', 'GDSAFRAD', 'GDSHLPLS', 'GDSHOME', 'GDSMEMRY', 'GDSWRTLS', 'GDSHOPLS', 'GDSBETER' ]
    
    df_GDSSHORT['GDS_TOT'] = (0 == df_GDSSHORT[negative_questions]).sum(axis=1) + (1 == df_GDSSHORT[positive_questions]).sum(axis=1)
    
    # Check for at least 12 non-missing answers, otherwise NaN
    na_rows = (np.logical_not(np.isnan(df_GDSSHORT[negative_questions + positive_questions])).sum(axis=1) < 12)
    df_GDSSHORT.set_value(na_rows,'GDS_TOT',np.nan)
    
    return ['GDS_TOT']

def calculate_scores_SCOPAAUT(df_SCOPAAUT):
    """
    Calculates SCPOPA-AUT subscores and total by summing the relevant questions:
    
    Gastrointestinal: sum 1-7
    Urinary: sum 8-13 (beware coding 9 for "uses a catheter", which scores as 3=="often")
    Cardiovascular: sum 14-16
    Thermoregulatory: sum 17,18,20-21
    Pupillomotor: 19
    Sexual: 22+23 (male) or 24+25 (female) (beware 9 coding for not-applicable - scores 0)
    Total: out of 23  (23A and 26 do not contribute to the score)
    
    Returns a list of added column names
    """
    # First check for 9 and change to 3 (urinary) or 0 (sexual)
    urinary_columns = df_SCOPAAUT[[ 'SCAU8', 'SCAU9', 'SCAU10', 'SCAU11', 'SCAU12', 'SCAU13' ]].copy().replace([9],[3])
    sexual_columns = df_SCOPAAUT[[ 'SCAU22', 'SCAU23', 'SCAU24', 'SCAU25']].copy().replace([9],[0])
    
    # Sum total for each SCOPA-AUT subscale
    df_SCOPAAUT['SCOPAAUT_gastrointestinal'] = df_SCOPAAUT[[ 'SCAU1', 'SCAU2', 'SCAU3', 'SCAU4', 'SCAU5', 'SCAU6', 'SCAU7' ]].sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_urinary'] = urinary_columns.sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_cardiovascular'] = df_SCOPAAUT[[ 'SCAU14', 'SCAU15', 'SCAU16' ]].sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_thermoregulatory'] = df_SCOPAAUT[[ 'SCAU17', 'SCAU18', 'SCAU20', 'SCAU21' ]].sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_pupillomotor'] = df_SCOPAAUT[[ 'SCAU19' ]].sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_sexual'] = sexual_columns.sum(axis=1)
    df_SCOPAAUT['SCOPAAUT_TOT'] = df_SCOPAAUT[['SCOPAAUT_gastrointestinal', 
                                               'SCOPAAUT_urinary', 
                                               'SCOPAAUT_cardiovascular', 
                                               'SCOPAAUT_thermoregulatory', 
                                               'SCOPAAUT_pupillomotor', 
                                               'SCOPAAUT_sexual']].sum(axis=1)
    return ['SCOPAAUT_gastrointestinal','SCOPAAUT_urinary','SCOPAAUT_cardiovascular','SCOPAAUT_thermoregulatory','SCOPAAUT_pupillomotor','SCOPAAUT_sexual','SCOPAAUT_TOT']

def calculate_scores_STAI(df_STAI):
    """
    Calculates STAI subscores (state and trait) and total by summing the relevant questions:
    
    Each question is scored between 1-4, so we simply sum the "affirmative" questions, 
    and invert the scale
    
    State anxiety:
      Positive anxiety (3,4,6,7,9,12,13,14,17,18)
      Negative anxiety  (1,2,5,8,10,11,15,16,19,20)
    Trait anxiety: 
      Positive anxiety (22,24,25,28,29,31,32,35,37,38,40)
      Negative anxiety (21,23,26,27,30,33,34,36,39)
    
    Returns a list of added column names
    """
    state_cols_positive = ['STAIAD3', 'STAIAD4', 'STAIAD6', 'STAIAD7', 'STAIAD9', 'STAIAD12', 'STAIAD13', 'STAIAD14', 'STAIAD17', 'STAIAD18']
    state_cols_negative = ['STAIAD1', 'STAIAD2', 'STAIAD5',  'STAIAD8', 'STAIAD10', 'STAIAD11', 'STAIAD15', 'STAIAD16', 'STAIAD19', 'STAIAD20']
    trait_cols_positive = ['STAIAD22', 'STAIAD24', 'STAIAD25', 'STAIAD28', 'STAIAD29', 'STAIAD31', 'STAIAD32', 'STAIAD35', 'STAIAD37', 'STAIAD38', 'STAIAD40']
    trait_cols_negative = ['STAIAD21', 'STAIAD23', 'STAIAD26', 'STAIAD27', 'STAIAD30', 'STAIAD33', 'STAIAD34', 'STAIAD36', 'STAIAD39']
    
    df_STAI['STAI_TOT_State'] = df_STAI[state_cols_positive].sum(axis=1) + (5-df_STAI[state_cols_negative]).sum(axis=1)
    df_STAI['STAI_TOT_Trait'] = df_STAI[trait_cols_positive].sum(axis=1) + (5-df_STAI[trait_cols_negative]).sum(axis=1)
    df_STAI['STAI_TOT'] = df_STAI[['STAI_TOT_State','STAI_TOT_Trait']].sum(axis=1)
    
    return ['STAI_TOT_State','STAI_TOT_Trait','STAI_TOT']

def calculate_scores_QUIPCS(df_QUIPCS):
    """
    Calculates total scores for QUIPCS:
    
    Raw total:
      QUIPCS_TOT = sum()
    
    Suggested calculation in Derived_Variable_Definitions_and_Score_Calculations.csv:
      QUIP_Number = any(cols_A==1) + any(cols_B==1) + any(cols_C==1) + any(cols_D==1) + sum(cols_E==1)
    
    Returns a list of added column names
    """
    # Impulsive-Compulsive behaviour: Gambling, Sex, Shopping, Eating, Other, Medication
    cols_A = ['TMGAMBLE','CNTRLGMB']
    cols_B = ['TMSEX','CNTRLSEX']
    cols_C = ['TMBUY','CNTRLBUY']
    cols_D = ['TMEAT','CNTRLEAT']
    cols_E = ['TMTORACT','TMTMTACT','TMTRWD']
    cols_F = ['TMDISMED','CNTRLDSM']
    df_QUIPCS['QUIPCS_TOT'] = df_QUIPCS[cols_A].sum(axis=1)
    
    #* Calculation from Derived_Variable_Definitions_and_Score_Calculations.csv
    df_QUIPCS['QUIP_Number'] = 1*df_QUIPCS[cols_A].any(axis=1) + 1*df_QUIPCS[cols_B].any(axis=1) + 1*df_QUIPCS[cols_C].any(axis=1) + 1*df_QUIPCS[cols_D].any(axis=1) + df_QUIPCS[cols_E].sum(axis=1)
    
    return ['QUIPCS_TOT','QUIP_Number']

def calculate_scores_EPWORTH(df_EPWORTH):
    """
    Calculates total score for EPWORTH:
    
      ESS_TOT = sum()
    
    Returns a list of added column names
    """
    cols = ['ESS1','ESS2','ESS3','ESS4','ESS5','ESS6','ESS7','ESS8']
    df_EPWORTH['ESS_TOT'] = df_EPWORTH[cols].sum(axis=1)
    df_EPWORTH['ESS_Sleepy'] = df_EPWORTH['ESS_TOT'] >= 10
    return ['ESS_TOT','ESS_Sleepy']

def calculate_scores_REMSLEEP(df_REMSLEEP):
    """
    Calculates total score for REMSLEEP, according to 
    Derived_Variable_Definitions_and_Score_Calculations.csv:
      1 point each for "Yes" to the following variables:
        DRMVIVID, DRMAGRAC, DRMNOCTB, SLPLMBMV, SLPINJUR, DRMVERBL, DRMFIGHT, DRMUMV, DRMOBJFL, MVAWAKEN, DRMREMEM, SLPDSTRB
      Add a further 1 point if any of these are "Yes":
        STROKE, HETRA, PARKISM, RLS, NARCLPSY, DEPRS, EPILEPSY, BRNINFM, CNSOTH
    - If any of the previous variables are missing, then RBD score is missing.
    - Subjects with score >=5 are RBD Positive.  Subjects with score <5 are RBD Negative.
    
    Returns a list of added column names
    """
    cols_sum = ['DRMVIVID','DRMAGRAC','DRMNOCTB','SLPLMBMV','SLPINJUR','DRMVERBL',
                'DRMFIGHT','DRMUMV','DRMOBJFL','MVAWAKEN','DRMREMEM','SLPDSTRB']
    cols_neuro = ['STROKE','HETRA','PARKISM','RLS','NARCLPSY','DEPRS','EPILEPSY','BRNINFM']
    
    df_REMSLEEP['RBD_TOT'] = df_REMSLEEP[cols_sum].sum(axis=1) + 1*(df_REMSLEEP[cols_neuro].any(axis=1))
    
    # Missing values cause the total score to be missing
    m = np.isnan(df_REMSLEEP[cols_sum + cols_neuro]).any(axis=1)
    df_REMSLEEP.set_value(m,'RBD_TOT',np.nan)
    
    return ['RBD_TOT']

def calculate_scores_UPSIT(df_UPSIT):
    """
    Calculates total score for UPSIT:
    
      UPSIT_TOT = sum()
    
    Also sets the total to NaN where any individual score is missing.
    
    Could also look into the COMM column to understand other low scores 
    that might've been entered as zero, rather than missing.
    
    Returns a list of added column names
    """
    UPSIT_cols = ['UPSITBK1','UPSITBK2','UPSITBK3','UPSITBK4']
    df_UPSIT['UPSIT_TOT'] = df_UPSIT[UPSIT_cols].sum(axis=1)
    
    #* Potentially corrupted data due to missing/expired smell test booklets
    rows_equal_to_zero = (df_UPSIT[UPSIT_cols]==0).any(axis=1)
    rows_missing = np.isnan(df_UPSIT[UPSIT_cols]).any(axis=1)
    #df_UPSIT.set_value(rows_equal_to_zero | rows_missing,'UPSIT_TOT',np.nan)
    df_UPSIT.set_value(rows_missing,'UPSIT_TOT',np.nan)
    
    return ['UPSIT_TOT']

def calculate_scores_NUPDRS1(df_NUPDRS1):
    """
    Calculates total score for part 1 of the UPDRS 
    
      NP1_TOT = sum()
    
    Returns a list of added column names
    """
    cols = ['NP1COG','NP1HALL','NP1DPRS','NP1ANXS','NP1APAT','NP1DDS']
    df_NUPDRS1['NP1_TOT'] = df_NUPDRS1[cols].sum(axis=1)
    return ['NP1_TOT']

def calculate_scores_NUPDRS1p(df_NUPDRS1p):
    """
    Calculates total score for part 1p of the UPDRS 
    
      NP1p_TOT = sum()
    
    Returns a list of added column names
    """
    cols = ['NP1SLPN','NP1SLPD','NP1PAIN','NP1URIN','NP1CNST','NP1LTHD','NP1FATG']
    df_NUPDRS1p['NP1p_TOT'] = df_NUPDRS1p[cols].sum(axis=1)
    return ['NP1p_TOT']

def calculate_scores_NUPDRS2p(df_NUPDRS2p):
    """
    Calculates total score for part 2 of the UPDRS 
    
      NP2_TOT = sum()
    
    Returns a list of added column names
    """
    cols = ['NP2SPCH','NP2SALV','NP2SWAL','NP2EAT',
            'NP2DRES','NP2HYGN','NP2HWRT','NP2HOBB','NP2TURN',
            'NP2TRMR','NP2RISE','NP2WALK','NP2FREZ']
    df_NUPDRS2p['NP2_TOT'] = df_NUPDRS2p[cols].sum(axis=1)
    return ['NP2_TOT']
        
def calculate_scores_NUPDRS3(df_NUPDRS3):
    """
    Calculates total score for part 3 of the UPDRS 
    
      NP3_TOT = sum()
    
    Returns a list of added column names
    
    Note the typo in the column named PN3RIGRL (as of March 2018).
    """
    cols = ['NP3SPCH','NP3FACXP','NP3RIGN','NP3RIGRU','NP3RIGLU',
            'PN3RIGRL',
            'NP3RIGLL','NP3FTAPR','NP3FTAPL','NP3HMOVR','NP3HMOVL','NP3PRSPR','NP3PRSPL',
            'NP3TTAPR','NP3TTAPL','NP3LGAGR','NP3LGAGL','NP3RISNG','NP3GAIT','NP3FRZGT',
            'NP3PSTBL','NP3POSTR','NP3BRADY','NP3PTRMR','NP3PTRML','NP3KTRMR','NP3KTRML',
            'NP3RTARU','NP3RTALU','NP3RTARL','NP3RTALL','NP3RTALJ','NP3RTCON']
    df_NUPDRS3['NP3_TOT'] = df_NUPDRS3[cols].sum(axis=1)
    return ['NP3_TOT']

def calculate_scores_NUPDRS4(df_NUPDRS4):
    """
    Calculates total score for part 4 of the UPDRS 
    
      NP4_TOT = sum()
    
    Returns a list of added column names
    """
    cols = ['NP4WDYSK','NP4DYSKI','NP4OFF','NP4FLCTI','NP4FLCTX','NP4DYSTN']
    df_NUPDRS4['NP4_TOT'] = df_NUPDRS4[cols].sum(axis=1)
    return ['NP4_TOT']

def calculate_scores_NUPDRS_TOT(df_NUPDRS1,df_NUPDRS1p,df_NUPDRS2p,df_NUPDRS3_or_3A,df_NUPDRS4):
    """
    Calculates total score UPDRS, pre/post dose: 
      NP1_TOT + NP1p_TOT + NP2_TOT + NP3_TOT
    
    Returns the total score
    """
    UPDRS_TOT = df_NUPDRS1['NP1_TOT'] + df_NUPDRS1p['NP1p_TOT'] + df_NUPDRS2p['NP2_TOT'] + df_NUPDRS3_or_3A['NP3_TOT'] 
    # UPDRS_TOT = UPDRS_TOT + df_NUPDRS4['NP4_TOT']
    return UPDRS_TOT

def calculate_scores_PENEURO(df_PENEURO):
    """
    Calculates subscores and total score for PENEURO:
    
      PENEURO_R = sum(test results on RHS)
      PENEURO_L = sum(test results on LHS)
      PENEURO_TOT = PENEURO_R + PENEURO_L
    
    Returns a list of added column names
    """
    # Muscle strength: R and L
    cols_MSR = ['MSRARSP','MSRLRSP']
    cols_MSL = ['MSLARSP','MSLLRSP']
    cols_MS = cols_MSR + cols_MSL
    
    # Coordination: R and L
    cols_COR = ['COFNRRSP','COHSRRSP']
    cols_COL = ['COFNLRSP','COHSLRSP']
    cols_CO = cols_COR + cols_COL
    
    # Sensory: R and L
    cols_SER = ['SENRARSP','SENRLRSP']
    cols_SEL = ['SENLARSP','SENLLRSP']
    cols_SE = cols_SER + cols_SEL
    
    # Reflex: R and L
    cols_RFR = ['RFLRARSP','RFLRLRSP']
    cols_RFL = ['RFLLARSP','RFLLLRSP']
    cols_RF = cols_RFR + cols_RFL
    
    # Plantar: R and L
    cols_PR = ['PLRRRSP']
    cols_PL = ['PLRLRSP']
    cols_P = cols_PR + cols_PL
    
    cols_R = cols_MSR + cols_COR + cols_SER + cols_RFR + cols_PR
    cols_L = cols_MSL + cols_COL + cols_SEL + cols_RFL + cols_PL
    df_PENEURO['PENEURO_R'] = df_PENEURO[cols_R].sum(axis=1)
    df_PENEURO['PENEURO_L'] = df_PENEURO[cols_L].sum(axis=1)
    df_PENEURO['PENEURO_TOT'] = df_PENEURO[cols_R + cols_L].sum(axis=1)
    
    return ['PENEURO_R','PENEURO_L','PENEURO_TOT']

def calculate_scores_PECN(df_PECN):
    """
    Calculates total score for PECN:
    
      CN_TOT = sum(cols)
    
    Returns a list of added column names
    """
    # 
    cols = ['CN1RSP','CN2RSP','CN346RSP','CN5RSP','CN7RSP','CN8RSP','CN910RSP','CN11RSP','CN12RSP']
    df_PECN['CN_TOT'] = df_PECN[cols].sum(axis=1)
    return ['CN_TOT']


#*** 3.0 Tables of interest
#* Created earlier by this script: BIOLAB
#* Suggest adding MRI volumes (perhaps using GIF, or FreeSurfer, for example)

#* PPMI raw tables
tables_of_interest = {
    'PDMEDUSE':'Use_of_PD_Medication',
    'CLINDX':'Clinical_Diagnosis_and_Management',
    'PDFEAT':'PD_Features',
    'VITAL':'Vital_Signs',
    'PECN':'Neurological_Exam___Cranial_Nerves',
    'PENEURO':'General_Neurological_Exam',
    'MOCA':'Montreal_Cognitive_Assessment__MoCA_',
    'SFT':'Semantic_Fluency',
    'LNSPD':'Letter___Number_Sequencing__PD_',
    'HVLT':'Hopkins_Verbal_Learning_Test',
    'LINEORNT':'Benton_Judgment_of_Line_Orientation',
    'SDM':'Symbol_Digit_Modalities',
    'GDSSHORT':'Geriatric_Depression_Scale__Short_',
    'STAI':'State_Trait_Anxiety_Inventory',
    'QUIPCS':'QUIP_Current_Short',
    'EPWORTH':'Epworth_Sleepiness_Scale', 
    'REMSLEEP':'REM_Sleep_Disorder_Questionnaire',
    'NUPDRS1':'MDS_UPDRS_Part_I',
    'NUPDRS1p':'MDS_UPDRS_Part_I__Patient_Questionnaire',
    'NUPDRS2p':'MDS_UPDRS_Part_II__Patient_Questionnaire',
    'NUPDRS3':'MDS_UPDRS_Part_III__Post_Dose_',
    'NUPDRS4':'MDS_UPDRS_Part_IV',
    'UPSIT':'University_of_Pennsylvania_Smell_ID_Test',
    'SCOPAAUT':'SCOPA_AUT',
    'DATSCAN':'DaTscan_Imaging',
    'DATSCAN_RESULTS':'DATScan_Analysis',
    #'DATSCAN_visual':'DaTSCAN_SPECT_Visual_Interpretation_Assessment',
    'AVIMAG_meta':'AV_133_Image_Metadata',
    'AVIMAG_RESULTS':'AV_133_SBR_Results',
    'DTIROI':'DTI_Regions_of_Interest',
    'MRI':'Magnetic_Resonance_Imaging'
}

valz = list(tables_of_interest.values())
keyz = list(tables_of_interest.keys())
columns_of_interest = {
    'PDMEDUSE':['PDMEDYN','ONLDOPA','ONDOPAG','ONOTHER'],
    'CLINDX':['INFODT', 'PRIMDIAG', 'DCRTREM', 'DCRIGID', 'DCBRADY'],
    'PDFEAT':['PDDXDT','PDDXEST','DXTREMOR','DXRIGID','DXBRADY','DOMSIDE'], 
    'VITAL':['WGTKG', 'HTCM', 'TEMPC', 'SYSSUP', 'DIASUP', 'HRSUP', 'SYSSTND', 'DIASTND', 'HRSTND'],
    'PECN':['CN1RSP','CN2RSP','CN346RSP','CN5RSP','CN7RSP','CN8RSP','CN910RSP','CN11RSP','CN12RSP'],
    'PENEURO':['MSRARSP','MSRACM','MSLARSP','MSLACM','MSRLRSP','MSRLCM','MSLLRSP','MSLLCM',
               'COFNRRSP','COFNRCM','COFNLRSP','COFNLCM','COHSRRSP','COHSRCM','COHSLRSP','COHSLCM',
               'SENRARSP','SENRACM','SENLARSP','SENLACM','SENRLRSP','SENRLCM','SENLLRSP','SENLLCM',
               'RFLRARSP','RFLRACM','RFLLARSP','RFLLACM','RFLRLRSP','RFLRLCM','RFLLLRSP','RFLLLCM',
               'PLRRRSP','PLRRCM','PLRLRSP','PLRLCM'],
    'MOCA':['MCATOT',
            'MCAALTTM','MCACUBE','MCACLCKC','MCACLCKN','MCACLCKH',
            'MCALION','MCARHINO','MCACAMEL',
            'MCAFDS','MCABDS','MCAVIGIL','MCASER7',
            'MCASNTNC','MCAVFNUM','MCAVF',
            'MCAABSTR',
            'MCAREC1','MCAREC2','MCAREC3','MCAREC4','MCAREC5',
            'MCADATE','MCAMONTH','MCAYR','MCADAY','MCAPLACE','MCACITY'],
    'SFT':['DVS_SFTANIM', 'DVT_SFTANIM', 'VLTANIM', 'VLTVEG', 'VLTFRUIT'],
    'LNSPD':['LNS_TOTRAW','DVS_LNS'],
    'HVLT':['HVLTRT1', 'HVLTRT2', 'HVLTRT3', 'HVLTRDLY', 'HVLTREC', 'HVLTFPRL', 'HVLTFPUN',
            'DVT_TOTAL_RECALL', 'DVT_DELAYED_RECALL', 'DVT_RETENTION', 'DVT_RECOG_DISC_INDEX'],
    'LINEORNT':['JLO_TOTRAW', 'JLO_TOTCALC', 'DVS_JLO_MSSA', 'DVS_JLO_MSSAE'],
    'SDM':['SDMTOTAL', 'DVSD_SDM', 'DVT_SDM'],
    'GDSSHORT':['GDSSATIS','GDSDROPD','GDSEMPTY','GDSBORED','GDSGSPIR','GDSAFRAD','GDSHAPPY','GDSHLPLS','GDSHOME','GDSMEMRY','GDSALIVE','GDSWRTLS','GDSENRGY','GDSHOPLS','GDSBETER'],
    'STAI':['STAIAD1','STAIAD2','STAIAD3','STAIAD4','STAIAD5','STAIAD6','STAIAD7','STAIAD8','STAIAD9','STAIAD10',
            'STAIAD11','STAIAD12','STAIAD13','STAIAD14','STAIAD15','STAIAD16','STAIAD17','STAIAD18','STAIAD19','STAIAD20',
            'STAIAD21','STAIAD22','STAIAD23','STAIAD24','STAIAD25','STAIAD26','STAIAD27','STAIAD28','STAIAD29','STAIAD30',
            'STAIAD31','STAIAD32','STAIAD33','STAIAD34','STAIAD35','STAIAD36','STAIAD37','STAIAD38','STAIAD39','STAIAD40'],
    'QUIPCS':['TMGAMBLE','CNTRLGMB','TMSEX','CNTRLSEX','TMBUY','CNTRLBUY','TMEAT','CNTRLEAT','TMTORACT','TMTMTACT','TMTRWD','TMDISMED','CNTRLDSM'],
    'EPWORTH':['ESS1','ESS2','ESS3','ESS4','ESS5','ESS6','ESS7','ESS8'], 
    'REMSLEEP':['DRMVIVID','DRMAGRAC','DRMNOCTB','SLPLMBMV','SLPINJUR','DRMVERBL','DRMFIGHT','DRMUMV','DRMOBJFL','MVAWAKEN','DRMREMEM','SLPDSTRB','STROKE','HETRA','PARKISM','RLS','NARCLPSY','DEPRS','EPILEPSY','BRNINFM','CNSOTH'],
    'NUPDRS1':['NP1COG','NP1HALL','NP1DPRS','NP1ANXS','NP1APAT','NP1DDS'],
    'NUPDRS1p':['NP1SLPN','NP1SLPD','NP1PAIN','NP1URIN','NP1CNST','NP1LTHD','NP1FATG'],
    'NUPDRS2p':['NP2SPCH','NP2SALV','NP2SWAL','NP2EAT','NP2DRES',
                'NP2HYGN','NP2HWRT','NP2HOBB','NP2TURN','NP2TRMR','NP2RISE','NP2WALK','NP2FREZ'],
    'NUPDRS3':['NHY','NP3SPCH','NP3FACXP','NP3RIGN','NP3RIGRU','NP3RIGLU',
               'PN3RIGRL',
               'NP3RIGLL','NP3FTAPR','NP3FTAPL','NP3HMOVR','NP3HMOVL','NP3PRSPR','NP3PRSPL',
               'NP3TTAPR','NP3TTAPL','NP3LGAGR','NP3LGAGL','NP3RISNG','NP3GAIT','NP3FRZGT',
               'NP3PSTBL','NP3POSTR','NP3BRADY','NP3PTRMR','NP3PTRML','NP3KTRMR','NP3KTRML',
               'NP3RTARU','NP3RTALU','NP3RTARL','NP3RTALL','NP3RTALJ','NP3RTCON'],
    'NUPDRS3A':['ANNUAL_TIME_BTW_DOSE_NUPDRS','PD_MED_USE'],
    'NUPDRS4':['NP4WDYSK','NP4DYSKI','NP4OFF','NP4FLCTI','NP4FLCTX','NP4DYSTN'],
    'UPSIT':['UPSITBK1','UPSITBK2','UPSITBK3','UPSITBK4'],
    'SCOPAAUT':['SCAU1','SCAU2','SCAU3','SCAU4','SCAU5','SCAU6','SCAU7',
                'SCAU8','SCAU9','SCAU10','SCAU11','SCAU12','SCAU13',
                'SCAU14','SCAU15','SCAU16', 
                'SCAU17','SCAU18', 
                'SCAU19', 
                'SCAU20','SCAU21',
                'SCAU22','SCAU23', 'SCAU24','SCAU25',
                'SCAU23A','SCAU23AT','SCAU26A','SCAU26AT','SCAU26B','SCAU26BT','SCAU26C','SCAU26CT','SCAU26D','SCAU26DT'],
    'DATSCAN':['DATSCAN','DATSCNDT','VSINTRPT'],
    'DATSCAN_RESULTS':['CAUDATE_R','CAUDATE_L','PUTAMEN_R','PUTAMEN_L'],
    'AVIMAG_META':['pass_qc'],
    'AVIMAG_RESULTS':['RCAUD-S','RPUTANT-S','RPUTPOST-S','LCAUD-S','LPUTANT-S','LPUTPOST-S'],
    'MRI':['MRICMPLT','MRIDT','MRIWDTI','MRIWRSS','MRIXFRYN','MRIRSLT',
           'PDMEDYN','ONLDOPA','ONDOPAG','ONOTHER','PDMEDDT'],
    'DTIROI':['Measure','Tissue','ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','REF1','REF2']
}
# Note that DATSCAN column VSINTRPT == 1/2 (3): visual interpretation is consistent/inconsistent (NA) with evidence

#*** 3.1 Subscores and total scores: MOCA, GDSSHORT, SCOPAAUT, STAI, QUIPCS, EPWORTH, REMSLEEP, UPSIT, PENEURO, PECN, UPDRS
df_MOCA = PPMI_df[np.where([1 if n.lower()==tables_of_interest['MOCA'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['MOCA'].extend(calculate_scores_MOCA(df_MOCA))

df_GDSSHORT = PPMI_df[np.where([1 if n.lower()==tables_of_interest['GDSSHORT'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['GDSSHORT'].extend(calculate_scores_GDSSHORT(df_GDSSHORT))

df_SCOPAAUT = PPMI_df[np.where([1 if n.lower()==tables_of_interest['SCOPAAUT'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['SCOPAAUT'].extend(calculate_scores_SCOPAAUT(df_SCOPAAUT))

df_STAI = PPMI_df[np.where([1 if n.lower()==tables_of_interest['STAI'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['STAI'].extend(calculate_scores_STAI(df_STAI))

df_QUIPCS = PPMI_df[np.where([1 if n.lower()==tables_of_interest['QUIPCS'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['QUIPCS'].extend(calculate_scores_QUIPCS(df_QUIPCS))

df_EPWORTH = PPMI_df[np.where([1 if n.lower()==tables_of_interest['EPWORTH'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['EPWORTH'].extend(calculate_scores_EPWORTH(df_EPWORTH))

df_REMSLEEP = PPMI_df[np.where([1 if n.lower()==tables_of_interest['REMSLEEP'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['REMSLEEP'].extend(calculate_scores_REMSLEEP(df_REMSLEEP))

df_UPSIT = PPMI_df[np.where([1 if n.lower()==tables_of_interest['UPSIT'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['UPSIT'].extend(calculate_scores_UPSIT(df_UPSIT))

df_PENEURO = PPMI_df[np.where([1 if n.lower()==tables_of_interest['PENEURO'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['PENEURO'].extend(calculate_scores_PENEURO(df_PENEURO))

df_PECN = PPMI_df[np.where([1 if n.lower()==tables_of_interest['PECN'].lower() else 0 for n in PPMI_df_names])[0][0]]
columns_of_interest['PECN'].extend(calculate_scores_PECN(df_PECN))

df_NUPDRS1 = PPMI_df[np.where([1 if n.lower()==tables_of_interest['NUPDRS1'].lower() else 0 for n in PPMI_df_names])[0][0]]
df_NUPDRS1p = PPMI_df[np.where([1 if n.lower()==tables_of_interest['NUPDRS1p'].lower() else 0 for n in PPMI_df_names])[0][0]]
df_NUPDRS2p = PPMI_df[np.where([1 if n.lower()==tables_of_interest['NUPDRS2p'].lower() else 0 for n in PPMI_df_names])[0][0]]
df_NUPDRS3_ = PPMI_df[np.where([1 if n.lower()==tables_of_interest['NUPDRS3'].lower() else 0 for n in PPMI_df_names])[0][0]]
df_NUPDRS3A = df_NUPDRS3_.copy()
df_NUPDRS3A = df_NUPDRS3A.loc[df_NUPDRS3A.PAG_NAME=='NUPDRS3A']
df_NUPDRS3 = df_NUPDRS3_.copy()
df_NUPDRS3 = df_NUPDRS3.loc[df_NUPDRS3.PAG_NAME=='NUPDRS3']
df_NUPDRS4 = PPMI_df[np.where([1 if n.lower()==tables_of_interest['NUPDRS4'].lower() else 0 for n in PPMI_df_names])[0][0]] 

columns_of_interest['NUPDRS1'].extend(calculate_scores_NUPDRS1(df_NUPDRS1))
columns_of_interest['NUPDRS1p'].extend(calculate_scores_NUPDRS1p(df_NUPDRS1p))
columns_of_interest['NUPDRS2p'].extend(calculate_scores_NUPDRS2p(df_NUPDRS2p))
calculate_scores_NUPDRS3(df_NUPDRS3A)
columns_of_interest['NUPDRS3'].extend(calculate_scores_NUPDRS3(df_NUPDRS3))
columns_of_interest['NUPDRS4'].extend(calculate_scores_NUPDRS4(df_NUPDRS4))

# ### Data Preparation step 4: PPMIMERGE
# - 4.1 Subjects table: demographics and covariates
#   - Demographics: `RANDOM`+`SCREEN`: gender, ethnicity (HISPLAT), race (Amer/Asian/Afro/Pacific/White/NotSpecified)
#   - Covariates: 
#     - `FAMHXPD`: family history of PD (perhaps convert to binary)
#     - `MUTRSLT`: genetics (GENECAT = 1,2,3: LRRK, SNCA, GBA)
# - 4.2 Visits table - selected only (main visits and unscheduled - not including telephone visits)
#   - `RS1` (rescreen) overrides `SC`
# - 4.3 PPMIMERGE base table using `df_crossjoin(df_subjects,df_visits)`
# - 4.4 Merge tables of interest using incremental left-joins
# - 4.5 Cleanup - remove empty rows

#*** 4.1 Subjects table
df_subjects = df_RANDOM_enrolled.copy()
df_subjects['APPRDX_enrol'] = df_subjects['APPRDX']
df_subjects.drop('APPRDX',axis=1,inplace=True)

#* Add covariates *#

#* Family History of PD
# N 1st-degree relatives with PD = sum(Mum, Dad, Full Siblings, Kids)
# N 2nd-degree relatives with PD = sum(Half Siblings, Grandparents, Aunts, Uncles)
df_FAMHXPD = PPMI_df[np.where([1 if n=="Family_History__PD_" else 0 for n in PPMI_df_names])[0][0]].copy()
df_tmp = df_FAMHXPD.copy()
df_tmp['FAMHXPD_N1stDegree'] = df_tmp[['BIOMOMPD','BIODADPD','FULSIBPD','KIDSPD']].fillna(0).sum(axis=1)
df_tmp['FAMHXPD_N2ndDegree'] = df_tmp[['HAFSIBPD','MAGPARPD','PAGPARPD','MATAUPD','PATAUPD']].fillna(0).sum(axis=1)
# Remove duplicates: keeping most recent LAST_UPDATE
df_tmp = df_tmp.sort_values(by=['PATNO','LAST_UPDATE'])
df_tmp = df_tmp[~df_tmp.duplicated(subset='PATNO',keep='last')]

#* Genetic mutation category: GENECAT = 1,2,3 (LRRK2, GBA, SNCA) from MUTRSLT 
df_MUTRSLT = PPMI_df[np.where([1 if n=="Genetic_Testing_Results" else 0 for n in PPMI_df_names])[0][0]].copy()
# Remove duplicates: keeping first REC_ID (there were fewer NaN in the LRRKCD column)
df_MUTRSLT = df_MUTRSLT.sort_values(by=['PATNO','REC_ID','GENECAT','LRRKCD','MUTRSLT','LAST_UPDATE'])
df_MUTRSLT = df_MUTRSLT[~df_MUTRSLT.duplicated(subset=['PATNO'],keep='first')]
df_MUTRSLT = df_MUTRSLT[['PATNO','GENECAT','MUTRSLT','LRRKCD']]

#* Add GENECAT,'LRRKCD','MUTRSLT' columns and new FAMHXPD columns
df_subjects = df_subjects.merge(df_MUTRSLT,on=['PATNO'],how='left')
df_subjects = df_subjects.merge(df_tmp[['PATNO','FAMHXPD_N1stDegree','FAMHXPD_N2ndDegree']],on=['PATNO'],how='left')

#* Add disease duration: ENROLLDT-PDFEAT[PDDXDT]
PDDURAT = 'Years_since_DX_bl' # label
df_PDFEAT_tmp = PPMI_df[np.where([1 if n==tables_of_interest['PDFEAT'] else 0 for n in PPMI_df_names])[0][0]].copy()
df_PDFEAT_tmp = df_PDFEAT_tmp[['PATNO']+columns_of_interest['PDFEAT']]
df_PDFEAT_tmp = df_PDFEAT_tmp.merge(df_subjects[['PATNO','ENROLLDT']],on=['PATNO'],how='left')
df_PDFEAT_tmp[PDDURAT] = (pd.to_datetime(df_PDFEAT_tmp['ENROLLDT']) - pd.to_datetime(df_PDFEAT_tmp['PDDXDT'])).apply(lambda x: round(x.days/daysInAYear,3))
df_subjects = df_subjects.merge(df_PDFEAT_tmp[['PATNO',PDDURAT]],on=['PATNO'],how='left')

df_subjects = df_subjects[['PATNO', 'APPRDX_enrol', PDDURAT, 'GENDER', 'BIRTHDT',
                           'GENECAT', 'MUTRSLT', 'LRRKCD', # genetics info
                           'FAMHXPD_N1stDegree', 'FAMHXPD_N2ndDegree', 
                           'EDUCYRS','HANDED',
                           'HISPLAT', 'RAINDALS', 'RAASIAN', 'RABLACK', 'RAHAWOPI', 'RAWHITE', 'RANOS',
                           'ENROLLDT', 'CONSNTDT'
                          ]]

#*** 4.2 Visits table (Selected visits only - NB: All visits can be obtained as EVENT_IDs_all = df_EVENT_ID.CODE)
EVENT_IDs = ['SC', 'SCBL', 'RS1', 'BL', 'V01', 'V02', 'V03', 'V04', 'V05', 'V06', 'V07', 'V08', 'V09', 'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'ST', 'PW', 'U01', 'U02', 'U03', 'U04', 'U05', 'U06', 'UT1']
df_visits = pd.DataFrame(data={'EVENT_ID':EVENT_IDs})

#*** 4.3 Create massive dataframe for all individuals and visits (can downsize it later)
df_base = df_crossjoin(df_subjects,df_visits)
df_base.reset_index(drop=True,inplace=True)


# ### Data Preparation step 4 (continue): PPMIMERGE
# - 4.4 Merge tables of interest using incremental left-joins
# - 4.5 Cleanup - remove empty rows
df_PPMIMERGE = df_base.copy()

key = 'PDMEDUSE'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

#* One way to check that values actually merged into PPMIMERGE:
# df_PPMIMERGE.loc[df_PPMIMERGE.PATNO.isin(df_.PATNO),['INFODT_'+key]+cols]

key = 'CLINDX'
cols = ['PATNO','EVENT_ID'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'VITAL'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

# CSF etc.
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_BIOLAB_SemiWide,on=Keys,how='left',suffixes=('','_BIOANALYS'))

key = 'MOCA'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'REMSLEEP'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'EPWORTH'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'UPSIT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'SCOPAAUT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'NUPDRS1'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'NUPDRS1p'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'NUPDRS2p'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'NUPDRS3'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
#******************** BEWARE!!!!!!!!!!
#* Here (and below for NUPDRS3A), I fixed a typo in the PPMI data: PN3RIGRL => NP3RIGRL
cols_renamed = [s.replace('PN3','NP3') for s in cols] 
df_ = df_NUPDRS3[cols].copy()
df_.rename(columns=dict(zip(cols,cols_renamed)),inplace=True)
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'NUPDRS3A'
cols = columns_of_interest['NUPDRS3'] + columns_of_interest[key]
df_ = df_NUPDRS3A[Keys + ['INFODT'] + cols].copy()
#* Fix typo (PN3=>NP3) & add suffix to 3A (post-drug) scores
cols_renamed = [s.replace('PN3','NP3') + '_'+key for s in cols] 
df_.rename(columns=dict(zip(cols,cols_renamed)),inplace=True)
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

#*** Fix typo in UPDRS3/3A column: PN should be NP (PN3RIGRL => NP3RIGRL)
columns_of_interest['NUPDRS3'] = [s.replace('PN3','NP3') for s in columns_of_interest['NUPDRS3']]

key = 'NUPDRS4'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

#*** Total UPDRS score: NP1_TOT + NP2_TOT + NP3_TOT/NP3A_TOT
#* Sum across columns
df_PPMIMERGE['UPDRS_TOT'] = df_PPMIMERGE[['NP1_TOT','NP2_TOT','NP3_TOT']].sum(axis=1)
df_PPMIMERGE['UPDRS_TOT_A'] = df_PPMIMERGE[['NP1_TOT','NP2_TOT','NP3_TOT_NUPDRS3A']].sum(axis=1)

key = 'GDSSHORT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'QUIPCS'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'SFT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'LNSPD'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'LINEORNT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'SDM'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'STAI'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'HVLT'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'PECN'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'PENEURO'
cols = ['PATNO','EVENT_ID','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][cols]
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_,on=Keys,how='left',suffixes=('','_'+key))

key = 'DATSCAN'
df_DATSCAN = PPMI_df[np.where([1 if n.lower()=="DATScan_Analysis".lower() else 0 for n in PPMI_df_names])[0][0]]  [['PATNO','EVENT_ID',
    'CAUDATE_R', 'CAUDATE_L',
    'PUTAMEN_R', 'PUTAMEN_L']].copy()
df_DATSCAN_meta = PPMI_df[np.where([1 if n.lower()=="DaTscan_Imaging".lower() else 0 for n in PPMI_df_names])[0][0]]  [['PATNO','EVENT_ID','INFODT',
    'DATSCAN', 'DATSCNDT', 'SCNLOC', 'SCNINJCT', 'DATXFRYN', 'VSINTRPT']].copy()
df_DATSCAN = pd.merge(df_DATSCAN,df_DATSCAN_meta)
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_DATSCAN,on=Keys,how='left',suffixes=('','_'+key))

key = 'AVIMAG'
df_AVIMAG_Imaging = PPMI_df[np.where([1 if n.lower()=="AV_133_Imaging".lower() else 0 for n in PPMI_df_names])[0][0]]  [['PATNO','EVENT_ID','INFODT']].copy()
df_AVIMAG_Meta = PPMI_df[np.where([1 if n.lower()=="AV_133_Image_Metadata".lower() else 0 for n in PPMI_df_names])[0][0]]  [['PATNO','EVENT_ID',
    'scan_date', 'pass_qc', 'scan_quality_rating_pet', 'ligand']].copy()
df_AVIMAG_RESULTS = PPMI_df[np.where([1 if n.lower()=="AV_133_SBR_Results".lower() else 0 for n in PPMI_df_names])[0][0]]  [['SUBJECT_NUMBER', 'VISIT_NUMBER', 'SCAN_DATE', 'TIMEPOINT',
    'RCAUD-S', 'RPUTANT-S', 'RPUTPOST-S', 'LCAUD-S', 'LPUTANT-S', 'LPUTPOST-S']].copy()
df_AVIMAG_RESULTS.rename(index=str, 
                         columns={"SUBJECT_NUMBER": "PATNO", "VISIT_NUMBER": "EVENT_ID"}, 
                         inplace=True)
df_AVIMAG_RESULTS['PATNO'] = df_AVIMAG_RESULTS['PATNO'].apply(lambda x: str(x))
#"SCAN_DATE": "INFODT"
df_AVIMAG = pd.merge(pd.merge(df_AVIMAG_Imaging,df_AVIMAG_RESULTS,on=['PATNO', 'EVENT_ID']),df_AVIMAG_Meta,on=['PATNO', 'EVENT_ID'])
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_AVIMAG,on=Keys,how='left',suffixes=('','_'+key))

key = 'MRI'
cols = ['PATNO','EVENT_ID','INFODT','MRIDT',
        'MRICMPLT', 'MRIWDTI', 'MRIWRSS', 'MRIXFRYN', 'MRIRSLT',
        'PDMEDDT', 'COMM']
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_MRI = PPMI_df[n][cols].copy()
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_MRI,on=Keys,how='left',suffixes=('','_'+key))

key = 'DTIROI'
cols = ['PATNO','INFODT'] + columns_of_interest[key]
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_DTIROI = PPMI_df[n][cols].copy()

#*** Unstack the Measure column: E1, E2, E3, FA
#* NOTE: occasionally there are multiple rows per visit
#        e.g., PATNO==3166 has three values per measure for INFODT=='Sep-11'
#  I take the mean value
Keys_DTIROI = ['PATNO','INFODT']
df_DTIROI_E1 = df_DTIROI[df_DTIROI.Measure=='E1']
df_DTIROI_E2 = df_DTIROI[df_DTIROI.Measure=='E2']
df_DTIROI_E3 = df_DTIROI[df_DTIROI.Measure=='E3']
df_DTIROI_FA = df_DTIROI[df_DTIROI.Measure=='FA']
df_DTIROI_E1_ = df_DTIROI_E1[Keys_DTIROI].drop_duplicates()
df_DTIROI_E2_ = df_DTIROI_E2[Keys_DTIROI].drop_duplicates()
df_DTIROI_E3_ = df_DTIROI_E3[Keys_DTIROI].drop_duplicates()
df_DTIROI_FA_ = df_DTIROI_FA[Keys_DTIROI].drop_duplicates()
df_DTIROI_E1_.reset_index(inplace=True,drop=True)
df_DTIROI_E2_.reset_index(inplace=True,drop=True)
df_DTIROI_E3_.reset_index(inplace=True,drop=True)
df_DTIROI_FA_.reset_index(inplace=True,drop=True)
#*** Calculate the means for each Measure, and rename the columns
cols0 = ['ROI1','ROI2','ROI3','ROI4','ROI5','ROI6','REF1','REF2']
prefix = 'DTI_E1_'
for k in range(0,df_DTIROI_E1_.shape[0]):
    p = df_DTIROI_E1_.loc[k,Keys_DTIROI]
    tmp = df_DTIROI_E1.loc[(df_DTIROI_E1[Keys_DTIROI]==p).sum(axis=1)==2]
    for c in cols0:
        df_DTIROI_E1_.loc[k,prefix+c] = tmp[c].mean(axis=0)
prefix = 'DTI_E2_'
for k in range(0,df_DTIROI_E2_.shape[0]):
    p = df_DTIROI_E2_.loc[k,Keys_DTIROI]
    tmp = df_DTIROI_E2.loc[(df_DTIROI_E2[Keys_DTIROI]==p).sum(axis=1)==2]
    for c in cols0:
        df_DTIROI_E2_.loc[k,prefix+c] = tmp[c].mean(axis=0)
prefix = 'DTI_E3_'
for k in range(0,df_DTIROI_E3_.shape[0]):
    p = df_DTIROI_E3_.loc[k,Keys_DTIROI]
    tmp = df_DTIROI_E3.loc[(df_DTIROI_E3[Keys_DTIROI]==p).sum(axis=1)==2]
    for c in cols0:
        df_DTIROI_E3_.loc[k,prefix+c] = tmp[c].mean(axis=0)
prefix = 'DTI_FA_'
for k in range(0,df_DTIROI_FA_.shape[0]):
    p = df_DTIROI_FA_.loc[k,Keys_DTIROI]
    tmp = df_DTIROI_FA.loc[(df_DTIROI_FA[Keys_DTIROI]==p).sum(axis=1)==2]
    for c in cols0:
        df_DTIROI_FA_.loc[k,prefix+c] = tmp[c].mean(axis=0)
df_DTIROI_unstacked = df_DTIROI_FA_.merge(df_DTIROI_E1_,on=Keys_DTIROI,how='left').merge(df_DTIROI_E2_,on=Keys_DTIROI,how='left').merge(df_DTIROI_E3_,on=Keys_DTIROI,how='left')
#* Merge with MRI to get EVENT_ID
df_DTIROI_ = pd.merge(df_MRI.loc[df_MRI.MRIWDTI==1,['PATNO','EVENT_ID','MRIDT']],df_DTIROI_unstacked,left_on=['PATNO','MRIDT'],right_on=['PATNO','INFODT'],how='right')

#* Calculate axial/radial/mean diffusivity and FA
# NOTE: This assumes that Measures E1/E2/E3 are DTI eigenvalues,
#       but they might not be, since FA doesn't match FA_calc
#       (see commented-out `unexplained_factor` below)
def calc_anisotropy(x,coll='ROI1'):
    l1 = x['DTI_E1_' + coll]
    l2 = x['DTI_E2_' + coll]
    l3 = x['DTI_E3_' + coll]
    AD = l1
    RD = (l2 + l3)/2
    MD = (l1 + l2 + l3)/3
    FA = ( (1.0/2.0) * ((l1-l2)**2.0 + (l1-l3)**2.0 + (l2-l3)**2.0) / (l1**2.0 + l2**2.0 + l3**2.0) ) ** (1/2.0)
    #unexplained_factor = (2/3.0)**(1/2.0)
    #FA = unexplained_factor * FA
    #FA = (2/3.0)**(1/2.0) * ( (l1**2.0 + l2**2.0 + l3**2.0) - (l1*l2 + l1*l3 + l2*l3) )**(1/2.0) 
    return (FA,AD,RD,MD)

for k in range(0,len(cols0)):
    coll = cols0[k]
    (FA,AD,RD,MD) = calc_anisotropy(x=df_DTIROI_,coll=coll)
    df_DTIROI_['DTI_FA_'+coll+'_calc'] = FA
    df_DTIROI_['DTI_AD_'+coll+'_calc'] = AD
    df_DTIROI_['DTI_RD_'+coll+'_calc'] = RD
    df_DTIROI_['DTI_MD_'+coll+'_calc'] = MD

plt.plot(df_DTIROI_.DTI_FA_REF1,df_DTIROI_.DTI_FA_REF1_calc,'.')
plt.plot(plt.gca().get_xlim(),plt.gca().get_xlim(),label='Reference')
plt.legend()
plt.xlabel('FA PPMI')
plt.ylabel('FA calculated (with fudge factor)')
plt.show()
#* Merge PPMIMERGE and DTI
df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_DTIROI_,on=Keys,how='left',suffixes=('','_DTI'))

# ## MRI volumes: you'll have to add this information manually
include_MRI_volumes = False
#* Here's my code, which used the parcellation results from the Geodesic Information Flows algorithm,
#* along with meta data (mostly image ID) from a couple of IDA searches
if include_MRI_volumes:
    ## * Preprocessed structural MRI: T1-Anatomical
    #******************** FIX ME *******************
    PPMI_T1_Anatomical_IDASEARCH_csv = os.path.join("/Users/noxtoby/Documents/Research/UCLPOND/Projects/201803-PPMI-EBM/Data/","PPMI-T1-Anatomical_4_06_2018.csv")
    #PPMI_T1_Anatomical_IDASEARCH_csv = os.path.join("path/to","PPMI-T1-Anatomical_IDASEARCH.csv")
    #****************** END FIX ME *****************
    df_MRI_ref1 = pd.read_csv(PPMI_T1_Anatomical_IDASEARCH_csv)
    df_MRI_ref1.rename(index=str, columns={"Image Data ID": "ImageID", "Subject": "PATNO"},inplace=True)
    df_MRI_ref1.drop(["Type","Modality","Format","Downloaded"],axis=1,inplace=True)
    #* Visit code as text, rather than as a number
    #******************** FIX ME *******************
    another_IDA_search_results_csv_for_some_reason = "/Users/noxtoby/Documents/Research/UCLPOND/Projects/201803-PPMI-EBM/Data/idaSearch_4_06_2018.csv"
    #****************** END FIX ME *****************
    df_MRI_ref0 = pd.read_csv(another_IDA_search_results_csv_for_some_reason)
    df_MRI_ref0.rename(index=str, columns={"Image ID": "ImageID", "Subject ID": "PATNO", "Visit": "Visit_desc"},inplace=True)
    df_MRI_ref0.drop(["Imaging Protocol","Modality"],axis=1,inplace=True)
    #* Tidy and merge
    df_MRI_ref0 = df_MRI_ref0[['PATNO', 'ImageID', 'Visit_desc', 'Study Date', 'Research Group', 'Description', ]]
    df_MRI_ref1 = df_MRI_ref1[['PATNO', 'ImageID', 'Visit', 'Acq Date' ]]
    df_MRI_ref = pd.merge(df_MRI_ref0,df_MRI_ref1,on=['PATNO','ImageID'])[['PATNO', 'ImageID', 'Visit', 'Visit_desc', 'Study Date', 'Research Group', 'Description']]
    df_MRI_ref['PATNO'] = df_MRI_ref['PATNO'].apply(lambda x: str(x))
    
    #* Generate EVENT_ID
    df_MRI_ref["EVENT_ID"] = df_MRI_ref.Visit_desc.map({
        'Screening':'SC', 'Baseline':"BL", 
        'Month 12':'V04', 'Month 24':'V06', 'Month 48':'V10',
        'Symptomatic Therapy':'ST', 'Premature Withdrawal':'PW',
        'Unscheduled Visit 01':'UT01', 'Unscheduled Visit 02':'UT02'})
    
    #* GIF volumes
    #******************** FIX ME *******************
    df_MRI_ROIvols = pd.read_csv("/Users/noxtoby/Documents/Research/UCLPOND/Projects/201803-PPMI-EBM/Data/PPMI-T1-Anatomical_GIF3Vols.csv")
    #****************** END FIX ME *****************
    ROI_labels = df_MRI_ROIvols.columns[2:]
    ext = '.nii'
    # Extract Image IDs from filenames in GIF CSV
    for k in range(len(df_MRI_ROIvols.filename)):
        SI = df_MRI_ROIvols.filename[k].strip(ext).split("_")
        series_id = SI[-2][1:]
        image_id = SI[-1][1:]
        df_MRI_ROIvols.set_value(k,"PATNO",SI[1])
        df_MRI_ROIvols.set_value(k,"SeriesID",series_id)
        df_MRI_ROIvols.set_value(k,"ImageID",image_id)
    df_MRI_ROIvols['PATNO'] = df_MRI_ROIvols['PATNO'].apply(lambda x: str(x))
    df_MRI_ROIvols.ImageID = pd.to_numeric(df_MRI_ROIvols.ImageID,errors='coerce')
    df_MRI_ROIvols.SeriesID = pd.to_numeric(df_MRI_ROIvols.SeriesID,errors='coerce')
    df_MRI_ROIvols.rename(columns={'filename':'MRI_filename'},inplace=True)
    # Optionally Reorder columns
    reorder = False
    ROI_labels_GIF2 = [
        '3rd Ventricle','3rd Ventricle (Posterior part)','4th Ventricle','5th Ventricle',
        'Right Inf Lat Vent','Left Inf Lat Vent','Right Lateral Ventricle','Left Lateral Ventricle',
        'Right Ventral DC','Left Ventral DC','Pons','Brain Stem',
        'Left Basal Forebrain','Right Basal Forebrain',
        'Right Accumbens Area','Left Accumbens Area','Right Amygdala','Left Amygdala',
        'Right Caudate','Left Caudate',
        'Right Hippocampus','Left Hippocampus',
        'Right Pallidum','Left Pallidum','Right Putamen','Left Putamen','Right Thalamus Proper','Left Thalamus Proper',
        'Right Cun cuneus','Left Cun cuneus',
        'Right Cerebellum Exterior','Left Cerebellum Exterior',
        'Cerebellar Vermal Lobules I-V','Cerebellar Vermal Lobules VI-VII','Cerebellar Vermal Lobules VIII-X',
        'Right Cerebellum White Matter','Left Cerebellum White Matter','Right Cerebral Exterior','Left Cerebral Exterior',
        'Right Cerebral White Matter','Left Cerebral White Matter',
        'Right ACgG anterior cingulate gyrus','Left ACgG anterior cingulate gyrus',
        'Right AIns anterior insula','Left AIns anterior insula',
        'Right AOrG anterior orbital gyrus','Left AOrG anterior orbital gyrus',
        'Right AnG angular gyrus','Left AnG angular gyrus','Right Calc calcarine cortex','Left Calc calcarine cortex',
        'Right CO central operculum','Left CO central operculum',
        'Right Ent entorhinal area','Left Ent entorhinal area','Right FO frontal operculum','Left FO frontal operculum',
        'Right FRP frontal pole','Left FRP frontal pole','Right FuG fusiform gyrus','Left FuG fusiform gyrus',
        'Right GRe gyrus rectus','Left GRe gyrus rectus',
        'Right IOG inferior occipital gyrus','Left IOG inferior occipital gyrus',
        'Right ITG inferior temporal gyrus','Left ITG inferior temporal gyrus',
        'Right LiG lingual gyrus','Left LiG lingual gyrus',
        'Right LOrG lateral orbital gyrus','Left LOrG lateral orbital gyrus',
        'Right MCgG middle cingulate gyrus','Left MCgG middle cingulate gyrus',
        'Right MFC medial frontal cortex','Left MFC medial frontal cortex',
        'Right MFG middle frontal gyrus','Left MFG middle frontal gyrus',
        'Right MOG middle occipital gyrus','Left MOG middle occipital gyrus',
        'Right MOrG medial orbital gyrus','Left MOrG medial orbital gyrus',
        'Right MPoG postcentral gyrus medial segment','Left MPoG postcentral gyrus medial segment',
        'Right MPrG precentral gyrus medial segment','Left MPrG precentral gyrus medial segment',
        'Right MSFG superior frontal gyrus medial segment','Left MSFG superior frontal gyrus medial segment',
        'Right MTG middle temporal gyrus','Left MTG middle temporal gyrus',
        'Right OCP occipital pole','Left OCP occipital pole',
        'Right OFuG occipital fusiform gyrus','Left OFuG occipital fusiform gyrus',
        'Right OpIFG opercular part of the inferior frontal gyrus','Left OpIFG opercular part of the inferior frontal gyrus',
        'Right OrIFG orbital part of the inferior frontal gyrus','Left OrIFG orbital part of the inferior frontal gyrus',
        'Right PCgG posterior cingulate gyrus','Left PCgG posterior cingulate gyrus',
        'Right PCu precuneus','Left PCu precuneus','Right PHG parahippocampal gyrus','Left PHG parahippocampal gyrus',
        'Right PIns posterior insula','Left PIns posterior insula','Right PO parietal operculum','Left PO parietal operculum',
        'Right PoG postcentral gyrus','Left PoG postcentral gyrus',
        'Right POrG posterior orbital gyrus','Left POrG posterior orbital gyrus',
        'Right PP planum polare','Left PP planum polare','Right PrG precentral gyrus','Left PrG precentral gyrus',
        'Right PT planum temporale','Left PT planum temporale','Right SCA subcallosal area','Left SCA subcallosal area',
        'Right SFG superior frontal gyrus','Left SFG superior frontal gyrus',
        'Right SMC supplementary motor cortex','Left SMC supplementary motor cortex',
        'Right SMG supramarginal gyrus','Left SMG supramarginal gyrus',
        'Right SOG superior occipital gyrus','Left SOG superior occipital gyrus',
        'Right SPL superior parietal lobule','Left SPL superior parietal lobule',
        'Right STG superior temporal gyrus','Left STG superior temporal gyrus',
        'Right TMP temporal pole','Left TMP temporal pole',
        'Right TrIFG triangular part of the inferior frontal gyrus',
        'Left TrIFG triangular part of the inferior frontal gyrus',
        'Right TTG transverse temporal gyrus','Left TTG transverse temporal gyrus',
        'TIV',
        'Right vessel','Left vessel','Optic Chiasm',
        'Right Lesion','Left Lesion',
        'Background and skull','Non-ventricular CSF']
    ROI_labels_GIF3 = [
        '3rd Ventricle','4th Ventricle','5th Ventricle',
        'Right Lateral Ventricle','Left Lateral Ventricle','Right Inf Lat Vent','Left Inf Lat Vent',
        'Right Accumbens Area','Left Accumbens Area','Right Amygdala','Left Amygdala',
        'Pons','Brain Stem','Right Ventral DC','Left Ventral DC','Left Basal Forebrain','Right Basal Forebrain',
        'Right Caudate','Left Caudate','Right Putamen','Left Putamen','Right Thalamus Proper','Left Thalamus Proper',
        'Right Cerebellum Exterior','Left Cerebellum Exterior','Right Cerebellum White Matter','Left Cerebellum White Matter',
        'Right Cerebral Exterior','Left Cerebral Exterior','3rd Ventricle (Posterior part)',
        'Right Hippocampus','Left Hippocampus',
        'Right Pallidum','Left Pallidum',
        'Cerebellar Vermal Lobules I-V','Cerebellar Vermal Lobules VI-VII','Cerebellar Vermal Lobules VIII-X',
        'Right Ventricular Lining','Left Ventricular Lining','Optic Chiasm',
        'Right Temporal White Matter','Right Insula White Matter','Right Cingulate White Matter','Right Frontal White Matter',
        'Right Occipital White Matter','Right Parietal White Matter',
        'Corpus Callosum',
        'Left Temporal White Matter','Left Insula White Matter','Left Cingulate White Matter','Left Frontal White Matter',
        'Left Occipital White Matter','Left Parietal White Matter',
        'Right Claustrum','Left Claustrum','Right ACgG anterior cingulate gyrus','Left ACgG anterior cingulate gyrus',
        'Right AIns anterior insula','Left AIns anterior insula','Right AOrG anterior orbital gyrus',
        'Left AOrG anterior orbital gyrus','Right AnG angular gyrus','Left AnG angular gyrus','Right Calc calcarine cortex',
        'Left Calc calcarine cortex','Right CO central operculum','Left CO central operculum','Right Cun cuneus','Left Cun cuneus',
        'Right Ent entorhinal area','Left Ent entorhinal area','Right FO frontal operculum','Left FO frontal operculum',
        'Right FRP frontal pole','Left FRP frontal pole','Right FuG fusiform gyrus','Left FuG fusiform gyrus',
        'Right GRe gyrus rectus','Left GRe gyrus rectus','Right IOG inferior occipital gyrus','Left IOG inferior occipital gyrus',
        'Right ITG inferior temporal gyrus','Left ITG inferior temporal gyrus','Right LiG lingual gyrus','Left LiG lingual gyrus',
        'Right LOrG lateral orbital gyrus','Left LOrG lateral orbital gyrus','Right MCgG middle cingulate gyrus',
        'Left MCgG middle cingulate gyrus','Right MFC medial frontal cortex','Left MFC medial frontal cortex',
        'Right MFG middle frontal gyrus','Left MFG middle frontal gyrus','Right MOG middle occipital gyrus',
        'Left MOG middle occipital gyrus','Right MOrG medial orbital gyrus','Left MOrG medial orbital gyrus',
        'Right MPoG postcentral gyrus medial segment','Left MPoG postcentral gyrus medial segment',
        'Right MPrG precentral gyrus medial segment','Left MPrG precentral gyrus medial segment',
        'Right MSFG superior frontal gyrus medial segment','Left MSFG superior frontal gyrus medial segment',
        'Right MTG middle temporal gyrus','Left MTG middle temporal gyrus','Right OCP occipital pole','Left OCP occipital pole',
        'Right OFuG occipital fusiform gyrus','Left OFuG occipital fusiform gyrus',
        'Right OpIFG opercular part of the inferior frontal gyrus','Left OpIFG opercular part of the inferior frontal gyrus',
        'Right OrIFG orbital part of the inferior frontal gyrus','Left OrIFG orbital part of the inferior frontal gyrus',
        'Right PCgG posterior cingulate gyrus','Left PCgG posterior cingulate gyrus','Right PCu precuneus','Left PCu precuneus',
        'Right PHG parahippocampal gyrus','Left PHG parahippocampal gyrus','Right PIns posterior insula',
        'Left PIns posterior insula','Right PO parietal operculum','Left PO parietal operculum','Right PoG postcentral gyrus',
        'Left PoG postcentral gyrus','Right POrG posterior orbital gyrus','Left POrG posterior orbital gyrus',
        'Right PP planum polare','Left PP planum polare','Right PrG precentral gyrus','Left PrG precentral gyrus',
        'Right PT planum temporale','Left PT planum temporale','Right SCA subcallosal area','Left SCA subcallosal area',
        'Right SFG superior frontal gyrus','Left SFG superior frontal gyrus','Right SMC supplementary motor cortex',
        'Left SMC supplementary motor cortex','Right SMG supramarginal gyrus','Left SMG supramarginal gyrus',
        'Right SOG superior occipital gyrus','Left SOG superior occipital gyrus','Right SPL superior parietal lobule',
        'Left SPL superior parietal lobule','Right STG superior temporal gyrus','Left STG superior temporal gyrus',
        'Right TMP temporal pole','Left TMP temporal pole','Right TrIFG triangular part of the inferior frontal gyrus',
        'Left TrIFG triangular part of the inferior frontal gyrus','Right TTG transverse temporal gyrus',
        'Left TTG transverse temporal gyrus',
        'Right vessel','Left vessel','Right Lesion','Left Lesion',
        'Non-Brain Outer Tissue','NonBrain Low','NonBrain Mid','NonBrain High','Non-ventricular CSF'
        ]
    
    if reorder:
        df_MRI_ROIvols = df_MRI_ROIvols[['PATNO','ImageID','SeriesID','filenames'] + ROI_labels_GIF3]
    
    #* Merge MRI and GIF: for EVENT_ID
    df_MRI_vols = pd.merge(df_MRI_ROIvols,df_MRI_ref[['EVENT_ID','PATNO','ImageID']],on=['PATNO','ImageID'],how='left')
    #* Merge PPMIMERGE and MRI
    df_PPMIMERGE = pd.merge(df_PPMIMERGE,df_MRI_vols,on=Keys,how='left',suffixes=('','_MRI'))
else:
    print('ALERT: No MRI. If you want to include regional brain volumes (e.g., from FreeSurfer), here''s where you could be doing it.')

# # Optional: Derived Variable Definitions and Score Calculations
# See the spreadsheet `Derived_Variable_Definitions_and_Score_Calculations.csv` 
# by Christopher S. Coffey, Chelsea J. Caspell-Garcia, and Eric D. Foster, University of Iowa [eric-foster@uiowa.edu](mailto:eric-foster@uiowa.edu)
# 
# I've done the ones that are ~~struckout~~
# 
# - ~~Disease duration at enrollment: ENROLLDT minus PDDXDT~~
# - ~~TD/PIGD classification (UPDRS):~~
#   - Tremor score = mean of NP2TRMR, NP3PTRMR, NP3PTRML, NP3KTRMR, NP3KTRML, NP3RTARU, NP3RTALU, NP3RTARL, NP3RTALL, NP3RTALJ, NP3RTCON.
#   - PIGD score = mean of NP2WALK, NP2FREZ, NP3GAIT, NP3FRZGT, NP3PSTBL.
#   - Ratio = Tremor score / PIGD score.
#     - Classification TD: if ratio >= 1.15, OR if PIGD score = 0 and Tremor score > 0
#     - Classification PIGD: if ratio <= 0.9
#     - Indeterminate: if 0.9 < ratio < 1.15, OR if Tremor score and PIGD score = 0
# 
# I haven't done these ones:
# 
# - Genotypes from BIOANALYS: pivot table on TESTNAME and TESTVALUE for
#   - SNCA rs3910105: TESTNAME = "rs3910105"
#   - SNCA rs356181:  TESTNAME = "rs356181"
# - MAPT H1/H2
#   1. Extract the following 5 SNPs from the genetics data and count the number of minor alleles in each SNP:
#     - rs17652121 (minor allele = C)
#     - rs8070723 (minor allele = G)
#     - rs1052587 (minor allele = C)
#     - rs16940799 (minor allele = C)
#     - rs17652748 (minor allele = T)
#   2. Genotype
#     - H1/H1: sum of minor alleles for all 5 SNPs is <5
#     - H1/H2: sum = 5
#     - H2/H2: sum = 10
df_PPMIMERGE['TD_score'] = df_PPMIMERGE[['NP2TRMR', 'NP3PTRMR', 'NP3PTRML', 'NP3KTRMR', 'NP3KTRML', 'NP3RTARU', 'NP3RTALU', 'NP3RTARL', 'NP3RTALL', 'NP3RTALJ', 'NP3RTCON']].mean(axis=1)
df_PPMIMERGE['PIGD_score'] = df_PPMIMERGE[['NP2WALK', 'NP2FREZ', 'NP3GAIT', 'NP3FRZGT', 'NP3PSTBL']].mean(axis=1)
df_PPMIMERGE['TD_PIGD_ratio'] = df_PPMIMERGE.TD_score / df_PPMIMERGE.PIGD_score.replace({ 0 : np.nan })
TD = ( df_PPMIMERGE.TD_PIGD_ratio >= 1.15 ) | ( (df_PPMIMERGE.TD_score > 0)&(df_PPMIMERGE.PIGD_score == 0) )
PIGD = (df_PPMIMERGE.TD_PIGD_ratio <= 0.90)
Indeterminate = ( df_PPMIMERGE.TD_PIGD_ratio.between(0.9,1.15,inclusive=False) ) | ( (df_PPMIMERGE.TD_score==0) & (df_PPMIMERGE.PIGD_score==0) )
TD_PIGD_class = pd.Series(np.where(TD, 'TD', ''))
TD_PIGD_class.loc[PIGD] = 'PIGD'
TD_PIGD_class.loc[Indeterminate] = 'Indeterminate'
df_PPMIMERGE['TD_PIGD_class'] = TD_PIGD_class

# ## Check for missing INFODT
# 1. Replace missing with available ones from other tables
#   - Take the earliest/latest/average/consensus(mode?)
#   - Define a precedence ranking and take the first available
# 2. Simultaneously convert new `INFODT_` to datetime

#*** This could be the ranking order
INFODT_cols = ['INFODT',
               'INFODT_CLINDX','INFODT_VITAL','INFODT_BIOANALYS','INFODT_MOCA','INFODT_REMSLEEP',
               'INFODT_EPWORTH','INFODT_UPSIT',
               'INFODT_NUPDRS1','INFODT_NUPDRS2p','INFODT_NUPDRS3','INFODT_NUPDRS3A','INFODT_NUPDRS4',
               'INFODT_GDSSHORT','INFODT_QUIPCS','INFODT_SFT','INFODT_LNSPD','INFODT_LINEORNT','INFODT_SDM',
               'INFODT_STAI','INFODT_HVLT','INFODT_PECN','INFODT_PENEURO','INFODT_DATSCAN','INFODT_AVIMAG',
               'INFODT_MRI']

def fill_missing_INFODT(df_in, INFODT_cols, INFODT_new = 'INFODT_'):
    """
    My ad hoc policy here was to use the earliest (oldest) available date.
    I decided this after hand-checking a few entries that had more than one 
    unique INFODT across tables to discover that more tests/exams were performed 
    at the earliest available date.
    
    For these cases, visits were up to 3 months apart.
    
    A better way to do this might be to explicitly search for missing INFODTs 
    before merging tables, and even using different INFODTs for different biomarkers.
    """
    approaches = ['consensus','precedence']
    approach = approaches[0]
    
    df = df_in.copy()
    df.reset_index(inplace=True)
    
    if approach == 'consensus':
        print('Taking a consensus from all INFODTs: df[INFODTs].mode(axis=1)')
        #* Find the mode/min/max/mean per row
        mode = df[INFODT_cols].mode(axis=1)
        mode = mode.apply(pd.to_datetime)
        #* Since the mode can return more than one value per row, take the first one
        mode_ = mode[0]
        #* DIAGNOSTIC: print out the rows with multiple modes
        #for k in mode.index:
        #    v = mode.loc[k]
        #    v = v[~v.isnull().values]
        #    if len(v)>1:
        #        print('index = {0}\n  {1}'.format(k,v))
        df[INFODT_new] = mode_
        
    elif approach == 'precedence':
        for c in INFODT_cols:
            if ~df[c].isnull():
                df[INFODT_new] = df[c]
                break
    else :
        print('ERROR in fill_missing_INFODT(): unknown approach.')
    
    # #* SLOW VERSION: Loop through each row.
    # for k in range(0,df.shape[0]):
    #     #* Extract and Sort dates (increasing, NaN at end)
    #     infodates = df.loc[k,INFODT_cols]
    #     infodates.sort_values(inplace=True)
    #     #* Set new
    #     df.loc[k,INFODT_new] = infodates[0]
    
    return df

df_PPMIMERGE = fill_missing_INFODT(df_PPMIMERGE,INFODT_cols)

# # Validate PPMIMERGE
# Find difference between biomarker values in original table and in PPMIMERGE
def validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys=['PATNO','EVENT_ID']):
    # Use tables_of_interest and columns_of_interest and 
    # check whether column-data in PPMIMERGE matches the original table,
    # using table keys: PATNO,EVENT_ID
    key = df[keys]
    key_PPMIMERGE = df_PPMIMERGE[keys]
    data = df[keys+cols]
    data = data.loc[~data[keys].isnull().any(axis=1)] #* Remove null keys, if they exist
    data_MERGE = df_PPMIMERGE[keys+cols]
    cols_MERGE = [c + '_MERGE' for c in cols]
    #* Use pandas.merge to extract corresponding data from PPMIMERGE (I realise this is kinda circular, since I used pd.merge() to generate PPMIMERGE)
    df_ = data.merge(data_MERGE,on=keys,how='left',suffixes=('','_MERGE'))
    
    #* Calculate difference: original data, minus PPMIMERGE
    tmp = df_[cols].apply(lambda x: pd.to_numeric(x))
    tmp_MERGE = df_[cols_MERGE].apply(lambda x: pd.to_numeric(x))
    dif = tmp.values - tmp_MERGE.values
    
    #* DIAGNOSTIC
    #print(tmp.head())
    #print(tmp_MERGE.head())
    # dif = df_[cols].values - df_[[c + '_MERGE' for c in cols]].values
    
    #* The number expected to match (handle null values) = len(cols), if all values are not null
    num_expected_to_match = np.sum(~df_[cols].isnull(),axis=1) # column of values, one row per key
    
    #* DIAGNOSTIC
    #print(len(cols))
    #print(np.sum(num_expected_to_match!=len(cols)))
    #print(num_expected_to_match.values[num_expected_to_match!=len(cols)])
    #print(dif)
    #print('Number of columns expected to match: {0} of {1}'.format(num_expected_to_match,len(cols)))
    
    #* Identify unmatched rows, where difference is not zero
    unmatched_rows = df_.loc[np.sum((dif)==0,axis=1) != num_expected_to_match].copy()
    n = unmatched_rows.shape[0]
    
    return (n,unmatched_rows)

# ## Manual: dataframes requiring intervention
keys = ['PATNO','EVENT_ID']
#* DTI biomarkers
cols = ['DTI_FA_ROI1', 'DTI_FA_ROI2', 'DTI_FA_ROI3', 'DTI_FA_ROI4', 'DTI_FA_ROI5', 'DTI_FA_ROI6', 'DTI_FA_REF1', 'DTI_FA_REF2', 
        'DTI_E1_ROI1', 'DTI_E1_ROI2', 'DTI_E1_ROI3', 'DTI_E1_ROI4', 'DTI_E1_ROI5', 'DTI_E1_ROI6', 'DTI_E1_REF1', 'DTI_E1_REF2', 
        'DTI_E2_ROI1', 'DTI_E2_ROI2', 'DTI_E2_ROI3', 'DTI_E2_ROI4', 'DTI_E2_ROI5', 'DTI_E2_ROI6', 'DTI_E2_REF1', 'DTI_E2_REF2', 
        'DTI_E3_ROI1', 'DTI_E3_ROI2', 'DTI_E3_ROI3', 'DTI_E3_ROI4', 'DTI_E3_ROI5', 'DTI_E3_ROI6', 'DTI_E3_REF1', 'DTI_E3_REF2']
c = keys + cols
df = df_DTIROI_unstacked.loc[df_DTIROI_unstacked.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
(DTIROI_validated,unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if DTIROI_validated==0:
    print('DTI biomarkers validated')
else:
    print('DTI biomarkers failed validation')
#* AV Imaging
cols = ['RCAUD-S', 'RPUTANT-S', 'RPUTPOST-S', 'LCAUD-S', 'LPUTANT-S', 'LPUTPOST-S']
c = Keys + cols
df = df_AVIMAG.loc[df_AVIMAG.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
(AVIMAG_validated,unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if AVIMAG_validated==0:
    print('AV Imaging biomarkers validated')
else:
    print('AV Imaging biomarkers failed validation')
#* BIO
cols = ['ABeta 1-42', 'Abeta 42',
       'Apolipoprotein A1', 'CSF Alpha-synuclein', 'CSF Hemoglobin',
       'EGF ELISA', 'HDL', 'IgG', 'IgG3', 'IgG3/IgG', 'LDL', 'MTDNA_DELETION',
       'MTDNA_ND1_CN', 'MTDNA_ND4_CN', 'NDNA_B2M_CN', 'NDNA_B2M_CN_v2', 'PD2',
       'PD2 Peptoid', 'PD2 peptoid', 'Serum IGF-1', 'Total Cholesterol',
       'Total tau', 'Triglycerides', 'p-Tau181P', 'pTau', 'tTau']
c = Keys + cols
df = df_BIOLAB_SemiWide.loc[df_BIOLAB_SemiWide.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
(BIOLAB_SemiWide_validated,unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if BIOLAB_SemiWide_validated==0:
    print('BIOLAB biomarkers validated')
else:
    print('BIOLAB failed validation')
#* DATSCAN
cols = ['CAUDATE_L','CAUDATE_R','PUTAMEN_L','PUTAMEN_R']
c = Keys + cols
df = df_DATSCAN.loc[df_DATSCAN.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
(DATSCAN_validated,unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if DATSCAN_validated==0:
    print('DATSCAN biomarkers validated')
else:
    print('DATSCAN failed validation')
#* UPDRS3: beware PN3RIGRL vs NP3RIGRL (my renamed version)
cols = ['NP3_TOT', 'NP3SPCH', 'NP3FACXP', 'NP3RIGN', 'NP3RIGRU',
       'NP3RIGLU', 'PN3RIGRL', 'NP3RIGLL', 'NP3FTAPR', 'NP3FTAPL', 'NP3HMOVR',
       'NP3HMOVL', 'NP3PRSPR', 'NP3PRSPL', 'NP3TTAPR', 'NP3TTAPL', 'NP3LGAGR',
       'NP3LGAGL', 'NP3RISNG', 'NP3GAIT', 'NP3FRZGT', 'NP3PSTBL', 'NP3POSTR',
       'NP3BRADY', 'NP3PTRMR', 'NP3PTRML', 'NP3KTRMR', 'NP3KTRML', 'NP3RTARU',
       'NP3RTALU', 'NP3RTARL', 'NP3RTALL', 'NP3RTALJ', 'NP3RTCON', 'NHY']
c = keys + cols
df = df_NUPDRS3A.loc[df_NUPDRS3A.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
cols_ = [col.replace('PN3','NP3') for col in cols]
df.rename(columns={'PN3RIGRL':'NP3RIGRL'},inplace=True)
#* Add the suffix I used above
cols__ = [col+'_NUPDRS3A' for col in cols_]
df.rename(columns=dict(zip(cols_,cols__)),inplace=True)
(NUPDRS3A_validated,NUPDRS3A_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols__,keys)
if NUPDRS3A_validated==0:
    print('NUPDRS3A biomarkers validated')
else:
    print('NUPDRS3A failed validation')

df = df_NUPDRS3.loc[df_NUPDRS3.PATNO.isin(df_RANDOM_enrolled.PATNO),c].copy()
df.rename(columns={'PN3RIGRL':'NP3RIGRL'},inplace=True)
(NUPDRS3_validated,NUPDRS3_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols_,keys)
if NUPDRS3_validated==0:
    print('NUPDRS3 biomarkers validated')
else:
    print('NUPDRS3 failed validation')

#* CLINDX
key = 'CLINDX'
cols = columns_of_interest[key][1:] # Removes INFODT
c = ['PATNO','EVENT_ID'] + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(CLINDX_validated,CLINDX_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if CLINDX_validated==0:
    print('CLINDX biomarkers validated')
else:
    print('CLINDX failed validation')

#* MRI
if include_MRI_volumes:
    key = 'MRI'
    cols = columns_of_interest[key]
    cols = [cols[0]] + cols[2:-5] # Removes MRIDT, PDMEDYN, ONLDOPA, ONDOPAG, ONOTHER
    c = ['PATNO','EVENT_ID'] + cols
    n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
    df_ = PPMI_df[n][c]
    df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
    (MRI_validated,MRI_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
    if MRI_validated==0:
        print('MRI biomarkers validated')
    else:
        print('MRI failed validation')

#* Use of PD medication
key = 'PDMEDUSE'
cols = columns_of_interest[key]
c = keys + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(PDMEDUSE_validated,PDMEDUSE_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if PDMEDUSE_validated==0:
    print('{0} biomarkers validated'.format(key))
else:
    print('{0} failed validation'.format(key))
    print(PDMEDUSE_unmatched_rows)

#* Neuro exam
key = 'PENEURO'
cols = columns_of_interest[key]
#* Remove non-numeric columns
cols_to_remove = ['MSRACM','MSLACM','MSRLCM','MSLLCM','COFNRCM','COFNLCM','COHSRCM','COHSLCM','SENRACM','SENLACM',
                  'SENRLCM','SENLLCM','RFLRACM','RFLLACM','RFLRLCM','RFLLLCM','PLRRCM','PLRLCM']
cols = [col for col in cols if col not in cols_to_remove]

c = keys + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(PENEURO_validated,PENEURO_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if PENEURO_validated==0:
    print('{0} biomarkers validated'.format(key))
else:
    print('{0} failed validation'.format(key))
    print(PENEURO_unmatched_rows)

#* QUIPCS
key = 'QUIPCS'
cols = columns_of_interest[key]
#* Remove non-numeric columns
cols_to_remove = ['TMDISMED','CNTRLDSM']
cols = [col for col in cols if col not in cols_to_remove]

c = keys + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(QUIPCS_validated,QUIPCS_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if QUIPCS_validated==0:
    print('{0} biomarkers validated'.format(key))
else:
    print('{0} failed validation'.format(key))
    print(QUIPCS_unmatched_rows)

#* SCOPA
key = 'SCOPAAUT'
cols = columns_of_interest[key]
#* Remove non-numeric columns
cols_to_remove = ['SCAU23AT','SCAU26AT','SCAU26BT','SCAU26CT','SCAU26DT']
cols = [col for col in cols if col not in cols_to_remove]

c = keys + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(SCOPAAUT_validated,SCOPAAUT_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if SCOPAAUT_validated==0:
    print('{0} biomarkers validated'.format(key))
else:
    print('{0} failed validation'.format(key))
    print(SCOPAAUT_unmatched_rows)

#* Benton
key = 'LINEORNT'
cols = columns_of_interest[key]
#* Remove non-numeric columns
cols_to_remove = []
cols = [col for col in cols if col not in cols_to_remove]

c = keys + cols
n = np.where([1 if n.lower()==tables_of_interest[key].lower() else 0 for n in PPMI_df_names])[0][0]
df_ = PPMI_df[n][c]
df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
(LINEORNT_validated,LINEORNT_unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
if LINEORNT_validated==0:
    print('{0} biomarkers validated'.format(key))
else:
    print('{0} failed validation'.format(key))
    print(LINEORNT_unmatched_rows)

# ## Validate the rest automatically using `tables_of_interest` and `columns_of_interest`
dfs_of_interest = ['NUPDRS1','NUPDRS1p','NUPDRS2p','NUPDRS4',
                   'EPWORTH',
                   'GDSSHORT',
                   'HVLT',
                   'LNSPD',
                   'MOCA',
                   'PECN',
                   'REMSLEEP',
                   'SDM',
                   'SFT',
                   'STAI',
                   'UPSIT',
                   'VITAL']

dfs_validated = []
unmatched_rows_ = []
for k in dfs_of_interest:
    n = np.where([1 if n.lower()==tables_of_interest[k].lower() else 0 for n in PPMI_df_names])[0][0]
    cols = columns_of_interest[k]
    df_ = PPMI_df[n][keys + cols]
    df = df_.loc[df_.PATNO.isin(df_RANDOM_enrolled.PATNO)].copy()
    (v,unmatched_rows) = validate_PPMIMERGE(df_PPMIMERGE,df,cols,keys)
    if v==0:
        print('{0} biomarkers validated'.format(k))
    else:
        print('{0} failed validation'.format(k))
        df_ = df.copy()
    dfs_validated.append(v)
    unmatched_rows_.append(unmatched_rows)

# ### PPMIMERGE: remove (nearly-)empty rows from the large cross-join
threshold = 0.9 * df_PPMIMERGE.shape[1]
n_null = np.sum(df_PPMIMERGE.isnull(),axis=1)
empty_rows = n_null > threshold
plt.hist(n_null)
a = plt.gca()
plt.plot(threshold*np.ones([2,1]),a.get_ylim())
plt.show()

df_PPMIMERGE_ = df_PPMIMERGE.loc[~empty_rows].copy()

# ## Merge `BL` and `SC`, where they both exist
# Tests performed at screening usually weren't repeated at baseline. Notably:
# - MOCA, GDSSHORT, DATSCAN, PECN, PENEURO, STAI
# Uses `pd.DataFrame.fillna()` (NOTE: since visit order is chronological in PPMIMERGE, `SC` is before `BL` and we can forward fill)
boo = (df_PPMIMERGE_.EVENT_ID=='SC') & (~df_PPMIMERGE_.PATNO.isin(['42418R','42415R']))
x = df_PPMIMERGE_.loc[boo,'MCATOT']
x_PATNO = df_PPMIMERGE_.loc[boo,'PATNO']
#* Identify SC and BL rows
scbl = (df_PPMIMERGE_.EVENT_ID=='BL') | (df_PPMIMERGE_.EVENT_ID=='SC')
#* Pad selected columns: copy, fillna, reinsert
c = columns_of_interest['MOCA'] + columns_of_interest['GDSSHORT'] + columns_of_interest['DATSCAN_RESULTS']     + columns_of_interest['PECN'] + columns_of_interest['PENEURO']     + columns_of_interest['STAI']
#* Add VSINTRPT, too
c = c + ['VSINTRPT']
dftmp = df_PPMIMERGE_.loc[scbl,Keys+c].copy()
dftmp.fillna(method='pad',axis=0,inplace=True)
df_PPMIMERGE_.loc[scbl,Keys+c] = dftmp
#* Drop SC visit
df_PPMIMERGE__ = df_PPMIMERGE_.loc[df_PPMIMERGE_.EVENT_ID!='SC']

#* DIAGNOSTIC
# boo = (df_PPMIMERGE__.EVENT_ID=='BL') & (~df_PPMIMERGE__.PATNO.isin(['42418R','42415R','42357R']))
# y = df_PPMIMERGE__.loc[boo,'MCATOT']
# y_PATNO = df_PPMIMERGE__.loc[boo,'PATNO']
# plt.plot(x_PATNO,x,'r+')
# plt.hold
# plt.plot(y_PATNO,y,'b.')
# plt.show()
# print(len(x_PATNO))
# print(len(y_PATNO))


# ### PPMIMERGE: add `INFODT_bl` (for `Years_bl`)
#*** Calculate Years_bl (NOTE: SC visits can have negative Years_bl)
#* First: INFODT_bl
def fill_baseline(df,col='INFODT_'):
    """
    Returns copy of DataFrame with baseline values filled
    
    Baseline rows identified by
       EVENT_ID == 'BL'
    Subject rows identified by
       PATNO
    """
    df_out = df.copy()
    rowz_bl = df.EVENT_ID=='BL'
    rowz_sc = df.EVENT_ID=='SC'
    for PATNO in df['PATNO'].unique():
        rowz = df.PATNO==PATNO
        if sum(rowz & rowz_bl)==1:
            INFODT_bl = df.loc[rowz & rowz_bl,col]
            df_out.loc[rowz,col+'bl'] = INFODT_bl.values
        else:
            print('Bugger! PATNO = {0} has {1} bl rows. Searching for SC visit instead.'.format(PATNO,sum(rowz & rowz_bl)))
            if sum(rowz & rowz_sc)==1:
                INFODT_sc = df.loc[rowz & rowz_sc,col]
                df_out.loc[rowz,col+'bl'] = INFODT_sc.values
            else:
                print('Bugger! PATNO = {0}; {1}/{2} bl/sc rows'.format(PATNO,sum(rowz & rowz_bl),sum(rowz & rowz_sc)))
    return df_out

df_PPMIMERGE___ = fill_baseline(df_PPMIMERGE__,col='INFODT_')

#*** Then Years_bl and Age
def calculate_years_bl(df,infodt_col='INFODT_',infodt_bl_col='INFODT_bl'):
    yrsbl = 'Years_bl'
    visit_dt = pd.to_datetime(df[infodt_col])
    visit_dt_bl = pd.to_datetime(df[infodt_bl_col])
    delt = visit_dt - visit_dt_bl
    return delt
#* Keep three decimal places: PPMI dates are only precise to within one month = 0.083 years
df_PPMIMERGE___['Years_bl'] = calculate_years_bl(df_PPMIMERGE___).apply(lambda x: round(x.days/daysInAYear,3))

#* Age
Age = df_PPMIMERGE___['INFODT_'] - pd.to_datetime(df_PPMIMERGE___['BIRTHDT'])
Age = Age.apply(lambda x: round(x.days/daysInAYear,3))
df_PPMIMERGE___['Age'] = Age


# ## Reorder the columns, for the fun of it
cols = df_PPMIMERGE___.columns
cols_list = cols.tolist()
#* Start with some demographics
cols_ordered = ['index','PATNO','APPRDX_enrol','EVENT_ID','Age','Years_bl','GENDER',
                'FAMHXPD_N1stDegree','FAMHXPD_N2ndDegree',
                'ENROLLDT','CONSNTDT','GENECAT', 'MUTRSLT', 'LRRKCD','BIRTHDT',
                'HISPLAT','RAINDALS','RAASIAN','RABLACK','RAHAWOPI','RAWHITE','RANOS'
               ]
#* Append the rest
for c in cols_ordered:
    k = np.where([d==c for d in cols_list])[0][0]
    cols_list.pop(k)
cols_ordered.extend(cols_list)
#* Reorder columns
df_PPMIMERGE___ = df_PPMIMERGE___[cols_ordered]

# # Write PPMIMERGE to CSV
df_PPMIMERGE___.to_csv(PPMIMERGE_csv,index=False)
#* Printout
print('PPMIMERGE.shape: {0},{1}'.format(df_PPMIMERGE.shape[0],df_PPMIMERGE.shape[1]))
print('PPMIMERGE_.shape: {0},{1}'.format(df_PPMIMERGE_.shape[0],df_PPMIMERGE_.shape[1]))
print('PPMIMERGE___.shape: {0},{1}'.format(df_PPMIMERGE___.shape[0],df_PPMIMERGE___.shape[1]))
print('Number of unique PATNOs: {0}'.format(len(df_PPMIMERGE___.PATNO.unique())))

