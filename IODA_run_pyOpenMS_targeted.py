# Author: Louis Felix Nothias, louisfelix.nothias@gmail.com, June 2020
import os
import subprocess
import sys
from logzero import logger, logfile
import datetime
import zipfile
from datetime import date
from IODA_exclusion_workflow import get_all_file_paths
from subprocess import call
from pyopenms import MSExperiment, MzMLFile, MassTraceDetection, ElutionPeakDetection, FeatureMap, FeatureFindingMetabo, FeatureXMLFile, FeatureGroupingAlgorithmKD, ConsensusMap, ColumnHeader, ConsensusXMLFile

import glob
import pandas as pd
import gc

def IODA_targeted_workflow(blank_mzML:str,sample_mzML:str,ppm_tolerance:float,noise_level:float, chrom_peak_snr:float, elements_alphabet:str):
    # Test samples
        #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_SAMPLE = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_SAMPLE = "https://drive.google.com/file/d/1_lOYEtsmEPAlfGVYbzJpLePPSitUp1yh/view?usp=sharing"
        #input_BLANK = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"
        #input_SAMPLE = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_F1-1_F2-1_NIST-1_To-1_20181227135238.mzML"

    
    os.system('rm OpenMS_workflow/logfile_IODA_OpenMS_from_mzML.txt')
    logfile('OpenMS_workflow/logfile_IODA_OpenMS_from_mzML.txt')
    OpenMS_output_folder = "OpenMS_output"
    OpenMS_folder = "OpenMS_workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r OpenMS_workflow/OpenMS_input/*')
    os.system('rm -r OpenMS_workflow/OpenMS_output/*')
    os.system('mkdir download_results')
        
        
        
    #large_noise = 5E5
    #narrow_noise = 1E5
    #ppm_error = 10

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-targeted WORKFLOW')
    logger.info('======')
    logger.info('Path to the input files: ')
    logger.info('Blank: '+blank_mzML)
    logger.info('Sample: '+sample_mzML)
    
    # Collect the mzML and copy
    def download_copy_mzML(input_file):
        # Test samples
            #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/OpenMS_input/Blank.mzML"
            #input_mzML = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
            #input_mzML = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
            #input_mzML = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"

        today = str(date.today())
        now = datetime.datetime.now()
        logger.info(now)
        logger.info('STARTING the IODA-targeted WORKFLOW with OpenMS')
        logger.info('======')
        logger.info('Getting the mzML, please wait ...')

        if input_file.startswith(('http','ftp')):
            if 'google' in input_file:
                logger.info('This is the Google Drive download link:'+str(input_file))
                logger.info('Downloading the mzML, please wait ...')
                url_id = input_file.split('/', 10)[5]
                prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
                input_file = prefixe_google_download+url_id
                bashCommand1 = "wget --no-check-certificate '"+input_file+"' -O "+os.path.join(OpenMS_folder+"/OpenMS_input/", os.path.basename(input_file)[:-4] + ".mzML")+" || rm -f "+os.path.join(OpenMS_folder+"/OpenMS_input/", os.path.basename(input_file))
                cp1 = subprocess.run(bashCommand1,shell=True)
                try:
                    cp1
                except subprocess.CalledProcessError:
                    raise
            if 'massive.ucsd.edu' in input_file:
                logger.info('This is the MassIVE repository link: '+str(input_file))
                logger.info('Downloading the mzML, please wait ... ')
                bashCommand4 = "wget -r "+input_file+" -O "+os.path.join(OpenMS_folder+"/OpenMS_input/", os.path.basename(input_file)[:-4] + ".mzML")+" || rm -f "+os.path.join(OpenMS_folder+"/OpenMS_input/", os.path.basename(input_file))
                cp4 = subprocess.run(bashCommand4,shell=True)
                try:
                    cp4
                except subprocess.CalledProcessError:
                    raise

        elif input_file.endswith(('.raw','.RAW')):
            logger.info('Thermo RAW file detected')
            logger.info('This is the input file path: '+str(input_file))
            bashCommand5 = "mono ThermoRawFileParser/ThermoRawFileParser.exe -i="+input_file+" --logging=1 --ignoreInstrumentErrors --excludeExceptionData --output_file "+os.path.join(OpenMS_folder+"/OpenMS_input/", os.path.basename(input_file)[:-3] + ".mzML")
            logger.info('The file is converting to mzML thanks ThermoRawFileParser v1.3.4, please wait few seconds ...: '+str(input_file))
            logger.info(str(bashCommand5))
            cp5 = subprocess.run(bashCommand5,shell=True)
            try:
                cp5
            except subprocess.CalledProcessError:
                raise

        else:
            #Check the file path is correct for local upload
            logger.info('This is the input file path: '+str(input_file))
            bashCommand3 = "cp "+input_file+" "+os.path.join("OpenMS_workflow/OpenMS_input/", os.path.basename(input_file))
            cp3 = subprocess.run(bashCommand3,shell=True)
            try:
                cp3
            except subprocess.CalledProcessError:
                raise
        # Error getting the file ! PLEASE VERY THE PATH TO THE FILE OR DOWNLOAD LINK ...

        if input_file.endswith(('.raw','.RAW')):
            try:
                f = open(os.path.join("OpenMS_workflow/OpenMS_input/", os.path.basename(input_file)[:-3] + ".mzML"))
                f.close()
            except subprocess.CalledProcessError:
                logger.info('There was an error getting the file !')
            logger.info('The mzML file was found !')

        elif input_file.endswith(('.mzML','.mzml')):
            try:
                f = open(os.path.join("OpenMS_workflow/OpenMS_input/", os.path.basename(input_file)[:-5] + ".mzML"))
                f.close()
            except subprocess.CalledProcessError:
                logger.info('There was an error getting the file !')
            logger.info('The mzML file was found !')

        logger.info('Copying the mzML to the OpenMS input folder')

    # Downloading the two input mzML
    logger.info('Copying the mzML files ...')
    download_copy_mzML(blank_mzML)
    download_copy_mzML(sample_mzML)

    # Feature Detection
    input_mzml_files = glob.glob("OpenMS_workflow/OpenMS_input/*.mzML") # introduce a set of mzML files directory

    # Mass trace detection

    #mtd
    mass_error_ppm = 10.0
    noise_threshold_int = 5e5
    chrom_peak_snr = 3.0
    chrom_fwhm = 10.0
    reestimate_mt_sd = 'true'
    quant_method = 'max_height'
    trace_termination_criterion = 'outlier'
    trace_termination_outliers = 3
    min_sample_rate = 0.3
    min_trace_length = 1.0
    max_trace_length = 100.0

        #epd
    enabled = 'true'
    width_filtering = 'fixed'
    min_fwhm = 0.0
    max_fwhm = 30.0
    masstrace_snr_filtering = 'false'

        #ffm
    local_rt_range = 7.0
    local_mz_range = 6.0
    charge_lower_bound = 1
    charge_upper_bound = 3
    report_summed_ints = 'false'
    enable_RT_filtering = 'false'
    isotope_filtering_model = 'none'
    mz_scoring_13C = 'false'
    use_smoothed_intensities = 'true'
    report_convex_hulls = 'false'
    remove_single_traces = 'false'
    mz_scoring_by_elements = 'false'
    elements = elements_alphabet #'CHNOPS'

        #width_filtering
    for filename in input_mzml_files: # for each file in the set of files
        print("Mass Trace Detection: ", filename) #print the filename
        exp = MSExperiment()    
        MzMLFile().load(filename, exp) # load each mzML file to an OpenMS file format (MSExperiment)

        mass_traces = [] # introduce an empty list where the mass traces will be loaded
        mtd = MassTraceDetection()
        mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
        mtd_par.setValue("mass_error_ppm", mass_error_ppm) # high-res instrument, orbitraps
        mtd_par.setValue("noise_threshold_int", noise_threshold_int) # data-dependent (usually works for orbitraps)
        mtd_par.setValue("chrom_peak_snr", chrom_peak_snr)
        mtd_par.setValue("reestimate_mt_sd", reestimate_mt_sd)
        mtd_par.setValue("quant_method", quant_method)
        mtd_par.setValue("trace_termination_criterion", trace_termination_criterion)
        mtd_par.setValue("trace_termination_outliers", trace_termination_outliers)
        mtd_par.setValue("min_sample_rate", min_sample_rate)
        mtd_par.setValue("min_trace_length", min_trace_length)
        mtd_par.setValue("max_trace_length", max_trace_length)
        mtd.setParameters(mtd_par) # set the new parameters
        mtd.run(exp, mass_traces, 0) # run mass trace detection

        # Elution peak detection

        print("Elution Peak Detection: ", filename)
        mass_traces_deconvol = []
        epd = ElutionPeakDetection()
        epd_par = epd.getDefaults()
        epd_par.setValue("width_filtering", width_filtering) # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
        epd_par.setValue("chrom_fwhm", chrom_fwhm)
        epd_par.setValue("min_fwhm", min_fwhm) 
        epd_par.setValue("max_fwhm", max_fwhm) 
        epd_par.setValue("enabled", enabled) 
        epd_par.setValue("masstrace_snr_filtering", masstrace_snr_filtering) 
        epd.setParameters(epd_par)
        epd.detectPeaks(mass_traces, mass_traces_deconvol)

        # Feature detection

        print("Feature Detection: ", filename)
        feature_map_FFM = FeatureMap() # output features 
        chrom_out = [] # output chromatograms 
        ffm = FeatureFindingMetabo()
        ffm_par = ffm.getDefaults() 
        ffm_par.setValue("remove_single_traces", remove_single_traces) # remove mass traces without satellite isotopic traces
        ffm_par.setValue("local_rt_range", local_rt_range)  
        ffm_par.setValue("charge_lower_bound", charge_lower_bound)      
        ffm_par.setValue("charge_upper_bound", charge_upper_bound)  
        ffm_par.setValue("report_summed_ints", report_summed_ints)      
        ffm_par.setValue("enable_RT_filtering", enable_RT_filtering)      
        ffm_par.setValue("isotope_filtering_model", isotope_filtering_model)      
        ffm_par.setValue("mz_scoring_13C", mz_scoring_13C)    
        ffm_par.setValue("use_smoothed_intensities", use_smoothed_intensities)      
        ffm_par.setValue("report_convex_hulls", report_convex_hulls)  
        ffm_par.setValue("mz_scoring_by_elements", mz_scoring_by_elements)      
        ffm_par.setValue("elements", elements)       
        ffm.setParameters(ffm_par)
        ffm.run(mass_traces_deconvol, feature_map_FFM, chrom_out)
        feature_map_FFM.setUniqueIds() # Assigns a new, valid unique id per feature
        feature_map_FFM.setPrimaryMSRunPath([filename.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
        FeatureXMLFile().store("OpenMS_workflow/OpenMS_output/"+os.path.basename(filename)[:-5] + ".featureXML", feature_map_FFM)

    print("Finished Feature Detection")

        # load feature files 

    input_feature_files = glob.glob('OpenMS_workflow/OpenMS_output/*.featureXML') # set of feature files

    feature_maps = [] # empty list to fill with FeatureMaps: the OpenMS file format for feature files
    for featurexml_file in input_feature_files:
        fmap = FeatureMap()
        FeatureXMLFile().load(featurexml_file, fmap) # load each file to a feature map
        feature_maps.append(fmap) # append all maps to the empty list 

    # Feature grouping

    feature_grouper = FeatureGroupingAlgorithmKD()

    consensus_map = ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, ColumnHeader())
        file_description.filename = os.path.basename(feature_map.getMetaValue('spectra_data')[0].decode())
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description

    keep_subelements = 'true'
    mz_unit = 'ppm'
    nr_partitions = 50
    warp_enabled = 'false'
    warp_rt_tol = 30.0
    warp_mz_tol = 10.0
    warp_max_pairwise_log_fc = 0.0
    warp_min_rel_cc_size = 0.0
    warp_max_nr_conflicts = -1
    link_rt_tol = 30
    link_mz_tol = 10
    link_charge_merging = 'With_charge_zero'
    link_adduct_merging = 'Any'
    distance_RT_exponent = 1.0
    distance_MZ_weight = 2.0
    distance_intensity_log_transform = 1.0
    distance_intensity_weight = 1.0
    LOWESS_span = 0.0
    LOWESS_num_iterations = 3
    LOWESS_delta = -1.0
    LOWESS_interpolation_type = 'cspline'
    LOWESS_extrapolation_type = 'four-point-linear'

    feature_grouper_par = feature_grouper.getDefaults()
    feature_grouper_par.setValue("keep_subelements", keep_subelements)
    feature_grouper_par.setValue("mz_unit", mz_unit) 
    feature_grouper_par.setValue("nr_partitions", nr_partitions) 
    feature_grouper_par.setValue("warp:enabled", warp_enabled) 
    feature_grouper_par.setValue("warp:rt_tol", warp_rt_tol) 
    feature_grouper_par.setValue("warp:mz_tol", warp_mz_tol) 
    feature_grouper_par.setValue('warp:max_pairwise_log_fc', warp_max_pairwise_log_fc) 
    feature_grouper_par.setValue("warp:min_rel_cc_size", warp_min_rel_cc_size) 
    feature_grouper_par.setValue("warp:max_nr_conflicts", warp_max_nr_conflicts) 
    feature_grouper_par.setValue("link:rt_tol", link_rt_tol)
    feature_grouper_par.setValue("link:mz_tol", link_mz_tol) 
    feature_grouper_par.setValue("link:adduct_merging", link_adduct_merging) 
    feature_grouper_par.setValue("link:charge_merging", link_charge_merging) 
    feature_grouper_par.setValue("distance_RT:exponent", distance_RT_exponent) 
    feature_grouper_par.setValue("distance_MZ:weight", distance_MZ_weight)
    feature_grouper_par.setValue("distance_intensity:log_transform", distance_intensity_log_transform) 
    feature_grouper_par.setValue("distance_intensity:weight", distance_intensity_weight) 
    feature_grouper_par.setValue("LOWESS:span", LOWESS_span) 
    feature_grouper_par.setValue("LOWESS:num_iterations", LOWESS_num_iterations)
    feature_grouper_par.setValue("LOWESS:delta", LOWESS_delta)
    feature_grouper_par.setValue("LOWESS:num_iterations", LOWESS_num_iterations)
    feature_grouper_par.setValue("LOWESS:interpolation_type", LOWESS_interpolation_type)
    feature_grouper_par.setValue("LOWESS:extrapolation_type", LOWESS_extrapolation_type)

        #print(feature_grouper_par)

    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setUniqueIds()
    consensus_map.setColumnHeaders(file_descriptions)

    ConsensusXMLFile().store('OpenMS_workflow/OpenMS_output/consensus.consensusXML', consensus_map)

    df = consensus_map.get_df()

    df = df.rename(columns={'mz':'Mass [m/z]', 
                                'RT':'retention_time'})
    df = df.drop(columns=['sequence'])
    df.to_csv('OpenMS_workflow/OpenMS_output/consensus.csv', index=False)

    # Error with the OpenMS workflow. No output files.
    try:
        f = open('OpenMS_workflow/OpenMS_output/consensus.csv')
        f.close()
    except:
        logger.info('There was an issue with the pyOpenMS workflow ! See the log below.')
        raise

    logger.info('======')
    logger.info('Completed the OpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ..')
    get_all_file_paths('OpenMS_workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('Completed zipping up the OpenMS workflow result files')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE IODA-targeted WORKFLOW')


def MS2Planner_Curve_OpenMS(sample_mzML:str):
    # Test samples
        #source_mzML1 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #source_mzML2 = "https://raw.githubusercontent.com/lfnothias/IODA_MS/master/tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_SAMPLE = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_latex_MS1_2uL.mzML"
        #input_BLANK = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_SAMPLE = "https://drive.google.com/file/d/1_lOYEtsmEPAlfGVYbzJpLePPSitUp1yh/view?usp=sharing"
        #input_BLANK = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"
        #input_SAMPLE = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_F1-1_F2-1_NIST-1_To-1_20181227135238.mzML"

    TOPPAS_Pipeline = "MS1_MS2Planner_Curve_mzTab.toppas"
    OpenMS_output_folder = "OpenMS_output"
    OpenMS_input_folder = "OpenMS_input"
    OpenMS_folder = "OpenMS_workflow"

    os.system('rm '+OpenMS_folder+'/logfile_MS2Planner_OpenMS.txt')
    logfile(OpenMS_folder+'/logfile_MS2Planner_OpenMS.txt')
    os.system('rm download_results/MS2Planner_OpenMS_results.zip')
    os.system('rm -r '+OpenMS_folder+'/'+OpenMS_output_folder+'/OpenMS_out/MS2Planner*')
    os.system('mkdir download_results')
    os.system('mv '+OpenMS_folder+'/'+OpenMS_output_folder+'/TOPPAS.log '+OpenMS_folder+'/'+OpenMS_output_folder+'/TOPPAS_FFM.log')
    os.system('mkdir '+OpenMS_folder+'/'+OpenMS_output_folder+'/OpenMS_out/MS2Planner_mzTab')


    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the MS2Planner Curve processing')
    logger.info('======')
    logger.info('Path to the input files: ')
    logger.info('    Sample: '+sample_mzML)

    logger.info('======')
    logger.info('Download the latest version of the workflow from the repository ...')

    try:
        bashCommand0 = "wget https://github.com/lfnothias/IODA_MS/raw/MS2Planner_merge_w_master/"+OpenMS_folder+'/'+TOPPAS_Pipeline+" -O "+OpenMS_folder+'/'+TOPPAS_Pipeline
        logger.info(bashCommand0)
        cp0 = subprocess.run(bashCommand0,shell=True)
        cp0
        a_file = open(OpenMS_folder+'/'+TOPPAS_Pipeline, "r")
        list_of_lines = a_file.readlines()
    except subprocess.CalledProcessError:
        logger.info('ERROR getting the reference workflow')

    # Preserve the original mzML file names in the OpenMS workflow for local files
    if sample_mzML.startswith('http'):
        if 'google' in sample_mzML:
            pass
    else:
        sample_filename = str(sample_mzML.split('/', 10)[-1])
        list_of_lines = [sub.replace('LISTITEM value="OpenMS_input/Sample.mzML', 'LISTITEM value="'+OpenMS_input_folder+'/'+sample_filename) for sub in list_of_lines]

    # Replace OpenMS workflow parameters
    # Write out the file
    a_file = open(OpenMS_folder+'/'+TOPPAS_Pipeline, "w")
    a_file.writelines(list_of_lines)
    a_file.close()

    logger.info('======')
    logger.info('Initializing the OpenMS workflow')

    try:
        vdisplay = Xvfb()
        vdisplay.start()
    except subprocess.CalledProcessError:
        raise

    logger.info('======')
    logger.info('Running the OpenMS workflow, this usually takes couple minutes, please wait ...')

    bashCommand4 = "cd "+OpenMS_folder+" && /openms-build/bin/ExecutePipeline -in "+TOPPAS_Pipeline+" -out_dir "+OpenMS_output_folder
    try:
        cp4 = subprocess.run(bashCommand4,shell=True)
        cp4
    except subprocess.CalledProcessError as e:
        logger.info("!!! There was an error with OpenMS workflow, See the log !")
        logger.info(e.output)
        f = open(OpenMS_folder+'/OpenMS_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger(file_contents)
        raise

    vdisplay.stop()

    # Error with the OpenMS workflow. No output files.
    try:
        mzTab_file = os.listdir(OpenMS_folder+"/OpenMS_output/OpenMS_out/MS2Planner_mzTab/")[0]
        f = open(OpenMS_folder+'/OpenMS_output/OpenMS_out/MS2Planner_mzTab/'+mzTab_file)
        f.close()
    except:
        logger.info('There was an issue with the OpenMS workflow ! See the log below. Not enough ressources might be available (memory, ...).')
        f = open(OpenMS_folder+'/OpenMS_output/TOPPAS.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    logger.info('======')
    logger.info('Completed the OpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ..')

    try:
        cmdA = 'cp '+OpenMS_folder+'/'+OpenMS_output_folder+'/TOPPAS.log '+OpenMS_folder+'/'+OpenMS_output_folder+'/MS2Planner_mzTab/TOPPAS_PathFinder.log'
        os.system(cmdA)
        cmdB ='cp -r '+OpenMS_folder+'/*MS2Planner*'+' '+OpenMS_folder+'/'+OpenMS_output_folder+'/MS2Planner_mzTab/'
        os.system(cmdB)
        get_all_file_paths(OpenMS_folder+'/'+OpenMS_output_folder+'/OpenMS_out/MS2Planner_mzTab','download_results/MS2Planner_OpenMS_results.zip')
    except:
        raise

    logger.info('======')
    logger.info('Completed zipping up the OpenMS workflow result files')

    logger.info('======')
    logger.info('NOW CONTINUE WITH THE REST OF THE WORKFLOW')



if __name__ == "__main__":
    IODA_targeted_workflow(str(sys.argv[1]),str(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
