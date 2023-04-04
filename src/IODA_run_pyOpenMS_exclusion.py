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
from pyopenms import MSExperiment, MzMLFile, MassTraceDetection, ElutionPeakDetection, FeatureMap, FeatureFindingMetabo, FeatureXMLFile
import glob
import pandas as pd
import gc

def IODA_exclusion_workflow(input_mzML,ppm_error,narrow_noise_threshold,large_noise_threshold):
    # Test samples
        #source_mzML = "https://raw.githubusercontent.com/lfnothias/IODA_MS/test2/tests/Euphorbia/exclusion/OpenMS_input/Blank.mzML"
        #input_mzML = "tests/Euphorbia/Targeted/OpenMS_input/Euphorbia_rogers_latex_Blank_MS1_2uL.mzML"
        #input_mzML = "https://drive.google.com/file/d/11p2Jau2T-gCQb9KZExWdC7dy8AQWV__l/view?usp=sharing"
        #input_mzML = "ftp://massive.ucsd.edu/MSV000083306/peak/QE_C18_mzML/QEC18_blank_SPE_20181227092326.mzML"

    os.system('rm OpenMS_workflow/logfile_IODA_OpenMS_from_mzML.txt')
    logfile('OpenMS_workflow/logfile_IODA_OpenMS_from_mzML.txt')
    OpenMS_output_folder = "OpenMS_output"
    OPENMS_folder = "OpenMS_workflow"
    os.system('rm download_results/IODA_OpenMS_results.zip')
    os.system('rm -r OpenMS_workflow/OpenMS_input/*')
    os.system('rm -r OpenMS_workflow/OpenMS_output/OPENMS_out/')
    os.system('mkdir download_results')

    today = str(date.today())
    now = datetime.datetime.now()
    logger.info(now)
    logger.info('STARTING the IODA-exclusion WORKFLOW with OpenMS')
    logger.info('======')
    logger.info('Getting the mzML, please wait ...')

    if input_mzML.startswith(('http','ftp')):
        if 'google' in input_mzML:
            logger.info('This is the Google Drive download link:'+str(input_mzML))
            logger.info('Downloading the mzML, please wait ...')
            url_id = input_mzML.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_mzML = prefixe_google_download+url_id
            bashCommand1 = "wget --no-check-certificate '"+input_mzML+"' -O "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)[:-4] + ".mzML")+" || rm -f "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML))
            cp1 = subprocess.run(bashCommand1,shell=True)
            try:
                cp1
            except subprocess.CalledProcessError:
                raise
        if 'massive.ucsd.edu' in input_mzML:
            logger.info('This is the MassIVE repository link: '+str(input_mzML))
            logger.info('Downloading the mzML, please wait ... ')
            bashCommand4 = "wget -r "+input_mzML+" -O "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)[:-4] + ".mzML")+" || rm -f "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML))
            cp4 = subprocess.run(bashCommand4,shell=True)
            try:
                cp4
            except subprocess.CalledProcessError:
                raise
    
    elif input_mzML.endswith(('.raw','.RAW')):
        logger.info('Thermo RAW file detected')
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand5 = "mono ThermoRawFileParser/ThermoRawFileParser.exe -i="+input_mzML+" --logging=1 --ignoreInstrumentErrors --excludeExceptionData --output_file "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)[:-4] + ".mzML")
        logger.info('The file is converting to mzML thanks ThermoRawFileParser v1.3.4, please wait few seconds ...: '+str(input_mzML))
        logger.info(str(bashCommand5))
        cp5 = subprocess.run(bashCommand5,shell=True)
        try:
            cp5
        except subprocess.CalledProcessError:
            raise
        
    else:
        #Check the file path is correct for local upload
        logger.info('This is the input file path: '+str(input_mzML))
        bashCommand3 = "cp "+input_mzML+" "+os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML))
        cp3 = subprocess.run(bashCommand3,shell=True)
        try:
            cp3
        except subprocess.CalledProcessError:
            raise
    # Error getting the file ! PLEASE VERY THE PATH TO THE FILE OR DOWNLOAD LINK ...
    
    if input_mzML.endswith(('.raw','.RAW')):
        try:
            f = open(os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)[:-4] + ".mzML"))
            f.close()
        except subprocess.CalledProcessError:
            logger.info('There was an error getting the file !')
        logger.info('The mzML file was found !')
        
    elif input_mzML.endswith(('.mzML','.mzml')):
        try:
            f = open(os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)[:-5] + ".mzML"))
            f.close()
        except subprocess.CalledProcessError:
            logger.info('There was an error getting the file !')
        logger.info('The mzML file was found !')

    logger.info('Copying the mzML to the OpenMS input folder')

    logger.info('======')
    logger.info('Changing variables of the OpenMS workflow')
    logger.info('   ppm error = '+str(ppm_error))
    logger.info('   narrow peak/feature noise threshold = '+str(narrow_noise_threshold))
    logger.info('   large peak/feature noise_threshold = '+str(large_noise_threshold))


    #"== The noise level must be a float or an integer, such as 6.0e05 =="
    try:
        float(large_noise_threshold)
    except subprocess.CalledProcessError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")

    #"== The noise level must be a float or an integer, such as 6.0e05 =="
    try:
        float(narrow_noise_threshold)
    except subprocess.CalledProcessError:
        logger.info("== The noise level must be a float or an integer, such as 6.0e05 ==")
    #"== The ppm error must be a float or an integer, such as 10 ppm =="
    try:
        float(ppm_error)
    except subprocess.CalledProcessError:
        logger.info("== The ppm error must be a float or an integer, such as 10 ppm ==")

    logger.info('======')
    logger.info('Initializing the pyOpenMS workflow')

    logger.info('======')
    logger.info('Running the pyOpenMS workflow, this usually takes less than a minute, please wait ...')

    
    # Feature Detection on large features

    def pyopenms_exclusion_large(filepath, mass_error_ppm, noise_threshold_int, chrom_peak_snr):
        # 1.1) Mass trace detection

        #mtd
        #mass_error_ppm = 10.0
        #noise_threshold_int = 4e5
        #chrom_peak_snr = 3.0

        chrom_fwhm = 30.0   #large=30  narrow = 10
        reestimate_mt_sd = 'true'
        quant_method = 'max_height'
        trace_termination_criterion = 'outlier'
        trace_termination_outliers =  3 # large = 3 narrow =2
        min_sample_rate = 0.3
        min_trace_length = 10.0   #large=10  narrow = 1
        max_trace_length = 200000.0  #large=200000  narrow = 60

        #epd
        #enabled = 'true'
        width_filtering = 'off'
        min_fwhm = 5.0
        max_fwhm = 1000.0
        masstrace_snr_filtering = 'false'

        #ffm
        local_rt_range = 7.0
        local_mz_range = 0.3
        charge_lower_bound = 1
        charge_upper_bound = 3
        report_summed_ints = 'true'
        enable_RT_filtering = 'false'
        isotope_filtering_model = 'none'
        mz_scoring_13C = 'false'
        use_smoothed_intensities = 'true'
        report_convex_hulls = 'true'
        remove_single_traces = 'false'
        mz_scoring_by_elements = 'false'
        elements = 'CHNOPS'

        print("Mass Trace Detection: ", filepath) #print the filename
        exp = MSExperiment()    
        MzMLFile().load(filepath, exp) # load each mzML file to an OpenMS file format (MSExperiment)
        
        mass_traces = [] # introduce an empty list where the mass traces will be loaded
        mtd = MassTraceDetection()
        mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
        mtd_par.setValue("mass_error_ppm", float(mass_error_ppm)) # high-res instrument, orbitraps
        mtd_par.setValue("noise_threshold_int", float(noise_threshold_int)) # data-dependent (usually works for orbitraps)
        mtd_par.setValue("reestimate_mt_sd", reestimate_mt_sd)
        mtd_par.setValue("quant_method", quant_method)
        mtd_par.setValue("trace_termination_criterion", trace_termination_criterion)
        mtd_par.setValue("trace_termination_outliers", trace_termination_outliers)
        mtd_par.setValue("min_sample_rate", min_sample_rate)
        mtd_par.setValue("min_trace_length", min_trace_length)
        mtd_par.setValue("max_trace_length", max_trace_length)
        mtd.setParameters(mtd_par) # set the new parameters
        mtd.run(exp, mass_traces, 0) # run mass trace detection
        #print(mtd_par)
        # Elution peak detection

        print("Elution Peak Detection: ", filepath)
        mass_traces_deconvol = []
        epd = ElutionPeakDetection()
        epd_par = epd.getDefaults()
        epd_par.setValue("width_filtering", width_filtering) # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
        epd_par.setValue("chrom_fwhm", chrom_fwhm)
        epd_par.setValue("min_fwhm", min_fwhm) 
        epd_par.setValue("max_fwhm", max_fwhm) 
        #epd_par.setValue("enabled", enabled) 
        epd_par.setValue("masstrace_snr_filtering", masstrace_snr_filtering) 
        epd.setParameters(epd_par)
        epd.detectPeaks(mass_traces, mass_traces_deconvol)
        #print(epd_par)
        
        # Feature detection

        print("Feature Detection: ", filepath)
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
        feature_map_FFM.setPrimaryMSRunPath([filepath.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
        FeatureXMLFile().store(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(filepath)[:-5] + "_large.featureXML"), feature_map_FFM)

        #Export a df
        df = feature_map_FFM.get_df()
        df = df.rename(columns={'intensity': os.path.basename(filepath), 
                                'mz':'Mass [m/z]', 
                               'RT':'retention_time',
                               'RTstart':'rt_start',
                               'RTend':'rt_end'})
        df = df[['Mass [m/z]','retention_time','charge',os.path.basename(filepath),'rt_start','rt_end','quality']]
        df.to_csv(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(filepath)[:-5] + "_large.csv"), index=False)

        logger.info('======')
        logger.info('"Finished feature detection of large features"')
        logger.info('======')
        
        del mtd
        del epd
        del exp
        del df
        del ffm
    
    
    # 1) Feature Detection on narrow features
    def pyopenms_exclusion_narrow(filepath, mass_error_ppm, noise_threshold_int, chrom_peak_snr):
        # 1.1) Mass trace detection

        #mtd
        #mass_error_ppm = 10.0
        #noise_threshold_int = 3e5
        #chrom_peak_snr = 3.0

        chrom_fwhm = 10.0   #large=30  narrow = 10
        reestimate_mt_sd = 'true'
        quant_method = 'max_height'
        trace_termination_criterion = 'outlier'
        trace_termination_outliers =  3 # large = 3 narrow =2
        min_sample_rate = 0.3
        min_trace_length = 1.0   #large=10  narrow = 1
        max_trace_length = 60.0  #large=200000  narrow = 60

        #epd
        #enabled = 'true'
        width_filtering = 'fixed'
        min_fwhm = 0.0
        max_fwhm = 30.0
        masstrace_snr_filtering = 'false'

        #ffm
        local_rt_range = 7.0
        local_mz_range = 0.3
        charge_lower_bound = 1
        charge_upper_bound = 3
        report_summed_ints = 'true'
        enable_RT_filtering = 'false'
        isotope_filtering_model = 'none'
        mz_scoring_13C = 'false'
        use_smoothed_intensities = 'true'
        report_convex_hulls = 'true'
        remove_single_traces = 'false'
        mz_scoring_by_elements = 'false'
        elements = 'CHNOPS'

        print("Mass Trace Detection: ", filepath) #print the filename
        exp = MSExperiment()    
        MzMLFile().load(filepath, exp) # load each mzML file to an OpenMS file format (MSExperiment)

        mass_traces = [] # introduce an empty list where the mass traces will be loaded
        mtd = MassTraceDetection()
        mtd_par = mtd.getDefaults() # get the default parameters in order to edit them
        mtd_par.setValue("mass_error_ppm", float(mass_error_ppm)) # high-res instrument, orbitraps
        mtd_par.setValue("noise_threshold_int", float(noise_threshold_int)) # data-dependent (usually works for orbitraps)
        mtd_par.setValue("reestimate_mt_sd", reestimate_mt_sd)
        mtd_par.setValue("quant_method", quant_method)
        mtd_par.setValue("trace_termination_criterion", trace_termination_criterion)
        mtd_par.setValue("trace_termination_outliers", trace_termination_outliers)
        mtd_par.setValue("min_sample_rate", min_sample_rate)
        mtd_par.setValue("min_trace_length", min_trace_length)
        mtd_par.setValue("max_trace_length", max_trace_length)
        mtd.setParameters(mtd_par) # set the new parameters
        mtd.run(exp, mass_traces, 0) # run mass trace detection
        #print(mtd_par)
        # Elution peak detection

        print("Elution Peak Detection: ", filepath)
        mass_traces_deconvol = []
        epd = ElutionPeakDetection()
        epd_par = epd.getDefaults()
        epd_par.setValue("width_filtering", width_filtering) # The fixed setting filters out mass traces outside the [min_fwhm: 1.0, max_fwhm: 60.0] interval
        epd_par.setValue("chrom_fwhm", chrom_fwhm)
        epd_par.setValue("min_fwhm", min_fwhm) 
        epd_par.setValue("max_fwhm", max_fwhm) 
        #epd_par.setValue("enabled", enabled) 
        epd_par.setValue("masstrace_snr_filtering", masstrace_snr_filtering) 
        epd.setParameters(epd_par)
        epd.detectPeaks(mass_traces, mass_traces_deconvol)
        #print(epd_par)
        # Feature detection

        print("Feature Detection: ", filepath)
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
        feature_map_FFM.setPrimaryMSRunPath([filepath.encode()]) # Sets the file path to the primary MS run (usually the mzML file)
        FeatureXMLFile().store(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(filepath)[:-5] + "_narrow.featureXML"), feature_map_FFM)
        
        #Export a df
        df = feature_map_FFM.get_df()
        df = df.rename(columns={'intensity': os.path.basename(filepath), 
                                'mz':'Mass [m/z]', 
                               'RT':'retention_time',
                               'RTstart':'rt_start',
                               'RTend':'rt_end'})
        df = df[['Mass [m/z]','retention_time','charge',os.path.basename(filepath),'rt_start','rt_end','quality']]
        df.to_csv(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(filepath)[:-5] + "_narrow.csv"),index=False)

        del mtd
        del epd
        del exp
        del df
        del ffm
        
        logger.info('======')
        logger.info('"Finished feature detection of narrow features"')
        logger.info('======')
    
    patterns = ['raw', 'mzML', 'RAW', 'mzml']
    for pattern in patterns:
        input_mzML = input_mzML.replace(pattern, 'mzML')
    
    try:
        pyopenms_exclusion_narrow(os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)), ppm_error, narrow_noise_threshold, 3)
    except CalledProcessError as e:
        logger.info("!!! There was an error with OpenMS workflow for narrow features, please check your input files and parameters !!!")
        logger.info(e.output)
        raise    
        
    try:
        pyopenms_exclusion_large(os.path.join(OPENMS_folder+"/OpenMS_input/", os.path.basename(input_mzML)), ppm_error, large_noise_threshold, 3)
    except CalledProcessError as e:
        logger.info("!!! There was an error with OpenMS workflow for large features, please check your input files and parameters !!!")
        logger.info(e.output)
        raise    
        

    # Error with the OpenMS workflow. No output files.
    try:
        f = open(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(input_mzML)[:-5] + "_narrow.csv"))
        f.close()
        f = open(os.path.join(OPENMS_folder+"/OpenMS_output/", os.path.basename(input_mzML)[:-5] + "_large.csv"))
        f.close()
    except:
        logger.info('There was an issue with the pyOpenMS workflow ! See the log below.')
        #f = open(OPENMS_folder+'/OpenMS_output/TOPPAS.log', 'r')
        #file_contents = f.read()
        #logger.info(file_contents)
        raise

    gc.collect()
    
    logger.info('======')
    logger.info('Completed the pyOpenMS workflow')
    logger.info('======')
    logger.info('Zipping up the OpenMS workflow results ...')
    get_all_file_paths('OpenMS_workflow/','download_results/IODA_OpenMS_results.zip')

    logger.info('======')
    logger.info('NOW YOU CAN CONTINUE WITH THE REST OF THE WORKFLOW')

if __name__ == "__main__":
    IODA_run_OPENMS_exclusion(str(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))