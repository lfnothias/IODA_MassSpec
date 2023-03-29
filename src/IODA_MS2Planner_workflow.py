
# coding: utf-8
# Author: Louis Felix Nothias, louisfelix.nothias@gmail.com, June 2020
import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from io import StringIO
import warnings
from format_to_qexactive_list import *
from zipfile import ZipFile
from logzero import logger, logfile
import datetime
import subprocess
from subprocess import call
import pathlib

def process_input_table(input_filename: str, output_filename: str):
    """Take an MZmine3 or an mzTab containing two samples, output a table with mz, charge, rt, intensities."""
    
    if input_filename.endswith('.csv'):
        print("Input file is a MZmine3 feature table")
        
        # Process the MZmine3 table
        table = pd.read_csv(input_filename, index_col=False, on_bad_lines='skip')
            
        cols_to_drop = ['row ID', 'row ion mobility','row ion mobility unit','row CCS','correlation group ID','annotation network number','partners','neutral M mass', 'auto MS2 verify','identified by n=']
        for col in cols_to_drop:
            if col in table.columns:
                table.drop(col, axis=1, inplace=True)
            else:
                print(f"Column {col} not found. This supports MZmine3 feature table")
        
        table.drop(table.columns[table.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
        # Replace MZmine adducts into charge
        def replace_best_ion(best_ion):
            if pd.isnull(best_ion) or best_ion != best_ion:
                return 0
            if ']2+' in best_ion:
                return '2'
            elif ']+' in best_ion:
                return '1'
            elif ']2-' in best-ion:
                return '-1'
            elif ']1-' in best-ion:
                return '-1'
            else:
                return '0'

        table['best ion'] = table['best ion'].apply(lambda x: replace_best_ion(x))
        table = table.rename(columns={'row m/z': 'Mass [m/z]'})
        table = table.rename(columns={'best ion': 'charge'})
        table = table.rename(columns={'row retention time': 'retention_time'})
        #Convert min into seconds for MS2Planner
        table['retention_time'] = table['retention_time']*60

        #Detection of blank
        #print('#Deducing the blank sample by comparing the sum of feature intensity between samples')

        # Get all columns containing 'Peak area' in their name
        peak_area_cols = table.filter(like='Peak area')

        # Loop through each column and calculate its sum
        sums = {}
        for col in peak_area_cols:
            sums[col] = table[col].sum()

        # Find the column with the highest sum
        highest_sum_col = max(sums, key=sums.get)
        smallest_sum_col = min(sums, key=sums.get)

        logger.info('- For sample ' + str(smallest_sum_col) + ', the sum of feature intensities = ' + '{:.2e}'.format(sums[smallest_sum_col]))
        logger.info('- For sample ' + str(highest_sum_col) + ', the sum of feature intensities = ' + '{:.2e}'.format(sums[highest_sum_col]))
        logger.info('- The blank sample is assumed to be '+str(smallest_sum_col)+' in the feature table')
        logger.info('- The sample is assumed to be '+str(highest_sum_col)+' in the feature table')
        table = table.fillna(0).sort_values('retention_time')

        # Define the new order of columns
        new_order = ['Mass [m/z]', 'retention_time', 'charge', smallest_sum_col, highest_sum_col]

        # Reorder the columns using reindex
        table = table.reindex(columns=new_order)

        table.to_csv(output_filename, sep=',', index=False)
        
        return highest_sum_col, smallest_sum_col
    
    elif input_filename.endswith('.mzTab'):
        print("Input file is an mzTab")

        df = pd.read_csv(input_filename, sep='\t', on_bad_lines='skip')
        # rest of your code goes here

        # Get the metadata
        metadata = []
        start_row_consensus = 1
        for row in df['1.0.0']:
            metadata.append(row)
            start_row_consensus += 1

        # Change the type of the list
        metadata = [str(item) for item in metadata]
        [type(item) for item in metadata]

        # Get the filenames
        Filenames = []
        for x in metadata:
            if x.startswith('file:/'):
                x = x.split('/')[-1]
                #Remove duplicates
                if x not in Filenames:
                    Filenames.append(x[:-5])

        Filename1 = Filenames[0]
        Filename2 = Filenames[1]

        # Display error message for additional samples
        for x in Filenames:
            if x == "ms_run[2]-location":
                Filename3 = Filenames[2]
                logger.info('Warning: There is more than two samples in that mzTab file. We support only two samples currently')
            if x == "ms_run[4]-location":
                Filename4 = Filenames[3]
                logger.info('Warning: There is more than three samples in that mzTab file. We support only two samples currently.')

        # Read and edit the table
        main_df = pd.read_csv(input_filename,sep='\t',index_col=0, skiprows=range(0,start_row_consensus))

        # Get the columns of interest
        feat_int = main_df.loc[:, main_df.columns.str.startswith('peptide_abundance_study_variable')]
        feat_mz = main_df.loc[:, main_df.columns.str.startswith('mass_to_charge')]
        feat_charge = main_df.loc[:, main_df.columns.str.startswith('charge')]
        feat_ret = main_df[['retention_time']]

        # Concat into a master table
        table = pd.concat([feat_mz,feat_ret,feat_charge,feat_int], axis=1)

        #Detection of blank
        #print('#Deducing the blank sample by comparing the sum of feature intensity between samples')
        column1_sum = table['peptide_abundance_study_variable[1]'].sum()
        logger.info('- For sample '+Filename1+' the sum of feature intensities is = '+str("{:.2e}".format(column1_sum)))
        column2_sum = table['peptide_abundance_study_variable[2]'].sum()
        logger.info('- For sample '+Filename2+' the sum of feature intensities = '+str("{:.2e}".format(column2_sum)))
        if column1_sum > column2_sum:
        #    logger.info('- The blank sample is assumed to be '+str(Filename2)+' in the mzTab-M')
        #    logger.info('- The samples is assumed to be '+str(Filename1)+' in the mzTab-M')
            table.rename(columns={'peptide_abundance_study_variable[1]':Filename2}, inplace=True)
            table.rename(columns={'peptide_abundance_study_variable[2]':Filename1}, inplace=True)
            highest_sum_col = Filename1
            smallest_sum_col = Filename2
        elif column1_sum <= column2_sum:
        #    logger.info('- The blank sample is assumed to be '+str(Filename1)+' in the mzTab-M')
        #    logger.info('- The samples is assumed to be '+str(Filename2)+' in the mzTab-M')
            table.rename(columns={'peptide_abundance_study_variable[1]':Filename1}, inplace=True)
            table.rename(columns={'peptide_abundance_study_variable[2]':Filename2}, inplace=True)
            highest_sum_col = Filename2
            smallest_sum_col = Filename1
        #Replace the sample headers for mandatory samples
        table.rename(columns={'mass_to_charge':"Mass [m/z]"}, inplace=True)

        table = table.fillna(0).sort_values('retention_time')
        table.to_csv(output_filename, sep=',', index=False)
        
        return highest_sum_col, smallest_sum_col

def get_all_file_paths(directory,output_zip_path):
    # initializing empty file paths list
    file_paths = []

    # crawling through directory and subdirectories
    for root, directories, files in os.walk(directory):
        for filename in files:
            # join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)

            # writing files to a zipfile
    with ZipFile(output_zip_path,'w') as zip:
        # writing each file one by one
        for file in file_paths:
            zip.write(file)

    logger.info('All files zipped successfully!')

# Here we limit the number of feature with the same RT to the value max_same_RT
def limit_number_of_RT_same_RT(feature_table, samplename, max_same_RT, output_filename):
    initial_nb_features = feature_table.shape[0]-1
    feature_table = feature_table.sort_values(by=['retention_time', samplename], ascending=[True, False])
    feature_table = feature_table.groupby('retention_time').head(max_same_RT)
    feature_table.to_csv(output_filename[:-4]+'_filtered.csv', sep=',', index=False)
    afterfiltering_nb_features = feature_table.shape[0]-1
    logger.info('   Remaining features = '+str(afterfiltering_nb_features)+' after same RT filtering with top '+str(max_same_RT))
        
    return feature_table


def filter_feature_table(df_master, sample, blank, min_intensity_value, ratio):
    df_master_targeted = df_master[(df_master[sample] != 0)]
    df_master_targeted_int = df_master[(df_master[sample] > min_intensity_value)]
    df_master_targeted_ratio = df_master[(df_master[sample]/df_master[blank] > ratio)]
    df_master_targeted_filtered = df_master[(df_master[sample] > min_intensity_value) & (df_master[sample]/df_master[blank] > ratio)]

    return df_master_targeted_filtered

def count_files_with_pattern(directory, pattern):
    count = 0
    for filename in os.listdir(directory):
        if pattern in filename:
            count += 1
    return count

def make_plot_MS2Planner(output_filename, directory, count, mode):
    table_list_MS2Planner = []
    for x in range(count):
        table_list_MS2Planner.append(
            output_filename[:-4]+'_filtered_MS2Planner_'+mode+'_path_'+str(x+1)+'.csv')

    try:
        # print(table_list_MS2Planner)
        make_plot_MS2Planner_RT_mz(table_list_MS2Planner, directory)
        make_plot_MS2Planner_mz_int(table_list_MS2Planner, directory)
        make_plot_MS2Planner_RT_int(table_list_MS2Planner, directory)
    except:
        raise
    
def generate_outputs(output_dir, mode):
    # Cleaning files first
    logger.info('Cleaning and zipping workflow results files ...')
    #mkdir Exactive
    os.system('mkdir '+output_dir+'/Exactive')
    # mv files Exactive
    os.system('mv '+output_dir+'/*_QE* '+output_dir+'/Exactive')

    #mkdir MQL
    os.system('mkdir '+output_dir+'/MaxQuantLive')

    # mv files MQL
    os.system('mv '+output_dir+'/*_MQL* '+output_dir+'/MaxQuantLive')

    #mkdir Exploris
    os.system('mkdir '+output_dir+'/Exploris_DDA_MS2')
    os.system('mkdir '+output_dir+'/Exploris_tMS2')
    # mv files  Exploris
    os.system('mv '+output_dir+'/*_Exploris_DDA_MS2* '+output_dir+'/Exploris_DDA_MS2')
    os.system('mv '+output_dir+'/*_Exploris_tMS2* '+output_dir+'/Exploris_tMS2')

    # mkdir intermediate files
    os.system('mkdir '+output_dir+'/intermediate_files')
    os.system('mkdir '+output_dir+'/plots')
    os.system('mkdir '+output_dir+'/log')

    # mv
    os.system('mv '+output_dir+'/*_plot_* '+output_dir+'/plots')
    os.system('mv '+output_dir+'/logfile.txt '+output_dir+'/log')
    os.system('cp MS2Planner.log '+output_dir+'/log/MS2Planner_'+mode+'.log')
    os.system('mv '+output_dir+'/*.csv '+output_dir+'/intermediate_files')

    get_all_file_paths(output_dir,'download_'+output_dir+'/IODA_MS2Planner_'+mode+'_results.zip')

    logger.info('======')
    logger.info('END OF THE MS2Planner WORKFLOW - '+mode+' mode')
    logger.info('======')
    print(' ')

# Run the MS2Planner workflow with baseline method

def MS2Planner_baseline(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, 
                        win_len:float, isolation:float, delay:float, rt_margin:float, 
                        max_same_RT:float, transient_time:float,
                        polarity:str, apex_int_percent:float=0.6,
                        pretarget_rt_margin:float=0, posttarget_rt_margin:float=0, 
                        RF_base_value:int=np.nan, CEs:str or List[str] = np.nan
                        ):

    output_dir = 'results_targeted_MS2Planner_baseline'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')

    logger.info('STARTING THE MS2Planner WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename.endswith('.csv'):
        logger.info('File path to a MZmine3 feature table was specified by the user')
        pass
    elif input_filename.endswith('.mzTab'):
        logger.info('File path to a mzTab was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
        print(input_filename)
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')

    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))
    
    elif input_filename.endswith('.csv'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1]
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))

    elif input_filename.endswith('.mzTab'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))


    # Convert to Table
    logger.info('======')
    logger.info('Converting to intermediate table format ...')
    samples = process_input_table(input_filename, output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = str(samples[0])
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = str(samples[1])
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('= PARAMETERS ==')
    logger.info(' ')
    logger.info('== MS2Planner parameters')
    ratio = intensity_ratio
    logger.info('    Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('    Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('    Retention time window (secs) for binning target ions = ' +str(win_len))
    logger.info('    Isolation window (m/z) = ' +str(isolation))
    logger.info('    Transient time (ms) = ' +str(transient_time))
    logger.info('    Retention time margin (sec.) = ' +str(delay))
    experiments = num_path
    logger.info('    Number of iterative experiment(s) = ' + str(experiments))
    logger.info(' ')
    logger.info('== Output target list parameters')
    logger.info('    Polarity = ' + str(polarity))
    logger.info('    Exactive serie retention time window for target ion list (sec. = ' + str(rt_margin*2)+')')
    logger.info('    MQL retention time window for target ion list = apex - ' + str(rt_margin)+' secs)')
    logger.info('    Exploris: Pretarget retention time margin for target ion list (sec. = ' + str(pretarget_rt_margin)+')')
    logger.info('    Exploris: Posttarget retention time margin for target ion list (sec. = ' + str(posttarget_rt_margin)+')')
    logger.info('    Exploris serie RF base value = ' + str(RF_base_value))
    logger.info('    Exploris serie CEs = ' + str(CEs))

    logger.info(' ')
    logger.info('== Pre-filtering input table ===')
    logger.info('   Number of features in the table = '+ str(feature_table.shape[0]-1))
    feature_table = filter_feature_table(feature_table, samplename, blank_samplename, min_intensity, ratio)
    logger.info('   Remaining features = '+str(feature_table.shape[0]-1)+ ' after the intensity/ratio filtering')
    feature_table = limit_number_of_RT_same_RT(feature_table, samplename, max_same_RT, output_filename)
    logger.info(' ')

    # Running the table processing
    logger.info('Running MS2Planner in Baseline mode ...')
    #Clean up the log
    try:
        f = open('MS2Planner.log', 'w')
        f.truncate(0)
        f.close()
    except:
        pass
    try:
        run_MS2Planner_baseline(output_filename[:-4]+'_filtered.csv', output_filename[:-4]+'_filtered_MS2Planner_baseline.csv', 
                                intensity_threshold, intensity_ratio, num_path, win_len, isolation, delay)
    except:
        logger.info('There was an issue with the MS2Planner ! See the log below.')
        f = open('MS2Planner.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise
    logger.info('======')

    f = open('MS2Planner.log', 'r')
    file_contents = f.read()
    logger.info(file_contents)

    Test_MS2Planner_Output = pathlib.Path(output_filename[:-4]+'_filtered_MS2Planner_baseline_path_1.csv')

    try:
        if Test_MS2Planner_Output.exists ():
            logger.info("MS2Planner Path output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running MS2Planner Path !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running MS2Planner Path !!!")
    except:
        raise

    logger.info('Preparing results ...')
    directory = os.path.dirname(output_filename)  # Replace with the path to your folder
    pattern = 'MS2Planner_baseline_path'  # Replace with the pattern you want to search for
    count = count_files_with_pattern(directory, pattern)

    min_scan = 0 #Hardcoded to keep the same def function with Curve mode. Parameter only used in Curve mode.

    ## GENERATE OUTPUT TABLES

    #Format for Exactive
    for i in range(count):
        generate_Exactive_table_from_MS2Planner(
        output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'.csv', 
        output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'_QE.csv',
        polarity=polarity, rt_margin=rt_margin)

    #Format for MaxQuant.Live
    for i in range(count):
        generate_MQL_tMS2_table_from_MS2Planner(
                    output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'.csv', 
                    output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'_MQL.txt',
                    rt_margin=rt_margin, delay=delay, 
                    transient_time=transient_time, 
                    polarity=polarity, 
                    apex_int_percent=apex_int_percent)
    
    #Format for Exploris
    for i in range(count):
        generate_Exploris_DDAMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'_Exploris_DDA_MS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, transient_time=transient_time,
            polarity=polarity, CEs=CEs, apex_int_percent=apex_int_percent)

    for i in range(count):
        generate_Exploris_tMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_baseline_path_'+str(i+1)+'_Exploris_tMS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, transient_time=transient_time,
            polarity=polarity, RF_base_value=RF_base_value, CEs=CEs)
    
    make_plot_MS2Planner(output_filename, directory, count, 'baseline')

    logger.info('======')
    generate_outputs(output_dir, 'baseline')

    print(' ')


# Run the MS2Planner workflow with apex method
def MS2Planner_apex(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, 
                    intensity_accu:float, isolation:float, delay:float, min_scan:float, max_scan:float, 
                    rt_margin:float, max_same_RT:int, transient_time:float,
                    polarity:str, apex_int_percent:float, 
                    pretarget_rt_margin:float=0, posttarget_rt_margin:float=0, 
                    RF_base_value:float=np.nan, CEs:str or List[str]=np.nan
):

    output_dir = 'results_targeted_MS2Planner_apex'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')


    logger.info('STARTING THE MS2Planner WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename.endswith('.csv'):
        logger.info('File path to a MZmine3 feature table was specified by the user')
        pass
    elif input_filename.endswith('.mzTab'):
        logger.info('File path to a mzTab was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
        print(input_filename)
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')

    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))
    
    elif input_filename.endswith('.csv'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1]
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))

    elif input_filename.endswith('.mzTab'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))


    # Convert to Table
    logger.info('======')
    logger.info('Converting to intermediate table format ...')
    samples = process_input_table(input_filename, output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = str(samples[0])
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = str(samples[1])
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('= PARAMETERS ==')
    logger.info(' ')
    logger.info('== MS2Planner parameters')
    logger.info(' ')
    ratio = intensity_ratio
    logger.info('    Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('    Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('    Precursor ion intensity to accumulate in the MS2 scan = ' +str("{:.2e}".format(intensity_accu)))
    logger.info('    Isolation window (m/z) = ' +str(isolation))
    logger.info('    Transient time (msec.) = ' +str(transient_time))
    logger.info('    Retention time margin (sec.) = ' +str(delay))
    logger.info('    Delay between targeted MS2 scans (sec.)= ' +str(delay))
    logger.info('    Minimum MS2 scan duty cycle (sec.)= ' +str(min_scan))
    logger.info('    Maximum MS2 scan duty cycle (sec.)= ' +str(max_scan))
    experiements = num_path
    logger.info('    Number of iterative experiment(s) = ' + str(experiements))

    logger.info('== Output target list parameters')
    logger.info('    Polarity = ' + str(polarity))
    logger.info('    Orbitrap transient + overhead time for MaxQuant.Live (msec.) = ' + str(transient_time))
    logger.info('    Exactive serie retention time window for target ion list (sec. = ' + str(rt_margin*2)+')')
    logger.info('    MQL retention time window for target ion list = apex - ' + str(rt_margin)+' secs)')
    logger.info('    Exploris: Pretarget retention time margin for target ion list (sec. = ' + str(pretarget_rt_margin)+')')
    logger.info('    Exploris: Posttarget retention time margin for target ion list (sec. = ' + str(posttarget_rt_margin)+')')
    logger.info('    Exploris serie RF base value = ' + str(RF_base_value))
    logger.info('    Exploris serie CEs = ' + str(CEs))

    logger.info(' ')
    logger.info('== Pre-filtering input table ===')
    logger.info('   Number of features in the table = '+ str(feature_table.shape[0]-1))
    feature_table = filter_feature_table(feature_table, samplename, blank_samplename, min_intensity, ratio)
    logger.info('   Remaining features = '+str(feature_table.shape[0]-1)+ ' after the intensity/ratio filtering')
    feature_table = limit_number_of_RT_same_RT(feature_table, samplename, max_same_RT, output_filename)
    logger.info(' ')
    
    # Running the table processing
    logger.info('Running MS2Planner in Apex mode ...')
    #Clean up the log
    try:
        f = open('MS2Planner.log', 'w')
        f.truncate(0)
        f.close()
    except:
        pass
    try:
        run_MS2Planner_apex(output_filename[:-4]+'_filtered.csv', output_filename[:-4]+'_filtered_MS2Planner.csv', intensity_threshold,
                             intensity_ratio, num_path, intensity_accu, isolation, delay, min_scan, max_scan)
    except:
        logger.info('There was an issue with the MS2Planner ! See the log below.')
        f = open('MS2Planner.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    logger.info('======')

    f = open('MS2Planner.log', 'r')
    file_contents = f.read()
    logger.info(file_contents)

    logger.info('======')

    Test_MS2Planner_Output = pathlib.Path(output_filename[:-4]+'_filtered_MS2Planner_apex_path_1.csv')
    try:
        if Test_MS2Planner_Output.exists ():
            logger.info("MS2Planner output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running MS2Planner !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running MS2Planner !!!")
    except:
        raise
    
    logger.info('Preparing results ...')
    directory = os.path.dirname(output_filename)  # Replace with the path to your folder
    pattern = 'MS2Planner_apex_path'  # Replace with the pattern you want to search for
    count = count_files_with_pattern(directory, pattern)

    ## GENERATE OUTPUT TABLES
    
    #Format for Exactive
    for i in range(count):
        generate_Exactive_table_from_MS2Planner(
        output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'.csv', 
        output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'_QE.csv',
        polarity=polarity, rt_margin=rt_margin)

    #Format for MaxQuant.Live
    for i in range(count):
        generate_MQL_tMS2_table_from_MS2Planner(
                    output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'.csv', 
                    output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'_MQL.txt',
                    rt_margin=rt_margin, delay=delay, transient_time=transient_time, polarity=polarity, apex_int_percent=apex_int_percent)
    
    #Format for Exploris
    for i in range(count):
        generate_Exploris_DDAMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'_Exploris_DDA_MS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, transient_time=transient_time,
            polarity=polarity, CEs=CEs, apex_int_percent=apex_int_percent)

    for i in range(count):
        generate_Exploris_tMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_apex_path_'+str(i+1)+'_Exploris_tMS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, 
            transient_time=transient_time, polarity=polarity,
            RF_base_value=RF_base_value, CEs=CEs)
    
    logger.info('======')
    make_plot_MS2Planner(output_filename,directory,count, 'apex')

    generate_outputs(output_dir, 'apex')
    print(' ')



# Run the MS2Planner workflow with CURVE method
def MS2Planner_curve(input_filename:int, num_path:int, intensity_ratio:float, intensity_threshold:float, 
                        input_filename_curve:int, intensity_accu:float, rt_tolerance_curve:float, mz_tolerance_curve:float, 
                        isolation:float, delay:float, min_scan:float, max_scan:float, cluster:str, 
                        rt_margin:float, transient_time:float, max_same_RT:float, 
                        apex_int_percent:float, polarity:str,
                        pretarget_rt_margin:float, posttarget_rt_margin:float, 
                        RF_base_value:float=np.nan, CEs:str or List[str]=np.nan
):

    output_dir = 'results_targeted_MS2Planner_curve'
    os.system('rm -r '+output_dir)
    os.system('rm -r download_'+output_dir)
    os.system('mkdir '+output_dir)
    os.system('mkdir download_'+output_dir)
    logfile(output_dir+'/logfile.txt')

    logger.info('STARTING THE MS2Planner WORKFLOW')
    if input_filename.startswith('http'):
        logger.info('File path was specified by the user')
        pass
    elif input_filename.endswith('.csv'):
        logger.info('File path to a MZmine3 feature table was specified by the user')
        pass
    elif input_filename.endswith('.mzTab'):
        logger.info('File path to a mzTab was specified by the user')
        pass
    elif input_filename == 'OpenMS_generated':
        logger.info('The mzTab was generated with the IODA-OpenMS workflow')
        path_input_folder = "TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/"
        mzTab_file = os.listdir("TOPPAS_Workflow/toppas_output/TOPPAS_out/Targeted_MzTab/")[0]
        input_filename = path_input_folder+mzTab_file
        print(input_filename)
    else:
        logger.info("the input_filename variable should be a valid path/download link or must be: 'OpenMS_generated', when using the OpenMS workflow online")

    now = datetime.datetime.now()
    logger.info(now)

    logger.info('======')
    logger.info('Getting the mzTab')

    if input_filename.startswith('http'):
        if 'google' in input_filename:
            logger.info('This is the Google Drive download link: '+str(input_filename))
            url_id = input_filename.split('/', 10)[5]
            prefixe_google_download = 'https://drive.google.com/uc?export=download&id='
            input_filename = prefixe_google_download+url_id
            output_filename = output_dir+'/Converted_mzTab.csv'

        else:
            output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
            logger.info('This is the input file path: '+str(input_filename))
            logger.info('This is the output file path: '+str(output_filename))
    
    elif input_filename.endswith('.csv'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1]
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))

    elif input_filename.endswith('.mzTab'):
        output_filename = output_dir+'/'+input_filename.split('/', 10)[-1][:-6]+'.csv'
        logger.info('This is the input file path: '+str(input_filename))
        logger.info('This is the output file path: '+str(output_filename))


    # Convert to Table
    logger.info('======')
    logger.info('Converting to intermediate table format ...')
    samples = process_input_table(input_filename, output_filename)
    logger.info('======')

    # Read the table to get the filenames
    feature_table = pd.read_csv(output_filename)
    samplename = str(samples[0])
    logger.info('Assumed sample filename: '+samplename)
    blank_samplename = str(samples[1])
    logger.info('Assumed blank filename: ' +blank_samplename)
    logger.info('======')

    # User-defined parameters
    logger.info('= PARAMETERS ==')
    logger.info(' ')
    logger.info('== MS2Planner parameters')
    ratio = intensity_ratio
    logger.info('    Ratio between sample/blank for ion filtering = ' + str(ratio))
    min_intensity = intensity_threshold
    logger.info('    Minimum intensity for ion filtering in sample = '+ str("{:.2e}".format(min_intensity)))
    logger.info('    Precursor ion intensity to accumulate in the MS2 scan = ' +str("{:.2e}".format(intensity_accu)))
    logger.info('    Input file for curve data : ' +str(input_filename_curve))
    logger.info('    Restriction parameter : ' +str(rt_tolerance_curve))
    logger.info('    Mass accuracy (m/z): ' +str(mz_tolerance_curve))
    logger.info('    Isolation window (m/z) = ' +str(isolation))
    logger.info('    Delay between targeted MS2 scans (sec.)= ' +str(delay))
    logger.info('    Transient time (msec.) = ' +str(transient_time))
    logger.info('    Retention time margin (sec.) = ' +str(delay))
    logger.info('    Minimum MS2 scan duty cycle (sec.)= ' +str(min_scan))
    logger.info('    Maximum MS2 scan duty cycle (sec.)= ' +str(max_scan))
    experiements = num_path
    logger.info('    Number of iterative experiment(s) = ' + str(experiements))
    logger.info('    Mode for the curve mode: '+str(cluster))
    logger.info(' ')
    logger.info('== Output target list parameters')
    logger.info('    Polarity = ' + str(polarity))
    logger.info('    Orbitrap transient + overhead time for MaxQuant.Live (msec.) = ' + str(transient_time))
    logger.info('    Exactive serie retention time window for target ion list (sec. = ' + str(rt_margin*2)+')')
    logger.info('    MQL retention time window for target ion list = apex - ' + str(rt_margin)+' secs)')
    logger.info('    Exploris: Pretarget retention time margin for target ion list (sec. = ' + str(pretarget_rt_margin)+')')
    logger.info('    Exploris: Posttarget retention time margin for target ion list (sec. = ' + str(posttarget_rt_margin)+')')
    logger.info('    Exploris serie RF base value = ' + str(RF_base_value))
    logger.info('    Exploris serie CEs = ' + str(CEs))
    logger.info(' ')
    logger.info('=== FILTERING THE INPUT FEATURE TABLE ===')
    logger.info('   Number of features in the table = '+ str(feature_table.shape[0]-1))
    feature_table = filter_feature_table(feature_table, samplename, blank_samplename, min_intensity, ratio)
    logger.info('   Remaining features = '+str(feature_table.shape[0]-1)+ ' after the intensity/ratio filtering')
    feature_table = limit_number_of_RT_same_RT(feature_table, samplename, max_same_RT, output_filename)
    logger.info(' ')

    #Running MS2Planner
    logger.info('Running MS2Planner in Curve mode ...')
    #Clean up the log
    #Clean up the log
    try:
        f = open('MS2Planner.log', 'w')
        f.truncate(0)
        f.close()
    except:
        pass

    try:
        run_MS2Planner_curve(output_filename[:-4]+'_filtered.csv', output_filename[:-4]+'_filtered_MS2Planner.csv', intensity_threshold, 
                             intensity_ratio, num_path, input_filename_curve, intensity_accu, rt_tolerance_curve, mz_tolerance_curve, isolation,
                               delay, min_scan, max_scan, cluster)
    except:
        logger.info('There was an issue with the MS2Planner ! See the log below.')
        f = open('MS2Planner.log', 'r')
        file_contents = f.read()
        logger.info(file_contents)
        raise

    logger.info('======')

    f = open('MS2Planner.log', 'r')
    file_contents = f.read()
    logger.info(file_contents)

    Test_MS2Planner_Output = pathlib.Path(output_filename[:-4]+'_filtered_MS2Planner_curve_path_1.csv')
    try:
        if Test_MS2Planner_Output.exists ():
            logger.info("MS2Planner output found")
        else:
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            print("Problem when running MS2Planner !!!")
            print("<---------- !!!!!!!!!!!!!!!!!! ---------->")
            logger.info("Problem when running MS2Planner !!!")
    except:
        raise

    logger.info('======')
    logger.info('Preparing results ...')
    
    directory = os.path.dirname(output_filename)  # Replace with the path to your folder
    pattern = 'MS2Planner_curve_path'  # Replace with the pattern you want to search for
    count = count_files_with_pattern(directory, pattern)

    ## GENERATE OUTPUT TABLES
    
    #Format for Exactive
    for i in range(count):
        generate_Exactive_table_from_MS2Planner(
        output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'.csv', 
        output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'_QE.csv',
        polarity=polarity, rt_margin=rt_margin)

    #Format for MaxQuant.Live
    for i in range(count):
        generate_MQL_tMS2_table_from_MS2Planner(
                    output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'.csv', 
                    output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'_MQL.txt',
                    rt_margin=rt_margin, delay=delay, transient_time=transient_time, polarity=polarity, apex_int_percent=apex_int_percent)
    
    #Format for Exploris
    for i in range(count):
        generate_Exploris_DDAMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'_Exploris_DDA_MS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, transient_time=transient_time,
            polarity=polarity, CEs=CEs, apex_int_percent=apex_int_percent)

    for i in range(count):
        generate_Exploris_tMS2_table_from_MS2Planner(
            output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'.csv', 
            output_filename[:-4]+'_filtered_MS2Planner_curve_path_'+str(i+1)+'_Exploris_tMS2.txt',
            pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin, 
            transient_time=transient_time, polarity=polarity,
            RF_base_value=RF_base_value, CEs=CEs)
    
    make_plot_MS2Planner(output_filename,directory,count, 'curve')

    generate_outputs(output_dir, 'curve')


### MS2Planner
# This parse one line from the MS2Planner output and create a output table per path. The rows to skip define which line/path is parsed.
def MS2Planner_format(input_filename: str, output_filename: str, rows_to_skip:int):
    df_path = pd.read_csv(input_filename, sep=' ', header=None, skiprows=rows_to_skip, on_bad_lines='skip')

    #Make a list for the first row
    df_path_list = df_path.iloc[0].values.tolist()
    df_path_list.pop(0)
    nfeatures = int(len(df_path_list)/8)

    # Convert the list into a nested list
    target_list = []
    for entries in range(nfeatures):
            try:
                while len(df_path_list) > 8:
                    target_list.append(df_path_list[:8])
                    indexes = [0,1,2,3,4,5,6,7]
                    for index in sorted(indexes, reverse=True):
                        del df_path_list[index]
            except:
                continue
    #Make a dataframe
    target_table = pd.DataFrame(target_list)
    target_table = target_table.rename(columns={0: 'Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end', 5: 'intensity', 6: 'rt_apex',7: 'charge'})
    logger.info('Valid target ions in path'+str(rows_to_skip+1)+' = '+str(target_table.shape[0]))
    target_table.to_csv(output_filename, sep=',', index=False)

# This parse MS2Planner output file and create output tables formatted for XCalibur and MaxQuant.live
def make_MS2Planner_targeted_lists_from_table(
        input_filename:str,
        pretarget_rt_margin:float = 0,
        posttarget_rt_margin:float = 0,
        rt_margin:float = 0, 
        delay:float = 0,
        transient_time:float = 0,
        polarity:str = 'Positive',
        apex_int_percent:float=0.6, 
        RF_base_value:int = ''
):
    #Format for Exactive
    generate_Exactive_table_from_MS2Planner(input_filename, input_filename[:-4]+'_QE.csv', polarity=polarity, rt_margin=rt_margin)

    #Format for MaxQuant.Live
    generate_MQL_tMS2_table_from_MS2Planner(input_filename, input_filename[:-4]+'_MQL.txt', rt_margin=rt_margin, delay=delay, transient_time=transient_time, 
                                            polarity=polarity, apex_int_percent=apex_int_percent)
    
    #Format for Exploris
    generate_Exploris_DDAMS2_table_from_MS2Planner(input_filename, input_filename[:-4]+'_Exploris_DDA_MS2.csv',  pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin,
                                                   delay=delay, polarity=polarity, CEs=CEs, apex_int_percent=apex_int_percent)

    generate_Exploris_tMS2_table_from_MS2Planner(input_filename, input_filename[:-4]+'_Exploris_tMS2.csv',  pretarget_rt_margin=pretarget_rt_margin,  posttarget_rt_margin=posttarget_rt_margin,
                                                 RF_base_value=RF_base_value, CEs=CEs, polarity=polarity)
    
def run_MS2Planner_baseline(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, win_len:float, isolation:float, delay:float):
    cmd_baseline = ('python3 MS2Planner/path_finder.py baseline '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -win_len '+str(win_len)+' -isolation '+str(isolation)+' -delay '+str(delay))
    logger.info('Command: '+cmd_baseline)
    os.system(cmd_baseline)

def run_MS2Planner_apex(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, intensity_accu:float, isolation:float, delay:float, min_scan:float, max_scan:float):
    cmd_apex = ('python3 MS2Planner/path_finder.py apex '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -intensity_accu '+str(intensity_accu)+' -isolation '+str(isolation)+' -delay '+str(delay)+' -min_scan '+str(min_scan)+' -max_scan '+str(max_scan))
    logger.info('Command: '+cmd_apex)
    os.system(cmd_apex)

def run_MS2Planner_curve(input_filename:str, output_filename:str, intensity_threshold:float, intensity_ratio:float, num_path:int, input_filename_curve:str, intensity_accu:float, rt_tolerance_curve:float, mz_tolerance_curve:float, isolation:float, delay:float, min_scan:float, max_scan:float, cluster:str):
    cmd_curve = ('python3 MS2Planner/path_finder.py curve '+input_filename+' '+output_filename+' '+str(intensity_threshold)+' '+str(intensity_ratio)+' '+str(num_path)+' -infile_raw '+str(input_filename_curve)+' -intensity_accu '+str(intensity_accu)+' -restriction '+str(rt_tolerance_curve)+' '+str(mz_tolerance_curve)+' -isolation '+str(isolation)+' -delay '+str(delay)+' -min_scan '+str(min_scan)+' -max_scan '+str(max_scan)+' -cluster '+str(cluster))
    logger.info('Command: '+cmd_curve)
    logger.info('MS2Planner in Curve mode can take up to 10 minutes to complete ... please wait')
    logger.info('Output of MS2Planner: '+output_filename)
    try:
        cp0 = subprocess.run(cmd_curve,shell=True)
        cp0
    except subprocess.CalledProcessError:
        logger.info('ERROR running MS2Planner ...')

#MS2Planner generate mz / rt figures



def make_plot_MS2Planner_RT_mz(table_list_MS2Planner, output_dir):
    colors = ['blue', 'violet', 'orange', 'red', 'yellow','green','purple','black','brown','pink','grey']
    labels = []
    scaling_factor = 50
    min_size = 2

    #labels.append('Intensity scaling factor: '+str(scaling_factor))

    for i, table_path in enumerate(table_list_MS2Planner):
        table = pd.read_csv(table_path, sep=',', header=0)
        #print(table_path)

        sizes = (table['intensity'] / table['intensity'].max()) * scaling_factor
        sizes = sizes.apply(lambda x: max(x, min_size))

        plt.scatter('rt_apex', 'Mass [m/z]', data=table, marker='o', color=colors[i % len(colors)],  edgecolors='k', linewidths=0.2, s=sizes, alpha=0.8)

        label = f"Inj. {i + 1}, n = {table.shape[0]}, median = {table['intensity'].median():.2e}, mean = {table['intensity'].mean():.2e}"
        labels.append(label)
        plt.ylabel('Mass [m/z]')
        plt.xlabel('Ret. time apex (s)')
        plt.legend(labels=labels, fontsize=4)
        plt.title('Distribution of MS2Planner targets. Node size shows intensity with a scale factor of '+str(scaling_factor), fontsize=8)
        plt.savefig(table_path[:-4]+'_plot_mz_rt.png', dpi=300)


    plt.savefig(table_list_MS2Planner[0][:-6]+'_ALL_plot_mz_rt.png', dpi=300)
    plt.close()

def make_plot_MS2Planner_mz_int(table_list_MS2Planner, output_dir):
    colors = ['blue', 'violet', 'orange', 'red', 'yellow','green','purple','black','brown','pink','grey']
    labels = []

    for i, table_path in enumerate(table_list_MS2Planner):
        table = pd.read_csv(table_path, sep=',', header=0)

        plt.scatter('Mass [m/z]', 'intensity', data=table, marker='o', color=colors[i % len(colors)],  edgecolors='k', linewidths=0.2, s=2, alpha=0.8)

        label = f"Inj. {i + 1}, n = {table.shape[0]}, median = {table['intensity'].median():.2e}, mean = {table['intensity'].mean():.2e}"
        labels.append(label)
        
        plt.yscale("log")
        plt.ylabel('Intensity')
        plt.xlabel('Mass [m/z]')
        plt.legend(labels=labels, fontsize=4)
        plt.title('Distribution of MS2Planner targets.', fontsize=8)
        plt.savefig(table_path[:-4]+'_plot_mz_int.png', dpi=300)

    plt.savefig(table_list_MS2Planner[0][:-6]+'_ALL_plot_mz_int.png', dpi=300)
    plt.close()

def make_plot_MS2Planner_RT_int(table_list_MS2Planner, output_dir):
    colors = ['blue', 'violet', 'orange', 'red', 'yellow','green','purple','black','brown','pink','grey']
    labels = []

    for i, table_path in enumerate(table_list_MS2Planner):
        table = pd.read_csv(table_path, sep=',', header=0)

        plt.scatter('rt_apex', 'intensity', data=table, marker='o', color=colors[i % len(colors)],  edgecolors='k', linewidths=0.2, s=2, alpha=0.8)

        label = f"Inj. {i + 1}, n = {table.shape[0]}, median = {table['intensity'].median():.2e}, mean = {table['intensity'].mean():.2e}"
        labels.append(label)

        plt.yscale("log")
        plt.xlabel('Ret. time apex (s)')
        plt.ylabel('Intensity')
        plt.legend(labels=labels, fontsize=4)
        plt.title('Distribution of MS2Planner targets', fontsize=8)
        plt.savefig(table_path[:-4]+'_plot_RT_int.png', dpi=300)

    plt.savefig(table_list_MS2Planner[0][:-6]+'_ALL_plot_rt_int.png', dpi=300)
    plt.close()


def main():
    # call your IODA_MS2Planner_workflow function here with command line arguments
    IODA_MS2Planner_workflow(str(sys.argv[1]),int(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))

if __name__ == "__main__":
    main()
