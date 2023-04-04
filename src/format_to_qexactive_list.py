import pandas as pd
import numpy as np
import sys
import os
from io import StringIO
import warnings
from typing import List

## ====== EXCLUSION LISTS =========
# EXACTIVE EXCLUSION
def generate_Exactive_exclusion_table(input_table: str, 
                                      blank_samplename:str, 
                                      output_filename:str,
                                      rt_margin: float = 0,           # retention time margin to include before and after the start and end times (in seconds)
                                      polarity: str = 'Positive',     # ionization mode ('Positive' or 'Negative')
                                      ):
    
    """Format a table with mz, charge, rt, into a standard Exactive exclusion list
    
    Args:
        input_table (str): Path to the input table file. Mandatory columns: #Mass [m/z]',1: 'retention_time',2: 'rt_start',3: 'rt_end'})

        blank_samplename (str): Name of the blank sample column
        output_filename (str): Desired name and path for the output MQL table file
        rt_margin (float, optional): Retention time margin to include before and after the start and end times (in seconds)
        polarity (str, optional): Ionization mode ('Positive' or 'Negative'). Defaults to 'Positive'.
        
    Returns:
        None        
    """

    # Prepare the columns
    df_master = pd.read_csv(input_table, sep=',', header=0)
    df_master['Start [min]']=(df_master['rt_start']-rt_margin)/60
    df_master['End [min]']=(df_master['rt_end']+rt_margin)/60
    #Format to scientific notation the intensity
    df_master[blank_samplename] = df_master[blank_samplename].astype(float).map(lambda n: '{:.2E}'.format(n))

    #Build a comment (optional)
    df_master['block1'] = round(df_master['retention_time']*1/60,3)
    df_master['block2'] = round(df_master['retention_time'],2)

    df_master.loc[df_master['Start [min]'] < df_master['block1'].min(), 'Start [min]'] = df_master['block1'].min()  ### This to prevent negative value
    df_master.loc[df_master['End [min]'] > df_master['block1'].max(),'End [min]'] = df_master['block1'].max()  ### This to prevent value higher than the acquisition method

    df_master['block1'] = df_master['block1'].astype(str)
    df_master['block2'] = df_master['block2'].astype(str)
    df_master['for_comments'] = 'Apex = '+df_master['block1']+' (min) or '+df_master['block2']+' (sec) and int. = '+ df_master[blank_samplename].astype(str)

    #Make the output table
    df = pd.DataFrame(data=None)
    df['Mass [m/z]'] = df_master["Mass [m/z]"].round(decimals=4)
    df['Formula [M]'] = ''
    df['Formula type'] = ''
    df['Species'] = '' # '=-H'

    # Polarity
    if not isinstance(polarity, str):
        raise TypeError("Polarity should be a string")

    polarity = polarity.capitalize()
    if polarity not in ('Positive', 'Negative'):
        raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
    if polarity == 'Negative':
        df['CS [z]'] = [-1 if x == 0 or pd.isna(x) else -x for x in df_master['charge']]
        df['Polarity'] = 'Negative'

    else:
        df['CS [z]'] = [1 if x == 0 or pd.isna(x) else x for x in df_master['charge']]
        df['Polarity'] = 'Positive'

    df["Start [min]"] = df_master["Start [min]"].round(decimals=3)
    df["End [min]"] = df_master["End [min]"].round(decimals=3)
    df['(N)CE'] = '' #Can be empty ONLY INCLUSION
    df['(N)CE type'] = '' #Can be empty ONLY INCLUSION
    df['MSX ID'] = '' #Can be empty ONLY INCLUSION
    df['Comment'] = df_master['for_comments']
    df.to_csv(output_filename, index = False, sep=',')

# EXPLORIS EXCLUSION
def generate_Exploris_exclusion_table(
        input_table: str, 
        blank_samplename:str,
        output_filename:str,
        rt_margin: float = 0,            # retention time margin to include before and after the start and end times (in seconds)
        polarity: str = 'Positive',     # ionization mode ('Positive' or 'Negative')
        apex_int_percent:float = 0.1        # Percent of the apex intensity - used to define the intensity threshold
)-> pd.DataFrame :
    
    """Format a table with mz, charge, rt_start, rt_end, intensities into a Exploris exclusion list"""

    # Prepare the columns
    df_master = pd.read_csv(input_table, sep=',', header=0)
    df_master['Start [min]']=((df_master['rt_start'])-rt_margin)/60
    df_master['End [min]']=((df_master['rt_end'])+rt_margin)/60
    #Format to scientific notation the intensity
    df_master[blank_samplename] = df_master[blank_samplename].astype(float).map(lambda n: '{:.2E}'.format(n))

    #Build a comment (optional)
    df_master['block1'] = round(df_master['retention_time']*1/60,3)
    df_master['block2'] = round(df_master['retention_time'],2)

    df_master.loc[df_master['Start [min]'] < df_master['block1'].min(), 'Start [min]'] = df_master['block1'].min()  ### This to prevent negative value
    df_master.loc[df_master['End [min]'] > df_master['block1'].max(),'End [min]'] = df_master['block1'].max()  ### This to prevent value higher than the acquisition method

    df_master['block1'] = df_master['block1'].astype(str)
    df_master['block2'] = df_master['block2'].astype(str)
    df_master['for_comments'] = 'Apex = '+df_master['block1']+' min or '+df_master['block2']+' sec and int. = '+ df_master[blank_samplename].astype(float).map(lambda n: '{:.2E}'.format(n)).astype(str)

    #Make the output table
    df = pd.DataFrame(data=None)
    df['Compound'] = df_master['for_comments']
    df['Formula'] = ''
    df['Adduct'] = str("(no adduct)")
    df['m/z'] = df_master["Mass [m/z]"].round(decimals=4)

    polarity = polarity.capitalize()
    if polarity not in ('Positive', 'Negative'):
        raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
    if polarity == 'Negative':
        df['z'] = [-1 if (x != 0 and not np.isnan(x)) else x for x in df_master['charge']] #if 0 sequence is buging
    else:
        df['z'] = [1 if (x != 0 and not np.isnan(x)) else x for x in df_master['charge']] #if 0 sequence is buging

    df["t start (min)"] = df_master["Start [min]"].round(decimals=3)
    df["t stop (min)"] = df_master["End [min]"].round(decimals=3)
    df["Intensity Threshold"] = round(pd.to_numeric(df_master[blank_samplename], errors='coerce'), 0) * apex_int_percent


    df.to_csv(output_filename, index = False, sep=',')  


# MaxQuant.Live EXCLUSION  # not tested
def generate_MQL_exclusion_table(input_table:str, 
                           blank_samplename:str, 
                           output_filename:str, 
                           polarity:str = 'Positive',      # Ionization mode ('Positive' or 'Negative'). Defaults to 'Positive'.
                           apex_int_percent:float = 0.1    # Percent of the apex intensity - used to define the intensity threshold
) ->pd.DataFrame:
    """Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list
    
    Args:
        input_table (str): Path to the input table file. Mandatory columns: #Mass [m/z]',1: 'retention_time',2: 'duration',3: 'rt_start',4: 'rt_end',5: 'intensity'})

        blank_samplename (str): Name of the blank sample column
        output_filename (str): Desired name and path for the output MQL table file
        polarity (str, optional): Ionization mode ('Positive' or 'Negative'). Defaults to 'Positive'.
        
    Returns:
    - df(pd.DataFrame): DataFrame with formatted data at the output_filename\      
    """

    #Mass [m/z]',1: 'mz_isolation',2: 'duration',3: 'rt_start',4: 'rt_end',5: 'intensity'})
    df = pd.read_csv(input_table)

    df2 = df[['Mass [m/z]','retention_time','charge',blank_samplename]]
    df2 = df2.copy()
    df2.loc[:, 'Mass'] = df2['Mass [m/z]'].round(decimals=5)
    df2.loc[:, 'Retention time'] = (df2['retention_time'] / 60)
    df2 = df2.copy().rename(columns={blank_samplename:'Apex intensity'})
    df2['Apex intensity']=df2['Apex intensity'].round(decimals=0)*apex_int_percent
    df2['placeholder'] = np.arange(len(df2)) + 1 #Mandatory for import
    df2['Modified sequence'] = np.arange(len(df2)) + 1 #Mandatory for import. Arbitrary string.
    #Polarity
    polarity = polarity.capitalize()
    if polarity not in ('Positive', 'Negative'):
        raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
    if polarity == 'Negative':
        df2['Charge'] = [-1 if (x != 0 and not np.isnan(x)) else x for x in df['charge']] #if 0 sequence is buging
    else:
        df2['Charge'] = [1 if (x != 0 and not np.isnan(x)) else x for x in df['charge']] #if 0 sequence is buging

    df2['Retention time'] = df2['Retention time'].round(decimals=4)
    df2['MaxIt'] = ''
    df2["Colission Energies"] = ''
    
    # Find the 75 percentile value of 'Apex Intensity' column
    q3 = df2['Apex intensity'].quantile(.75)
    # Update the 'RealtimeCorrection' column with 'TRUE' value where the 'Apex Intensity' value is greater than or equal to q3.
    df2.loc[df2['Apex intensity'] >= q3, 'RealtimeCorrection'] = 'TRUE'

    df2['TargetedMs2'] = 'FALSE'
    df2['Targetedlabeled'] = 'FALSE'
    df2['TargetedMultiInjection'] = 'FALSE'
    df2['TopNExclusion'] = 'TRUE'
    df2['Fragments mz'] = ''
    df2['NCE Factors'] = ''

    df_out = df2[['placeholder','Modified sequence','Mass', 'Charge',\
                  'Retention time','Apex intensity','Fragments mz', 'MaxIt',\
                  'NCE Factors', 'Colission Energies','RealtimeCorrection','TargetedMs2',\
                  'Targetedlabeled','TargetedMultiInjection','TopNExclusion']]
    df_out = df_out.copy()
    df_out.rename(columns={'placeholder': ''}, inplace=True)

    df_out.to_csv(output_filename, index=None, sep='\t')

## ====== TARGETED LISTS - MS2Planner =========

# EXACTIVE TARGETED
def generate_Exactive_table_from_MS2Planner(
    input_table: str,               # path to the input table file
    output_filename: str,           # desired name and path for the output QExactive list
    polarity: str = 'Positive',     # ionization mode ('Positive' or 'Negative')
    rt_margin: float = 0            # retention time margin to include before and after the start and end times (in seconds)
) -> pd.DataFrame:

    """
    Formats a table with m/z, charge, rt_start, rt_end, intensities into a standard QExactive inclusion/exclusion list.

    Args:
        input_table (str): Path to the input table file.
        output_filename (str): Desired name and path for the output QExactive list.
        polarity (str, optional): Ionization mode. Allowed values are 'Positive' or 'Negative'. Defaults to 'Positive'.
        rt_margin (float, optional): Retention time margin to include before and after the start and end times (in seconds). Defaults to 0.
    
    Returns:
    - df(pd.DataFrame): DataFrame with formatted data at the output_filename
    """

    # Prepare the columns
    df_master = pd.read_csv(input_table, sep=',', header=0)   # read the input table into a Pandas dataframe
    df_master['rt_apex'] = pd.to_numeric(df_master['rt_apex'], errors='coerce')

    threshold = df_master['rt_end'].max() + rt_margin
    df_master.loc[df_master['rt_start'] > (rt_margin*2), 'Start [min]'] = (df_master['rt_start'] - rt_margin)/60
    df_master.loc[df_master['rt_end'] < threshold, 'End [min]'] = (df_master['rt_end'] + rt_margin)/60

    # Find the minimum value in the 'rt_start' column that is greater than 0
    # and assign it to the 'min_size' variable
    min_size = df_master.loc[df_master['rt_start'] > 0, 'rt_start'].min()

    # Apply a lambda function to the 'rt_start' column in the dataframe (df_master)
    # This function replaces any value in the 'rt_start' column smaller than 'min_size' with 'min_size'
    df_master['rt_start'] = df_master['rt_start'].apply(lambda x: max(x, min_size))

    # Format to scientific notation the intensity
    df_master['intensity'] = df_master['intensity'].astype(float).map(lambda n: '{:.2E}'.format(n))

    # Build a comment (optional)
    df_master['block1'] = round(df_master['rt_apex']*1/60,3)
    df_master['block2'] = round(df_master['rt_apex'],2)
    df_master['block1'] = df_master['block1'].astype(str)
    df_master['block2'] = df_master['block2'].astype(str)
    df_master['for_comments'] = 'Apex = '+df_master['block1']+' (min) or '+df_master['block2']+' (sec) and int. = '+ df_master['intensity'].astype(str)

    # Make the output table
    df = pd.DataFrame(data=None)
    df['Mass [m/z]'] = df_master["Mass [m/z]"].round(decimals=4)
    df['Formula [M]'] = ''
    df['Formula type'] = ''
    df['Species'] = ''

    # Polarity

    if not isinstance(polarity, str):
        raise TypeError("Polarity should be a string")

    polarity = polarity.capitalize()

    if polarity not in ('Positive', 'Negative'):
        raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
    if polarity == 'Negative':
        df['CS [z]'] = [-1 if x == 0 or pd.isna(x) else -x for x in df_master['charge']]
    else:
        df['CS [z]'] = [1 if x == 0 or pd.isna(x) else x for x in df_master['charge']]

    df['Polarity'] = polarity
    df["Start [min]"] = df_master["Start [min]"].round(decimals=3)
    df["End [min]"] = df_master["End [min]"].round(decimals=3)
    df['(N)CE'] = ''
    df['(N)CE type'] = ''
    df['MSX ID'] = ''
    df['Comment'] = df_master['for_comments']

    df.to_csv(output_filename, index=False, sep=',')   # write the output table to a CSV file


# MaxQuant.Live TARGETED
def generate_MQL_tMS2_table_from_MS2Planner(
        input_table:str,                # Path to the input table file
        output_filename:str,            # Desired name and path for the output MQL file
        rt_margin:float = 0,            # Retention time margin to include before and after the start and end times (in seconds)
        delay:float = 0,                # Delay time between inclusion and fragmentation events (in seconds)
        transient_time:float = 0,       # Transient time of the instrument (in seconds)
        polarity:str = 'Positive',      # Ionization mode ('Positive' or 'Negative')
        apex_int_percent:float = 0.6    # Percent of the apex intensity - used to define the intensity threshold
) -> pd.DataFrame:
    """
    Format a table with mz, charge, rt, intensities into a standard MaxQuantLive list
    
    Args:
        input_table (str): Path to the input table file
        output_filename (str): Desired name and path for the output MQL file
        rt_margin (float, optional): Retention time margin to include before and after the start and end times (in seconds). Defaults to 0.
        delay (float, optional): Delay time between inclusion and fragmentation events (in seconds). Defaults to 0.
        transient_time (float, optional): Transient time of the instrument (in seconds). Defaults to 0.
        polarity (str, optional): Ionization mode ('Positive' or 'Negative'). Defaults to 'Positive'.
        apex_int_percent (float, optional): Percent of the apex intensity - used to define the intensity threshold. Defaults to 0.6.  
    """

    df = pd.read_csv(input_table)
    df2 = df[['Mass [m/z]','duration','rt_start','rt_end','intensity','rt_apex','charge']]
    df2 = df2.copy()  # create a new copy of the DataFrame
    df2.rename(columns={'Mass [m/z]':'Mass'}, inplace=True)  # modify the copy using the .rename() method
    # This is to set the RT acquisition window
    df2['Retention time']= round((df2['rt_apex']-rt_margin)/60, 3)  #MQL requires in minutes
    df2['Apex intensity']= (df2['intensity']*apex_int_percent).round(decimals=0)
    df2['placeholder'] = np.arange(len(df2)) + 1 #Mandatory for import
    df2['Modified sequence'] = np.arange(len(df2)) + 1 #Mandatory for import. Arbitrary number.

    polarity = polarity.capitalize()
    if polarity not in ('Positive', 'Negative'):
        raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
    if polarity == 'Negative':
        df2['Charge'] = [-1 if (x != 0 and not np.isnan(x)) else x for x in df2['charge']] #if 0 sequence is buging
    else:
        df2['Charge'] = [1 if (x != 0 and not np.isnan(x)) else x for x in df2['charge']] #if 0 sequence is buging

    df2['Retention time'] = df2['Retention time'].round(decimals=4)
    #The MaxIT is derived from the 'duration' column. We take only 95% of it. We substract the transient time (MS2 and overhead time).
    if transient_time == 0:
        df2['MaxIt'] = (df2['duration']*1000)-transient_time-delay #MQL MaxIT in ms
    else:
        df2['MaxIt'] = ''

    df2["Colission Energies"] = ''

    # Find the 75 percentile value of 'Apex Intensity' column
    q3 = df2['Apex intensity'].quantile(.75)
    # Update the 'RealtimeCorrection' column with 'TRUE' value where the 'Apex intensity' value is greater than or equal to q3.
    df2.loc[df2['Apex intensity'] >= q3, 'RealtimeCorrection'] = 'TRUE'
    
    df2['TargetedMs2'] = 'TRUE'
    df2['TargetedLabeled'] = 'FALSE'
    df2['TargetedMultiInjection'] = 'FALSE'
    df2['TopNExclusion'] = 'TRUE'
    df2['Fragments mz'] = ''
    df2['NCE Factors'] = ''

    df_out = df2[['placeholder','Modified sequence','Mass', 'Charge',\
                  'Retention time','Apex intensity','Fragments mz', 'MaxIt',\
                  'NCE Factors', 'Colission Energies','RealtimeCorrection','TargetedMs2',\
                  'TargetedLabeled','TargetedMultiInjection','TopNExclusion']]
    df_out = df_out.copy()  # create a new copy of the DataFrame
    df_out.rename(columns={'placeholder':''}, inplace=True)
    
    df_out.to_csv(output_filename, index=None, sep='\t')


# Exploris TARGETED DDA-MS2
def generate_Exploris_DDAMS2_table_from_MS2Planner(
    input_table: str,                   # File path to the input table containing mz, charge, rt, and intensities
    output_filename:str,                # File name for the output MS2Planner target list in Exploris format
    pretarget_rt_margin:float = 0,      # Time (in secs) to subtract to the rt_start
    posttarget_rt_margin:float = 0,     # Time (in secs) to add to the rt_end
    transient_time:float = '',          # Transient time (in msecs)
    polarity:str = 'Positive',          # Polarity either 'Positive' or 'Negative'
    CEs:str or List[str] = '',          # Colission energy(ies), stepping is possible like '25,35,45' or ['25,35,45','55,75']
    apex_int_percent:float = 0.6        # Percent of the apex intensity - used to define the intensity threshold
    ) -> pd.DataFrame:
    
        """
    Format a table with mz, charge, rt, and intensities into an Exploris series inclusion/exclusion list. 
    
    Args:
    - input_table (str): inputfile path to the input table containing mz, charge, rt, and intensities.
    - output_filename (str): outputfile path for the output MS2Planner target list in Exploris format.
    - pretarget_rt_margin (float): time (in secs) to subtract to the rt_start. Default = 0
    - posttarget_rt_margin (float): time (in secs) to add to the rt_end. Default = 0
    - CEs (list[str]/str): collision energy(ies), stepping is possible like '25,35,45' or ['25,35,45'],['55;75']. Can be set to '' for default value.
    - Percent (float): of the apex intensity - used to define the intensity threshold
   
   Returns:
    - df(pd.DataFrame): DataFrame with formatted data at the output_filename
    """
        
    # Prepare the columns
        df_master = pd.read_csv(input_table)

    # Calculate the minimum value of rt_start column by adding pretarget_rt_margin to the maximum value of rt_start column 
        min_rt_value = df_master['rt_start'].max() + pretarget_rt_margin
        # Subtract pretarget_rt_margin from rt_start column for all rows that have values greater than min_rt_value
        df_master.loc[df_master['rt_start'] > min_rt_value, 'rt_start'] = (df_master['rt_start'] - pretarget_rt_margin)
        
        # Calculate the maximum value of rt_end column by subtracting posttarget_rt_margin from the maximum value of rt_end column 
        max_rt_value = df_master['rt_end'].max() - posttarget_rt_margin
        # Add posttarget_rt_margin to rt_end column for all rows that have values less than max_rt_value
        df_master.loc[df_master['rt_end'] < max_rt_value, 'rt_end'] = (df_master['rt_end'] + posttarget_rt_margin)
        
        df_master['for_comments'] = 'Apex = ' + df_master['rt_apex'].round(3).astype(str) + ' (sec) ' + (df_master['rt_apex']/60).round(2).astype(str) + ' (min) and int. = '+ df_master['intensity'].astype(float).map(lambda n: '{:.2E}'.format(n)).astype(str)

        #Make the output table
        df = pd.DataFrame(data=None)
        df['Compound'] = df_master['for_comments']
        df['Formula'] = ''
        df['Adduct'] = str("(no adduct)")
        df['m/z'] = df_master["Mass [m/z]"].round(decimals=4)
        
        #Polarity
        if not isinstance(polarity, str):
            raise TypeError("Polarity should be a string")
            
        polarity = polarity.capitalize()
        if polarity not in ('Positive', 'Negative'):
            raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
        
        if polarity == 'Negative':
            df['z'] = [-1 if (x != 0 and not np.isnan(x)) else x for x in df_master['charge']]
        else:
            df['z'] = [1 if (x != 0 and not np.isnan(x)) else x for x in df_master['charge']]

        #RT
        df["rt start (min)"] = round(df_master["rt_start"]/60, 3)
        df["rt stop (min)"] = round(df_master["rt_end"]/60, 3)
        df["Intensity Threshold"] = round(df_master["intensity"]*apex_int_percent,0)

        # If multiple collision energies are provided
        if CEs != '':
            if type(CEs) == list:
                for x in CEs:
                    df_new = df.copy()
                    df_new["HCD Collision Energies (%)"] = str(x)
                    df = pd.concat([df, df_new], axis=0)
                    df["Maximum Injection Time (ms)"] = df["Maximum Injection Time (ms)"]*len(CEs)*1.25   # reducing the max IT to account for multiple MS2 scan 

            else:
                df["HCD Collision Energies (%)"] = str(CEs)
                df["Maximum Injection Time (ms)"] = round(((df_master["duration"])/1000)-transient_time,0)

        df.to_csv(output_filename, index = False, sep=',')

# EXPLORIS TARGETED tMS2
def generate_Exploris_tMS2_table_from_MS2Planner(
    input_table: str,                   # File path to the input table containing mz, charge, rt, and intensities
    output_filename:str,                # File name for the output MS2Planner target list in Exploris format
    pretarget_rt_margin:float = 0,      # Time (in secs) to subtract to the rt_start
    posttarget_rt_margin:float = 0,     # Time (in secs) to add to the rt_end
    transient_time:float = 15,          # Transient time (in msecs)
    polarity:str = 'Positive',          # Polarity either 'Positive' or 'Negative'
    RF_base_value:float = np.nan,           # Base RF lens value (%)
    CEs:str or List[str] = np.nan         # Normalized Colission energy(ies), stepping is possible like '25,35,45' or ['25,35,45','55,75']
) -> pd.DataFrame:
    
        """
    Format a table with mz, charge, rt, and intensities into an Exploris series inclusion/exclusion list. 
    
    Args:
    - input_table (str): inputfile path to the input table containing mz, charge, rt, and intensities.
    - output_filename (str): outputfile path for the output MS2Planner target list in Exploris format.
    - pretarget_rt_margin (float): time (in secs) to subtract to the rt_start. Default = 0
    - posttarget_rt_margin (float): time (in secs) to add to the rt_end. Default = 0
    - polarity (str): polarity of the ions. Either 'Positive' or 'Negative'. 
    - RF_base_value (float): base RF lens value (%). Value must be >= 0. Default is 70%
    - CEs (list[str]/str): collision energy(ies), stepping is possible like '25,35,45' or ['25,35,45'],['55;75']. Can be set to '' for default value.
    
    Returns:
    - df(pd.DataFrame): DataFrame with formatted data at the output_filename
    """
        
    # Prepare the columns
        df_master = pd.read_csv(input_table)

        # Check if polarity is valid
        if not isinstance(polarity, str):
            raise TypeError("Polarity should be a string")
            
        polarity = polarity.capitalize()
        if polarity not in ('Positive', 'Negative'):
            raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
        
        
        # Calculate the minimum value of rt_start column by adding pretarget_rt_margin to the maximum value of rt_start column 
        min_rt_value = df_master['rt_start'].max() + pretarget_rt_margin
        # Subtract pretarget_rt_margin from rt_start column for all rows that have values greater than min_rt_value
        df_master.loc[df_master['rt_start'] > min_rt_value, 'rt_start'] = (df_master['rt_start'] - pretarget_rt_margin)
        
        # Calculate the maximum value of rt_end column by subtracting posttarget_rt_margin from the maximum value of rt_end column 
        max_rt_value = df_master['rt_end'].max() - posttarget_rt_margin
        # Add posttarget_rt_margin to rt_end column for all rows that have values less than max_rt_value
        df_master.loc[df_master['rt_end'] < max_rt_value, 'rt_end'] = (df_master['rt_end'] + posttarget_rt_margin)
        
        df_master['for_comments'] = 'Apex = ' + df_master['rt_apex'].round(3).astype(str) + ' (sec) ' + (df_master['rt_apex']/60).round(2).astype(str)+' (min) and int. = '+ df_master['intensity'].astype(float).map(lambda n: '{:.2E}'.format(n)).astype(str)

        #Make the output table
        df = pd.DataFrame(data=None)
        df['Compound'] = df_master['for_comments']
        df['Formula'] = ''
        df['Adduct'] = str("(no adduct)")
        df['Precursor (m/z)'] = df_master["Mass [m/z]"].round(decimals=4)

        #Polarity
        if not isinstance(polarity, str):
            raise TypeError("Polarity should be a string")
            
        polarity = polarity.capitalize()
        if polarity not in ('Positive', 'Negative'):
            raise ValueError("Invalid polarity. Allowed values are 'Positive' or 'Negative'")
                
        if polarity == 'Negative':
            df['Precursor Charge (z)'] = [-1 if x != x or x == 0 else -x for x in df_master['charge']] 
        else:
            df['Precursor Charge (z)'] = [1 if x != x or x == 0 else x for x in df_master['charge']] 

        df["rt start (min)"] = round(df_master["rt_start"]/60, 3)
        df["rt stop (min)"] = round(df_master["rt_end"]/60, 3)
        df["Isolation Window (m/z)"] = round(df_master["mz_isolation"],0)
        df["Orbitrap Resolution"] = '15000'
        
        # Normalize RF Lens based on Precursor column
        if not np.isnan(RF_base_value):
            precursor_min = df['Precursor (m/z)'].min()
            precursor_max = df['Precursor (m/z)'].max()
            df['RF Lens (%)'] = (df['Precursor (m/z)'] - precursor_min)/(precursor_max - precursor_min) * RF_base_value
            df["RF Lens (%)"] = np.clip(df["RF Lens (%)"], 40, 100)             # Clip the values for min/max

        # Normalize Normalized AGC Target based on intensity column
        df['Normalized AGC Target (%)'] = 50 * (df_master['intensity'] - df_master['intensity'].min()) / (df_master['intensity'].max() - df_master['intensity'].min()) + 75

        df["Maximum Injection Time (ms)"] = round(((df_master["duration"])/1000)-transient_time,0)

        df["Microscans"] = '1'
        df["Data Type"] = 'Profile'
        df["Polarity"] = str(polarity)
        
        # If multiple collision energies are provided
        if not np.isnan(CEs):
            if type(CEs) == list:
                for x in CEs:
                    df_new = df.copy()
                    df_new["HCD Collision Energies (%)"] = str(x)
                    df = pd.concat([df, df_new], axis=0)
                    df["Maximum Injection Time (ms)"] = df["Maximum Injection Time (ms)"]*len(CEs)*1.25   # reducing the max IT to account for multiple MS2 scan 

            else:
                df["HCD Collision Energies (%)"] = str(CEs)

        df.to_csv(output_filename, index = False, sep=',')
