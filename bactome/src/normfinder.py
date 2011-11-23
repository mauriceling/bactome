"""
Python implementation of NormFinder.

NormFinder reference: Andersen CL, Jensen JL, Orntoft TF. 2004.
Normalization of real-time quantitative reverse transcription-PCR 
data: a model-based variance estimation approach to identify genes 
suited for normalization, applied to bladder and colon cancer data 
sets. Cancer Research 64(15):5245-5250.

Date created: 23rd November 2011

Licence: Python Software Foundation License version 2
"""

def normfinder(data, header=True):
    """
    NormFinder implementation.
    
    @param data: List of values in the format of 
        [[<row identifier>, value, value, ...]]
    @param header: Flag if first row of data is sample header
        (default = True)
    @return: Dictionary of result - 
        {<data row number>: [<row identifier>, <NormFinder Index>]
    """
    return result
    
