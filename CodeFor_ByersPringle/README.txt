This code retrieves data from the EZfate collection of precomputed
Lagrangian pathways (https://github.com/JamiePringle/EZfate) and
calculates the quantities plotted in the Byers and Pringle manuscript
for a number of different larval depths and vertical behaviors, as
described in the manuscript.

NOTE WELL: The file "connectivityUtilities.R" is part of the EZfate
project, version 1.02. If EZfate changes how or where data is stored,
this file may fail to retrieve or process the data correctly. In this
case, you should replace it with the most recent version of
"connectivityUtilities.R" from the EZfate github page.

The files are:

    README.txt: This file
    
    01_Get_and_Trim_data.R: The file that downloads EZfate Lagrangian
    pathways and computes statistics.  If you want to change the
    definition of habitat, the vertical behavior, or the depth of
    larvae, This is the file to alter. Note well, by default this makes data
    for the entire globe, and for many larval durations and seasons. This 
    takes a while to calculate. If you don't need all of this calculated, you 
    can alter it to only calculate what you need! This can save significant
    time.

    averagePlaceOnSphere.R: A library that computes the average
    location of a set of points given in latitude and longitude in a
    least-squares minimum distance sense, using a Haversine formula
    for distance on a spherical Earth. Also returns the standard
    deviation of the distance from all of the input points to the mean
    location.

    connectivityUtilities.R: The library to manipulate connectivity
    from the EZfate project, and to download the data for these
    calculations. Please go to EZfate for documentation of this code,
    and how to use it.

    connectivityData: a directory which will have EZfate connectivity
    data placed in it when the code in 01_Get_and_Trim_data.R or
    connectivtyUtilities is run.

    dataFiles: a directory where the results of the calculations of
    01_Get_and_Trim_data.R will be placed.

    plot_results_on_map.R: This is a little example of how to process
    and display the quantities calculated by 01_Get_and_Trim_data.R or
    derived from that data using equations in the paper. See below. 

The results of the calculations in 01_Get_and_Trim_data.R will be
stored in CSV files in dataFiles. These, when loaded into R, will
result in a data frame with the following fields:

       lonFrom: the longitude the Lagrangian particle started from
       
       latFrom: the latitude the Lagrangian particle started from

       Ladv: the distance in kilometers from the starting point to the
       mean ending point of the Lagrangian pathways

       Ldiff: the standard deviation of the distance in kilometers of
       the ending points of the Lagrangian pathways from the mean
       ending point

       fracReturn: the fraction of the particles launched from that
       point which returned to ANY habitat

       lonMeanTo: the longitude of the average ending point of the
       Lagrangian pathways

       latMeanTo: the latitude of the average ending point of the
       Lagrangian pathways

All other quantities in this paper can be straightforwardly calculated
from this data. To walk you through this process, plot_results_on_map.R
reads in the data that is calculated 01_Get_and_Trim_data.R, subsets it
to only display a portion of the Carribbean and Florida, calculates
the physical adversity (equation 3 in the paper), and then plots Ladv, Ldiff,
the fraction retained in the domain, and the physical adversity. The data is
plotted with plotly, so the resulting plots are interactive.


