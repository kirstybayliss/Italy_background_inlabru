This README file contains information about the downloadable files of the DISS3.2.1, this file ends with an <EOF> tag.
Last updated on 05/04/2018

--- Citation ---
DISS Working Group (2018). Database of Individual Seismogenic Sources (DISS), Version 3.2.1: A compilation of potential sources for earthquakes larger than M 5.5 in Italy and surrounding areas. http://diss.rm.ingv.it/diss/, © INGV 2018 - Istituto Nazionale di Geofisica e Vulcanologia - DOI:10.6092/INGV.IT-DISS3.2.1.

Except where otherwise noted, content of this file and the files here described are licensed under a Creative Commons Attribution 4.0 International (CC BY-SA 4.0) licence.

DISS 3.2.1 is made available for desktop PCs in the following file formats:
MapInfo (TAB), ESRI Shapefile (SHP), Google Earth (KMZ).
For each format, there are several suites of files with the following extensions:
1) MapInfo (TAB): .tab, .map, .id, .dat, .ind, and one .wor file.
2) ESRI Shapefile (SHP): .shp, .shx, .prj, .dbf, .lyr, and one .mxd file.
3) Google Earth (KMZ): .kmz, .jpg.

For formats 1) and 2) there are the following files:
ISS321 - Individual Seismogenic Sources;
CSSTOP321 - Top traces of Composite Seismogenic Sources;
CSSPLN321 - Dipping planes of Composite Seismogenic Sources;
DSS321 - Debated Seismogenic Sources; 
SUBDCNT321 - Subduction depth contours;
SUBDZON321 - Subduction zone;
DISS321 - Container file with layering information.

All geographic information is provided in latitude/longitude, datum WGS84.

For format 3) the kmz file is a container itself with embedded layering information.

For additional information on the database structure not indicated here please refer to
Basili, R., Kastelic, V., Valensise, G., and DISS Working Group 2009 (2009), DISS3 tutorial series: Guidelines for compiling records of the Database of Individual Seismogenic Sources, version 3, Rapporti Tecnici INGV, no. 108, 20 p.,http://portale.ingv.it/produzione-scientifica/rapporti-tecnici-ingv/archivio/rapporti-tecnici-2009/

or otherwise send an email message to sorgenti@ingv.it.

Parametric attributes for each category of seismogenic sources are as follows:
--- ISS321 ---
IDSource	The DISS identifier assigned to the record (text)
SourceName	Name of the seismogenic source (text)
Length		Fault length measured along strike (km)
Width		Fault width measured along dip (km)
MinDepth	Depth of the fault upper edge from sea level (km)
MaxDepth	Depth of the fault lower edge from sea level (km)
Strike		Fault strike right-hand rule (deg, 0-360 CW from North)
Dip			Fault dip (deg, 0-90 from horizontal) 
Rake		Fault rake (deg, 0-360 CCW from horizontal)
AvgDispl	Average displacement (m)
SRMin		Minimum slip rate value (mm/y)
SRMax		Maximum slip rate value (mm/y)
RecIntMin	Minimum recurrence interval (y)
RecIntMax	Maximum recurrence interval (y)
LatestEq	Date or age of the most recent earthquake associated with the source, if known (text)
MaxMag		Maximum moment magnitude (Mw)
LatUL		Latitude of upper-left corner
LonUL		Longitude of upper-left corner
LatUR		Latitude of upper-right corner
LonUR		Longitude of upper-right corner
LatLR		Latitude of lower-right corner
LonLR		Longitude of lower-right corner
LatLL		Latitude of lower-left corner
LonLL		Longitude of lower-left corner
Created		Date of creation of the record (text)
Updated		Date of latest update of the record (text)
LengthQ		Qualifier of the length   
WidthQ		Qualifier of the width    
MinDepthQ	Qualifier of the minimum depth 
MaxDepthQ	Qualifier of the maximum depth 
StrikeQ		Qualifier of the strike angle   
DipQ		Qualifier of the dip angle      
RakeQ		Qualifier of the rake angle     
AvgDisplQ	Qualifier of the Average displacement 
SlipRateQ	Qualifier of the slip rate values
RecIntQ		Qualifier of the recurrence interval values  
MaxMagQ		Qualifier of the maximum magnitude  
LocationQ	Qualifier of the fault location
LinkToInfo	WWW link to full information of the record

--- CSSTOP321 - CSSPLN320 ---
IDSource	The DISS identifier assigned to the record (text)
SourceName	Name of the seismogenic source (text)
MinDepth	Depth of the fault upper edge from sea level (km)
MaxDepth	Depth of the fault lower edge from sea level (km)
StrikeMin	Minimum fault strike right-hand rule (deg, 0-360 CW from North)
StrikeMax	Maximum fault strike right-hand rule (deg, 0-360 CW from North)
DipMin		Minimum fault dip (deg, 0-90 from horizontal)
DipMax		Maximum fault dip (deg, 0-90 from horizontal)
RakeMin		Minimum fault rake (deg, 0-360 CCW from horizontal)
RakeMax		Maximum fault rake (deg, 0-360 CCW from horizontal)
SRMin		Minimum slip rate value (mm/y)
SRMax		Maximum slip rate value (mm/y)
MaxMag		Maximum moment magnitude (Mw)
Created		Date of creation of the record (text)
Updated		Date of latest update of the record (text)
MinDepthQ	Qualifier of the minimum depth 
MaxDepthQ	Qualifier of the maximum depth 
StrikeQ		Qualifier of the strike angle   
DipQ		Qualifier of the dip angle      
RakeQ		Qualifier of the rake angle     
SlipRateQ	Qualifier of the slip rate values
MaxMagQ		Qualifier of the maximum magnitude  
LinkToInfo	WWW link to full information of the record

--- DSS321 ---
IDSource	The DISS identifier assigned to the record (text)
SourceName	Name of the seismogenic source (text)
Created		Date of creation of the record (text)
Updated		Date of latest update of the record (text)
LinkToInfo	WWW link to full information of the record

--- SUBDCNT321 ---
IDSource	The DISS identifier assigned to the record (text)
Depth_km	Depth of the contour line (km)

--- SUBDZON321 ---
IDSource	The DISS identifier assigned to the record (text)
SourceName	Name of the seismogenic source (text)
MinDepth	Depth of the seismogenic interface upper limit from sea level (km)
MaxDepth	Depth of the seismogenic interface lower limit from sea level (km)
DipDir		Average dip direction of the slab (N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW)
ConvAzMin	Minimum direction of plate convergence (deg, 0-360 CW from North)
ConvAzMax	Maximum direction of plate convergence (deg, 0-360 CW from North)
ConvRMin	Minimum convergence rate value (mm/y)
ConvRMax	Maximum convergence rate value (mm/y)
MaxMag		Maximum moment magnitude (Mw)
Created		Date of creation of the record (text)
Updated		Date of latest update of the record (text)
MinDepthQ	Qualifier of the minimum seismogenic interface depth 
MaxDepthQ	Qualifier of the maximum seismogenic interface depth 
DipDirQ		Qualifier of the dip direction  
ConvAzQ		Qualifier of the direction values of plate convergence  
ConvRQ		Qualifier of the convergence rate values  
MaxMagQ		Qualifier of the maximum magnitude  
LinkToInfo	WWW link to full information of the record

--- Common values for "qualifier" fields in all tables ---
1 = LD for Literature Data;
2 = OD for Original Data;
3 = ER for Empirical Relationship;
4 = AR for Analytical Relationship;
5 = EJ for Expert Judgment.

<EOF>