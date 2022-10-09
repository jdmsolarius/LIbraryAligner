Buchser Lab / GEiC 'Library Aligner' 2017

11/5/2020

Three parts are needed to run LA (Library Aligner)

1: Reference sequences (should be provided as a tab delimited text file. The two columns are gRNA_name <tab> gRNA_sequence)
2: FastQ sequences (can be in multiple sub folders)
3: (optional) Settings file

The settings file is an XML file, use the default one provided in the folder to get you started.  The most important parts are:

PathReferenceList	<Path to the Reference File(s)>
PathFastQToSearch	<Path to the FastQ File(s)>
AttemptJoin	TRUE/FALSE
JoinRevType	0 (0 or 1 changes changes whether the R2 is reverse complemented or the R1 is before the Join, usually 0)

Run from command line >

dotnet run --project "..\Library Aligner\Library Aligner.csproj" "..\PathToSettings\LASettings.xml"

After the run, open the ref_f.txt file in Excel or another Spreadsheet to get the results for Read Counts.  Open the res_f.txt file to get the results per READ (warning, this file can be quite large, usually not of interest).