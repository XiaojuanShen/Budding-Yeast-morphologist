# Budding-Yeast-morphologist

Image analysis tool that partitions microscopy images into single-cell areas or mother-bud pairs areas with identified bud neck positions
=========================================================

This software contains various primitive c++ structures, templates and functions,
which are to facilitate modeling image feature measurents, as capture similarity
between populations of circular cells.

It includes a number of processes that are designed to allow the analysis of
collections of images by designing a pipeline of sequential tasks that manipulates
and analyzes each given image. As processes have no external dependencies, it can
easily be used on computer clusters or LSF architectures.

This code is available under the following licence:

 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.

# Availability
============

The original version is also available at 
http://www.moseslab.csb.utoronto.ca/louis-f/unsupervised/

# Installation
============

Download the source code, under "Processes and Code"

clone the repository

git clone https://github.com/jormungant/Budding-Yeast-morphologist.git

The installation only require a make command:
	./make -C ./PMCode/
(or)	cd ./PMCode ; ./make

The command produces the following processes:

PMExtractFeatures
PMFindMultiCover
PMGTRpvalues
PMHiddenMapDirect
PMMPDpvalues
PMMakeDisplay
PMProfiles
PMSegmentation
PMTiffManip
PMTiffOper
PMWaterShed
PMMakeConfidenceMatrices

# Pipeline Example
================

The provided code (see Installation) also includes an example pipeline in the
/example/ directory. In tha directory, there is a shell script "example.sh" which
will analyze a example image "example.tif".

After the shellscript is run, two output files will be produce:
	example_preview.tif	(Displayable representation of the segmented cells)
	example_data.txt	(Tab delimited text file of feature measurements)

Other intermediates file will be produced, some are tiff files that uses floating
point representation, so that some image software may not be able to recognize them.

# Help
====

Any processes has an embedded help, which is called using a -h flag, for example:

./PMMakeDisplay -h
