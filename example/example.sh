echo Parse Tiff Files for "red" channels:
../bin/PMTiffManip -f 2,4,6,8 example.tif example_red.tif 

echo Separate Foreground and Background in Images:
../bin/PMSegmentation -B 1.0f -b example_seg.scp -d example_dist.tif example_red.tif example_seg.tif

echo Identify Cell centers:
../bin/PMFindMultiCover -o example_ellfit.mcv -R2 1.0 6.0 example_red.tif example_seg.tif example_dist.tif

echo Find Cell sub-segments:
../bin/PMWaterShed -M -b 1.0 example_red.tif example_watershed.tif

echo Cell idenfication from circle coordinated directed sub-segment agglomeration: 
../bin/PMHiddenMapDirect -M -G example_watershed.tif -C example_cellseg.tif example_seg.tif example_dist.tif example_ellfit.mcv

echo Parse Tiff Files for "green" channel:
../bin/PMTiffManip -f 1,3,5,7 example.tif example_gre.tif 

echo Compile Quality measuers into a condifence matrix binary file
../bin/PMMakeConfidenceMatrices Conf_Matrices.scp Quality_measures.txt

echo Measure GFP Intensity and spatial spread:
../bin/PMExtractFeatures -c Conf_Matrices.scp -m example_ellfit.mcv -a 0.1 10 0.9 0.1 50 -t example_data.txt -S example_seg.tif example_red.tif example_gre.tif example_cellseg.tif

echo Produce a Displayable RGB Tiff file showing the cell boundary and bud necks:
../bin/PMMakeDisplay example_data.txt example_cellseg.tif example_red.tif example_gre.tif example_preview.tif

