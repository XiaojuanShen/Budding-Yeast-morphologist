.PHONY: clean

COMPILER_FLAGS = -Wno-format -Wno-write-strings
all: ./bin/PMMakeConfidenceMatrices ./bin/PMHiddenMapDirect ./bin/PMSegmentation ./bin/PMTiffManip ./bin/PMTiffOper ./bin/PMWaterShed ./bin/PMExtractFeatures ./bin/PMProfiles ./bin/PMGTRpvalues ./bin/PMMPDpvalues ./bin/PMFindMultiCover ./bin/PMMakeDisplay ./bin/PMFlexToTiff

./src/primitive.o: ./src/primitive.cpp ./src/primitive.h ./src/primitive_tem.h ./src/primitive_stats.h ./src/DataGrid.hpp ./src/primitive_stats.h ./src/Vector.hpp ./src/primstats_tem.h
	g++ ${COMPILER_FLAGS} -o ./src/primitive.o -c ./src/primitive.cpp

./bin/PMProfiles: ./src/Modeling.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_PROFILES ./src/main_caller.cpp ./src/Modeling.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMProfiles

./src/Modeling.o: ./src/Modeling.cpp
	g++ ${COMPILER_FLAGS} -c ./src/Modeling.cpp -o ./src/Modeling.o

./bin/PMTiffManip : ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_TIFF_FILE_MANIPULATION ./src/main_caller.cpp ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMTiffManip

./bin/PMTiffOper: ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_TIFF_FILE_OPERATION ./src/main_caller.cpp ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMTiffOper

./bin/PMFlexToTiff: ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_REMOVE_FACTORS ./src/main_caller.cpp ./src/TiffManip.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMFlexToTiff

./src/TiffManip.o: ./src/TiffManip.cpp
	g++ ${COMPILER_FLAGS} -c ./src/TiffManip.cpp -o ./src/TiffManip.o

./bin/PMHiddenMapDirect: ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_EXTRACT_HIDDENMAP_DIRECT ./src/main_caller.cpp ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMHiddenMapDirect

./bin/PMSegmentation: ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_HMM_SEGMENTATION_AND_DISTANCE ./src/main_caller.cpp ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMSegmentation

./bin/PMFindMultiCover: ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_FIND_CELL_COVER_UPGRADE ./src/main_caller.cpp ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMFindMultiCover

./bin/PMWaterShed : ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_GEN_META_PIXELS ./src/main_caller.cpp ./src/Segmentation.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMWaterShed

./src/Segmentation.o : ./src/Segmentation.cpp
	g++ ${COMPILER_FLAGS} -c ./src/Segmentation.cpp -o ./src/Segmentation.o

./bin/PMExtractFeatures: ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_HIDDENMAP_BASE_DATAEXTRACTION ./src/main_caller.cpp ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMExtractFeatures

./bin/PMMakeDisplay : ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_MAKE_DISPLAYABLE ./src/main_caller.cpp ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMMakeDisplay

./bin/PMMakeConfidenceMatrices: ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_GETLABEL_CONFIDENCE ./src/main_caller.cpp ./src/Extraction.o ./src/Madstructs.o ./src/primitive.o -o ./bin/PMMakeConfidenceMatrices

./src/Extraction.o : ./src/Extraction.cpp
	g++ ${COMPILER_FLAGS} -c ./src/Extraction.cpp -o ./src/Extraction.o

./bin/PMMPDpvalues : ./src/Pvalue.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_MULTI_PAIR_DISTANCE_PVAL ./src/main_caller.cpp ./src/Madstructs.o ./src/primitive.o ./src/Pvalue.o -o ./bin/PMMPDpvalues

./bin/PMGTRpvalues : ./src/Pvalue.o ./src/Madstructs.o ./src/primitive.o
	g++ ${COMPILER_FLAGS} -DTaSk_Id=TASK_CDTGTR_PVALUES ./src/main_caller.cpp ./src/Madstructs.o ./src/primitive.o ./src/Pvalue.o -o ./bin/PMGTRpvalues

./src/Pvalue.o : ./src/Pvalue.cpp
	g++ ${COMPILER_FLAGS} -c ./src/Pvalue.cpp -o ./src/Pvalue.o

./src/Madstructs.o: ./src/Madstructs.cpp ./src/Madstructs.h
	g++ ${COMPILER_FLAGS} -c ./src/Madstructs.cpp -o ./src/Madstructs.o

clean:
	rm ./src/*\.o
