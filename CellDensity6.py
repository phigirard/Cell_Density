#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	CellDensity.py
#	Release v6.0
#
#	Copyright 2025 - BSD-3-Clause license
#                                                                             
#******************************************************************************/


#@ File nucleiFile (label="Select the nuclei image:", style="file") 
#@ String msg1 (visibility=MESSAGE, value="Choosing the best bin size for a cell density map is crucial because it directly affects the resolution and interpretability of your results.", required=False) 
#@ String msg2 (visibility=MESSAGE, value="The best value is based on cell size : bin size = 2 x cell size", required=False) 
#@ String msg3 (visibility=MESSAGE, value=" Automatic is calculating the cell size average with bin size = 2 x cell size average", required=False) 
#@ String AutomaticBin (label="Bin Size Selection", choices={"Automatic", "Manual (value below)"}, style="radioButtonHorizontal", value = choices[0],persist=True) 
#@ Integer BinSize (label="Manual bin Size in pixels:", value=15, persist=True) 

#@ DatasetIOService io
#@ UIService uiService
#@ LogService log
#@ CommandService command

#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

#import class ImageJ1
from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import WindowManager
from ij import Prefs
from ij.io import Opener
from ij.gui import  GenericDialog
from ij.gui import  Roi
from ij.gui import  ShapeRoi
from ij.gui import  YesNoCancelDialog
from ij.gui import  WaitForUserDialog
from ij.plugin import RoiEnlarger
from ij.plugin.frame import RoiManager
from ij.plugin.frame import ThresholdAdjuster
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.measure import Calibration
from ij.plugin.filter import Analyzer
from ij.plugin.filter import MaximumFinder
from ij.plugin.filter import ThresholdToSelection
from ij.process import ImageProcessor
from ij.process import ImageConverter
from ij.process import ImageStatistics


from mpicbg.ij.clahe import Flat
from net.imglib2.img.display.imagej import ImageJFunctions as IJF #wrapper to convert Dataset from Stardist to ImagePlus
from de.csbdresden.stardist import StarDist2D #import StarDist plugin
from de.csbdresden.stardist import Opt #import StarDist plugin
from net.haesleinhuepf.clij2 import CLIJ2 #import CLIJ2 plugin




import os
import sys
import csv
import math 

from java.io import File
from java.lang import Double
from java.awt import Color
from java.awt import Polygon
from java.awt.event import AdjustmentListener  

from datetime import datetime as dt



#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------
#---------------          CONSTANTS            -----------------
#---------------------------------------------------------------

Prefs.blackBackground = True

#Constants for CSV Files
IJ.run("Conversions...", "scale")
IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column save_row")

outputType = ["ROI Manager","Label Image","Both","Polygons"] #output for Stardist plugin

MAXSIZE = Double.POSITIVE_INFINITY # for ParticleAnalyzer


weigth  = 0.9 # bin size weightweightweightweight based on mean area calculation


# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')



#---------------------------------------------------------------
#----------------- All Classes for analysis  -----------------
#---------------------------------------------------------------

## threshold class
class ThresholdPreviewer(AdjustmentListener):  
	def __init__(self, imp, sliders):  
		self.imp = imp   
		self.sliders = sliders  
		self.threshold()
		
	def adjustmentValueChanged(self, event):  
		self.threshold()   
		
	def threshold(self):
		IJ.setThreshold(self.imp, 0 ,self.sliders.get(0).getValue())
		
	def getMaxThreshold(self) :
		return self.sliders.get(0).getValue()

#---------------------------------------------------------------
#----------------- All Functions for analysis  -----------------
#---------------------------------------------------------------


# Adjust the string size with 0
def adjustSizeNum(S_ , length_):
	newS_ = str(S_)
	while len(newS_)<length_ :
		newS_ = "0"+newS_
	return newS_

def changeValues(imp_,valIN_, ValOUT_):
	ip_ = imp_.getProcessor()
	pixels_ = ip_.getPixels()
	for  i in range (len (pixels_)):
		if (pixels_[i] == valIN_):
			pixels_[i] = ValOUT_
	imp_.updateAndDraw()
	return

def changeValue2NAN(imp_,val_):
	changeValues(imp_,val_,  Double.NaN)
	return

	
# Generic Dialog for threshold method
def thresholdImageUI(ip_, min_ , max_): 
	imp_=ImagePlus("AreaNuclei", ip_)	 
	imp_.show()	
	gd = GenericDialog("Nuclei Threshold")
	gd.addSlider("Nuclei area ", min_, max_, min_)
	sliders = gd.getSliders()
	previewer = ThresholdPreviewer(imp_, sliders)  
	for slider in sliders :
		slider.addAdjustmentListener(previewer)  

	gd.showDialog()  
	
	maxThreshold = min_
	if gd.wasOKed():  
		maxThreshold = previewer.getMaxThreshold()
	imp_.hide()
	imp_.close()
	return maxThreshold		

#---------------------------------------------------------------
#     -----------------       Start       -----------------
#---------------------------------------------------------------  


# clear the console automatically when not in headless mode
uiService.getDefaultUI().getConsolePane().clear()

#close Result Table if opened
if IJ.isResultsWindow() :
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	tw.close()
rt= ResultsTable()

#reset RoiManager or open one
rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()
rm.reset()


#close threshold window 
windowThreshold = WindowManager.getWindow("Threshold")
if (windowThreshold != None) :
	windowThreshold.setVisible(False)
			

# initialize ClearCL context and convenience layer
clij2 = CLIJ2.getInstance();

			
#convert Files from #@ parameters to String and extract the main directory of the data
nucleiPath = nucleiFile.getCanonicalPath()


#create folder for analysis 
srcDir = os.path.dirname(nucleiPath)
basename = os.path.basename(os.path.splitext(nucleiPath)[0]).lower()
imageDir = os.path.join(srcDir, basename) 
if not os.path.exists(imageDir):
	os.makedirs(imageDir)
	
	

#open the files and remove LUT
impNuclei = Opener().openImage(nucleiPath)
IJ.run(impNuclei, "Grays", "")
IJ.run(impNuclei, "Smooth", "stack")


# Get calibration and set time & pixel units 
cal = impNuclei.getCalibration()
pix2phys = cal.getX(1)	
physUnit = cal.getUnit()
if physUnit == "micron" : # convert to mm
	pix2phys = pix2phys/1000
# Remove physical calibration for processing steps
impNuclei.setCalibration(Calibration())


# Get width, height and frame 
nbSlice = impNuclei.getStackSize()
width = impNuclei.width
height = impNuclei.height
depth = impNuclei.getBitDepth()


# Create Nucleai and Heatmap cell density stacks
stackLabels = ImageStack(width, height)
stackCellDensity = stackLabels.duplicate()
stackNeighbors = stackLabels.duplicate()


impNAN = IJ.createImage("impNAN", "32-bit black" , width, height, 1)
changeValue2NAN(impNAN,0)
ipNAN = impNAN.getProcessor()

areaThres =1
thres_min = 0
thres_max = int(math.pow(2, depth))

for slic in range(nbSlice) : 
	rm.reset()
	rt.reset()
	#Add suffix to the name if stack
	suffix = ""
	if (nbSlice > 1) :
		print "Process image "+str(slic+1)+"/"+str(nbSlice)
		suffix = adjustSizeNum(slic+1, 3)
	
	## ## ## ## ## ## ## ## ## ## ## 
	## Processing onto Nucleus image
	## ## ## ## ## ## ## ## ## ## ## 
	impNuclei.setSlice(slic+1)
	Nuclei_slice = impNuclei.crop("whole-slice")
	
	#Segmentation of nuclei with Stardist
	res = command.run(StarDist2D, False,
			 "input", Nuclei_slice, 
			 "modelChoice", "Versatile (fluorescent nuclei)",
			 "normalizeInput",True, 
			 "percentileBottom",1, "percentileTop",100.0,
			 "probThresh",0.5, "nmsThresh", 0.4, 
			 "outputType",outputType[2],
			 "nTiles",1, 
			 "excludeBoundary",2, 
			 "roiPosition", "Automatic",
			 "verbose",False, 
			 "showCsbdeepProgress",False, 
			 "showProbAndDist",False).get()	
	label = res.getOutput("label")
	labelname = "Label_Image_Nuclei"+suffix+".tif"
	impLabel = IJF.wrap( label, "Label"+suffix )
	IJ.run(impLabel, "Enhance Contrast", "saturated=0.35")
	#uiService.show(label)
	rm = RoiManager.getRoiManager()
	rm.runCommand(Nuclei_slice,"Deselect")
	#rm.runCommand("Save", os.path.join(imageDir,"RoiSet-Nuclei"+suffix+".zip"))
	
	
	#Manual threshold of nuclei based on area (only on the first image)
	rmAllNuclei = rm.getRoisAsArray()
	ip_AreaNuclei = ipNAN.duplicate()
	areaNuclei = []
	centroidNuclei = []
	for roi in rmAllNuclei :
		Nuclei_slice.setRoi(roi)
		stats = Nuclei_slice.getStatistics(Measurements.CENTROID+Measurements.AREA)
		areaN = stats.area
		ip_AreaNuclei.setColor(areaN)
		ip_AreaNuclei.fill(roi)
		areaNuclei.append(areaN)
		centroidNuclei.append([stats.xCentroid, stats.yCentroid])
	if slic == 0 :		
		areaThres = thresholdImageUI(ip_AreaNuclei, areaThres, max(areaNuclei))
	roi2Delete =[]
	
	#Create image with dot at the centroid of each nucleus (use for MaximumFinder)
	impResultNuclei = IJ.createImage("Result Nucleai", "16-bit black", width, height, 1)
	ipResultNuclei = impResultNuclei.getProcessor()
	for i, roi in enumerate(rmAllNuclei) :
		if areaNuclei[i] <= areaThres :
			roi2Delete.append(i)
		else:
			impResultNuclei.setRoi(roi)
			ipResultNuclei.setColor(Color.WHITE)
			centroidN = centroidNuclei[i]
			ipResultNuclei.fillOval(int(centroidN[0])-2,int(centroidN[1])-2,4,4)
	if roi2Delete :
		rm.setSelectedIndexes(roi2Delete)
		rm.runCommand(impLabel,"Delete")
		rm.runCommand(impLabel,"Deselect")
		rmAllNuclei = rm.getRoisAsArray() # all the nuclei with area above the thresholded value
		rm.setSelectedIndexes(range(rm.getCount()))
		rm.runCommand(impLabel,"Combine")
		IJ.setBackgroundColor(0, 0, 0)
		IJ.run(impLabel, "Clear Outside", "")	
	stackLabels.addSlice(impLabel.getProcessor())
	rm.reset()
	impResultNuclei.killRoi()
	impResultNuclei.updateAndDraw()
	
	# push image to GPU
	IJ.run(impResultNuclei, "32-bit", "")
	input = clij2.push(impResultNuclei)
	
	
	# reserve memory for output, same size and type as input
	thresholded = clij2.create(input)
	thresholded2 = clij2.create(input)
	voronoi_diagram = clij2.create(input)
	touch_matrix = clij2.create(input)
	inverted_voronoi = clij2.create(input)
	labelled = clij2.create(input)
	labelled_extended = clij2.create(input)
	intensity_values = clij2.create(input)
	intensity_map = clij2.create(input)
	neighbors_map = clij2.create(input)
	
	
	
	# Make Voronoi diagram
	clij2.thresholdOtsu(input, thresholded)
	clij2.copy(thresholded, thresholded2)
	clij2.release(input)
	clij2.voronoiOctagon(thresholded, voronoi_diagram)
	
	
	# convert to Float
	connected = clij2.create(thresholded2)
	centroid = clij2.create(thresholded2)
	density_values = clij2.create(thresholded2)
	clij2.connectedComponentsLabelingBox(thresholded2, connected)	
	clij2.reduceLabelsToCentroids(connected, centroid)
	
	
	# invert
	clij2.binaryNot(voronoi_diagram, inverted_voronoi)
	
	#Generate a label map and extend it to make labels touch
	clij2.connectedComponentsLabelingBox(inverted_voronoi, labelled)
	
	# Extend labels so that they touch
	clij2.maximum2DBox(labelled, labelled_extended, 1, 1)
		
	# Determine touch matrix
	clij2.generateTouchMatrix(labelled_extended, touch_matrix)

	# count neighbors and make a parametric image
	clij2.countTouchingNeighbors(touch_matrix, intensity_values)
	clij2.replaceIntensities(labelled_extended, intensity_values, intensity_map)

	# put black frames between cells
	clij2.mask(intensity_map, inverted_voronoi, neighbors_map)
	

	# convert the result back to imglib2 and show it
	impResultNuclei = clij2.pull(labelled_extended)#inverted_voronoi)
	impNeighborsMap = clij2.pull(neighbors_map)
	impNeighborsMap.setTitle("Neighbors Map "+suffix)
	IJ.run(impNeighborsMap, "Enhance Contrast", "saturated=0.35")
	stackNeighbors.addSlice(impNeighborsMap.getProcessor())
	
	#cleanup memory on GPU
	clij2.release(thresholded)
	clij2.release(voronoi_diagram)
	clij2.release(touch_matrix)
	clij2.release(inverted_voronoi)
	clij2.release(labelled)
	clij2.release(labelled_extended)
	clij2.release(intensity_values)
	clij2.release(intensity_map)
	clij2.release(neighbors_map)
	clij2.release(connected)
	
	#IJ.run(impResultNuclei, "Multiply...", "value=255")
	imageRoi = Roi(0,0,width,height)
	borderRoi = ShapeRoi(imageRoi).xor(ShapeRoi(RoiEnlarger.enlarge(imageRoi, -2)))
	rmCell = []
	areaCells = [] # array of cell area (in pixels^2)
	for roi in rmAllNuclei :
		impResultNuclei.setRoi(roi)
		stats = impResultNuclei.getStatistics(Measurements.CENTROID)
		impResultNuclei.killRoi()
		IJ.doWand(impResultNuclei, int(stats.xCentroid), int(stats.yCentroid), 0.0, "8-connected") #associate nuclei roi to cell roi
		roiCell = impResultNuclei.getRoi()
		interRoi = (borderRoi.clone()).and(ShapeRoi(roiCell))
		if interRoi.getLength()==0 :
			impResultNuclei.setRoi(roiCell)
			areaCells.append(impResultNuclei.getStatistics(Measurements.AREA).area)
			rmCell.append(roiCell)
			rm.addRoi(roiCell)
	
	rm.setSelectedIndexes(range(rm.getCount()))
	rm.runCommand(impResultNuclei,"Combine")
	IJ.setBackgroundColor(0, 0, 0);
	IJ.run(impResultNuclei, "Clear Outside", "")
	IJ.setThreshold(impResultNuclei, 0, 1)
	roiNAN = ThresholdToSelection.run(impResultNuclei)
	
	
	sumAreaCells = sum(areaCells) # in pixels^2
	binSize = int( round( math.sqrt( weigth * 2 * sumAreaCells/len(areaCells)) ) )	
	clij2.countNonZeroPixels2DSphere(centroid, density_values,binSize,binSize)
	impCellDensityMap = clij2.pull(density_values)
	clij2.release(centroid)
	clij2.release(density_values)
	ipCellDensityMap = impCellDensityMap.getProcessor()
	ipCellDensityMap.setColor(float('nan'))
	ipCellDensityMap.fill(roiNAN)	
		
		
	binSizePhys = binSize * pix2phys # conversion to mm
	binAreaInvert = 1/(math.pi * binSizePhys * binSizePhys) # in mm^-2
	IJ.run(impCellDensityMap, "Multiply...", "value="+str(binAreaInvert))
	impCellDensityMap.setTitle("Cell Density Map "+suffix)
	IJ.run(impCellDensityMap, "Enhance Contrast", "saturated=0.35")
	stackCellDensity.addSlice(impCellDensityMap.getProcessor())
	
	meanCellDensity =  len(rmCell)/ sumAreaCells / (pix2phys* pix2phys)
	print "Image "+str(slic+1)+"- Mean Cell Density (per mm2) = "+str(meanCellDensity)

impLabels=ImagePlus("Labels", stackLabels)
#IJ.run(impLabels, "glasbey", "")
IJ.saveAs(impLabels, "TIFF",os.path.join(imageDir,"impLabels"))
impLabels.show()
impCellDensity=ImagePlus("CellDensity", stackCellDensity)
IJ.run(impCellDensity, "Fire", "")
IJ.saveAs(impCellDensity, "TIFF",os.path.join(imageDir,"impCellDensity"))
impCellDensity.show()
impNeighbors=ImagePlus("Neighbors", stackNeighbors)
IJ.run(impNeighbors, "glasbey", "")
IJ.saveAs(impNeighbors, "TIFF",os.path.join(imageDir,"impNeighbors"))
impNeighbors.show()

print 'End of Processing'	