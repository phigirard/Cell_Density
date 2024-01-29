#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	CellDensity.py
#	Release v2.0
#
#	Copyright 2024 - BSD-3-Clause license
#                                                                             
#******************************************************************************/


#@ File nucleiFile (label="Select the nuclei image:", style="file") 
#@ File impFile (label="Select the 2nd image (for thresholding step):", style="file")  

#@ DatasetIOService io
#@ UIService uiService
#@ LogService log
#@ CommandService command

#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#import class ImageJ1
from ij import IJ, ImagePlus,ImageStack, WindowManager, Prefs
from ij.io import Opener
from ij.gui import  GenericDialog, Roi, ShapeRoi,WaitForUserDialog
from ij.plugin import RoiEnlarger
from ij.plugin import ImageCalculator as IC
from ij.plugin.frame import RoiManager, ThresholdAdjuster
from ij.measure import ResultsTable , Measurements, Calibration
from ij.plugin.filter import Analyzer, MaximumFinder, ThresholdToSelection
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.process import ImageProcessor,ImageConverter, ImageStatistics

from net.imglib2.img.display.imagej import ImageJFunctions as IJF #wrapper to convert Dataset from Stardist to ImagePlus


import os
import sys
import csv
import math 

from java.lang import Double
from java.awt import Color
from java.awt.event import AdjustmentListener  
#---------------------------------------------------------------------------------------------------------
#import StarDist plugin
from de.csbdresden.stardist import StarDist2D, Opt
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------
#---------------          CONSTANTS            -----------------
#---------------------------------------------------------------

Prefs.blackBackground = True

#Constants for MaximaFinder
excludeOnEdges = False
tolerance = 1

#Constants for CSV Files
IJ.run("Conversions...", "scale")
IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column save_row")

outputType = ["ROI Manager","Label Image","Both","Polygons"] #output for Stardist plugin

MAXSIZE = Double.POSITIVE_INFINITY # for ParticleAnalyzer

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

			
			
#convert Files from #@ parameters to String and extract the main directory of the data
nucleiPath = nucleiFile.getCanonicalPath()
impPath = impFile.getCanonicalPath()


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
impImg = Opener().openImage(impPath)
IJ.run(impImg, "Grays", "")
IJ.run(impImg, "Smooth", "stack")

# Get calibration and set time & pixel units 
cal = impNuclei.getCalibration()
pix2phys = cal.getX(1)	
physUnit = cal.getXUnit()


# Remove physical calibration for processing steps
impNuclei.setCalibration(Calibration())


# Get width, height and frame 
nbSlice = impNuclei.getStackSize()
width = impNuclei.width
height = impNuclei.height
depth = impNuclei.getBitDepth()


# Create Ratio and area (nuclei and cell) Images  
stack_Nuclei = ImageStack(width, height)

impNAN = IJ.createImage("impNAN", "32-bit black" , width, height, 1)
changeValue2NAN(impNAN,0)
ipNAN = impNAN.getProcessor()
areaThres =0

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
	implabel = IJF.wrap( label, "Label"+suffix )
	#uiService.show(label)
	rm = RoiManager.getRoiManager()
	rm.runCommand(Nuclei_slice,"Deselect")
	#rm.runCommand("Save", os.path.join(imageDir,"RoiSet-Nuclei"+suffix+".zip"))
	
	
	#Manual threshold of nuclei based on area (only on the first image)
	rmAllNuclei = rm.getRoisAsArray()
	ip_AreaNuclei = ipNAN.duplicate()
	areaNuclei = []
	for roi in rmAllNuclei :
		Nuclei_slice.setRoi(roi)
		areaN = Nuclei_slice.getStatistics(Measurements.AREA).area
		ip_AreaNuclei.setColor(areaN)
		ip_AreaNuclei.fill(roi)
		areaNuclei.append(areaN)
	if slic == 0 :		
		areaThres = thresholdImageUI(ip_AreaNuclei, min(areaNuclei) -1, max(areaNuclei))
	roi2Delete =[]
	roiReduced = False
	for i, areaN in enumerate(areaNuclei) :
		if areaN <= areaThres :
			roiReduced =True
			roi2Delete.append(i)
	if roiReduced :
		rm.setSelectedIndexes(roi2Delete)
		rm.runCommand(Nuclei_slice,"Delete")
		rm.runCommand(Nuclei_slice,"Deselect")
		rmAllNuclei = rm.getRoisAsArray()	
	rm.reset()
	
	#Create image with dot at the centroid of each nucleus (use for MaximumFinder)
	impResultNuclei = IJ.createImage("Result Nuclei", "16-bit black", width, height, 1)
	ipResultNuclei = impResultNuclei.getProcessor()
	for roi in rmAllNuclei :
		impResultNuclei.setRoi(roi)
		stats = impResultNuclei.getStatistics(Measurements.CENTROID)
		ipResultNuclei.setColor(Color.WHITE)
		ipResultNuclei.fillOval(int(stats.xCentroid)-2,int(stats.yCentroid)-2,4,4)
	impResultNuclei.killRoi()
	impResultNuclei.updateAndDraw()
	#rm.reset()
	

	## ## ## ## ## ## ## ## ## ## ## 
	## Processing onto 2nde image
	## ## ## ## ## ## ## ## ## ## ## 
	
	# Extract slice from imp for the thesholding step
	impImg.setSlice(slic+1)
	imp_slice = impImg.crop("whole-slice")
	imp_slice.setCalibration(Calibration())
	
	# Threshold slice to get the signal/background ROIs
	if (slic == 0):
		IJ.run(imp_slice, "Find Edges", "")
		IJ.run(imp_slice, "Mean...", "radius=5")
		imp_slice.show()
		ta = ThresholdAdjuster()
		ta.setMethod("Triangle")
		ta.show()
		ta.update()
		IJ.setAutoThreshold(imp_slice, "Triangle dark")
		waitDialog = WaitForUserDialog("Manual threshold", "Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')") # human thresholding
		waitDialog.show()
		thres_min = imp_slice.getProcessor().getMinThreshold()
		thres_max = imp_slice.getProcessor().getMaxThreshold()
		ta.close()
		imp_slice.hide()
	IJ.setThreshold(imp_slice, thres_min, thres_max)
	cellsRoi = ShapeRoi(ThresholdToSelection.run(imp_slice))
	
	
	#Extract Rois of each cell (except at the edge of the image) based on MaximumFinder 
	mf = MaximumFinder()
	maxip = mf.findMaxima(ipResultNuclei, tolerance, MaximumFinder.SEGMENTED, excludeOnEdges)
	impResultNuclei=ImagePlus("Result Nuclei",maxip)
	rmCell = []
	rmNuclei = [] #not necessary
	imageRoi = Roi(0,0,width,height)
	borderRoi = ShapeRoi(imageRoi).xor(ShapeRoi(RoiEnlarger.enlarge(imageRoi, -2)))
	for roi in rmAllNuclei:
		impResultNuclei.setRoi(roi)
		stats = impResultNuclei.getStatistics(Measurements.CENTROID)
		impResultNuclei.killRoi()
		IJ.doWand(impResultNuclei, int(stats.xCentroid), int(stats.yCentroid), 0.0, "8-connected") #associate nuclei roi to cell roi
		cellRoi = impResultNuclei.getRoi()
		interRoi = (borderRoi.clone()).and(ShapeRoi(cellRoi))
		interRoiCells = (cellsRoi.clone()).and(ShapeRoi(roi))
		if interRoi.getLength()==0 and interRoiCells.getLength()!=0 :
			rmCell.append(cellRoi) 
			rmNuclei.append(roi) #not necessary
			rm.addRoi(cellRoi)
	rm.runCommand("Save",os.path.join(imageDir,"RoiSet-Cells"+suffix+".zip"))
	rm.setSelectedIndexes(range(rm.getCount()))
	rm.runCommand(impResultNuclei,"Combine")
	IJ.setBackgroundColor(0, 0, 0);
	IJ.run(impResultNuclei, "Clear Outside", "")
	#roiAll= impResultNuclei.getRoi()
	rm.reset()
	for roiN in rmNuclei:
		rm.addRoi(roiN)
	rm.runCommand("Save",os.path.join(imageDir,"RoiSet-Nucleus"+suffix+".zip"))
	rm.reset()
	#rm.addRoi(roiAll)
	areaCells = impResultNuclei.getStatistics(Measurements.AREA).area
	densityPix = len(rmCell)/ areaCells
	print "The cell area (in pixels^2) = ",areaCells
	print "Number of cells = ", len(rmCell)
	print "The cell density (per pixels^2)  = ",densityPix
	if physUnit == "micron" :
		pix2phys = pix2phys/1000
		densityPhys=densityPix/pix2phys/pix2phys
		print "The cell density (per mm^2)  = ",densityPhys

	stack_Nuclei.addSlice(impResultNuclei.getProcessor())
	
	
	


impCell=ImagePlus("impCell", stack_Nuclei)
IJ.saveAs(impCell, "TIFF",os.path.join(imageDir,"impCell"))
impCell.show()

print 'End of Processing'	