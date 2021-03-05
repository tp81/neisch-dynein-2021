#@ File inputFile
#@ File(style='directory') outputDir
#@ Integer(value=1, label="HRP channel (0:first, 1:second, ...)") HRP_channel
#@ Integer(value=0, label="CH2 channel (0:first, 1:second, ...)") CH2_channel
#@ OpService ops
#@ UIService ui

# AUTHOR Thomas Pengo, tpengo@umn.edu, 2021

# PARAMETER
regionSize = 47 # pixels (~10um)

from ij import IJ
from loci.plugins import BF
from net.imglib2.img.display.imagej import ImageJFunctions
from net.imglib2.view import Views
from net.imglib2.roi import Regions, Masks, BinaryMaskRegionOfInterest
from net.imglib2.type.logic import BitType
from net.imglib2.type.numeric.real import DoubleType
from ij.plugin import Duplicator

image, = BF.openImagePlus(inputFile.getAbsolutePath())
t = image.getTitle();

### SELECT REGION OF INTEREST

# Z Projection
from ij.plugin import ZProjector
image_mip = ZProjector.run(image,"max")

# ROI CREATION
image_mip.show()
image_mip.setRoi(image_mip.width/2,image_mip.height/2,regionSize,regionSize)

from ij.gui import WaitForUserDialog
myWait = WaitForUserDialog ("Choose ROI", "Please position the ROI.. Press 'OK' when done.")
myWait.show()

# CROP
image.setRoi(image_mip.getRoi())
from ij import IJ
image = Duplicator().run(image)
image.show()

image_mip.close()

image = ImageJFunctions.wrap(image)

# 
def roi(mask, image):
	# Convert ROI from R^n to Z^n.
	#discreteROI = Views.raster(Masks.toRealRandomAccessible(mask))
	# Apply finite bounds to the discrete ROI.
	boundedDiscreteROI = Views.interval(mask, image)
	# Create an iterable version of the finite discrete ROI.
	iterableROI = Regions.iterable(boundedDiscreteROI)
	return Regions.sample(iterableROI, image)


from ij import IJ
from loci.plugins import BF
from net.imglib2 import FinalInterval

w = image.dimension(0)
h = image.dimension(1)
ch = image.dimension(2)
d = image.dimension(3)

# CROP EACH CHANNEL
HRP, CH2 = [ ops.run(
    "crop", 
    image, 
    FinalInterval.createMinSize(0,0,i,0,w,h,1,d),
    True) 
                 for i in [HRP_channel, CH2_channel] ]

### BACKGROUND SUBTRACTION
                 
# Threshold HRP
HRP_m = ops.threshold().otsu(HRP)
HRP_m_f = ops.convert().float32(HRP_m)

# Subtract background
def subBkg(im):
	im_m_f = ops.convert().float32(ops.threshold().otsu(im))
	out = ops.create().img(im_m_f)
	im_m_f = ops.threshold().otsu(ops.run('invert',out,im_m_f,DoubleType(0),DoubleType(1)))
	bkg = ops.stats().mean(roi(im_m_f,im)).getRealFloat()

	im_bkg_sub = ops.run('subtract',ops.convert().int16(im),bkg)

	# Clip image to the input type
	clipped = ops.create().img(im_bkg_sub, im.firstElement())
	clip_op = ops.op("convert.clip", im.firstElement(), im_bkg_sub.firstElement())
	ops.convert().imageType(clipped, im_bkg_sub, clip_op)

	return bkg, clipped
	
HRP_bkg, HRP = subBkg(HRP)
CH2_bkg, CH2 = subBkg(CH2)

#### MEMBRANE MASK
from net.imglib2.algorithm.neighborhood import HyperSphereShape
shape = HyperSphereShape(2)

# Opening operation
HRP_m = ops.run("erode",None,HRP_m,shape,False)
HRP_m = ops.run("dilate",None,HRP_m,shape,False)

HRP_m_dilated = ops.run("dilate",None,HRP_m,shape,False)
HRP_m_eroded = ops.run("erode",None,HRP_m,shape,False)
membrane = ops.eval('dil-ero',{'dil':HRP_m_dilated,'ero':HRP_m_eroded})

ui.show(membrane)

#### RATIO CALCULATION
# Calculate mean(CH2(membrane))/mean(HRP(membrane))
CH2_mean = ops.stats().mean(roi(membrane,CH2)).getRealFloat()
HRP_mean = ops.stats().mean(roi(membrane,HRP)).getRealFloat()
ratio = CH2_mean / HRP_mean


#### RATIO IMAGE
ratio_image = ops.eval('(membrane*CH2)/(membrane*HRP)',{
	'membrane' : ops.convert().float32(membrane),
	'CH2' : ops.convert().float32(CH2),
	'HRP'  : ops.convert().float32(HRP) })

ui.show(ratio_image)

from ij import IJ
IJ.run('physics')
imp = IJ.getImage()
imp.setDisplayRange(0, 5);
IJ.saveAs('TIFF','{}/{}_ratioImage.tiff'.format(outputDir.getAbsolutePath(),t.replace('.czi','')))


#### SAVE TO TABLE

import os
outputFilePath = os.path.join(outputDir.getAbsolutePath(),'RESULTS.xls')
if not os.path.exists(outputFilePath):
	with open(outputFilePath,'w') as f:
		f.write('Filename\tHRP_bkg\tCH2_bkg\tHRP_mean_membrane\tCH2_mean_membrane\tratio_membrane\n')

with open(outputFilePath,'a') as f:
	f.write('\t'.join([t, str(HRP_bkg), str(CH2_bkg), str(HRP_mean), str(CH2_mean), str(ratio)])+'\n')

	
	
