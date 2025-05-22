#@ OpService ops
#@ File inputFile
#@ File(style='directory') outputDir

# AUTHOR Thomas Pengo, tpengo@umn.edu, 2021

from ij import IJ
from loci.plugins import BF
from net.imglib2.img.display.imagej import ImageJFunctions
from net.imglib2.algorithm.neighborhood import RectangleShape
from net.imglib2.img.array import ArrayImgs
from ij.plugin import ImageCalculator

imp, = BF.openImagePlus(inputFile.getAbsolutePath())
name = inputFile.getName()

from ij.plugin import ChannelSplitter
c1, c2 = ChannelSplitter.split(imp)
imp = c2

import re
pattern = r'^(?P<Animal>a[0-9][0-9]*) (?P<Genotype>.*) (?P<Wavelength1>[0-9]{3}) (?P<Protein1>.*) (?P<Wavelength2>[0-9]{3}) (?P<Protein2>.*)-(?P<SeqNo>[0-9][0-9]+)'
p = re.compile(pattern)

m = p.search(name)

# CREATE MASK FROM CH 1 (MOMENTS STACK APPLY)
# DILATION 1
#mask_1 = ImageJFunctions.wrap(ops.morphology().dilate(ops.threshold().moments(ImageJFunctions.wrap(c1)), RectangleShape(1, False)),c1.title+" ch1 MASK")
IJ.setAutoThreshold(c1, "Moments dark stack no-reset")
IJ.run(c1, "Convert to Mask", "black")
mask_1 = c1

# BKG SUB RAD=2
c2_0 = c2.duplicate()
c2_0.setTitle('TOQUANT')
IJ.run(c2, "Subtract Background...", "rolling=2 stack");

# CONTRAST STRETCH 0-65535 -> 8-bit
c2.getProcessor().setMinAndMax(0,65535)
IJ.run(c2, "8-bit", "")

# AUTO LOCAL THRESHOLD "BERNSEN" RAD 10
IJ.run(c2, "Auto Local Threshold", "method=Bernsen radius=10 parameter_1=0 parameter_2=0 white stack");
mask_2 = c2

# APPLY MASK FROM CH 1
ic = ImageCalculator();
imp = ic.run("Multiply create stack", mask_1, mask_2);

c2_0.show()
c2.show()
c1.show()
imp.show()

# 3D OBJECT COUNT, MIN 10 voxels
IJ.run("3D OC Options", "volume nb_of_obj._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid centre_of_mass dots_size=10 font_size=9 show_numbers white_numbers redirect_to=TOQUANT");
IJ.run(imp, "3D Objects Counter", "threshold=128 slice=10 min.=10 max.=5505024 statistics objects");

# EQ RADIUS CALCULATION
import math
def getRadius(vol):
	return math.pow(vol*3/4/3.14159265,1.0/3)

# SAVE RESULTS
from ij.measure import ResultsTable
rt = ResultsTable.getResultsTable()
for i in range(rt.size()):
	vol = rt.getValueAsDouble(0,i)
	rt.setValue("EqSphereRadius",i,getRadius(vol))

	for k in m.groupdict():
		rt.setValue(k,i,m.groupdict()[k])

IJ.saveAsTiff(IJ.getImage(), "{}/{}".format(outputDir.getAbsolutePath(),imp.getTitle().replace('.czi','_objects.tiff')))
IJ.getImage().close()

IJ.saveAs("Results", "{}/{}".format(outputDir.getAbsolutePath(),imp.getTitle().replace('.czi','_results.csv')))

IJ.run("Close All")
