
def measurements = ["Name", "Class", "Num Detections", "Num CTB", "Num Rabies", "Num TVA", "Num CTB: Rabies", "Num CTB: TVA", "Num CTB: Rabies: TVA", "Num Rabies: TVA", "DAPI: DAPI area "+GeneralTools.micrometerSymbol()+"^2"]

def annotations = getAnnotationObjects()

def resultsfolder = buildFilePath(PROJECT_BASE_DIR, "results")
mkdirs( resultsfolder )

// def imageName = getCurrentServer().getMetadata().getName()
def imageName = getProjectEntry().getImageName() // adjusted by Lukas, s.t. the image name includes 'LEFT' and/or 'RIGHT'

def resultsfile = new File(resultsfolder, imageName+"_regions.txt")

Utils.sendResultsToFile( measurements, annotations, resultsfile)

println("Completed")

import ch.epfl.biop.qupath.utils.*