// This script lists all parent annotations within annotations that are part of the 'Exclude' class.
// Created 02/02/2022 by Lukas van den Heuvel

import qupath.lib.regions.*
resetSelection()

def imageName = getProjectEntry().getImageName()
def regionsToExclude = []
def excludeAnnotations = getAnnotationObjects().findAll {it.getPathClass() == getPathClass("Exclude")}
def brainRegions = getAnnotationObjects().findAll {it.getPathClass() != getPathClass("Exclude")}

// Loop over annotations that contain the annotations to be excluded
excludeAnnotations.each{ann ->
    def roi = ann.getROI()
    
    brainRegions.each{region ->
        // Get the name of each brain region and the points of its boundary
        def allPoints = region.getROI().getAllPoints()
        def regionInROI = true
        // Loop over all points in the boundary of the brain region
        allPoints.each{pt ->
            // If there is at least 1 point not contained in the exclude annotation,
            // we should not exclude this region
            if (!roi.contains(pt.getX(), pt.getY())){
                regionInROI = false
            }
        }
        // If no points outside the exlude annotations were found,
        // add this region to regionsToExclude.
        if (regionInROI == true){
            regionsToExclude << region
        }
    }
}

// Remove 'child' annotations that are descendents of a parent
// We do this because regions to be excluded may not overlap!
hierarchy = getCurrentHierarchy()
parentRegionsToExclude = regionsToExclude.collect()

regionsToExclude.each{r1 ->
     annotation1 = getAnnotationObjects().findAll {it.getPathClass() == r1.getPathClass()}
     descendants = r1.getChildObjects()
     regionsToExclude.each{r2 ->
         if (descendants.contains(r2)){
             parentRegionsToExclude -= r2
         }
     }
}

// Print and select regions to exclude
print(parentRegionsToExclude)
hierarchy.getSelectionModel().selectObjects(parentRegionsToExclude)

// Save results to txt
def excludeFolder = buildFilePath(PROJECT_BASE_DIR, "regions_to_exclude")
mkdirs( excludeFolder )

def path = buildFilePath(excludeFolder, imageName+'_regions_to_exclude.txt')
outFile = new File(path)

outFile.withPrintWriter{pw ->
    parentRegionsToExclude.each{region ->
        pw.println(region.getPathClass())
    }
}
