// This script sets the correct channel names of an image.
// ch0 = DAPI
// ch1 = Rabies
// ch2 = TVA
// ch3 = CTB
//
// Then, it removes all existing annotations and detections from the image,
// and loads the ABBA annotations into it.
//
// Created 14/09/2021 by Lukas van den Heuvel.

setChannelNames(
     'DAPI',
     'Rabies',
     'TVA',
     'CTB'
)

def splitLeftRight = true

// Remove all anotations
def annotations = getAnnotationObjects()
removeObjects(annotations, false)

// Remove all detections
def detections = getDetectionObjects()
removeObjects(detections, false)

// Programmatically load all annotations
AtlasTools.loadWarpedAtlasAnnotations( getCurrentImageData(), splitLeftRight )
fireHierarchyUpdate() // Updates GUI, updates hierarchy (counting of detections within annotations is handled automatically by QuPath.)

import ch.epfl.biop.qupath.atlas.allen.api.*
