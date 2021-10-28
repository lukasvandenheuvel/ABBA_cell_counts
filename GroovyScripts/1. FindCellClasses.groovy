/* Detect And Cluster Cells, Count in All Brain Regions
    Rationale: To try and get all cells detected despite the background changes, we use cell detection on all channels
    Then we cluster the cells to get multiple positive cells when the detections overlap
    
    Finally, the clustered results are superimposed to the ABBA - Allen Brain BIOP Aligner
    
    This implies that the ABBA - Allen Brain BIOP Aligner has been run on this project from Fiji before starting
    https://biop.github.io/ijp-imagetoatlas/
    
    REQUIREMENTS
    ------------
    For QuPath: QuPath version 0.2.3 and the BIOP Tools extension
    https://c4science.ch/w/bioimaging_and_optics_platform_biop/image-processing/qupath2/
    
    FOLLOW-UP
    ---------
    The results should be exported using the `Export ABBA Results` script
    
    
    Author: Olivier Burri, EPFL - SV - PTECH - BIOP
    For Dr. bianca Ambrogina
    Date: 14.05/2021
*/ 

// User parameters, the settings for each channel can be copy pasted from a cell detection run
def categories = [ [ name: "TVA", parameters: {runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "TVA",  "requestedPixelSizeMicrons": 2.0,  "backgroundRadiusMicrons": 15.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.0,  "minAreaMicrons": 20.0,  "maxAreaMicrons": 700.0,  "threshold": 500.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}')} ],
                   [ name: "CTB", parameters: {runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "CTB",  "requestedPixelSizeMicrons": 2.0,  "backgroundRadiusMicrons": 15.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 3.0,  "minAreaMicrons": 50.0,  "maxAreaMicrons": 700.0,  "threshold": 27.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}')} ],
                   [ name: "Rabies", parameters: {runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "Rabies",  "requestedPixelSizeMicrons": 2.0,  "backgroundRadiusMicrons": 15.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 3.0,  "minAreaMicrons": 70.0,  "maxAreaMicrons": 1000.0,  "threshold": 360.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 0.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}')} ] 
                 ]

// Name of the DAPI threshold you want to set
def dapiThresholdName = "DAPI"

// Whether or not you want to remove RABIES+ cells
def removeRabies = true
def rabiesClassifier = "RabiesClassifier"

// ==== Script Start === 

clearAllObjects()
createSelectAllObject(true)
// Keep reference to rectangle to delete it later
def all = getSelectedObject()

// Run the cell detection for each channel
def cells = categories.collect{ cat ->
    println cat
    cat.parameters.run()
    // Retreive cells and set class.
    def cells = getDetectionObjects()
    cells.each{ it.setPathClass( getPathClass( cat.name ) ) }
    return cells
}.flatten()

if( removeRabies ) {
    // From the Rabnies+ cells, we are going to remove some detections (based on a classifier)
    def rabies = cells.findAll{ it.getPathClass().equals( getPathClass('Rabies') ) }
    
    println "Removing invalid Rabies+ cells based on an object classifier (added by Lukas)"
    def toRemove = filterDetectionsBasedOnClassfier(rabiesClassifier, rabies, 'Remove')
    cells.removeAll(toRemove) // Remove from list
    println toRemove.size() + " Rabies+ cells were removed."
    
    // The Classifier assigned to the Rabies+ cells that should be kept the class "Rabies: Keep".
    // We are now going to change the name back to "Rabies".
    rabies = cells.findAll{ it.getPathClass().toString().contains("Rabies") }
    rabies.each{ it.setPathClass( getPathClass('Rabies') ) }
}

// Make a copy of viable candidates, so we can remove them once they are chosen
def candidates = cells.collect()
println "Total Cells: " + cells.size()

println "Clustering..."

// This will contain the resulting classified cells
def allCells = []

cells.each{ c1 ->
    def pc = c1.getPathClass()
    
    if( candidates.contains( c1 ) ) {
        // Duplicate cell
        cell2 = PathObjects.createDetectionObject( c1.getROI(), c1.getPathClass(), c1.getMeasurementList() )
        
        // Remove current cell from candidates
        candidates.remove( c1 )
        
        // Find nearest cells contained within
        def nearest = findNearest( c1, candidates )
        
        // If there are cells contained in c1:
        if( nearest.size() > 0 ) {
            candidates.removeAll( nearest ) // Remove the nearest cells from candidates
            
            // Collect the classes that must be assigned to this detection.
            // That is, the classes of the nearest neighbours and the class of the cell c1 (is called pc).
            // Example: original class was Rabies, and nearest neighbours have classes TVA and CTB.
            // Then the new class should be Rabies : TVA : CTB.
            def classes = nearest.collect{ it.getPathClass() }
            classes.add( pc )
            // Sort classes so we do not end up with classes A:B and B:A
            classes.sort()
            
            // Add the new class
            cell2.setPathClass( getPathClass( classes.join( ":" ) ) )
        }
        
        allCells.add( cell2 )
       
    }
}

println "After clustering: "+allCells.size()
clearDetections()

addObjects(allCells)

//Add ABBA annotations
main = AtlasTools.getWarpedAtlasRegions( getCurrentImageData(), true )

// Remove rectangle
removeObject( all, true )

addObject( main )

selectAnnotations();
// Add DAPI Detection Area
addPixelClassifierMeasurements("DAPI", dapiThresholdName)

fireHierarchyUpdate() // Updates GUI, updates hierarchy (counting of detections within annotations is handled automatically by QuPath.)

// FUNCTIONS ----------------------------------------------

// Remove cells based on an object classifier (added by Lukas)
def filterDetectionsBasedOnClassfier(def classifierName, def cells, def removeLabel){
    // ----------
    // PARAMETERS
    // ----------
    // classifierName (string)
    // Name of the classifier to use. Train this classifier with Classify > Object Classification > Train Object Classifier.
    //
    // cells (list)
    // List of detection objects to classify.
    //
    // removeLabel (string) 
    // Label which the classifier gives to cells that should be removed (e.g. "Remove", if you trained your classifier that way).
    // --------- 

    // Load object classifier
    def imageData = getCurrentImageData()
    def classifier = loadObjectClassifier(classifierName)
    
    // Get a list of all Rabies+ cells, and apply the object classifier.
    // Rabies+ cells to keep will get the class label "Keep", all others will get the label "Remove"
    classifier.classifyObjects(imageData, cells, false)
    
    // Filter the objects with label "Remove"
    def toRemove = cells.findAll{ it.getPathClass().toString().contains(removeLabel) }
    return toRemove
}

// Convenience method to find nearest cells
def findNearest(def cell, def candidates ) {
    
    // It should be inside and not be of the same class as the cell being queried
    def inside = candidates.findAll{ it.getROI().contains( cell.getROI().getCentroidX() , cell.getROI().getCentroidY() ) 
                                     && !it.getPathClass().equals( cell.getPathClass() )
                                   }
 
    // Make sure that there is only one copy of each category   
    inside.unique{ it.getPathClass() }
    return inside
}
import ch.epfl.biop.qupath.atlas.allen.api.*