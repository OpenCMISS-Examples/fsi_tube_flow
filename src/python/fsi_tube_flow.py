#> This example program solves a FSI problem in a multi-block tube using OpenCMISS.
#>
#> By Chris Bradley
#>
#>

#================================================================================================================================
#  Start Program
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import math,numpy,csv,time,sys,os,pdb
from opencmiss.iron import iron

# Ensure output directories exist
if not os.path.exists('./output'):
    os.makedirs('./output')
if not os.path.exists('./output/Fluid'):
    os.makedirs('./output/Fluid')
if not os.path.exists('./output/Solid'):
    os.makedirs('./output/Solid')
if not os.path.exists('./output/Interface'):
    os.makedirs('./output/Interface')

LINEAR = 1
QUADRATIC = 2

numberOfSquareElements = 2
numberOfArmElements = 2
numberOfWallElements = 1
numberOfLengthElements = 2

pipeRadius = 0.015
lengthSize = 0.3
squareSizeRatio = 0.500
wallThickness = 0.002

uInterpolation = QUADRATIC
pInterpolation = LINEAR

# Time stepping parameters
startTime = 0.0
stopTime  = 15.0
timeStep  = 0.1

## Inlet velocity parameters
A = 0.25
B = 0.25
C = 5.0

# Material parmeters
# Fluid
Re = 1000

fluidDensity = 1060.0
#fluidDynamicViscosity = fluidDensity*(A+B)*(2.0*pipeRadius)/Re
fluidDynamicViscosity = 0.0035

fluidPInit = 0.0

# Solid
solidDensity = 1000.0
mooneyRivlin1 = 0.1                       # N / m^2
mooneyRivlin2 = 0.2

#solidPInit=-mooneyRivlin1
solidPInit=0.0

# Moving mesh
movingMeshKParameter   = 1.0       #default

RBS = False
#RBS = True

outputFrequency = 1 # Result output frequency

# Output flags
fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = iron.EquationsOutputTypes.NONE
#fluidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#fluidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#fluidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
fluidDynamicSolverOutputType = iron.SolverOutputTypes.NONE
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidLinearSolverOutputType = iron.SolverOutputTypes.NONE
#fluidLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidLinearSolverOutputType = iron.SolverOutputTypes.MATRIX
solidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#solidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
solidEquationsOutputType = iron.EquationsOutputTypes.NONE
#solidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#solidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#solidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#movingMeshEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
movingMeshEquationsOutputType = iron.EquationsOutputTypes.NONE
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.TIMING
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#movingMeshEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.NONE
#interfaceConditionOutputType = iron.InterfaceConditionOutputTypes.PROGRESS
interfaceEquationsOutputType = iron.EquationsOutputTypes.NONE
#interfaceEquationsOutputType = iron.EquationsOutputTypes.TIMING
#interfaceEquationsOutputType = iron.EquationsOutputTypes.PROGRESS
#interfaceEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#interfaceEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
# (NoOutput/ProgressOutput/TimingOutput/SolverOutput/SolverMatrixOutput)
movingMeshLinearSolverOutputType = iron.SolverOutputTypes.NONE
#movingMeshLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#movingMeshLinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fsiDynamicSolverOutputType = iron.SolverOutputTypes.NONE
#fsiDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fsiDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fsiNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fsiLinearSolverOutputType = iron.SolverOutputTypes.NONE
#fsiLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fsiLinearSolverOutputType = iron.SolverOutputTypes.MATRIX

# Set solver parameters
fsiDynamicSolverTheta    = [0.5]
nonlinearMaximumIterations      = 100000000 #default: 100000
nonlinearRelativeTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearSolutionTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearAbsoluteTolerance      = 1.0E-8    #default: 1.0E-10
nonlinearMaxFunctionEvaluations = 10000
nonlinearLinesearchAlpha        = 1.0
linearMaximumIterations      = 100000000 #default: 100000
linearRelativeTolerance      = 1.0E-6    #default: 1.0E-05
linearAbsoluteTolerance      = 1.0E-6    #default: 1.0E-10
linearDivergenceTolerance    = 1.0E5     #default: 1.0E5
linearRestartValue           = 30        #default: 30

progressDiagnostics = True
debug = True

#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================

if numberOfLengthElements == 0:
    numberOfDimensions = 2
else:
    numberOfDimensions = 3

fluidCoordinateSystemUserNumber     = 1
solidCoordinateSystemUserNumber     = 2
interfaceCoordinateSystemUserNumber = 3
  
solidRegionUserNumber = 1
fluidRegionUserNumber = 2
interfaceUserNumber   = 3

uBasisUserNumber = 1
pBasisUserNumber = 2
interfaceBasisUserNumber = 3

fluidMeshUserNumber     = 1
solidMeshUserNumber     = 2
interfaceMeshUserNumber = 3
movingMeshUserNumber    = 4
  
fluidDecompositionUserNumber     = 1
solidDecompositionUserNumber     = 2
interfaceDecompositionUserNumber = 3
  
fluidGeometricFieldUserNumber = 11
fluidEquationsSetFieldUserNumber = 12
fluidDependentFieldUserNumber = 13
fluidMaterialsFieldUserNumber = 14
fluidIndependentFieldUserNumber = 15
bcCellMLModelsFieldUserNumber = 16
bcCellMLStateFieldUserNumber = 17
bcCellMLParametersFieldUserNumber = 18
bcCellMLIntermediateFieldUserNumber = 19

solidGeometricFieldUserNumber     = 21
solidFibreFieldUserNumber     = 22
solidEquationsSetFieldUserNumber = 23
solidDependentFieldUserNumber = 24
solidMaterialsFieldUserNumber = 25
solidSourceFieldUserNumber = 26
  
movingMeshEquationsSetFieldUserNumber = 31
movingMeshDependentFieldUserNumber    = 32
movingMeshMaterialsFieldUserNumber    = 33
movingMeshIndependentFieldUserNumber  = 34

interfaceGeometricFieldUserNumber = 41
interfaceLagrangeFieldUserNumber  = 42

fluidEquationsSetUserNumber  = 1
solidEquationsSetUserNumber  = 2
movingMeshEquationsSetUserNumber = 3

bcCellMLUserNumber = 1
  
interfaceConditionUserNumber = 1
  
fsiProblemUserNumber = 1
 
#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

# Get the computational nodes info
computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

#iron.OutputSetOn("Testing")

#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

# Create a RC coordinate system for the fluid region
fluidCoordinateSystem = iron.CoordinateSystem()
fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber)
fluidCoordinateSystem.DimensionSet(3)
fluidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the solid region
solidCoordinateSystem = iron.CoordinateSystem()
solidCoordinateSystem.CreateStart(solidCoordinateSystemUserNumber)
solidCoordinateSystem.DimensionSet(3)
solidCoordinateSystem.CreateFinish()
# Create a RC coordinate system for the interface region
interfaceCoordinateSystem = iron.CoordinateSystem()
interfaceCoordinateSystem.CreateStart(interfaceCoordinateSystemUserNumber)
interfaceCoordinateSystem.DimensionSet(3)
interfaceCoordinateSystem.CreateFinish()

if (progressDiagnostics):
    print('Coordinate systems ... Done')
  
#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

# Create a fluid region
fluidRegion = iron.Region()
fluidRegion.CreateStart(fluidRegionUserNumber,iron.WorldRegion)
fluidRegion.label = 'FluidRegion'
fluidRegion.coordinateSystem = fluidCoordinateSystem
fluidRegion.CreateFinish()
# Create a solid region
solidRegion = iron.Region()
solidRegion.CreateStart(solidRegionUserNumber,iron.WorldRegion)
solidRegion.label = 'SolidRegion'
solidRegion.coordinateSystem = solidCoordinateSystem
solidRegion.CreateFinish()

if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')
    
numberOfNodesXi = uInterpolation+1
numberOfGaussXi = uInterpolation+1

uBasis = iron.Basis()
uBasis.CreateStart(uBasisUserNumber)
uBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
uBasis.numberOfXi = 3
if (uInterpolation == LINEAR):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
elif (uInterpolation == QUADRATIC):
    uBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*3
else:
    print('Invalid u interpolation')
    exit()
uBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
uBasis.CreateFinish()

pBasis = iron.Basis()
pBasis.CreateStart(pBasisUserNumber)
pBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
pBasis.numberOfXi = 3
pBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
pBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
pBasis.CreateFinish()

interfaceBasis = iron.Basis()
interfaceBasis.CreateStart(interfaceBasisUserNumber)
interfaceBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
interfaceBasis.numberOfXi = 2
if (uInterpolation == LINEAR):
    interfaceBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*2
elif (uInterpolation == QUADRATIC):
    interfaceBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*2
else:
    print('Invalid u interpolation')
    exit()
interfaceBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*2
interfaceBasis.CreateFinish()

numberOfLocalNodes = numberOfNodesXi*numberOfNodesXi*numberOfNodesXi
numberOfLocalInterfaceNodes = numberOfNodesXi*numberOfNodesXi
localNodeIdx000 = 0
localNodeIdx100 = numberOfNodesXi-1
localNodeIdx010 = numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx110 = numberOfNodesXi*numberOfNodesXi-1
localNodeIdx001 = numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx101 = numberOfNodesXi-1+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx011 = numberOfNodesXi*(numberOfNodesXi-1)+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx111 = numberOfLocalNodes-1

if (progressDiagnostics):
    print('Basis functions ... Done')
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

numberOfFluidNodesPerBlock = numberOfSquareElements*(numberOfNodesXi-1)*(numberOfArmElements*(numberOfNodesXi-1)+1)
numberOfFluidElementsPerBlock = numberOfSquareElements*numberOfArmElements
numberOfFluidNodesPerLength = 4*numberOfFluidNodesPerBlock+ \
                     (numberOfSquareElements*(numberOfNodesXi-1)-1)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
numberOfFluidElementsPerLength = 4*numberOfFluidElementsPerBlock+numberOfSquareElements*numberOfSquareElements
numberOfFluidNodes = numberOfFluidNodesPerLength*(numberOfLengthElements*(numberOfNodesXi-1)+1)
numberOfFluidElements = numberOfFluidElementsPerLength*numberOfLengthElements
numberOfSolidCircumfrentialElements = 4*numberOfSquareElements
numberOfSolidCircumfrentialNodes = numberOfSolidCircumfrentialElements*(numberOfNodesXi-1)
numberOfSolidLengthNodes = numberOfLengthElements*(numberOfNodesXi-1)+1
numberOfSolidWallNodes = numberOfWallElements*(numberOfNodesXi-1)+1
numberOfSolidNodesPerWall = numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
numberOfSolidNodes = numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes*numberOfSolidWallNodes
numberOfSolidElements = numberOfSolidCircumfrentialElements*numberOfLengthElements*numberOfWallElements
numberOfInterfaceNodes = numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
numberOfInterfaceElements = numberOfSolidCircumfrentialElements*numberOfLengthElements

if (debug):
    print('  Mesh Parameters:')
    print('    numberOfSquareElements: %d' % (numberOfSquareElements))
    print('    numberOfArmElements: %d' % (numberOfArmElements))
    print('    numberOfLengthElements: %d' % (numberOfLengthElements))
    print('    numberOfWallElements: %d' % (numberOfWallElements))
    print('    numberOfNodesXi: %d' % (numberOfNodesXi))
    print('    numberOfFluidNodesPerBlock: %d' % (numberOfFluidNodesPerBlock))
    print('    numberOfElementPerBlock: %d' % (numberOfFluidElementsPerBlock))
    print('    numberOfFluidNodesPerLength: %d' % (numberOfFluidNodesPerLength))
    print('    numberOfFluidElementsPerLength: %d' % (numberOfFluidElementsPerLength))
    print('    numberOfFluidNodes: %d' % (numberOfFluidNodes))
    print('    numberOfFluidElements: %d' % (numberOfFluidElements))
    print('    numberOfSolidCircumfrentialNodes: %d' % (numberOfSolidCircumfrentialNodes))
    print('    numberOfSolidCircumfrentialElements: %d' % (numberOfSolidCircumfrentialElements))
    print('    numberOfSolidLengthNodes: %d' % (numberOfSolidLengthNodes))
    print('    numberOfSolidWallNodes: %d' % (numberOfSolidWallNodes))
    print('    numberOfSolidNodesPerWall: %d' % (numberOfSolidNodesPerWall))
    print('    numberOfSolidNodes: %d' % (numberOfSolidNodes))
    print('    numberOfSolidElements: %d' % (numberOfSolidElements))
    print('    numberOfInterfaceNodes: %d' % (numberOfInterfaceNodes))
    print('    numberOfInterfaceElements: %d' % (numberOfInterfaceElements))
        
fluidNodes = iron.Nodes()
fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
fluidNodes.CreateFinish()

solidNodes = iron.Nodes()
solidNodes.CreateStart(solidRegion,numberOfSolidNodes)
solidNodes.CreateFinish()

fluidMesh = iron.Mesh()
fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,3)
fluidMesh.NumberOfElementsSet(numberOfFluidElements)
fluidMesh.NumberOfComponentsSet(2)

solidMesh = iron.Mesh()
solidMesh.CreateStart(solidMeshUserNumber,solidRegion,3)
solidMesh.NumberOfElementsSet(numberOfSolidElements)
solidMesh.NumberOfComponentsSet(2)

if (debug):
    print('  Fluid Elements:')

fluidUElements = iron.MeshElements()
fluidUElements.CreateStart(fluidMesh,1,uBasis)
fluidPElements = iron.MeshElements()
fluidPElements.CreateStart(fluidMesh,2,pBasis)

for zElementIdx in range(1,max(numberOfLengthElements+1,2)):
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yElementIdx in range(1,numberOfArmElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = xElementIdx+(yElementIdx-1)*numberOfSquareElements+(blockIdx-1)*numberOfSquareElements*numberOfArmElements+\
                                (zElementIdx-1)*numberOfFluidElementsPerLength
                if (xElementIdx == 1):
                    localNodes[localNodeIdx000] = (previousBlock-1)*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)+ \
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = (blockIdx-1)*numberOfFluidNodesPerBlock+numberOfNodesXi-1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                else:
                    localNodes[localNodeIdx000] = (blockIdx-1)*numberOfFluidNodesPerBlock+(xElementIdx-2)*(numberOfNodesXi-1)+(numberOfNodesXi-2)+1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSquareElements*(numberOfNodesXi-1))+\
                                                 (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+numberOfNodesXi-1
                localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                if(uInterpolation == QUADRATIC):
                    localNodes[1] = localNodes[localNodeIdx100] - 1
                    localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[5] = localNodes[4] + 1
                    localNodes[7] = localNodes[localNodeIdx110] - 1
                    localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
                    localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
                    localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
                    localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
                    localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
                    localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
                    localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
                    localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
                    localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
                    localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
                    localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
                    localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
                    localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
                    localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
                linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                               localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                if (debug):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                    if (uInterpolation == QUADRATIC):
                        print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                fluidPElements.NodesSet(elementNumber,linearNodes)
                fluidUElements.NodesSet(elementNumber,localNodes)
        previousBlock = blockIdx
    #Handle the square block
    if (numberOfSquareElements==1):
        elementNumber = elementNumber + 1
        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx110] = numberOfFluidNodesPerBlock+\
                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        if(uInterpolation == QUADRATIC):
            localNodes[1] = localNodes[localNodeIdx100] - 1
            localNodes[3] = localNodes[localNodeIdx000] - 1
            localNodes[4] = localNodes[localNodeIdx100] + 1
            localNodes[5] = localNodes[localNodeIdx110] - 1
            localNodes[7] = localNodes[localNodeIdx010] - 1
        localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
        linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                       localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
        if (uInterpolation == QUADRATIC):
            localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
            localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
            localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
            localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
            localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
            localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
            localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
            localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
            localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
            localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
            localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
            localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
            localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
            localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
        if (debug):
            print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                  (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
            if (uInterpolation == QUADRATIC):
                print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                                
        fluidPElements.NodesSet(elementNumber,linearNodes)
        fluidUElements.NodesSet(elementNumber,localNodes)
    else:
        for yElementIdx in range(1,numberOfSquareElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = 4*numberOfFluidElementsPerBlock+xElementIdx+(yElementIdx-1)*numberOfSquareElements+\
                                (zElementIdx-1)*numberOfFluidElementsPerLength
                if (yElementIdx == 1):
                    if (xElementIdx == 1):
                        #Bottom-left
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 3*numberOfFluidNodesPerBlock+numberOfArmElements*(numberOfNodesXi-1)*\
                                                      numberOfSquareElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 3*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1)
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Bottom-right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)-(numberOfNodesXi-2)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                        elif(uInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                    else:
                        #Bottom
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                        elif(uInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                elif (yElementIdx == numberOfSquareElements):
                    if (xElementIdx == 1):
                        #Top-left
                        localNodes[localNodeIdx000] = 2*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 2*numberOfFluidNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] + 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Top-right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfFluidNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                    else:
                        #Top
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfFluidNodesPerBlock-(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]-(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                else:
                    if (xElementIdx == 1):
                        #Left
                        localNodes[localNodeIdx000] = 3*numberOfFluidNodesPerBlock-(yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]-(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1) 
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Right
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx100] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                    else:
                        #Middle
                        localNodes[localNodeIdx000] = 4*numberOfFluidNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(uInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                               localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                if (uInterpolation == QUADRATIC):
                    localNodes[9] = localNodes[0]+numberOfFluidNodesPerLength
                    localNodes[10] = localNodes[1]+numberOfFluidNodesPerLength
                    localNodes[11] = localNodes[2]+numberOfFluidNodesPerLength
                    localNodes[12] = localNodes[3]+numberOfFluidNodesPerLength
                    localNodes[13] = localNodes[4]+numberOfFluidNodesPerLength
                    localNodes[14] = localNodes[5]+numberOfFluidNodesPerLength
                    localNodes[15] = localNodes[6]+numberOfFluidNodesPerLength
                    localNodes[16] = localNodes[7]+numberOfFluidNodesPerLength
                    localNodes[17] = localNodes[8]+numberOfFluidNodesPerLength
                    localNodes[19] = localNodes[10]+numberOfFluidNodesPerLength
                    localNodes[21] = localNodes[12]+numberOfFluidNodesPerLength
                    localNodes[22] = localNodes[13]+numberOfFluidNodesPerLength
                    localNodes[23] = localNodes[14]+numberOfFluidNodesPerLength
                    localNodes[25] = localNodes[16]+numberOfFluidNodesPerLength
                if (debug):
                    print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                    if (uInterpolation == QUADRATIC):
                        print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                        print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                
                fluidPElements.NodesSet(elementNumber,linearNodes)
                fluidUElements.NodesSet(elementNumber,localNodes)

fluidUElements.CreateFinish()
fluidPElements.CreateFinish()

if (debug):
    print('  Solid Elements:')

solidUElements = iron.MeshElements()
solidUElements.CreateStart(solidMesh,1,uBasis)
solidPElements = iron.MeshElements()
solidPElements.CreateStart(solidMesh,2,pBasis)

elementNumber = 0
for wallElementIdx in range(1,numberOfWallElements+1):
    for lengthElementIdx in range(1,numberOfLengthElements+1):
        for circumfrentialElementIdx in range(1,numberOfSolidCircumfrentialElements+1):
            localNodes = [0]*numberOfLocalNodes
            elementNumber = elementNumber + 1
            localNodes[localNodeIdx000] = (circumfrentialElementIdx-1)*(numberOfNodesXi-1) + 1 + \
                                          (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1) + \
                                          (wallElementIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
            if circumfrentialElementIdx == numberOfSolidCircumfrentialElements:
                localNodes[localNodeIdx100] = 1 + (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1) + \
                                              (wallElementIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
            else:
                localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + (numberOfNodesXi-1)
            localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localNodes[localNodeIdx001] = localNodes[localNodeIdx000] + numberOfSolidNodesPerWall*(numberOfNodesXi-1)
            localNodes[localNodeIdx101] = localNodes[localNodeIdx100] + numberOfSolidNodesPerWall*(numberOfNodesXi-1)
            localNodes[localNodeIdx011] = localNodes[localNodeIdx010] + numberOfSolidNodesPerWall*(numberOfNodesXi-1)
            localNodes[localNodeIdx111] = localNodes[localNodeIdx110] + numberOfSolidNodesPerWall*(numberOfNodesXi-1)
            if (uInterpolation == QUADRATIC):
                localNodes[1] = localNodes[localNodeIdx000] + 1
                localNodes[3] = localNodes[localNodeIdx000] + numberOfSolidCircumfrentialNodes
                localNodes[4] = localNodes[3] + 1
                localNodes[5] = localNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes
                localNodes[7] = localNodes[localNodeIdx010] + 1
                localNodes[9] = localNodes[0]+numberOfSolidNodesPerWall
                localNodes[10] = localNodes[1]+numberOfSolidNodesPerWall
                localNodes[11] = localNodes[2]+numberOfSolidNodesPerWall
                localNodes[12] = localNodes[3]+numberOfSolidNodesPerWall
                localNodes[13] = localNodes[4]+numberOfSolidNodesPerWall
                localNodes[14] = localNodes[5]+numberOfSolidNodesPerWall
                localNodes[15] = localNodes[6]+numberOfSolidNodesPerWall
                localNodes[16] = localNodes[7]+numberOfSolidNodesPerWall
                localNodes[17] = localNodes[8]+numberOfSolidNodesPerWall
                localNodes[19] = localNodes[10]+numberOfSolidNodesPerWall
                localNodes[21] = localNodes[12]+numberOfSolidNodesPerWall
                localNodes[22] = localNodes[13]+numberOfSolidNodesPerWall
                localNodes[23] = localNodes[14]+numberOfSolidNodesPerWall
                localNodes[25] = localNodes[16]+numberOfSolidNodesPerWall
            linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                           localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
            if (debug):
                print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                if (uInterpolation == QUADRATIC):
                    print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                        (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
            solidPElements.NodesSet(elementNumber,linearNodes)
            solidUElements.NodesSet(elementNumber,localNodes)

solidUElements.CreateFinish()
solidPElements.CreateFinish()

fluidMesh.CreateFinish()
solidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')    

#================================================================================================================================
#  Interface
#================================================================================================================================

if (progressDiagnostics):
    print('Interface ...')
    
# Create an interface between the two meshes
interface = iron.Interface()
interface.CreateStart(interfaceUserNumber,iron.WorldRegion)
interface.LabelSet('Interface')
# Add in the two meshes
solidMeshIndex = interface.MeshAdd(solidMesh)
fluidMeshIndex = interface.MeshAdd(fluidMesh)
interface.CoordinateSystemSet(interfaceCoordinateSystem)
interface.CreateFinish()

if (progressDiagnostics):
    print('Interface ... Done')
 
#================================================================================================================================
#  Interface Mesh
#================================================================================================================================

if (progressDiagnostics):
    print('Interface Mesh ...')
    
# Create an interface mesh
InterfaceNodes = iron.Nodes()
InterfaceNodes.CreateStartInterface(interface,numberOfInterfaceNodes)
InterfaceNodes.CreateFinish()

interfaceMesh = iron.Mesh()
interfaceMesh.CreateStartInterface(interfaceMeshUserNumber,interface,2)
interfaceMesh.NumberOfElementsSet(numberOfInterfaceElements)
interfaceMesh.NumberOfComponentsSet(1)

interfaceElements = iron.MeshElements()
interfaceElements.CreateStart(interfaceMesh,1,interfaceBasis)
        
if (debug):
    print('  Interface Elements:')
elementNumber = 0
for lengthElementIdx in range(1,numberOfLengthElements+1):
    for circumfrentialElementIdx in range(1,numberOfSolidCircumfrentialElements + 1):
        localNodes = [0]*numberOfLocalInterfaceNodes
        elementNumber = elementNumber + 1
        localNodes[localNodeIdx000] = (circumfrentialElementIdx-1)*(numberOfNodesXi-1) + 1 + \
                                      (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
        if circumfrentialElementIdx == numberOfSolidCircumfrentialElements:
            localNodes[localNodeIdx100] = 1 + (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
        else:
            localNodes[localNodeIdx100] = localNodes[localNodeIdx000] + (numberOfNodesXi-1)
        localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
        localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
        if (uInterpolation == QUADRATIC):
            localNodes[1] = localNodes[localNodeIdx000] + 1
            localNodes[3] = localNodes[localNodeIdx000] + numberOfSolidCircumfrentialNodes
            localNodes[4] = localNodes[3] + 1
            localNodes[5] = localNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes
            localNodes[7] = localNodes[localNodeIdx010] + 1
        linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110]]
        if (debug):
            print('    Element %8d; Nodes: %8d, %8d, %8d, %8d' % \
                  (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3]))
        if (uInterpolation == QUADRATIC):
            print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                  (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
        interfaceElements.NodesSet(elementNumber,localNodes)

interfaceElements.CreateFinish()

interfaceMesh.CreateFinish()

if (progressDiagnostics):
    print('Interface Mesh ... Done')
    

#================================================================================================================================
#  Mesh Connectivity
#================================================================================================================================

if (progressDiagnostics):
    print('Interface Mesh Connectivity ...')

# Couple the interface meshes
interfaceMeshConnectivity = iron.InterfaceMeshConnectivity()
interfaceMeshConnectivity.CreateStart(interface,interfaceMesh)
interfaceMeshConnectivity.BasisSet(interfaceBasis)
    
interfaceElementNumber = 0
interfaceNodes = [0]*(numberOfInterfaceNodes)
solidNodes = [0]*(numberOfInterfaceNodes)
fluidNodes = [0]*(numberOfInterfaceNodes)
localInterfaceNodes = [0]*numberOfLocalInterfaceNodes
localSolidNodes = [0]*numberOfLocalInterfaceNodes
localFluidNodes = [0]*numberOfLocalInterfaceNodes
previousBlock=4
for lengthElementIdx in range(1,numberOfLengthElements+1):
    for blockIdx in range(1,5):
        for circumfrentialElementIdx in range(1,numberOfSquareElements+1):
            interfaceElementNumber = interfaceElementNumber + 1
            if (debug):
                print('  Interface Element %8d:' % (interfaceElementNumber))        
            solidElementNumber = circumfrentialElementIdx+(blockIdx-1)*numberOfSquareElements+\
                                 (lengthElementIdx-1)*numberOfSolidCircumfrentialElements
            fluidElementNumber = circumfrentialElementIdx+(blockIdx-1)*numberOfFluidElementsPerBlock+\
                                 (lengthElementIdx-1)*numberOfFluidElementsPerLength
            # Map interface elements
            interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,solidMeshIndex,solidElementNumber)
            interfaceMeshConnectivity.ElementNumberSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber)
            if (debug):
                print('    Solid Element %8d; Fluid Element %8d' % (solidElementNumber,fluidElementNumber))        
            localInterfaceNodes[localNodeIdx000] = (circumfrentialElementIdx-1)*(numberOfNodesXi-1) + 1 + \
                                                   (blockIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                   (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localSolidNodes[localNodeIdx000] = (circumfrentialElementIdx-1)*(numberOfNodesXi-1) + 1 + \
                                               (blockIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                               (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes *(numberOfNodesXi-1)
            if (circumfrentialElementIdx==numberOfSquareElements and blockIdx==4):
                localInterfaceNodes[localNodeIdx100] = 1 + (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
                localSolidNodes[localNodeIdx100] = 1 + (lengthElementIdx-1)*numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            else:
                localInterfaceNodes[localNodeIdx100] = localInterfaceNodes[localNodeIdx000] + (numberOfNodesXi-1)
                localSolidNodes[localNodeIdx100] = localSolidNodes[localNodeIdx000] + (numberOfNodesXi-1)
            if (circumfrentialElementIdx == 1):
                localFluidNodes[localNodeIdx000] = (previousBlock-1)*numberOfFluidNodesPerBlock+\
                                                   numberOfSquareElements*(numberOfNodesXi-1)+ \
                                                   (lengthElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localFluidNodes[localNodeIdx100] = (blockIdx-1)*numberOfFluidNodesPerBlock+numberOfNodesXi-1+\
                                                   (lengthElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
            else:
                localFluidNodes[localNodeIdx000] = (blockIdx-1)*numberOfFluidNodesPerBlock+\
                                                   (circumfrentialElementIdx-2)*(numberOfNodesXi-1)+(numberOfNodesXi-2)+1+\
                                                   (lengthElementIdx-1)*numberOfFluidNodesPerLength*(numberOfNodesXi-1)
                localFluidNodes[localNodeIdx100] = localFluidNodes[localNodeIdx000]+numberOfNodesXi-1
            localInterfaceNodes[localNodeIdx010] = localInterfaceNodes[localNodeIdx000]+numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localSolidNodes[localNodeIdx010] = localSolidNodes[localNodeIdx000]+numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localFluidNodes[localNodeIdx010] = localFluidNodes[localNodeIdx000]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
            localInterfaceNodes[localNodeIdx110] = localInterfaceNodes[localNodeIdx100]+numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localSolidNodes[localNodeIdx110] = localSolidNodes[localNodeIdx100]+numberOfSolidCircumfrentialNodes*(numberOfNodesXi-1)
            localFluidNodes[localNodeIdx110] = localFluidNodes[localNodeIdx100]+numberOfFluidNodesPerLength*(numberOfNodesXi-1)
            if (uInterpolation == QUADRATIC):
                localInterfaceNodes[1] = localInterfaceNodes[localNodeIdx000]+1
                localInterfaceNodes[3] = localInterfaceNodes[localNodeIdx000]+numberOfSolidCircumfrentialNodes
                localInterfaceNodes[4] = localInterfaceNodes[3]+1
                localInterfaceNodes[5] = localInterfaceNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes
                localInterfaceNodes[7] = localInterfaceNodes[localNodeIdx010]+1
                localSolidNodes[1] = localSolidNodes[localNodeIdx000]+1
                localSolidNodes[3] = localSolidNodes[localNodeIdx000]+numberOfSolidCircumfrentialNodes
                localSolidNodes[4] = localSolidNodes[3]+1
                localSolidNodes[5] = localSolidNodes[localNodeIdx100] + numberOfSolidCircumfrentialNodes
                localSolidNodes[7] = localSolidNodes[localNodeIdx010]+1
                localFluidNodes[1] = localFluidNodes[localNodeIdx100]-1
                localFluidNodes[3] = localFluidNodes[localNodeIdx000]+numberOfFluidNodesPerLength
                localFluidNodes[4] = localFluidNodes[1]+numberOfFluidNodesPerLength
                localFluidNodes[5] = localFluidNodes[4]+1
                localFluidNodes[7] = localFluidNodes[localNodeIdx110]-1
            # Map interface xi
            for localNodeIdx2 in range(0,numberOfNodesXi):
                for localNodeIdx1 in range(0,numberOfNodesXi):
                    localNodeIdx = localNodeIdx1+localNodeIdx2*numberOfNodesXi
                    xi1=float(localNodeIdx1)/float(numberOfNodesXi-1)
                    xi2=float(localNodeIdx2)/float(numberOfNodesXi-1)
                    solidXi = [xi1,xi2,0.0]
                    fluidXi = [xi1,0.0,xi2]
                    interfaceNodes[localInterfaceNodes[localNodeIdx]-1]=localInterfaceNodes[localNodeIdx]
                    solidNodes[localInterfaceNodes[localNodeIdx]-1]=localSolidNodes[localNodeIdx]
                    fluidNodes[localInterfaceNodes[localNodeIdx]-1]=localFluidNodes[localNodeIdx]
                    interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,solidMeshIndex,solidElementNumber,localNodeIdx+1,1,solidXi)
                    interfaceMeshConnectivity.ElementXiSet(interfaceElementNumber,fluidMeshIndex,fluidElementNumber,localNodeIdx+1,1,fluidXi)
                    if (debug):
                        print('    Local node    %8d:' % (localNodeIdx+1))        
                        print('      Interface node    %8d:' % (localInterfaceNodes[localNodeIdx]))        
                        print('      Solid node        %8d; Solid xi = [ %.2f, %.2f, %.2f ]' % \
                              (localSolidNodes[localNodeIdx],solidXi[0],solidXi[1],solidXi[2]))
                        print('      Fluid node        %8d; Fluid xi = [ %.2f, %.2f, %.2f ]' % \
                              (localFluidNodes[localNodeIdx],fluidXi[0],fluidXi[1],fluidXi[2]))
        previousBlock=blockIdx
# Map interface nodes
interfaceMeshConnectivity.NodeNumberSet(interfaceNodes,solidMeshIndex,solidNodes,fluidMeshIndex,fluidNodes)        

interfaceMeshConnectivity.CreateFinish()

if (progressDiagnostics):
    print('Interface Mesh Connectivity ... Done')

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')
    
# Create a decomposition for the fluid mesh
fluidDecomposition = iron.Decomposition()
fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
fluidDecomposition.CalculateFacesSet(True)
fluidDecomposition.CreateFinish()

# Create a decomposition for the solid mesh
solidDecomposition = iron.Decomposition()
solidDecomposition.CreateStart(solidDecompositionUserNumber,solidMesh)
solidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
solidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
solidDecomposition.CalculateFacesSet(True)
solidDecomposition.CreateFinish()

# Create a decomposition for the interface mesh
interfaceDecomposition = iron.Decomposition()
interfaceDecomposition.CreateStart(interfaceDecompositionUserNumber,interfaceMesh)
interfaceDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
interfaceDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
interfaceDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')
    
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')
    
# Start to create a default (geometric) field on the fluid region
fluidGeometricField = iron.Field()
fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
# Set the decomposition to use
fluidGeometricField.MeshDecompositionSet(fluidDecomposition)
# Set the scaling to use
fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGeometry')
# Set the domain to be used by the field components.
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Finish creating the second field
fluidGeometricField.CreateFinish()

# Start to create a default (geometric) field on the solid region
solidGeometricField = iron.Field()
solidGeometricField.CreateStart(solidGeometricFieldUserNumber,solidRegion)
# Set the decomposition to use
solidGeometricField.MeshDecompositionSet(solidDecomposition)
solidGeometricField.meshDecomposition = solidDecomposition
# Set the scaling to use
solidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
solidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidGeometry')
# Set the domain to be used by the field components.
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
solidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Finish creating the first field
solidGeometricField.CreateFinish()

# Start to create a default (geometric) field on the Interface
interfaceGeometricField = iron.Field()
interfaceGeometricField.CreateStartInterface(interfaceGeometricFieldUserNumber,interface)
# Set the decomposition to use
interfaceGeometricField.MeshDecompositionSet(interfaceDecomposition)
# Set the scaling to use
interfaceGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
interfaceGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'InterfaceGeometry')
# Set the domain to be used by the field components.
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
interfaceGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Finish creating the first field
interfaceGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')
    
if (progressDiagnostics):
    print('Geometric Parameters ...')

armSize = (1.0-squareSizeRatio)*pipeRadius
squareSize = pipeRadius-armSize

if (debug):
    print('  Fluid Nodes:')
for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yNodeIdx in range(1,numberOfArmElements*(numberOfNodesXi-1)+2):
            for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
                nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                             (zNodeIdx-1)*numberOfFluidNodesPerLength
                nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
                if (nodeDomain == computationalNodeNumber):
                    if (yNodeIdx == numberOfArmElements*(numberOfNodesXi-1)+1):
                        #On the square
                        if (blockIdx == 1):
                            xPosition = squareSize - xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                        elif (blockIdx == 2):
                            xPosition = -xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = (numberOfSquareElements*(numberOfNodesXi-1)-xNodeIdx)*squareSize/ \
                                        (numberOfSquareElements*(numberOfNodesXi-1))
                        elif (blockIdx == 3):
                            xPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize
                            yPosition = -xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                        elif (blockIdx == 4):
                            xPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = -(numberOfSquareElements*(numberOfNodesXi-1)-xNodeIdx)*squareSize/ \
                                        (numberOfSquareElements*(numberOfNodesXi-1))
                    else:
                        #In the arm
                        #Work out the r, theta position
                        theta = xNodeIdx*math.pi/(2*numberOfSquareElements*(numberOfNodesXi-1))+(blockIdx-1)*math.pi/2.0
                        fraction = 1.0/(abs(math.sin(theta))+abs(math.cos(theta)))
                        armRadius = armSize+squareSize*(1.0-fraction)
                        radius = (numberOfArmElements*(numberOfNodesXi-1)-yNodeIdx+1)*armRadius/ \
                                 (numberOfArmElements*(numberOfNodesXi-1))+squareSize*fraction
                        xPosition = radius*math.cos(theta)
                        yPosition = radius*math.sin(theta)
                    zPosition = float(zNodeIdx-1)*lengthSize/float(numberOfLengthElements*(numberOfNodesXi-1))
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition,yPosition,zPosition))

    #Now handle square
    for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = 4*numberOfFluidNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)+\
                         (zNodeIdx-1)*numberOfFluidNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = squareSize*((xNodeIdx-1)-(yNodeIdx-1))/(numberOfSquareElements*(numberOfNodesXi-1))
                yPosition = squareSize*((xNodeIdx-1)+(yNodeIdx-1))/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize
                zPosition = float(zNodeIdx-1)*lengthSize/float(numberOfLengthElements*(numberOfNodesXi-1))
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition,yPosition,zPosition))
                        
# Update fields            
fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

innerRadius = pipeRadius
outerRadius = pipeRadius+wallThickness
if (debug):
    print('  Solid Nodes:')
for wallNodeIdx in range(1,numberOfWallElements*(numberOfNodesXi-1)+2):
    for lengthNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
        for circumfrentialNodeIdx in range(1,numberOfSolidCircumfrentialElements*(numberOfNodesXi-1)+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfSolidCircumfrentialNodes + \
                (wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
            nodeDomain = solidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                radius = innerRadius + (outerRadius - innerRadius)*float(wallNodeIdx-1)/float(numberOfSolidWallNodes)
                theta = float(circumfrentialNodeIdx-1)/float(numberOfSolidCircumfrentialNodes)*2.0*math.pi
                x = radius*math.cos(theta)
                y = radius*math.sin(theta)
                z = float(lengthNodeIdx-1)*lengthSize/float(numberOfSolidLengthNodes-1)
                solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
                solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
                solidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Position         = [ %.2f, %.2f, %.2f ]' % (x,y,z))

# Update fields            
solidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
solidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (debug):
    print('  Interface Nodes:')
for lengthNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    for circumfrentialNodeIdx in range(1,numberOfSolidCircumfrentialElements*(numberOfNodesXi-1)+1):
        nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfSolidCircumfrentialNodes 
        nodeDomain = computationalNodeNumber
        if (nodeDomain == computationalNodeNumber):
            radius = innerRadius
            theta = float(circumfrentialNodeIdx-1)/float(numberOfSolidCircumfrentialNodes)*2.0*math.pi
            x = radius*math.cos(theta)
            y = radius*math.sin(theta)
            z = float(lengthNodeIdx-1)*lengthSize/float(numberOfSolidLengthNodes-1)
            interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,x)
            interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,y)
            interfaceGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,z)
            if (debug):
                print('      Node        %d:' % (nodeNumber))
                print('         Position         = [ %.2f, %.2f, %.2f ]' % (x,y,z))


# Update fields            
interfaceGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
interfaceGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

# Export results
fluidFields = iron.Fields()
fluidFields.CreateRegion(fluidRegion)
fluidFields.NodesExport("FSITubeFluid","FORTRAN")
fluidFields.ElementsExport("FSITubeFluid","FORTRAN")
fluidFields.Finalise()
solidFields = iron.Fields()
solidFields.CreateRegion(solidRegion)
solidFields.NodesExport("FSITubeSolid","FORTRAN")
solidFields.ElementsExport("FSITubeSolid","FORTRAN")
solidFields.Finalise()
interfaceFields = iron.Fields()
interfaceFields.CreateInterface(interface)
interfaceFields.NodesExport("FSITubeInterface","FORTRAN")
interfaceFields.ElementsExport("FSITubeInterface","FORTRAN")
interfaceFields.Finalise()

#================================================================================================================================
#  Equations Set
#================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

# Create the equations set for the fluid region - Navier-Stokes
fluidEquationsSetField = iron.Field()
fluidEquationsSet = iron.EquationsSet()
if RBS:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.ALE_RBS_NAVIER_STOKES]
else:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.ALE_NAVIER_STOKES]
fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                              fluidEquationsSetSpecification,fluidEquationsSetFieldUserNumber,
                              fluidEquationsSetField)
fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
fluidEquationsSet.CreateFinish()

if RBS:
    # Set max CFL number (default 1.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0E20)
    # Set time increment (default 0.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,3,timeStep)
    # Set stabilisation type (default 1.0 = RBS)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,4,1.0)


# Create the equations set for the solid region 
solidEquationsSetField = iron.Field()
solidEquationsSet = iron.EquationsSet()
solidEquationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
                                  iron.EquationsSetTypes.FINITE_ELASTICITY,
                                  iron.EquationsSetSubtypes.MOONEY_RIVLIN]
solidEquationsSet.CreateStart(solidEquationsSetUserNumber,solidRegion,solidGeometricField,
                              solidEquationsSetSpecification,solidEquationsSetFieldUserNumber,
                              solidEquationsSetField)
solidEquationsSet.OutputTypeSet(solidEquationsSetOutputType)
solidEquationsSet.CreateFinish()

# Create the equations set for the moving mesh
movingMeshEquationsSetField = iron.Field()
movingMeshEquationsSet = iron.EquationsSet()
movingMeshEquationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
                                       iron.EquationsSetTypes.LAPLACE_EQUATION,
                                       iron.EquationsSetSubtypes.MOVING_MESH_LAPLACE]
movingMeshEquationsSet.CreateStart(movingMeshEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                                   movingMeshEquationsSetSpecification,movingMeshEquationsSetFieldUserNumber,
                                   movingMeshEquationsSetField)
movingMeshEquationsSet.OutputTypeSet(movingMeshEquationsSetOutputType)
movingMeshEquationsSet.CreateFinish()

if (progressDiagnostics):
    print('Equations Sets ... Done')


#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

# Create the equations set dependent field variables for dynamic Navier-Stokes
fluidDependentField = iron.Field()
fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber,fluidDependentField)
fluidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,4):
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
# Finish the equations set dependent field variables
fluidEquationsSet.DependentCreateFinish()

# Initialise the fluid dependent field
for componentIdx in range(1,4):
    fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)
# Initialise pressure component
fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,fluidPInit)

# Create the equations set dependent field variables for the solid equations set
solidDependentField = iron.Field()
solidEquationsSet.DependentCreateStart(solidDependentFieldUserNumber,solidDependentField)
solidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidDependent')
solidDependentField.VariableLabelSet(iron.FieldVariableTypes.DELUDELN,'SolidTraction')
for componentIdx in range(1,4):
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
solidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
solidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
solidDependentField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
solidEquationsSet.DependentCreateFinish()

# Initialise the solid dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
for componentIdx in range(1,4):
    solidGeometricField.ParametersToFieldParametersComponentCopy(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,\
                                                                 componentIdx,solidDependentField,iron.FieldVariableTypes.U,
                                                                 iron.FieldParameterSetTypes.VALUES,componentIdx)
solidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,solidPInit)

# Create the equations set dependent field variables for moving mesh
movingMeshDependentField = iron.Field()
movingMeshEquationsSet.DependentCreateStart(movingMeshDependentFieldUserNumber,movingMeshDependentField)
movingMeshDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,numberOfDimensions+1):
    movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    movingMeshDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
# Finish the equations set dependent field variables
movingMeshEquationsSet.DependentCreateFinish()

# Initialise dependent field moving mesh
for componentIdx in range(1,numberOfDimensions+1):
    movingMeshDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                         componentIdx,0.0)

fluidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
solidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
movingMeshDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
solidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
movingMeshDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')
     
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

# Create the equations set materials field variables for dynamic Navier-Stokes
fluidMaterialsField = iron.Field()
fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber,fluidMaterialsField)
# Finish the equations set materials field variables
fluidEquationsSet.MaterialsCreateFinish()
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,fluidDynamicViscosity)
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,fluidDensity)

# Create the solid materials field
solidMaterialsField = iron.Field()
solidEquationsSet.MaterialsCreateStart(solidMaterialsFieldUserNumber,solidMaterialsField)
solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.U,'SolidMaterials')
solidMaterialsField.VariableLabelSet(iron.FieldVariableTypes.V,'SolidDensity')
solidEquationsSet.MaterialsCreateFinish()
# Set Mooney-Rivlin constants c10 and c01 respectively
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,mooneyRivlin1)
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,mooneyRivlin2)
solidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES,1,solidDensity)

# Create the equations set materials field variables for moving mesh
movingMeshMaterialsField = iron.Field()
movingMeshEquationsSet.MaterialsCreateStart(movingMeshMaterialsFieldUserNumber,movingMeshMaterialsField)
# Finish the equations set materials field variables
movingMeshEquationsSet.MaterialsCreateFinish()

movingMeshMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,\
                                                     movingMeshKParameter)

if (progressDiagnostics):
    print('Materials Fields ... Done')
    
#================================================================================================================================
# Independent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Independent Fields ...')

# Create fluid mesh velocity independent field 
fluidIndependentField = iron.Field()
fluidEquationsSet.IndependentCreateStart(fluidIndependentFieldUserNumber,fluidIndependentField)
fluidIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidIndependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,numberOfDimensions+1):
    fluidIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
# Finish the equations set independent field variables
fluidEquationsSet.IndependentCreateFinish()
  
# Create the moving mesh independent field 
movingMeshIndependentField = iron.Field()
movingMeshEquationsSet.IndependentCreateStart(movingMeshIndependentFieldUserNumber,movingMeshIndependentField)
movingMeshIndependentField.VariableLabelSet(iron.FieldVariableTypes.U,'MovingMeshIndependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,numberOfDimensions+1):
    movingMeshIndependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)    
# Finish the equations set independent field variables
movingMeshEquationsSet.IndependentCreateFinish()

# Initialise independent field moving mesh
movingMeshIndependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,movingMeshKParameter)

if (progressDiagnostics):
    print('Independent Fields ... Done')

#================================================================================================================================
#  Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

# Fluid equations 
fluidEquations = iron.Equations()
fluidEquationsSet.EquationsCreateStart(fluidEquations)
fluidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
fluidEquations.outputType = fluidEquationsOutputType
fluidEquationsSet.EquationsCreateFinish()

# Solid equations
solidEquations = iron.Equations()
solidEquationsSet.EquationsCreateStart(solidEquations)
solidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
solidEquations.outputType = solidEquationsOutputType
solidEquationsSet.EquationsCreateFinish()

# Moving mesh equations
movingMeshEquations = iron.Equations()
movingMeshEquationsSet.EquationsCreateStart(movingMeshEquations)
movingMeshEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
movingMeshEquations.outputType = movingMeshEquationsOutputType
movingMeshEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

#================================================================================================================================
#  CellML
#================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

# Create CellML equations for the temporal boundary conditions
bcCellML = iron.CellML()
bcCellML.CreateStart(bcCellMLUserNumber,fluidRegion)
bcCellMLIdx = bcCellML.ModelImport("input/poiseuilleinlet.cellml")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/pipeRadius")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/length")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/dynamicViscosity")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/A")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/B")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/C")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/x")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/y")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletx")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inlety")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletz")
bcCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
bcCellML.FieldMapsCreateStart()
# Map geometric field to x and y
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/x",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/y",iron.FieldParameterSetTypes.VALUES)
# Map fluid velocity to ensure dependent field isn't cleared when the velocities are copied back
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/inletz",iron.FieldParameterSetTypes.VALUES)
# Map inletx, inlety and inletz to dependent field
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletz",iron.FieldParameterSetTypes.VALUES,
                                fluidDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES)
bcCellML.FieldMapsCreateFinish()

# Create the CellML models field
bcCellMLModelsField = iron.Field()
bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"BCModelMap")
bcCellML.ModelsFieldCreateFinish()

# Only evaluate BC on inlet nodes
bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0)
if (debug):
    print('  CellML Boundary Conditions:')
    print('    Inlet Model Set:')
for blockIdx in range(1,5):
    for yNodeIdx in range(2,numberOfArmElements*(numberOfNodesXi-1)+2):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                               1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
    for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        nodeNumber = 4*numberOfFluidNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debug):
                print('      Node        %d:' % (nodeNumber))

bcCellMLModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
bcCellMLModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                
# Create the CellML state field
bcCellMLStateField = iron.Field()
bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCState")
bcCellML.StateFieldCreateFinish()

# Create the CellML parameters field
bcCellMLParametersField = iron.Field()
bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"BCParameters")
bcCellML.ParametersFieldCreateFinish()

# Get the component numbers
pipeRadiusComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/pipeRadius")
lengthComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/length")
dynamicViscosityComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/dynamicViscosity")
AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/A")
BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/B")
CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/C")
# Set up the parameters field
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    pipeRadiusComponentNumber,pipeRadius)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    lengthComponentNumber,lengthSize)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    dynamicViscosityComponentNumber,fluidDynamicViscosity)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    AComponentNumber,A)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    BComponentNumber,B)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    CComponentNumber,C)

# Create the CELL intermediate field
bcCellMLIntermediateField = iron.Field()
bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCIntermediate")
bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

#================================================================================================================================
#  Interface Condition
#================================================================================================================================

if (progressDiagnostics):
    print('Interface Conditions ...')

# Create an interface condition between the two meshes
interfaceCondition = iron.InterfaceCondition()
interfaceCondition.CreateStart(interfaceConditionUserNumber,interface,interfaceGeometricField)
# Specify the method for the interface condition
interfaceCondition.MethodSet(iron.InterfaceConditionMethods.LAGRANGE_MULTIPLIERS)
# Specify the type of interface condition operator
interfaceCondition.OperatorSet(iron.InterfaceConditionOperators.SOLID_FLUID)
# Add in the dependent variables from the equations sets
interfaceCondition.DependentVariableAdd(solidMeshIndex,solidEquationsSet,iron.FieldVariableTypes.U)
interfaceCondition.DependentVariableAdd(fluidMeshIndex,fluidEquationsSet,iron.FieldVariableTypes.U)
# Set the label
interfaceCondition.LabelSet("FSI Interface Condition")
# Set the output type
interfaceCondition.OutputTypeSet(interfaceConditionOutputType)
# Finish creating the interface condition
interfaceCondition.CreateFinish()

if (progressDiagnostics):
    print('Interface Conditions ... Done')

if (progressDiagnostics):
    print('Interface Lagrange Field ...')
    
# Create the Lagrange multipliers field
interfaceLagrangeField = iron.Field()
interfaceCondition.LagrangeFieldCreateStart(interfaceLagrangeFieldUserNumber,interfaceLagrangeField)
interfaceLagrangeField.VariableLabelSet(iron.FieldVariableTypes.U,'InterfaceLagrange')
# Finish the Lagrange multipliers field
interfaceCondition.LagrangeFieldCreateFinish()

for componentIdx in range(1,numberOfDimensions+1):
    interfaceLagrangeField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)

interfaceLagrangeField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
interfaceLagrangeField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Interface Lagrange Field ... Done')

if (progressDiagnostics):
    print('Interface Equations ...')

# Create the interface condition equations
interfaceEquations = iron.InterfaceEquations()
interfaceCondition.EquationsCreateStart(interfaceEquations)
# Set the interface equations sparsity
interfaceEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
# Set the interface equations output
interfaceEquations.outputType = interfaceEquationsOutputType
# Finish creating the interface equations
interfaceCondition.EquationsCreateFinish()

if (progressDiagnostics):
    print('Interface Equations ... Done')

#================================================================================================================================
#  Problem
#================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a FSI problem
fsiProblem = iron.Problem()
if RBS:
    fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                               iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                               iron.ProblemSubtypes.FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE]
else:
    fsiProblemSpecification = [iron.ProblemClasses.MULTI_PHYSICS,
                               iron.ProblemTypes.FINITE_ELASTICITY_NAVIER_STOKES,
                               iron.ProblemSubtypes.FINITE_ELASTICITY_NAVIER_STOKES_ALE]
fsiProblem.CreateStart(fsiProblemUserNumber,fsiProblemSpecification)
fsiProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

#================================================================================================================================
#  Control Loop
#================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fsi problem control loop
fsiControlLoop = iron.ControlLoop()
fsiProblem.ControlLoopCreateStart()
fsiProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],fsiControlLoop)
fsiControlLoop.LabelSet('TimeLoop')
fsiControlLoop.TimesSet(startTime,stopTime,timeStep)
fsiControlLoop.TimeOutputSet(outputFrequency)
fsiProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = iron.Solver()
fsiDynamicSolver = iron.Solver()
fsiNonlinearSolver = iron.Solver()
fsiLinearSolver = iron.Solver()
movingMeshLinearSolver = iron.Solver()

fsiProblem.SolversCreateStart()
# Solvers for coupled FiniteElasticity NavierStokes problem
# Get the BC CellML solver
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
# Get the dynamic ALE solver
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fsiDynamicSolver)
fsiDynamicSolver.OutputTypeSet(fsiDynamicSolverOutputType)
fsiDynamicSolver.DynamicThetaSet(fsiDynamicSolverTheta)
# Get the dynamic nonlinear solver
fsiDynamicSolver.DynamicNonlinearSolverGet(fsiNonlinearSolver)
fsiNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
#fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
fsiNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
fsiNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
fsiNonlinearSolver.OutputTypeSet(fsiNonlinearSolverOutputType)
fsiNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
fsiNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
fsiNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
fsiNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
# Get the dynamic nonlinear linear solver
fsiNonlinearSolver.NewtonLinearSolverGet(fsiLinearSolver)
#fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
#fsiLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
#fsiLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
#fsiLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
#fsiLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
fsiLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
fsiLinearSolver.OutputTypeSet(fsiLinearSolverOutputType)
# Linear solver for moving mesh
fsiProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],3,movingMeshLinearSolver)
movingMeshLinearSolver.OutputTypeSet(movingMeshLinearSolverOutputType)
# Finish the creation of the problem solver
fsiProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

#================================================================================================================================
#  CellML Equations
#================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

# Create CellML equations and add BC equations to the solver
bcEquations = iron.CellMLEquations()
fsiProblem.CellMLEquationsCreateStart()
bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
fsiProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateStart()
# Get the fsi dynamic solver equations
fsiSolverEquations = iron.SolverEquations()
fsiDynamicSolver.SolverEquationsGet(fsiSolverEquations)
fsiSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
#fsiSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.FULL
fsiSolidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(solidEquationsSet)
fsiFluidEquationsSetIndex = fsiSolverEquations.EquationsSetAdd(fluidEquationsSet)
fsiInterfaceConditionIndex = fsiSolverEquations.InterfaceConditionAdd(interfaceCondition)
# Set the time dependence of the interface matrix to determine the interface matrix coefficient in the solver matrix
# (basicly position in big coupled matrix system)
interfaceEquations.MatrixTimeDependenceTypeSet(fsiSolidEquationsSetIndex,True, \
                                               [iron.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                iron.InterfaceMatricesTimeDependenceTypes.FIRST_ORDER_DYNAMIC])
interfaceEquations.MatrixTimeDependenceTypeSet(fsiFluidEquationsSetIndex,True, \
                                               [iron.InterfaceMatricesTimeDependenceTypes.STATIC,\
                                                iron.InterfaceMatricesTimeDependenceTypes.STATIC])
    
# Create the moving mesh solver equations
movingMeshSolverEquations = iron.SolverEquations()
# Get the linear moving mesh solver equations
movingMeshLinearSolver.SolverEquationsGet(movingMeshSolverEquations)
movingMeshSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
# Add in the equations set
movingMeshEquationsSetIndex = movingMeshSolverEquations.EquationsSetAdd(movingMeshEquationsSet)

# Finish the creation of the fsi problem solver equations
fsiProblem.SolverEquationsCreateFinish()


if (progressDiagnostics):
    print('Solver Equations ...')

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fluid boundary conditions
fsiBoundaryConditions = iron.BoundaryConditions()
fsiSolverEquations.BoundaryConditionsCreateStart(fsiBoundaryConditions)

if (debug):
    print('  Fluid Boundary Conditions:')
    print('    Inlet Boundary conditions:')
# Set inlet boundary conditions on the left hand edge
for blockIdx in range(1,5):
    for yNodeIdx in range(2,numberOfArmElements*(numberOfNodesXi-1)+2):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,3,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                 
for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
    for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        nodeNumber = 4*numberOfFluidNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,3,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (debug):
                print('      Node        %d:' % (nodeNumber))
                print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                 
if (debug):
    print('    Wall Boundary conditions:')
# Set no slip boundary conditions on the wall
for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    for blockIdx in range(1,5):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(zNodeIdx-1)*numberOfFluidNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Velocity         = [ ????, ????, %.2f ]' % (0.0))                 
if (debug):
    print('    No Pressure Boundary conditions:')

nodeNumber = 3*numberOfFluidNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)+ \
             numberOfLengthElements*(numberOfNodesXi-1)*numberOfFluidNodesPerLength
fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                              iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                              nodeNumber,4,iron.BoundaryConditionsTypes.PRESSURE,0.0)

# Set no pressure boundary conditions on the outlet
for blockIdx in range(1,5):
    for yElementIdx in range(2,numberOfArmElements+2):
        for xElementIdx in range(1,numberOfSquareElements+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xElementIdx*(numberOfNodesXi-1)+\
                         (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+ \
                         numberOfLengthElements*(numberOfNodesXi-1)*numberOfFluidNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
            if (nodeDomain == computationalNodeNumber):
                #fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                #                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                #                                nodeNumber,4,iron.BoundaryConditionsTypes.FIXED_OUTLET,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Pressure         = %.2f' % (0.0))                 
for yElementIdx in range(1,numberOfSquareElements):
    for xElementIdx in range(1,numberOfSquareElements):
        nodeNumber = 4*numberOfFluidNodesPerBlock+xElementIdx*(numberOfNodesXi-1)+numberOfSquareElements*(numberOfNodesXi-1)-1+\
                     (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSquareElements*(numberOfNodesXi-1)-1)+ \
                     numberOfLengthElements*(numberOfNodesXi-1)*numberOfFluidNodesPerLength
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            #fsiBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
            #                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
            #                                nodeNumber,4,iron.BoundaryConditionsTypes.FIXED_OUTLET,0.0)
            if (debug):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         = %.2f' % (0.0))                 

if (debug): 
    print('  Solid Boundary conditions:')
# Set solid boundary conditions
for wallNodeIdx in range (1,numberOfWallElements*(numberOfNodesXi-1)+2):
    #Set nodes on the axis to only slide on the axis.
    #Top y-axis node - fix in the x-direction
    nodeNumber = 1*numberOfSquareElements*(numberOfNodesXi-1)+\
                 (wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    if (nodeDomain == computationalNodeNumber):
        fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Displacement     = [ %.2f, ????, ???? ]' % (0.0))                         
    #Left x-axis node - fix in the y-direction
    nodeNumber = 2*numberOfSquareElements*(numberOfNodesXi-1)+\
                 (wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    if (nodeDomain == computationalNodeNumber):
        fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Displacement     = [ ????, %.2f, ???? ]' % (0.0))                         
    #Bottom y-axis node - fix in the x-direction
    nodeNumber = 3*numberOfSquareElements*(numberOfNodesXi-1)+\
                 (wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    if (nodeDomain == computationalNodeNumber):
        fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Displacement     = [ %.2f, ????, ???? ]' % (0.0))                         
    #Left x-axis node - fix in the y-direction
    nodeNumber = 4*numberOfSquareElements*(numberOfNodesXi-1)+\
                 (wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
    nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
    if (nodeDomain == computationalNodeNumber):
        fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Displacement     = [ ????, %.2f, ???? ]' % (0.0))                         
    #Fix all solid nodes at the begining of the tube in the z-direction    
    for circumfrentialNodeIdx in range (1,numberOfSolidCircumfrentialNodes+1):
        nodeNumber = circumfrentialNodeIdx+(wallNodeIdx-1)*numberOfSolidCircumfrentialNodes*numberOfSolidLengthNodes
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fsiBoundaryConditions.AddNode(solidDependentField,iron.FieldVariableTypes.U,1, \
                                          iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                          nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debug):
                print('      Node        %d:' % (nodeNumber))
                print('         Displacement     = [ ????, ????, %.2f ]' % (0.0))                 
                
if (debug):
    print('  Lagrange Boundary conditions:')
# Remove Lagrange multipliers where solid displacement and fluid velocity is zero
nodeNumber = 1*numberOfSquareElements*(numberOfNodesXi-1)
nodeDomain = computationalNodeNumber
if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Lagrange         = [ %.2f, ????, ???? ]' % (0.0))                 
nodeNumber = 2*numberOfSquareElements*(numberOfNodesXi-1)
nodeDomain = computationalNodeNumber
if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Lagrange         = [ ????, %.2f, ???? ]' % (0.0))                 
nodeNumber = 3*numberOfSquareElements*(numberOfNodesXi-1)
nodeDomain = computationalNodeNumber
if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Lagrange         = [ %.2f, ????, ???? ]' % (0.0))                 
nodeNumber = 4*numberOfSquareElements*(numberOfNodesXi-1)
nodeDomain = computationalNodeNumber
if (nodeDomain == computationalNodeNumber):
    fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                  iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                  nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (debug):
        print('      Node        %d:' % (nodeNumber))
        print('         Lagrange         = [ ????, %.2f, ???? ]' % (0.0))                 
for circumfrentialNodeIdx in range (1,numberOfSolidCircumfrentialNodes+1):
    nodeNumber = circumfrentialNodeIdx
    nodeDomain = computationalNodeNumber
    if (nodeDomain == computationalNodeNumber):
        fsiBoundaryConditions.SetNode(interfaceLagrangeField,iron.FieldVariableTypes.U,1, \
                                      iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                      nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (debug):
            print('      Node        %d:' % (nodeNumber))
            print('         Lagrange         = [ ????, ????, %.2f ]' % (0.0))                 

# Finish fsi boundary conditions
fsiSolverEquations.BoundaryConditionsCreateFinish()
        
# Start the creation of the moving mesh boundary conditions
movingMeshBoundaryConditions = iron.BoundaryConditions()
movingMeshSolverEquations.BoundaryConditionsCreateStart(movingMeshBoundaryConditions)
if (debug):
    print('  Moving Mesh Boundary Conditions:')
    print('    Fixed Wall Boundary conditions:')
# Inlet nodes
for blockIdx in range(1,5):
    for yNodeIdx in range(2,numberOfArmElements*(numberOfNodesXi-1)+2):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                     iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,
                                                     iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
    for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        nodeNumber = 4*numberOfFluidNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1,
                                                 iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,
                                                 iron.BoundaryConditionsTypes.FIXED_WALL,0.0)
            if (debug):
                print('      Node        %d:' % (nodeNumber))                
if (debug):
    print('    Moving Wall Boundary conditions:')
for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    for blockIdx in range(1,5):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfFluidNodesPerBlock+xNodeIdx+(zNodeIdx-1)*numberOfFluidNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                movingMeshBoundaryConditions.AddNode(movingMeshDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))

# Finish moving mesh boundary conditions
movingMeshSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#quit()

# Export results
fields = iron.Fields()
fields.CreateRegion(fluidRegion)
fields.NodesExport("FSITube","FORTRAN")
fields.ElementsExport("FSITube","FORTRAN")
fields.Finalise()

fsiLinearSolver.MumpsSetIcntl(14,20000)

# Create output directories
if not os.path.exists("output/Fluid"):
    os.makedirs("output/Fluid")
if not os.path.exists("output/Solid"):
    os.makedirs("output/Solid")
if not os.path.exists("output/Interface"):
    os.makedirs("output/Interface")

# Solve the problem
print('Solving problem...')
start = time.time()
fsiProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' %elapsed)
print('Problem solved!')
print('#')

