%% About defineModel.mlx
% This file defines the MATLAB interface to the library |Model|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for &lt;SHAPE&gt;, &lt;DIRECTION&gt;, etc. For more
% information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') Define MATLAB Interface for C++ Library>.



%% Setup. Do not edit this section.
function libDef = defineModel()
libDef = clibgen.LibraryDefinition("ModelData.xml");
%% OutputFolder and Libraries 
libDef.OutputFolder = "D:\VCBM_drug_model-github";
libDef.Libraries = "";

%% C++ enumeration |CellType| with MATLAB name |clib.Model.CellType| 
addEnumeration(libDef, "CellType", "int32",...
    [...
      "Cancer",...  % 1
      "Dead",...  % 3
      "Healthy",...  % 4
      "Empty",...  % 5
      "PSC",...  % 51
    ],...
    "MATLABName", "clib.Model.CellType", ...
    "Description", "clib.Model.CellType    Representation of C++ enumeration CellType."); % Modify help description values as needed.

%% C++ class |CellState| with MATLAB name |clib.Model.CellState| 
CellStateDefinition = addClass(libDef, "CellState", "MATLABName", "clib.Model.CellState", ...
    "Description", "clib.Model.CellState    Representation of C++ class CellState."); % Modify help description values as needed.

%% C++ class constructor for C++ class |CellState| 
% C++ Signature: CellState::CellState()
CellStateConstructor1Definition = addConstructor(CellStateDefinition, ...
    "CellState::CellState()", ...
    "Description", "clib.Model.CellState.CellState    Constructor of C++ class CellState."); % Modify help description values as needed.
validate(CellStateConstructor1Definition);

%% C++ class constructor for C++ class |CellState| 
% C++ Signature: CellState::CellState(CellState const & input1)
CellStateConstructor2Definition = addConstructor(CellStateDefinition, ...
    "CellState::CellState(CellState const & input1)", ...
    "Description", "clib.Model.CellState.CellState    Constructor of C++ class CellState."); % Modify help description values as needed.
defineArgument(CellStateConstructor2Definition, "input1", "clib.Model.CellState", "input");
validate(CellStateConstructor2Definition);

%% C++ class public data member |type| for C++ class |CellState| 
% C++ Signature: CellType CellState::type
addProperty(CellStateDefinition, "type", "clib.Model.CellType", ...
    "Description", "clib.Model.CellType    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class public data member |age| for C++ class |CellState| 
% C++ Signature: int CellState::age
addProperty(CellStateDefinition, "age", "int32", ...
    "Description", "int32    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class public data member |spring_length| for C++ class |CellState| 
% C++ Signature: double CellState::spring_length
addProperty(CellStateDefinition, "spring_length", "double", ...
    "Description", "double    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class public data member |sibling| for C++ class |CellState| 
% C++ Signature: Cell * CellState::sibling
addProperty(CellStateDefinition, "sibling", "clib.Model.Cell", 1, ... % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
   "Description", "clib.Model.Cell    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class public data member |X| for C++ class |CellState| 
% C++ Signature: double CellState::X
addProperty(CellStateDefinition, "X", "double", ...
    "Description", "double    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class public data member |Y| for C++ class |CellState| 
% C++ Signature: double CellState::Y
addProperty(CellStateDefinition, "Y", "double", ...
    "Description", "double    Data member of C++ class CellState."); % Modify help description values as needed.

%% C++ class |Params| with MATLAB name |clib.Model.Params| 
ParamsDefinition = addClass(libDef, "Params", "MATLABName", "clib.Model.Params", ...
    "Description", "clib.Model.Params    Representation of C++ class Params."); % Modify help description values as needed.

%% C++ class constructor for C++ class |Params| 
% C++ Signature: Params::Params(double p_0,double p_psc,double dmax,int gage,int page,double EC50)
ParamsConstructor1Definition = addConstructor(ParamsDefinition, ...
    "Params::Params(double p_0,double p_psc,double dmax,int gage,int page,double EC50)", ...
    "Description", "clib.Model.Params.Params    Constructor of C++ class Params."); % Modify help description values as needed.
defineArgument(ParamsConstructor1Definition, "p_0", "double");
defineArgument(ParamsConstructor1Definition, "p_psc", "double");
defineArgument(ParamsConstructor1Definition, "dmax", "double");
defineArgument(ParamsConstructor1Definition, "gage", "int32");
defineArgument(ParamsConstructor1Definition, "page", "int32");
defineArgument(ParamsConstructor1Definition, "EC50", "double");
validate(ParamsConstructor1Definition);

%% C++ class method |RandomDouble| for C++ class |Params| 
% C++ Signature: double Params::RandomDouble()
RandomDoubleDefinition = addMethod(ParamsDefinition, ...
    "double Params::RandomDouble()", ...
    "Description", "clib.Model.Params.RandomDouble    Method of C++ class Params."); % Modify help description values as needed.
defineOutput(RandomDoubleDefinition, "RetVal", "double");
validate(RandomDoubleDefinition);

%% C++ class method |WithProbability| for C++ class |Params| 
% C++ Signature: bool Params::WithProbability(double prob)
WithProbabilityDefinition = addMethod(ParamsDefinition, ...
    "bool Params::WithProbability(double prob)", ...
    "Description", "clib.Model.Params.WithProbability    Method of C++ class Params."); % Modify help description values as needed.
defineArgument(WithProbabilityDefinition, "prob", "double");
defineOutput(WithProbabilityDefinition, "RetVal", "logical");
validate(WithProbabilityDefinition);

%% C++ class public data member |gage| for C++ class |Params| 
% C++ Signature: int Params::gage
addProperty(ParamsDefinition, "gage", "int32", ...
    "Description", "int32    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class public data member |page| for C++ class |Params| 
% C++ Signature: int Params::page
addProperty(ParamsDefinition, "page", "int32", ...
    "Description", "int32    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class public data member |p_0| for C++ class |Params| 
% C++ Signature: double Params::p_0
addProperty(ParamsDefinition, "p_0", "double", ...
    "Description", "double    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class public data member |dmax| for C++ class |Params| 
% C++ Signature: double Params::dmax
addProperty(ParamsDefinition, "dmax", "double", ...
    "Description", "double    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class public data member |p_psc| for C++ class |Params| 
% C++ Signature: double Params::p_psc
addProperty(ParamsDefinition, "p_psc", "double", ...
    "Description", "double    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class public data member |EC50| for C++ class |Params| 
% C++ Signature: double Params::EC50
addProperty(ParamsDefinition, "EC50", "double", ...
    "Description", "double    Data member of C++ class Params."); % Modify help description values as needed.

%% C++ class |Cell| with MATLAB name |clib.Model.Cell| 
CellDefinition = addClass(libDef, "Cell", "MATLABName", "clib.Model.Cell", ...
    "Description", "clib.Model.Cell    Representation of C++ class Cell."); % Modify help description values as needed.

%% C++ class constructor for C++ class |Cell| 
% C++ Signature: Cell::Cell(double X,double Y,double spring_length,Cell * sibling,CellType type,int age)
CellConstructor1Definition = addConstructor(CellDefinition, ...
   "Cell::Cell(double X,double Y,double spring_length,Cell * sibling,CellType type,int age)", ...
   "Description", "clib.Model.Cell.Cell    Constructor of C++ class Cell."); % Modify help description values as needed.
defineArgument(CellConstructor1Definition, "X", "double");
defineArgument(CellConstructor1Definition, "Y", "double");
defineArgument(CellConstructor1Definition, "spring_length", "double");
defineArgument(CellConstructor1Definition, "sibling", "clib.Model.Cell", "input", 1); % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
defineArgument(CellConstructor1Definition, "type", "clib.Model.CellType");
defineArgument(CellConstructor1Definition, "age", "int32");
validate(CellConstructor1Definition);

%% C++ class constructor for C++ class |Cell| 
% C++ Signature: Cell::Cell(Cell * cell)
CellConstructor2Definition = addConstructor(CellDefinition, ...
   "Cell::Cell(Cell * cell)", ...
   "Description", "clib.Model.Cell.Cell    Constructor of C++ class Cell."); % Modify help description values as needed.
defineArgument(CellConstructor2Definition, "cell", "clib.Model.Cell", "input", 1); % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
validate(CellConstructor2Definition);

%% C++ class method |Renew| for C++ class |Cell| 
% C++ Signature: void Cell::Renew()
RenewDefinition = addMethod(CellDefinition, ...
    "void Cell::Renew()", ...
    "Description", "clib.Model.Cell.Renew    Method of C++ class Cell."); % Modify help description values as needed.
validate(RenewDefinition);

%% C++ class method |UpdateState| for C++ class |Cell| 
% C++ Signature: void Cell::UpdateState()
UpdateStateDefinition = addMethod(CellDefinition, ...
    "void Cell::UpdateState()", ...
    "Description", "clib.Model.Cell.UpdateState    Method of C++ class Cell."); % Modify help description values as needed.
validate(UpdateStateDefinition);

%% C++ class method |Infect| for C++ class |Cell| 
% C++ Signature: void Cell::Infect()
InfectDefinition = addMethod(CellDefinition, ...
    "void Cell::Infect()", ...
    "Description", "clib.Model.Cell.Infect    Method of C++ class Cell."); % Modify help description values as needed.
validate(InfectDefinition);

%% C++ class method |clearNeighbours| for C++ class |Cell| 
% C++ Signature: void Cell::clearNeighbours()
clearNeighboursDefinition = addMethod(CellDefinition, ...
    "void Cell::clearNeighbours()", ...
    "Description", "clib.Model.Cell.clearNeighbours    Method of C++ class Cell."); % Modify help description values as needed.
validate(clearNeighboursDefinition);

%% C++ class method |OnBoundary| for C++ class |Cell| 
% C++ Signature: bool Cell::OnBoundary()
OnBoundaryDefinition = addMethod(CellDefinition, ...
    "bool Cell::OnBoundary()", ...
    "Description", "clib.Model.Cell.OnBoundary    Method of C++ class Cell."); % Modify help description values as needed.
defineOutput(OnBoundaryDefinition, "RetVal", "logical");
validate(OnBoundaryDefinition);

%% C++ class method |TooCrowded| for C++ class |Cell| 
% C++ Signature: bool Cell::TooCrowded()
TooCrowdedDefinition = addMethod(CellDefinition, ...
    "bool Cell::TooCrowded()", ...
    "Description", "clib.Model.Cell.TooCrowded    Method of C++ class Cell."); % Modify help description values as needed.
defineOutput(TooCrowdedDefinition, "RetVal", "logical");
validate(TooCrowdedDefinition);

%% C++ class method |DistanceSquaredTo| for C++ class |Cell| 
% C++ Signature: double Cell::DistanceSquaredTo(Cell * cell)
DistanceSquaredToDefinition = addMethod(CellDefinition, ...
   "double Cell::DistanceSquaredTo(Cell * cell)", ...
   "Description", "clib.Model.Cell.DistanceSquaredTo    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(DistanceSquaredToDefinition, "cell", "clib.Model.Cell", "input", 1); % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
defineOutput(DistanceSquaredToDefinition, "RetVal", "double");
validate(DistanceSquaredToDefinition);

%% C++ class method |DistanceSquaredFromCentre| for C++ class |Cell| 
% C++ Signature: double Cell::DistanceSquaredFromCentre()
DistanceSquaredFromCentreDefinition = addMethod(CellDefinition, ...
    "double Cell::DistanceSquaredFromCentre()", ...
    "Description", "clib.Model.Cell.DistanceSquaredFromCentre    Method of C++ class Cell."); % Modify help description values as needed.
defineOutput(DistanceSquaredFromCentreDefinition, "RetVal", "double");
validate(DistanceSquaredFromCentreDefinition);

%% C++ class method |Necrotic| for C++ class |Cell| 
% C++ Signature: bool Cell::Necrotic(double distanceToBoundary,Params * parameters)
NecroticDefinition = addMethod(CellDefinition, ...
   "bool Cell::Necrotic(double distanceToBoundary,Params * parameters)", ...
   "Description", "clib.Model.Cell.Necrotic    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(NecroticDefinition, "distanceToBoundary", "double");
defineArgument(NecroticDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
defineOutput(NecroticDefinition, "RetVal", "logical");
validate(NecroticDefinition);

%% C++ class method |TooYoung| for C++ class |Cell| 
% C++ Signature: bool Cell::TooYoung(Params * parameters)
TooYoungDefinition = addMethod(CellDefinition, ...
   "bool Cell::TooYoung(Params * parameters)", ...
   "Description", "clib.Model.Cell.TooYoung    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(TooYoungDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
defineOutput(TooYoungDefinition, "RetVal", "logical");
validate(TooYoungDefinition);

%% C++ class method |LengthenSpring| for C++ class |Cell| 
% C++ Signature: void Cell::LengthenSpring(Params * parameters)
LengthenSpringDefinition = addMethod(CellDefinition, ...
   "void Cell::LengthenSpring(Params * parameters)", ...
   "Description", "clib.Model.Cell.LengthenSpring    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(LengthenSpringDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
validate(LengthenSpringDefinition);

%% C++ class method |DrugInducedDeath| for C++ class |Cell| 
% C++ Signature: bool Cell::DrugInducedDeath(Params * parameters,double * drugConcentration,int gridRadius)
DrugInducedDeathDefinition = addMethod(CellDefinition, ...
   "bool Cell::DrugInducedDeath(Params * parameters,double * drugConcentration,int gridRadius)", ...
   "Description", "clib.Model.Cell.DrugInducedDeath    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(DrugInducedDeathDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
defineArgument(DrugInducedDeathDefinition, "drugConcentration", "clib.array.Model.Double", "input", 1); % '<MLTYPE>' can be clib.array.Model.Double, or double
defineArgument(DrugInducedDeathDefinition, "gridRadius", "int32");
defineOutput(DrugInducedDeathDefinition, "RetVal", "logical");
validate(DrugInducedDeathDefinition);

%% C++ class method |Die| for C++ class |Cell| 
% C++ Signature: void Cell::Die()
DieDefinition = addMethod(CellDefinition, ...
    "void Cell::Die()", ...
    "Description", "clib.Model.Cell.Die    Method of C++ class Cell."); % Modify help description values as needed.
validate(DieDefinition);

%% C++ class method |Disintegrate| for C++ class |Cell| 
% C++ Signature: void Cell::Disintegrate()
DisintegrateDefinition = addMethod(CellDefinition, ...
    "void Cell::Disintegrate()", ...
    "Description", "clib.Model.Cell.Disintegrate    Method of C++ class Cell."); % Modify help description values as needed.
validate(DisintegrateDefinition);

%% C++ class method |PossiblyPSCInfectNeighbour| for C++ class |Cell| 
% C++ Signature: void Cell::PossiblyPSCInfectNeighbour(Params * parameters)
PossiblyPSCInfectNeighbourDefinition = addMethod(CellDefinition, ...
   "void Cell::PossiblyPSCInfectNeighbour(Params * parameters)", ...
   "Description", "clib.Model.Cell.PossiblyPSCInfectNeighbour    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(PossiblyPSCInfectNeighbourDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
validate(PossiblyPSCInfectNeighbourDefinition);

%% C++ class method |Move| for C++ class |Cell| 
% C++ Signature: void Cell::Move()
MoveDefinition = addMethod(CellDefinition, ...
    "void Cell::Move()", ...
    "Description", "clib.Model.Cell.Move    Method of C++ class Cell."); % Modify help description values as needed.
validate(MoveDefinition);

%% C++ class method |Proliferate| for C++ class |Cell| 
% C++ Signature: Cell * Cell::Proliferate(Params * parameters)
ProliferateDefinition = addMethod(CellDefinition, ...
   "Cell * Cell::Proliferate(Params * parameters)", ...
   "Description", "clib.Model.Cell.Proliferate    Method of C++ class Cell."); % Modify help description values as needed.
defineArgument(ProliferateDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
defineOutput(ProliferateDefinition, "RetVal", "clib.Model.Cell", 1);
validate(ProliferateDefinition);

%% C++ class public data member |currentState| for C++ class |Cell| 
% C++ Signature: CellState Cell::currentState
addProperty(CellDefinition, "currentState", "clib.Model.CellState", ...
    "Description", "clib.Model.CellState    Data member of C++ class Cell."); % Modify help description values as needed.

%% C++ class public data member |newState| for C++ class |Cell| 
% C++ Signature: CellState Cell::newState
addProperty(CellDefinition, "newState", "clib.Model.CellState", ...
    "Description", "clib.Model.CellState    Data member of C++ class Cell."); % Modify help description values as needed.

%% C++ class |Pancreas| with MATLAB name |clib.Model.Pancreas| 
PancreasDefinition = addClass(libDef, "Pancreas", "MATLABName", "clib.Model.Pancreas", ...
    "Description", "clib.Model.Pancreas    Representation of C++ class Pancreas."); % Modify help description values as needed.

%% C++ class constructor for C++ class |Pancreas| 
% C++ Signature: Pancreas::Pancreas(Pancreas * existing,Params * parameters)
PancreasConstructor1Definition = addConstructor(PancreasDefinition, ...
   "Pancreas::Pancreas(Pancreas * existing,Params * parameters)", ...
   "Description", "clib.Model.Pancreas.Pancreas    Constructor of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(PancreasConstructor1Definition, "existing", "clib.Model.Pancreas", "input", 1); % '<MLTYPE>' can be clib.Model.Pancreas, or clib.array.Model.Pancreas
defineArgument(PancreasConstructor1Definition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
validate(PancreasConstructor1Definition);

%% C++ class constructor for C++ class |Pancreas| 
% C++ Signature: Pancreas::Pancreas(Params * parameters)
PancreasConstructor2Definition = addConstructor(PancreasDefinition, ...
   "Pancreas::Pancreas(Params * parameters)", ...
   "Description", "clib.Model.Pancreas.Pancreas    Constructor of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(PancreasConstructor2Definition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
validate(PancreasConstructor2Definition);

%% C++ class method |CreateNewParticle| for C++ class |Pancreas| 
% C++ Signature: Pancreas * Pancreas::CreateNewParticle(Params * parameters)
CreateNewParticleDefinition = addMethod(PancreasDefinition, ...
   "Pancreas * Pancreas::CreateNewParticle(Params * parameters)", ...
   "Description", "clib.Model.Pancreas.CreateNewParticle    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(CreateNewParticleDefinition, "parameters", "clib.Model.Params", "input", 1); % '<MLTYPE>' can be clib.Model.Params, or clib.array.Model.Params
defineOutput(CreateNewParticleDefinition, "RetVal", "clib.Model.Pancreas", 1);
validate(CreateNewParticleDefinition);

%% C++ class method |InjectPoint| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::InjectPoint(int x,int y,double amount)
InjectPointDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::InjectPoint(int x,int y,double amount)", ...
    "Description", "clib.Model.Pancreas.InjectPoint    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(InjectPointDefinition, "x", "int32");
defineArgument(InjectPointDefinition, "y", "int32");
defineArgument(InjectPointDefinition, "amount", "double");
validate(InjectPointDefinition);

%% C++ class method |InjectFibre| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::InjectFibre(int x,int y,double amount)
InjectFibreDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::InjectFibre(int x,int y,double amount)", ...
    "Description", "clib.Model.Pancreas.InjectFibre    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(InjectFibreDefinition, "x", "int32");
defineArgument(InjectFibreDefinition, "y", "int32");
defineArgument(InjectFibreDefinition, "amount", "double");
validate(InjectFibreDefinition);

%% C++ class method |LoadCellsCoordinates| for C++ class |Pancreas| 
% C++ Signature: double * Pancreas::LoadCellsCoordinates()
LoadCellsCoordinatesDefinition = addMethod(PancreasDefinition, ...
   "double * Pancreas::LoadCellsCoordinates()", ...
   "Description", "clib.Model.Pancreas.LoadCellsCoordinates    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(LoadCellsCoordinatesDefinition, "RetVal", "double", 1);
validate(LoadCellsCoordinatesDefinition);

%% C++ class method |DetermineNeighbours| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::DetermineNeighbours()
DetermineNeighboursDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::DetermineNeighbours()", ...
    "Description", "clib.Model.Pancreas.DetermineNeighbours    Method of C++ class Pancreas."); % Modify help description values as needed.
validate(DetermineNeighboursDefinition);

%% C++ class method |HealthyCellsBeyondRadius| for C++ class |Pancreas| 
% C++ Signature: bool Pancreas::HealthyCellsBeyondRadius(double radius)
HealthyCellsBeyondRadiusDefinition = addMethod(PancreasDefinition, ...
    "bool Pancreas::HealthyCellsBeyondRadius(double radius)", ...
    "Description", "clib.Model.Pancreas.HealthyCellsBeyondRadius    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(HealthyCellsBeyondRadiusDefinition, "radius", "double");
defineOutput(HealthyCellsBeyondRadiusDefinition, "RetVal", "logical");
validate(HealthyCellsBeyondRadiusDefinition);

%% C++ class method |AddNewCell| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::AddNewCell(Cell * new_cell)
AddNewCellDefinition = addMethod(PancreasDefinition, ...
   "void Pancreas::AddNewCell(Cell * new_cell)", ...
   "Description", "clib.Model.Pancreas.AddNewCell    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(AddNewCellDefinition, "new_cell", "clib.Model.Cell", "input", 1); % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
validate(AddNewCellDefinition);

%% C++ class method |AddMoreTissue| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::AddMoreTissue(double moving_rim,double max_tumour_radius)
AddMoreTissueDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::AddMoreTissue(double moving_rim,double max_tumour_radius)", ...
    "Description", "clib.Model.Pancreas.AddMoreTissue    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(AddMoreTissueDefinition, "moving_rim", "double");
defineArgument(AddMoreTissueDefinition, "max_tumour_radius", "double");
validate(AddMoreTissueDefinition);

%% C++ class method |MoreTissueAddedIfNecessary| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::MoreTissueAddedIfNecessary()
MoreTissueAddedIfNecessaryDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::MoreTissueAddedIfNecessary()", ...
    "Description", "clib.Model.Pancreas.MoreTissueAddedIfNecessary    Method of C++ class Pancreas."); % Modify help description values as needed.
validate(MoreTissueAddedIfNecessaryDefinition);

%% C++ class method |DetermineBoundaryCells| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::DetermineBoundaryCells()
DetermineBoundaryCellsDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::DetermineBoundaryCells()", ...
    "Description", "clib.Model.Pancreas.DetermineBoundaryCells    Method of C++ class Pancreas."); % Modify help description values as needed.
validate(DetermineBoundaryCellsDefinition);

%% C++ class method |TumourRadius| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::TumourRadius()
TumourRadiusDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::TumourRadius()", ...
    "Description", "clib.Model.Pancreas.TumourRadius    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(TumourRadiusDefinition, "RetVal", "double");
validate(TumourRadiusDefinition);

%% C++ class method |DistanceToLine| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::DistanceToLine(Cell * cell)
DistanceToLineDefinition = addMethod(PancreasDefinition, ...
   "double Pancreas::DistanceToLine(Cell * cell)", ...
   "Description", "clib.Model.Pancreas.DistanceToLine    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(DistanceToLineDefinition, "cell", "clib.Model.Cell", "input", 1); % '<MLTYPE>' can be clib.Model.Cell, or clib.array.Model.Cell
defineOutput(DistanceToLineDefinition, "RetVal", "double");
validate(DistanceToLineDefinition);

%% C++ class method |TumourVolume| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::TumourVolume()
TumourVolumeDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::TumourVolume()", ...
    "Description", "clib.Model.Pancreas.TumourVolume    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(TumourVolumeDefinition, "RetVal", "double");
validate(TumourVolumeDefinition);

%% C++ class method |getPscRatio| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::getPscRatio()
getPscRatioDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::getPscRatio()", ...
    "Description", "clib.Model.Pancreas.getPscRatio    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(getPscRatioDefinition, "RetVal", "double");
validate(getPscRatioDefinition);

%% C++ class method |CreateInitialTumour| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::CreateInitialTumour()
CreateInitialTumourDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::CreateInitialTumour()", ...
    "Description", "clib.Model.Pancreas.CreateInitialTumour    Method of C++ class Pancreas."); % Modify help description values as needed.
validate(CreateInitialTumourDefinition);

%% C++ class method |SimulateOneHour| for C++ class |Pancreas| 
% C++ Signature: void Pancreas::SimulateOneHour()
SimulateOneHourDefinition = addMethod(PancreasDefinition, ...
    "void Pancreas::SimulateOneHour()", ...
    "Description", "clib.Model.Pancreas.SimulateOneHour    Method of C++ class Pancreas."); % Modify help description values as needed.
validate(SimulateOneHourDefinition);

%% C++ class method |SimulateOneDay| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::SimulateOneDay(int day)
SimulateOneDayDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::SimulateOneDay(int day)", ...
    "Description", "clib.Model.Pancreas.SimulateOneDay    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(SimulateOneDayDefinition, "day", "int32");
defineOutput(SimulateOneDayDefinition, "RetVal", "double");
validate(SimulateOneDayDefinition);

%% C++ class method |ReturnTotalNumberTumourCells| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnTotalNumberTumourCells()
ReturnTotalNumberTumourCellsDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnTotalNumberTumourCells()", ...
    "Description", "clib.Model.Pancreas.ReturnTotalNumberTumourCells    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnTotalNumberTumourCellsDefinition, "RetVal", "int32");
validate(ReturnTotalNumberTumourCellsDefinition);

%% C++ class method |ReturnTotalNumberDeadCells| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnTotalNumberDeadCells()
ReturnTotalNumberDeadCellsDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnTotalNumberDeadCells()", ...
    "Description", "clib.Model.Pancreas.ReturnTotalNumberDeadCells    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnTotalNumberDeadCellsDefinition, "RetVal", "int32");
validate(ReturnTotalNumberDeadCellsDefinition);

%% C++ class method |ReturnTotalNumberPSCCells| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnTotalNumberPSCCells()
ReturnTotalNumberPSCCellsDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnTotalNumberPSCCells()", ...
    "Description", "clib.Model.Pancreas.ReturnTotalNumberPSCCells    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnTotalNumberPSCCellsDefinition, "RetVal", "int32");
validate(ReturnTotalNumberPSCCellsDefinition);

%% C++ class method |ReturnTotalNumberHealthyCells| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnTotalNumberHealthyCells()
ReturnTotalNumberHealthyCellsDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnTotalNumberHealthyCells()", ...
    "Description", "clib.Model.Pancreas.ReturnTotalNumberHealthyCells    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnTotalNumberHealthyCellsDefinition, "RetVal", "int32");
validate(ReturnTotalNumberHealthyCellsDefinition);

%% C++ class method |TestingPoissonDist| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::TestingPoissonDist()
TestingPoissonDistDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::TestingPoissonDist()", ...
    "Description", "clib.Model.Pancreas.TestingPoissonDist    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(TestingPoissonDistDefinition, "RetVal", "int32");
validate(TestingPoissonDistDefinition);

%% C++ class method |ReturnDrugConcentrationDomain| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::ReturnDrugConcentrationDomain()
ReturnDrugConcentrationDomainDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::ReturnDrugConcentrationDomain()", ...
    "Description", "clib.Model.Pancreas.ReturnDrugConcentrationDomain    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnDrugConcentrationDomainDefinition, "RetVal", "double");
validate(ReturnDrugConcentrationDomainDefinition);

%% C++ class method |ReturnDrugConcentrationinFibre| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::ReturnDrugConcentrationinFibre()
ReturnDrugConcentrationinFibreDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::ReturnDrugConcentrationinFibre()", ...
    "Description", "clib.Model.Pancreas.ReturnDrugConcentrationinFibre    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnDrugConcentrationinFibreDefinition, "RetVal", "double");
validate(ReturnDrugConcentrationinFibreDefinition);

%% C++ class method |ReturnDrugConcentrationAout| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::ReturnDrugConcentrationAout()
ReturnDrugConcentrationAoutDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::ReturnDrugConcentrationAout()", ...
    "Description", "clib.Model.Pancreas.ReturnDrugConcentrationAout    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnDrugConcentrationAoutDefinition, "RetVal", "double");
validate(ReturnDrugConcentrationAoutDefinition);

%% C++ class method |ReturnCellPositions| for C++ class |Pancreas| 
% C++ Signature: double Pancreas::ReturnCellPositions(int index)
ReturnCellPositionsDefinition = addMethod(PancreasDefinition, ...
    "double Pancreas::ReturnCellPositions(int index)", ...
    "Description", "clib.Model.Pancreas.ReturnCellPositions    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(ReturnCellPositionsDefinition, "index", "int32");
defineOutput(ReturnCellPositionsDefinition, "RetVal", "double");
validate(ReturnCellPositionsDefinition);

%% C++ class method |ReturnNumberCells| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnNumberCells()
ReturnNumberCellsDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnNumberCells()", ...
    "Description", "clib.Model.Pancreas.ReturnNumberCells    Method of C++ class Pancreas."); % Modify help description values as needed.
defineOutput(ReturnNumberCellsDefinition, "RetVal", "int32");
validate(ReturnNumberCellsDefinition);

%% C++ class method |ReturnCellType| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::ReturnCellType(int index)
ReturnCellTypeDefinition = addMethod(PancreasDefinition, ...
    "int Pancreas::ReturnCellType(int index)", ...
    "Description", "clib.Model.Pancreas.ReturnCellType    Method of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(ReturnCellTypeDefinition, "index", "int32");
defineOutput(ReturnCellTypeDefinition, "RetVal", "int32");
validate(ReturnCellTypeDefinition);

%% C++ class constructor for C++ class |Pancreas| 
% C++ Signature: Pancreas::Pancreas(Pancreas const & input1)
PancreasConstructor3Definition = addConstructor(PancreasDefinition, ...
    "Pancreas::Pancreas(Pancreas const & input1)", ...
    "Description", "clib.Model.Pancreas.Pancreas    Constructor of C++ class Pancreas."); % Modify help description values as needed.
defineArgument(PancreasConstructor3Definition, "input1", "clib.Model.Pancreas", "input");
validate(PancreasConstructor3Definition);

%% C++ class public data member |fibreX| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::fibreX
addProperty(PancreasDefinition, "fibreX", "int32", ...
    "Description", "int32    Data member of C++ class Pancreas."); % Modify help description values as needed.

%% C++ class public data member |fibreY| for C++ class |Pancreas| 
% C++ Signature: int Pancreas::fibreY
addProperty(PancreasDefinition, "fibreY", "int32", ...
    "Description", "int32    Data member of C++ class Pancreas."); % Modify help description values as needed.

%% C++ class public data member |timeSinceInjection| for C++ class |Pancreas| 
% C++ Signature: long Pancreas::timeSinceInjection
addProperty(PancreasDefinition, "timeSinceInjection", "int32", ...
    "Description", "int32    Data member of C++ class Pancreas."); % Modify help description values as needed.

%% C++ class public data member |drugConcentration| for C++ class |Pancreas| 
% C++ Signature: double [19891] Pancreas::drugConcentration
addProperty(PancreasDefinition, "drugConcentration", "clib.array.Model.Double", [19891], ... % '<MLTYPE>' can be clib.array.Model.Double, or double
    "Description", "clib.array.Model.Double    Data member of C++ class Pancreas."); % Modify help description values as needed.

%% C++ class public data member |fibreConcentration| for C++ class |Pancreas| 
% C++ Signature: double * Pancreas::fibreConcentration
addProperty(PancreasDefinition, "fibreConcentration", "clib.array.Model.Double", 1, ... % '<MLTYPE>' can be clib.array.Model.Double, or double
   "Description", "clib.array.Model.Double    Data member of C++ class Pancreas."); % Modify help description values as needed.

%% C++ function |SeedAndGrowToStartVolumeM| with MATLAB name |clib.Model.SeedAndGrowToStartVolumeM|
% C++ Signature: Pancreas * SeedAndGrowToStartVolumeM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume)
SeedAndGrowToStartVolumeMDefinition = addFunction(libDef, ...
   "Pancreas * SeedAndGrowToStartVolumeM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume)", ...
   "MATLABName", "clib.Model.SeedAndGrowToStartVolumeM", ...
   "Description", "clib.Model.SeedAndGrowToStartVolumeM    Representation of C++ function SeedAndGrowToStartVolumeM."); % Modify help description values as needed.
defineArgument(SeedAndGrowToStartVolumeMDefinition, "p0", "double");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "psc", "double");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "dmax", "int32");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "gage", "int32");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "page", "int32");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "EC50", "double");
defineArgument(SeedAndGrowToStartVolumeMDefinition, "startVolume", "double");
defineOutput(SeedAndGrowToStartVolumeMDefinition, "RetVal", "clib.Model.Pancreas", 1);
validate(SeedAndGrowToStartVolumeMDefinition);

%% C++ function |SimulateWholeExperimentM| with MATLAB name |clib.Model.SimulateWholeExperimentM|
% C++ Signature: void SimulateWholeExperimentM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume,int timeSteps,double [] volumes)
SimulateWholeExperimentMDefinition = addFunction(libDef, ...
   "void SimulateWholeExperimentM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume,int timeSteps,double [] volumes)", ...
   "MATLABName", "clib.Model.SimulateWholeExperimentM", ...
   "Description", "clib.Model.SimulateWholeExperimentM    Representation of C++ function SimulateWholeExperimentM."); % Modify help description values as needed.
defineArgument(SimulateWholeExperimentMDefinition, "p0", "double");
defineArgument(SimulateWholeExperimentMDefinition, "psc", "double");
defineArgument(SimulateWholeExperimentMDefinition, "dmax", "int32");
defineArgument(SimulateWholeExperimentMDefinition, "gage", "int32");
defineArgument(SimulateWholeExperimentMDefinition, "page", "int32");
defineArgument(SimulateWholeExperimentMDefinition, "EC50", "double");
defineArgument(SimulateWholeExperimentMDefinition, "startVolume", "double");
defineArgument(SimulateWholeExperimentMDefinition, "timeSteps", "int32");
defineArgument(SimulateWholeExperimentMDefinition, "volumes", "clib.array.Model.Double", "input", 1); % '<MLTYPE>' can be clib.array.Model.Double, or double
validate(SimulateWholeExperimentMDefinition);

%% C++ function |PerformMultipleRunsM| with MATLAB name |clib.Model.PerformMultipleRunsM|
% C++ Signature: void PerformMultipleRunsM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume,int timeSteps,int iterations,double [] volumes)
PerformMultipleRunsMDefinition = addFunction(libDef, ...
   "void PerformMultipleRunsM(double p0,double psc,int dmax,int gage,int page,double EC50,double startVolume,int timeSteps,int iterations,double [] volumes)", ...
   "MATLABName", "clib.Model.PerformMultipleRunsM", ...
   "Description", "clib.Model.PerformMultipleRunsM    Representation of C++ function PerformMultipleRunsM."); % Modify help description values as needed.
defineArgument(PerformMultipleRunsMDefinition, "p0", "double");
defineArgument(PerformMultipleRunsMDefinition, "psc", "double");
defineArgument(PerformMultipleRunsMDefinition, "dmax", "int32");
defineArgument(PerformMultipleRunsMDefinition, "gage", "int32");
defineArgument(PerformMultipleRunsMDefinition, "page", "int32");
defineArgument(PerformMultipleRunsMDefinition, "EC50", "double");
defineArgument(PerformMultipleRunsMDefinition, "startVolume", "double");
defineArgument(PerformMultipleRunsMDefinition, "timeSteps", "int32");
defineArgument(PerformMultipleRunsMDefinition, "iterations", "int32");
defineArgument(PerformMultipleRunsMDefinition, "volumes", "clib.array.Model.Double", "input", 1); % '<MLTYPE>' can be clib.array.Model.Double, or double
validate(PerformMultipleRunsMDefinition);

%% C++ function |CreateNewParticle| with MATLAB name |clib.Model.CreateNewParticle|
% C++ Signature: Pancreas * CreateNewParticle(double p0,double psc,int dmax,int gage,int page,double EC50,Pancreas * pancreas)
CreateNewParticleDefinition = addFunction(libDef, ...
   "Pancreas * CreateNewParticle(double p0,double psc,int dmax,int gage,int page,double EC50,Pancreas * pancreas)", ...
   "MATLABName", "clib.Model.CreateNewParticle", ...
   "Description", "clib.Model.CreateNewParticle    Representation of C++ function CreateNewParticle."); % Modify help description values as needed.
defineArgument(CreateNewParticleDefinition, "p0", "double");
defineArgument(CreateNewParticleDefinition, "psc", "double");
defineArgument(CreateNewParticleDefinition, "dmax", "int32");
defineArgument(CreateNewParticleDefinition, "gage", "int32");
defineArgument(CreateNewParticleDefinition, "page", "int32");
defineArgument(CreateNewParticleDefinition, "EC50", "double");
defineArgument(CreateNewParticleDefinition, "pancreas", "clib.Model.Pancreas", "input", 1); % '<MLTYPE>' can be clib.Model.Pancreas, or clib.array.Model.Pancreas
defineOutput(CreateNewParticleDefinition, "RetVal", "clib.Model.Pancreas", 1);
validate(CreateNewParticleDefinition);

%% Validate the library definition
validate(libDef);

end
