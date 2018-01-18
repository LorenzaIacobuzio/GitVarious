#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "FilterMuonElectronThreeTrackVertex.cc.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class FilterMuonElectronThreeTrackVertex.cc
/// \Brief
/// Short description of your Analyzer
/// \EndBrief
///
/// \Detailed
/// Detailed description of your Analyzer\n\n
/// For examples of working Analyzer you can have a look at the examples in the Examples/ directory:\n
/// LKrPhotonMC \n
/// Pi0Reconstruction \n
/// SkimmingNoStrawMuonsTracks \n
/// Or the framework analyzers that can be found in the Analyzers/ directories: \n
/// CedarMCTester \n
/// VertexCDA \n
/// \n
/// All the following classes are available for you to use. Check their documentation for more information:\n
/// NA62Analysis::manip::TermManip \n
/// NA62Analysis::Analyzer \n
/// NA62Analysis::CounterHandler \n
/// NA62Analysis::DetectorAcceptance \n
/// NA62Analysis::EventFraction \n
/// NA62Analysis::MCSimple \n
/// NA62Analysis::NeuralNetwork \n
/// NA62Analysis::ParticleInterface \n
/// NA62Analysis::ParticleTree \n
/// NA62Analysis::StringBalancedTable \n
/// NA62Analysis::StringTable \n
/// NA62Analysis::UserMethods \n
/// NA62Analysis::Verbose \n
///
/// You might also be interested in checking the documentation for the following classes. However you should not
/// in general have to create instances of these. If necessary a pointer to the existing instance is usually
/// available or provided by specific methods.\n
/// NA62Analysis::Core::IOHisto \n
/// NA62Analysis::Core::IOTree \n
/// NA62Analysis::Core::IOHandler \n
/// NA62Analysis::Core::HistoHandler \n
/// NA62Analysis::Core::HistoHandler::Iterator \n
/// NA62Analysis::Core::PrimitiveReader \n
///
/// \EndDetailed

FilterMuonElectronThreeTrackVertex.cc::FilterMuonElectronThreeTrackVertex.cc(Core::BaseAnalysis *ba) : Analyzer(ba, "FilterMuonElectronThreeTrackVertex.cc")
{
	/// \MemberDescr
	/// \param ba : parent BaseAnalysis
	///
	/// Specify the trees you want to use and the event class corresponding\n
	/// Don't try to load MCTruth tree (RUN_0 or Event). Use the MCTruthEvent in Process function instead. Problems when opening twice the same tree.\n
	/// Example with RecoEvent\n
	///	\code
	/// 	RequestTree(new TRecoGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
	///		RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Digis");
	/// \endcode
	/// The first form can be used if the detector name is present in the class
	/// name.\n
	/// Example with MC Event\n
	///	\code
	/// 	RequestTree(new TGigaTrackerEvent);
	///		RequestTree("GigaTracker", new TGigaTrackerEvent);
	/// \endcode
	/// Example with generic tree\n
	///	\code
	///		RequestTree<MyClass>("MyTree", "BranchName", "MyClass", new MyClass);
	///		RequestTree("MyTree", "BranchName", "MyClass", new MyClass);
	/// 	RequestTree("MyTree", "BranchName", "int", new int);
	/// \endcode
	/// Requesting Trigger data\n
	///	\code
	///		RequestL0Data();
	///		RequestL1Data();
	///		RequestL2Data();
	/// \endcode
	//// \n
	/// Initialize DetectorAcceptance if needed\n
	/// use of global instance\n
	///	\code
	///		fDetectorAcceptanceInstance = GetDetectorAcceptanceInstance();
	/// \endcode
	/// use of local instance\n
	///	\code
	///		fDetectorAcceptanceInstance = new DetectorAcceptance("./NA62.root");
	/// \endcode
	/// \n
	/// Call one of: \n
	///	\code
	/// 	AddParam("paramName", &variableName, defaultValue);
	/// \endcode
	/// for each parameter of the analyzer. These parameters can be set when starting the FW from the command line with the -p option.\n
	/// paramName is the name of the parameter in the command line\n
	/// variableName is the name of the variable that should be declared in the definition of the class\n
	/// defaultValue is the default value if not specified in the command line\n
	/// The allowed types for parameters are the following: bool, int, long, float, double, char, string, TString\n
	/// \n
	/// A primitive file can be read in parallel to the event file. You can request to read primitives for a sub-detector
	/// with\n
	/// \code
	/// AddPrimitiveReader("detName", true/false);
	/// \endcode
	/// The boolean flag indicates if the primitives should be time sorted or kept in the original order they arrived from
	/// the detector.\n
	/// You can set the L0 matching window with \n
	/// \code
	/// SetL0MatchingWindowWidth("detName", val);
	/// SetL0MatchingWindowWidth("detName", timestamp, finetime);
	/// \endcode
	/// where val is a value in nanoseconds.
	/// \EndMemberDescr



}

void FilterMuonElectronThreeTrackVertex.cc::InitOutput(){
	/// \MemberDescr
	/// Register the output variables of the analyzer.\n
	/// Call: \n
	///	\code
	/// 	RegisterOutput("outputName", &variableName)
	/// \endcode
	/// for each variable that should be in the output of the Analyzer\n
	/// The name of the analyzer will be prepended to the outputName (to avoid collisions with other analyzers)\n
	/// variableName should be the name of a variable declared in the definition of the class\n
	/// \n
	/// To create a new TTree in the output file, call: \n
	///	\code
	/// 	void OpenNewTree("TTreeName", "TTreeTitle");
	/// \endcode
	/// TTreeName is the name of the TTree (will be used to refer to this TTree later)\n
	/// TTreeTitle is the title of the TTree\n
	/// \n
	/// To add a branch to the newly created TTree, call: \n
	///	\code
	/// 	void AddBranch<VariableType>("TTreeName", "BranchName", &pointer);
	/// \endcode
	/// VariableType is the type of the variable for this branch (fundamental data type or class)\n
	/// TTreeName is the name of the TTree to add this branch\n
	/// BranchName is the name of the branch\n
	/// pointer is a pointer to the variable (should be declared in the header file)\n
	/// \n
	/// To create a standard TTree containing KineParts (for candidates) use\n
	///	\code
	///     CreateStandardTree("TTreeName", "TTreeTitle");
	/// \endcode
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::InitHist(){
	/// \MemberDescr
	/// Book and Initialize histograms in this function.\n
	/// Same function to Book TH1, TH2, TGraph and TGraphAsymmErrors (anything derived from TH1 or TGraph)
	///	\code
	/// 	BookHisto(histogram*)
	/// \endcode
	/// If isAutotUpdate is true, this histogram will be drawn and updated regularly during the processing (default=false).\n
	/// The refresh interval can be set with (default=10):
	/// \code
	/// 	SetUpdateInterval(interval)
	/// \endcode
	/// Defining plots as AutoUpdate and setting the interval can also be done at runtime with a configuration file.\n
	/// \n
	/// Example of booking an histogram: \n
	///	\code
	/// 	BookHisto(new TH2I("PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
	/// \endcode
	/// THe histogram can be saved in a subdirectory of the output files in several ways:\n
	/// \code
	/// 	BookHisto(new TH2I("PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy), refresh, "dir1/dir2");
	/// 	BookHisto(new TH2I("dir1/dir2/PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
	/// 	BookHisto("dir1/dir2/PartEnergy", new TH2I("hPartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
	/// \endcode
	/// In the first case, the ID of the histo in the framework (to use with FillHisto, ...) and its name in the output file is "PartEnergy".
	/// In the second and third case, the ID is "dir1/dir2/PartEnergy" and the name in the output file is respectively "PartEnergy" and "hPartEnergy".
	/// Booking of counters and creation of EventFraction can be done here with\n
	///	\code
	///		BookCounter(name)
	///		NewEventFraction(name)
	/// \endcode
	/// Example\n
	///	\code
	/// 	BookCounter("Total");
	/// 	BookCounter("PassCuts");
	/// 	NewEventFraction("Cuts");
	/// \endcode
	/// Add the counters to the EventFraction\n
	///	\code
	/// 	AddCounterToEventFraction(EventFractionName, CounterName)
	/// \endcode
	/// Example\n
	///	\code
	/// 	AddCounterToEventFraction("Cuts", "Total");
	/// 	AddCounterToEventFraction("Cuts", "PassCuts");
	/// \endcode
	/// Then define which counter represents the sample size\n
	///	\code
	/// 	DefineSampleSizeCounter(EventFractionName, CounterName)
	/// \endcode
	/// Example\n
	///	\code
	/// 	DefineSampleSizeCounter("Cuts", "Total");
	/// \endcode
	/// You can check the content of the input file (directory, histograms, generic key) in any directory with the following methods
	/// \code
	/// 	GetListOfKeys("CedarMonitoring"); // CedarMonitoring directory
	/// 	GetListOfHisto(""); // Top directory
	/// 	GetListOfTH1("CedarMonitoring");
	/// 	GetListOfTH2("CedarMonitoring");
	/// 	GetListOfTGraph("CedarMonitoring");
	/// 	GetListOfDirs("");
	/// \endcode
	/// The first one returns a vector of IOHandler::keyPair containing the name and class name of the object
	/// \code
	/// 	vector<IOHandler::keyPair> keys = GetListOfKeys("CedarMonitoring");
	/// 	cout << keys[0].name << " " << keys[0].className << endl;
	/// \endcode
	/// The others returns a vector of TString containing the name of the objects.\n
	/// You can retrieve histograms from the input ROOT file under the directory "dir1/dir2" (Anything derived from TH1) with\n
	///	\code
	/// 	RequestHistogram("dir1/dir2", "HistogramName", appendOnNewFile);
	/// 	RequestHistogram("dir1", "dir2/HistogramName", appendOnNewFile);
	/// \endcode
	/// In the first case the ID in the framework is "HistogramName", while it is "dir2/HistogramName" in the second case.
	/// appendOnNewFile is a boolean. If set to true, each time a new file is opened the content
	/// of the histogram will be appended to the content of the previous one. If set to false, the content
	/// of the histogram is replaced each time a new file is opened.
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::DefineMCSimple(){
	/// \MemberDescr
	/// Setup of fMCSimple. You must specify the generated MC particles you want.\n
	/// Add particles you want to recover from fMCSimple\n
	///	\code
	/// 	int particleID = fMCSimple.AddParticle(parentID, pdgCode)
	///	\endcode
	/// parentID : 	0=no parent (=beam particle)\n
	/// 	...\n
	/// Example : you want to retrieve the kaon from the beam, the pi0 an pi+ from the beam kaon and the 2 photons coming from the previous pi0 decay :\n
	///	\code
	/// 	int kaonID = fMCSimple.AddParticle(0, 321) //Ask beam kaon (sequence ID=1)
	/// 	fMCSimple.AddParticle(kaonID, 211) //Ask pi+ from previous kaon (sequence ID=2)
	/// 	int pi0ID = fMCSimple.AddParticle(kaonID, 111) //Ask pi0 from previous kaon (sequence ID=3)
	/// 	fMCSimple.AddParticle(pi0ID, 22) //Ask first gamma from previous pi0 (sequence ID=4)
	/// 	fMCSimple.AddParticle(pi0ID, 22) //Ask second gamma from previous pi0 (sequence ID=4)
	///	\endcode
	///
	/// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::StartOfRunUser(){
	/// \MemberDescr
	/// This method is called at the beginning of the processing (corresponding to a start of run in the normal NA62 data taking)\n
	/// Do here your start of run processing if any
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::StartOfBurstUser(){
	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the beginning of the first file\n
	/// Do here your start/end of burst processing if any
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::ProcessSpecialTriggerUser(int iEvent, unsigned int triggerType){
	/// \MemberDescr
	/// \param iEvent : Special event number
	/// \param triggerType : Special trigger type (-1 if not known). For this
	/// variable to be filled, RequestL0SpecialTrigger must be called in the constructor.
	///
	/// Process method for special triggers. Called on each special trigger event after each start of burst.\n
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::Process(int iEvent){
	/// \MemberDescr
	/// \param iEvent : Event number
	///
	/// Main process method. Called on each event. Write you analysis here.\n
	/// You can retrieve MC particles from the fMCSimple set with (returns a vector<KinePart*>)\n
	/// \code
	/// 	fMCSimple["particleName"]
	/// 	fMCSimple[pdgID]
	/// \endcode
	/// Example\n
	/// \code
	/// 	fMCSimple["K+"][index]; //for the kaon
	/// 	fMCSimple["pi+"][index]; //for the positive pion
	/// 	fMCSimple["gamma"][index]; //for the photon
	/// \endcode
	/// The number in the brackets is the index of the particle (if you asked for two photons in the set, you can ask fMCSimple["gamma"][0] for the first one and fMCSimple["gamma"][1] for the second)\n
	/// \n
	/// If you need a property of a particle, you can make a call to fParticleInterface (instance of the ParticleInterface class).\n
	///	This class has two methods FindParticle that will return a TParticlePDG with the required particle. You can search by pdgID or by name.\n
	///	This class also provide two methods to switch between particle name and pdgID if necessary.\n
	///	Example\n
	/// \code
	/// 	double kaonMass = fParticleInterface->FindParticle(321).Mass();
	/// 	double pi0Lifetime = fParticleInterface->FindParticle("pi0").Lifetime();
	/// \endcode
	/// You can retrieve the events from the trees with\n
	/// \code
	/// 	GetEvent<eventClass>();
	/// 	GetEvent<eventClass>("Digis");
	/// 	(eventClass*)GetEvent("detectorName");
	/// 	(eventClass*)GetEvent("detectorName", "Digis");
	/// \endcode
	/// The first two forms can be used only if the class is of type
	/// T<i>DetectorName</i>Event or TReco<i>DetectorName</i>Event.\n
	/// You can retrieve data from generic TTrees with\n
	/// \code
	/// 	GetObject<MyClass>("treeName");			// For class objects
	/// 	GetPrimitiveObject<int>("treeName");	// For primitive datatype (int, double, Int_t, ...)
	/// \endcode
	/// You can retrieve full MC events if available ( GetWithMC() ) with\n
	/// \code
	/// 	GetMCEvent();
	/// 	GetMCEvent("Digis");
	/// \endcode
	/// You can retrieve EventHeader if available ( GetWithEventHeader() ) with\n
	/// \code
	/// 	GetEventHeader();
	/// 	GetEventHeader("Digis");
	/// \endcode
	/// You can retrieve Trigger data if requested with\n
	/// \code
	/// 	GetL0Data();
	/// 	GetL1Data();
	/// 	GetL2Data();
	/// \endcode
	/// You can retrieve the histograms you booked (for drawing, changing, filling, ...) with\n
	/// \code
	/// 	fHisto.GetTH1("histoName");// for TH1
	/// 	fHisto.GetTH2("histoName");// for TH2
	/// 	fHisto.GetTGraph("graphName");// for TGraph and TGraphAsymmErrors
	/// 	fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
	/// \endcode
	/// To fill the histograms you can use\n
	/// \code
	/// 	FillHisto("histoName", values)
	/// \endcode
	/// where values are the same parameters as if you call histogram->Fill(values) (x,y,weight,...)\n
	/// If the histogram is not found, an error message is printed\n
	/// \n
	/// Modify a counter with one of the following methods\n
	/// \code
	/// 	IncrementCounter(name)
	/// 	IncrementCounter(name, delta)
	/// 	DecrementCounter(name)
	/// 	DecrementCounter(name, delta)
	/// 	SetCounterValue(name, value)
	/// \endcode
	/// \n
	/// For use of fGeom, read DetectorAcceptance class.\n
	///	WARNING: this class provides "exact" results, there is not tolerance. If the particle\n
	///	passes in the sensitive volume of a detector it's considered as detected, wether it's close\n
	///	to the edge or not. But as the class gives you the position of the passage point and the \n
	///	estimated path length in the sensitive volume, you can apply further cuts from this\n
	///	information at your convenience.\n
	/// \n
	/// To use the output of a different analyzer, use\n
	/// \code
	/// 	outputType *var = GetOutput<outputType>("analyzerName.outputName", state);
	/// \endcode
	/// Where outputType is the variable type and state is of type outputState\n
	/// State is set with the state of the variable (kOUninit, kOInvalid ,kOValid). The value of the output should only be trusted if state == kOValid\n
	/// example :
	/// \code
	/// 	TLorentzVector vertex = *(TLorentzVector*)GetOutput("simpleVertexAnalyzer.vertex", state);
	/// \endcode
	/// Before starting the processing of an event, the state flag of each output variable is reset to kOUninit\n
	/// When setting the value of an output variable, don't forget to set appropriately the state flag to either kOValid or kOInvalid\n
	/// to indicate if the value can/can't be used in other analyzer\n
	/// \code
	/// 	SetOutputState("outputName", kOValid);
	/// \endcode
	/// If you want to append a candidate in one of your standard output Tree, use\n
	/// \code
	/// 	KinePart *candidate = CreateStandardCandidate("treeName");
	/// \endcode
	/// and fill the properties of your candidate. It will be automatically written in the output tree.\n
	///	\n
	/// If you want to save this event in your custom and standard TTrees (not the input tree replication), call\n
	/// \code
	/// 	FillTrees();
	/// 	FillTrees("treeName");
	/// \endcode
	/// This will call the Fill method of every TTree created in this analyzer.\n
	///	\n
	/// If you want to replicate this event in the output file, call\n
	/// \code
	/// 	FilterAccept();
	/// \endcode
	/// The structure of all the trees that have been opened (by all Analyzer) will be copied in the output file\n
	/// and the events for which at least one analyzer called FilterAccept() will be replicated in the output trees.\n
	/// For more information, please refer to the \ref eventFiltering "Event Filtering" page
	/// \n
	/// The primitives associated to the current event can be retrieved in two ways:\n
	/// \code
	/// FindAllPrimitiveInMatchingWindow("detName");
	/// \endcode
	/// allows to get all the primitives in a certain window around the event timestamp (L0MatchingWindowWidth)
	/// while
	/// \code
	/// FindMatchingPrimitive("detName");
	/// \endcode
	/// allows to get the primitive closest the the event timestamp (but nevertheless withing the
	/// L0MatchingWindowWidth).
	/// @see ROOT TParticlePDG for the particle properties
	/// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
	/// \EndMemberDescr
	if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
	if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}



}

void FilterMuonElectronThreeTrackVertex.cc::PostProcess(){
	/// \MemberDescr
	/// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
	/// during the Process.
	/// \EndMemberDescr

}

void FilterMuonElectronThreeTrackVertex.cc::EndOfBurstUser(){
	/// \MemberDescr
	/// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
	/// Do here your start/end of burst processing if any.
	/// Be careful: this is called after the event/file has changed.
	/// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::EndOfRunUser(){
	/// \MemberDescr
	/// This method is called at the end of the processing (corresponding to a end of run in the normal NA62 data taking)\n
	/// Do here your end of run processing if any\n
	/// \EndMemberDescr

}

void FilterMuonElectronThreeTrackVertex.cc::EndOfJobUser(){
    /// \MemberDescr
    /// Called at the end of job, just before exiting NA62Analysis.\n
    /// Do here all output operation needed (write plots, objects)
    /// \n
	/// If you want to save all plots, just call\n
	/// \code
	/// 	SaveAllPlots();
	/// \endcode
	/// Or you can just save the ones you want with\n
	/// \code
	/// 	histogram->Write()\n
	///		fHisto.Get...("histoname")->Write();
	/// \endcode
	/// \n
	/// To run over a set of histograms you can use Iterators (HistoHandler::IteratorTH1,
	/// HistoHandler::IteratorTH2, HistoHandler::IteratorTGraph). You can use it to run over
	/// all the histograms or only a subset of histogram by using one of the two forms of
	/// GetIterator...  (replace ... by TH1, TH2 or TGraph)\n
	/// \code
	/// 	GetIterator...()
	/// \endcode
	/// will get an Iterator running over all histograms of this type while
	/// \code
	/// 	GetIterator...("baseName")
	/// \endcode
	/// will get an Iterator running only over the histograms of this type whose name starts
	/// with baseName.\n
	/// For more details and examples on how to use the Iterator after getting it, please refer
	/// to the HistoHandler::Iterator documentation.\n
	/// Although this is described here, Iterators can be used anywhere after the
	/// histograms have been booked.
    /// \EndMemberDescr
}

void FilterMuonElectronThreeTrackVertex.cc::DrawPlot(){
	/// \MemberDescr
	/// This method is called at the end of processing to draw plots when the -g option is used.\n
	/// If you want to draw all the plots, just call\n
	/// \code
	/// 	DrawAllPlots();
	/// \endcode
	/// Or get the pointer to the histogram with\n
	/// \code
	/// 	fHisto.GetTH1("histoName");// for TH1
	/// 	fHisto.GetTH2("histoName");// for TH2
	/// 	fHisto.GetGraph("graphName");// for TGraph and TGraphAsymmErrors
	///     fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
	/// \endcode
	/// and manipulate it as usual (TCanvas, Draw, ...)\n
	/// \EndMemberDescr
}

FilterMuonElectronThreeTrackVertex.cc::~FilterMuonElectronThreeTrackVertex.cc(){
	/// \MemberDescr
	/// Destructor of the Analyzer. If you allocated any memory for class
	/// members, delete them here.
	/// \EndMemberDescr
}
