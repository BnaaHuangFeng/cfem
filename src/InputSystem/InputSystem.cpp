#include"InputSystem/InputSystem.h"
#include"Utils/MessagePrinter.h"
#include"Utils/Timer.h"
#include"fstream"
#include"petsc.h"
InputSystem::InputSystem(Timer *t_timer){
    m_timer=t_timer;
    m_readonly=false;
    m_completed=false;
};

InputSystem::InputSystem(int argc,char ** argv){
    init(argc,argv);
};

InputSystem::~InputSystem(){
};

void InputSystem::init(int argc,char *argv[]){
    m_readonly=false;
    m_completed=false;
    if(argc==1){
        // ./cfem or cfem
        MessagePrinter::printErrorTxt("invalid command line argc, the second argc must be '-i'");
        MessagePrinter::exitcfem();
    }
    else if(argc==3){
        // ./cfem -i input.json or cfem -i input.json
        if(string(argv[1]).find("-i")!=string::npos){
            if(string(argv[2]).size()<5){
                MessagePrinter::printErrorTxt("invalid input file name after '-i', it must be xxx.json");
                MessagePrinter::exitcfem();
            }
            m_inputfile_name=argv[2];
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid input file name, your input file must have the extension '.json'");
                MessagePrinter::exitcfem();
            }
        }
        else{
            MessagePrinter::printErrorTxt("invalid command line argc, the second argc must be '-i'");
            MessagePrinter::exitcfem();
        }
    }
    else{
        if(string(argv[1]).find("-i")!=string::npos){
            m_inputfile_name=argv[2];
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid command line argc, the second argc must be '-i'");
                MessagePrinter::exitcfem();
            }
        }
        else{
            MessagePrinter::printErrorTxt("invalid command line argc, the second argc must be '-i'");
            MessagePrinter::exitcfem();
        }
        for(int i=4-1;i<argc;i++){
            if(string(argv[i]).find("--read-only")!=string::npos){
                m_readonly=true;
            }
        }
    }
}

void InputSystem::readFile(){
    ifstream in;
    nlohmann::json meshJson,elementJson,stepJson,outputJson,bcJson;
    m_timer->startTimer();
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Start to read the input file",MessageColor::BLUE);
    in.open(m_inputfile_name.c_str(),ios::in);
    while(!in.is_open()){
        MessagePrinter::printWarningTxt("can\'t open the input file(name="+m_inputfile_name+")");
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
        cin>>m_inputfile_name;
        in.open(m_inputfile_name.c_str(),ios::in);
    }
    string str;
    stringstream strStream;
    strStream<<in.rdbuf();
    str=strStream.str();
    in.close();
    m_json=nlohmann::json::parse(str);
    /*******************************************
     * mesh reading                 ************
    ********************************************/
    if(m_json.contains("mesh")){
        readMeshBlock(m_json.at("mesh"));
        MessagePrinter::printNormalTxt("mesh reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("Mesh block is missing!");
        MessagePrinter::exitcfem();
    }
    /*******************************************
     * element reading              ************
    ********************************************/
    if(m_json.contains("element")){
        readElementBlock(m_json.at("element"));
        MessagePrinter::printNormalTxt("element reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("element block is missing!");
        MessagePrinter::exitcfem();
    }
    /*******************************************
     * material reading             ************
    ********************************************/
    if(m_json.contains("material")){
        readMaterialBlock(m_json.at("material"));
        MessagePrinter::printNormalTxt("material reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("material block is missing!");
        MessagePrinter::exitcfem();
    }
    /*******************************************
     * step reading                 ************
    ********************************************/
    if(m_json.contains("step")){
        readStepBlock(m_json.at("step"));
        MessagePrinter::printNormalTxt("step reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("Mesh block is missing!");
        MessagePrinter::exitcfem();
    }
    /*******************************************
     * output reading               ************
    ********************************************/
    if(m_json.contains("output")){
        readOutputBlock(m_json.at("output"));
        MessagePrinter::printNormalTxt("output reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("output block is missing!");
        MessagePrinter::exitcfem();
    }
    /*******************************************
     * BCs reading                 *************
    ********************************************/
    if(m_json.contains("bcs")){
        readBcBlock(m_json.at("bcs"));
        MessagePrinter::printNormalTxt("bcs reading completed!");
    }
    else{
        MessagePrinter::printErrorTxt("bcs block is missing!");
        MessagePrinter::exitcfem();
    }
    m_timer->endTimer();
    m_timer->printElapseTime("Input file reading is done",false);
    MessagePrinter::printNormalTxt("input file's reading is done",MessageColor::BLUE);
    MessagePrinter::printDashLine(MessageColor::BLUE);
}

bool InputSystem::readMeshBlock(nlohmann::json &t_json){
    if(t_json.contains("outputfile"))
        m_meshDes.s_outputMeshFile_Name=t_json.at("outputfile");
    if(t_json.contains("inputfile"))
        m_meshDes.s_inputMeshFile_Name=t_json.at("inputfile");
    m_meshDes.s_nx=0;
    m_meshDes.s_ny=0;
    m_meshDes.s_nz=0;
    m_meshDes.s_ifSaveMesh=true;
    m_meshDes.s_shape=MeshShape::COMPLEX;
    // read mesh mode
    if(t_json.contains("mode")){
        if(t_json.at("mode").is_string()){
            string mode_str=t_json.at("mode");
            if(mode_str=="structured"){
                m_meshDes.s_mode=MeshMode::STRUCTURED;
            }
            else if(mode_str=="unstructured"){
                m_meshDes.s_mode=MeshMode::UNSTRUCTURED;
            }
            else{
                MessagePrinter::printErrorTxt(mode_str+" is not a supported mesh type");
                MessagePrinter::exitcfem();
            }
        }
        else{
            MessagePrinter::printErrorTxt("input mesh->mode is not a string type!");
            MessagePrinter::exitcfem();
        }
    }
    // read problem dimension
    if(t_json.contains("dim")){
        if(t_json.at("dim").is_number_integer()){
            int dim=t_json.at("dim");
            switch (dim)
            {
            case 1:
                m_meshDes.s_dim=Dimension::ONE;
                break;
            case 2:
                m_meshDes.s_dim=Dimension::TWO;
                break;
            case 3:
                m_meshDes.s_dim=Dimension::THREE;
                break;
            default:
                MessagePrinter::printErrorTxt("imput->mesh->dim = "+to_string(dim));
                MessagePrinter::exitcfem();
                break;
            }
        }
        else{
            MessagePrinter::printErrorTxt("input mesh->dim is not a interger!");
            MessagePrinter::exitcfem();
        }
    }
    // read the complete mesh outer shape
    string shapeName=t_json.at("shape");
    if(shapeName=="rectangular"){
        m_meshDes.s_shape=MeshShape::RECTANGULAR;
    }
    else if(shapeName=="sin"){
        m_meshDes.s_shape=MeshShape::SIN;
    }
    else if(shapeName=="half-sin"){
        m_meshDes.s_shape=MeshShape::HALFSIN;
    }
    else{
        MessagePrinter::printErrorTxt(shapeName+" is not a supported mesh shape");
        MessagePrinter::exitcfem();
    }
    string meshTypeName=t_json.at("meshtype");
    if(meshTypeName=="quad4"){
        m_meshDes.s_type=MeshType::QUAD4;
    }
    else if(meshTypeName=="quad8"){
        m_meshDes.s_type=MeshType::QUAD8;
    }
    m_meshDes.s_ifSaveMesh=t_json.at("savemesh");
    m_meshDes.s_nx=t_json.at("nx");
    m_meshDes.s_ny=t_json.at("ny");
    m_meshDes.s_size_json=t_json.at("size");
    return true;
};

bool InputSystem::readElementBlock(nlohmann::json &t_json){
    for(auto it=t_json.begin();it!=t_json.end();it++){
        string elementName=it.key();
        m_ElDes.s_names.push_back(elementName);
        nlohmann::json elmt_json=t_json.at(elementName);/**< a json block of a kind of element type*/
        string elementTypeName=elmt_json.at("type");
        // read element type
        if(elementTypeName=="CPE4R"){
            m_ElDes.s_elmtTypes.push_back(ElementType::CPE4R);
        }
        // read element assignmen set
        m_ElDes.s_setNames.push_back(elmt_json.at("set"));
    }
    return true;
}

bool InputSystem::readMaterialBlock(nlohmann::json &t_json){
    for(auto it=t_json.begin();it!=t_json.end();it++){
        string materialName=it.key();
        m_MatDes.s_names.push_back(materialName);
        nlohmann::json mat_json=t_json.at(materialName);/**< a json block of a kind of element type*/
        string matTypeName=mat_json.at("type");
        if(matTypeName=="linearelastic"){
            m_MatDes.s_matType.push_back(MaterialType::LINEARELASTIC);
        }
        else if(matTypeName=="neohookean"){
            m_MatDes.s_matType.push_back(MaterialType::NEOHOOKEAN);
        }
        else if(matTypeName=="vonmises"){
            m_MatDes.s_matType.push_back(MaterialType::VONMISESPLAS);
        }
        nlohmann::json prop_json=mat_json.at("parameters");
        m_MatDes.s_properties.push_back(prop_json);
        m_MatDes.s_setName.push_back(mat_json.at("set"));
    }
    return true;
}

bool InputSystem::readStepBlock(nlohmann::json &t_json){
    // read solution algorithm method
    string algorithm=t_json.at("method");
    if(algorithm=="standard"){
        m_stepDes.s_algorithm=AlgorithmType::STANDARD;
    }
    else if(algorithm=="arclength_cylender"){
        m_stepDes.s_algorithm=AlgorithmType::ARCLENGTH_CYLENDER;
    }
    else{
        MessagePrinter::printErrorTxt(algorithm+" is not a supported solution algorithm method");
        MessagePrinter::exitcfem();
    }
    string nlsolver=t_json.at("nlsolver");
    if(nlsolver=="newtonls"){
        m_stepDes.s_SNESType=SNESNEWTONLS;
    }
    else if(nlsolver=="newtontr"){
        m_stepDes.s_SNESType=SNESNEWTONTR;
    }
    else if(nlsolver=="nrichardson"){
        m_stepDes.s_SNESType=SNESNRICHARDSON;
    }
    else if(nlsolver=="ksponly"){
        m_stepDes.s_SNESType=SNESKSPONLY;
    }
    else if(nlsolver=="ngmres"){
        m_stepDes.s_SNESType=SNESNGMRES;
    }
    else{
        MessagePrinter::printErrorTxt(nlsolver+" is not a supported SNESType.");
        MessagePrinter::exitcfem();
    }
    // read KSP solver method
    string kspType=t_json.at("kspsolver");
    if(kspType=="richardson"){
        m_stepDes.s_KSPType=KSPRICHARDSON;
    }
    else if(kspType=="chebyshev"){
        m_stepDes.s_KSPType=KSPCHEBYSHEV;
    }
    else if(kspType=="cg"){
        m_stepDes.s_KSPType=KSPCG;
    }
    else if(kspType=="gmres"){
        m_stepDes.s_KSPType=KSPGMRES;
    }
    else if(kspType=="bcgs"){
        m_stepDes.s_KSPType=KSPBCGS;
    }
    else{
        MessagePrinter::printErrorTxt(kspType+" is not a supported KSPType.");
        MessagePrinter::exitcfem();
    }
    // read preconditioner
    string pcType=t_json.at("preconditioner");
    if(pcType=="none"){
        m_stepDes.s_PCType=PCNONE;   
    }
    else if(pcType=="jacobi"){
        m_stepDes.s_PCType=PCJACOBI;   
    }
    else if(pcType=="sor"){
        m_stepDes.s_PCType=PCSOR;   
    }
    else if(pcType=="lu"){
        m_stepDes.s_PCType=PCLU;
    }
    else if(pcType=="qr"){
        m_stepDes.s_PCType=PCQR;
    }
    else if(pcType=="shell"){
        m_stepDes.s_PCType=PCSHELL;
    }
    else if(pcType=="bjacobi"){
        m_stepDes.s_PCType=PCBJACOBI;
    }
    else if(pcType=="mg"){
        m_stepDes.s_PCType=PCMG;
    }
    else if(pcType=="eisenstat"){
        m_stepDes.s_PCType=PCEISENSTAT;
    }
    else if(pcType=="ilu"){
        m_stepDes.s_PCType=PCILU;
    }
    else if(pcType=="icc"){
        m_stepDes.s_PCType=PCICC;
    }
    else if(pcType=="asm"){
        m_stepDes.s_PCType=PCASM;
    }
    else{
        MessagePrinter::printErrorTxt(pcType+" is not a supported PCType.");
        MessagePrinter::exitcfem();
    }
    // read total t
    if(t_json.at("t").is_number_float()){
        m_stepDes.s_t=t_json.at("t");
        if(m_stepDes.s_t<=0){
            MessagePrinter::printErrorTxt("t must be positive.");
            MessagePrinter::exitcfem();             
        }
    }
    else{
        MessagePrinter::printErrorTxt("t must be a float number.");
        MessagePrinter::exitcfem();              
    }
    // read dtmin
    if(t_json.at("dtmin").is_number_float()){
        m_stepDes.s_dtmin=t_json.at("dtmin");
        if(m_stepDes.s_dtmin<=0||m_stepDes.s_dtmin>m_stepDes.s_t){
            MessagePrinter::printErrorTxt("t must be in (0.0, t]");
            MessagePrinter::exitcfem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("dtmin must be a float number.");
        MessagePrinter::exitcfem();        
    }
    // read dtmax
    if(t_json.at("dtmax").is_number_float()){
        m_stepDes.s_dtmax=t_json.at("dtmax");
        if(m_stepDes.s_dtmax<m_stepDes.s_dtmin||m_stepDes.s_dtmax>m_stepDes.s_t){
            MessagePrinter::printErrorTxt("dtmax must be in [dtmin, t]");
            MessagePrinter::exitcfem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("dtmax must be a float number.");
        MessagePrinter::exitcfem();              
    }
    // read dt0
    if(t_json.at("dt0").is_number_float()){
        m_stepDes.s_dt0=t_json.at("dt0");
        if(m_stepDes.s_dt0<m_stepDes.s_dtmin||m_stepDes.s_dt0>m_stepDes.s_dtmax){
            MessagePrinter::printErrorTxt("dto must be in [dtmin, dtmax].");
            MessagePrinter::exitcfem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("dt0 must be a float number.");
        MessagePrinter::exitcfem();               
    }
    // read growth factor
    if(t_json.at("growth-factor").is_number_float()){
        m_stepDes.s_growFactor=t_json.at("growth-factor");
        if(m_stepDes.s_growFactor<1.0){
            MessagePrinter::printErrorTxt("growth-factor must be greater than 1.0.");
            MessagePrinter::exitcfem();            
        }
    }
    else{
        MessagePrinter::printErrorTxt("growth-factor must be a float number.");
        MessagePrinter::exitcfem();               
    }
    // read cutback factor
    if(t_json.at("cutback-factor").is_number_float()){
        m_stepDes.s_cutbackFactor=t_json.at("growth-factor");
        if(m_stepDes.s_cutbackFactor<=0.0||m_stepDes.s_cutbackFactor>=1.0){
            MessagePrinter::printErrorTxt("cutback-factor must be in (0.0, 1.0).");
            MessagePrinter::exitcfem();            
        }
    }
    else{
        MessagePrinter::printErrorTxt("cutback-factor must be a float number.");
        MessagePrinter::exitcfem();               
    }    
    // read converage tolerance
    m_stepDes.s_relTol=t_json.at("rel-tolerance");
    m_stepDes.s_absTol=t_json.at("abs-tolerance");
    m_stepDes.s_duTol=t_json.at("du-tolerance");
    m_stepDes.s_maxIterNum=t_json.at("maxiters");
    return true;
}
bool InputSystem::readOutputBlock(nlohmann::json &t_json){
    nlohmann::json field_json=t_json.at("field");
    nlohmann::json his_json_vec=t_json.at("history");
    string format;
    int variableNum;
    /** read Field Output Description*********************/
    // read Filed Output Format
    format=field_json.at("format");
    if(format=="vtu"){
        m_outDes.s_FD.s_format=FiledOutputFormat::VTU;
    }
    else{
        MessagePrinter::printErrorTxt(format+" is not a supported field output format.");
        MessagePrinter::exitcfem();
    }
    // read output interval
    m_outDes.s_FD.s_interval=field_json.at("interval");
    // read filed variable
    variableNum=static_cast<int>(field_json.at("variable").size());
    for(int i=0;i<variableNum;i++){
        string variableName=field_json.at("variable").at(i);
        if(variableName=="stress"){
            m_outDes.s_FD.s_varTypes.push_back(FieldVariableType::STRESS);
        }
        else if(variableName=="vonMises-stress"){
            m_outDes.s_FD.s_varTypes.push_back(FieldVariableType::VONMISES);
        }
        else if(variableName=="u"){
            m_outDes.s_FD.s_varTypes.push_back(FieldVariableType::U);
        }
        else if(variableName=="log-strain"){
            m_outDes.s_FD.s_varTypes.push_back(FieldVariableType::LOGSTRAIN);
        }
        else{
            MessagePrinter::printErrorTxt(variableName+" is not a supported field variable.");
            MessagePrinter::exitcfem();
        }
    }

    /** read historical Output Description*********************/
    for(auto it=his_json_vec.begin();it!=his_json_vec.end();it++){
        string his_json_name=it.key();
        nlohmann::json his_json=his_json_vec.at(his_json_name);
        HistoryOutputDescription HD; /***< History Output Description*/
        format=his_json.at("format");
        if(format=="csv"){
            HD.s_format=HistoryOutputFormat::CSV;
        }
        else{
            MessagePrinter::printErrorTxt(format+" is not a supported historical output format.");
            MessagePrinter::exitcfem();
        }
        // read output interval
        HD.s_interval=his_json.at("interval");
        // read historical variable
        variableNum=static_cast<int>(his_json.at("variable").size());
        for(int i=0;i<variableNum;i++){
            string variableName=his_json.at("variable").at(i);
            if(variableName=="U1"){
                HD.s_varTypes.push_back(HistoryVariableType::U1);
            }
            else if(variableName=="U2"){
                HD.s_varTypes.push_back(HistoryVariableType::U2);
            }
            else if(variableName=="U3"){
                HD.s_varTypes.push_back(HistoryVariableType::U3);
            }
            else if(variableName=="RF1"){
                HD.s_varTypes.push_back(HistoryVariableType::RF1);
            }
            else if(variableName=="RF2"){
                HD.s_varTypes.push_back(HistoryVariableType::RF2);
            }
            else if(variableName=="RF3"){
                HD.s_varTypes.push_back(HistoryVariableType::RF3);
            }
            else{
                MessagePrinter::printErrorTxt(variableName+" is not a supported historical variable.");
                MessagePrinter::exitcfem();
            }
        }
        HD.s_setName=his_json.at("set");
        m_outDes.s_HD.push_back(HD);
    }
    return true;
}
bool InputSystem::readBcBlock(nlohmann::json &t_json){
    for(auto it=t_json.begin();it!=t_json.end();it++){
        string bc_json_name=it.key();
        nlohmann::json bc_json=t_json.at(bc_json_name);
        SingleBCDes singleBCDes;
        singleBCDes.s_BCName=bc_json_name;
        nlohmann::json dof_json=bc_json.at("dofs");
        int dofNum=static_cast<int>(dof_json.size());
        for(int i=0;i<dofNum;i++){
            singleBCDes.s_presetDofIds.push_back(static_cast<int>(dof_json.at(i)));
        }
        singleBCDes.s_bcVals=bc_json.at("bcvalue");
        string setName=bc_json.at("set");
        singleBCDes.s_setName=setName;
        m_bcDes.push_back(singleBCDes);
    }
    return true;
}