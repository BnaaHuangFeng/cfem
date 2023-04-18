#pragma once
#include"InputSystem/DescriptionInfo.h"
#include<string>
#include "nlohmann/json.hpp"
class Timer;
/**
 * this class read the input file, and store the most fundamental information about the problem
 * it also implement the init of any system including and preallocation of variable.
*/
class InputSystem{
public:
    InputSystem():m_timer(NULL){};
    InputSystem(Timer *t_timer);
    InputSystem(int argc, char ** argv);
    ~InputSystem();
    /**
     * initialize the input file reading system
     * @param args integer number of total argv
     * @param argv char vector taken from command line
     */
    void init(int args,char *argv[]);
    /**
     * Read input file: xxx.json
    */
    void readFile();
    /**
     * Read the Mesh description block of input json file,
     * and set the Mesh description
    */
    bool readMeshBlock(nlohmann::json &t_json);
    /**
     * Read the Element description block of input json file,
     * and set the Mesh description
    */
    bool readElementBlock(nlohmann::json &t_json);
    /**
     * Read the Material description block of input json file,
     * and set the Material description
    */
    bool readMaterialBlock(nlohmann::json &t_json);
    /**
     * Read the Step description block of input json file,
     * and set the Step description
    */
    bool readStepBlock(nlohmann::json &t_json);
    /**
     * Read the Output description block of input json file,
     * and set the Output description
    */
    bool readOutputBlock(nlohmann::json &t_json);
    /**
     * Read the Boundaty condition description block of input json file,
     * and set the Boundaty condition description
    */
    bool readBcBlock(nlohmann::json &t_json);   
public:
    bool m_readonly;    /**< if it's true, it will only read the mesh block*/
    bool m_completed;   /**< if ture mean the input file's read is completed*/
    std::string m_inputfile_name;
    nlohmann::json m_json;          /**< input file json*/
    MeshDescription m_meshDes;      /**< Mesh Description*/
    MaterialDescription m_MatDes;   /**< Matrial description*/
    ElementDescription m_ElDes;     /**< Element Description*/
    StepDescriptiom m_stepDes;      /**< Step Descriptiom*/
    OutputDescription m_outDes;     /**< Output Description*/
    BCDescription m_bcDes;          /**< Boundaty condition Description*/
    Timer *m_timer;                 /**< a timer per processor*/
};