#pragma once

/**
 * this class read the input file, and store the most fundamental information about the problem
 * it also implement the init of any system including and preallocation of variable.
*/
class InputSystem{
public:
    InputSystem();
    InputSystem(int argc, char ** argv);
    ~InputSystem();
    /**
     * initialize the input file reading system
     * @param args integer number of total argv
     * @param argv char vector taken from command line
     */
    void init(int args,char *argv[]);
    
public:
    bool m_readonly;/**< if it's true, it will only read the mesh block*/
    bool m_completed;/**< if ture mean the input file's read is completed*/

};