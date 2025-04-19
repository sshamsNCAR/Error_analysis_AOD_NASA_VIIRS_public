/*
 * 
 */

#ifndef _ERROR_HANDLER_H_
#define _ERROR_HANDLER_H_

#include <iostream>
using namespace std;

class Error_Handler {
     public:
        static Error_Handler* getInstance();
        void warningMessage(const char* msg); 	
	void fatalMessage(const char* msg);
     private:
         Error_Handler(){};
         static Error_Handler *instance_; 
};

#endif /*_ERROR_HANDLER_H_*/
