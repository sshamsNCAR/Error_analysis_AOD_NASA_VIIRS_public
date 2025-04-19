#include "Error_Handler.h"
#include <stdlib.h> 

Error_Handler *Error_Handler::instance_ = NULL;

void Error_Handler::fatalMessage(const char* msg) 
{
   cout << msg << endl;
   abort();
}

void Error_Handler::warningMessage(const char* msg) 
{
   cout << msg << endl;
}

Error_Handler* Error_Handler::getInstance() 
{
   if (!instance_) instance_ = new Error_Handler();
   return instance_; 
}
