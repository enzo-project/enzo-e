#include <string.h>

#define UNIT_ASSERT(RESULT)		\
  printf ("%s %s:%d  %s::%s() %d\n",	\
	  (RESULT)?"pass":"FAIL",__FILE__,__LINE__,class_name,func_name,test_num++);


#define UNIT_MAX_NAME_LEN 40

char class_name[UNIT_MAX_NAME_LEN];
char func_name[UNIT_MAX_NAME_LEN];

int test_num = 1;

void UNIT_CLASS (char * c)
{
  strncpy (class_name,c,UNIT_MAX_NAME_LEN);
};

void UNIT_FUNC (char * f)
{
  strncpy (func_name,f,UNIT_MAX_NAME_LEN);
};

