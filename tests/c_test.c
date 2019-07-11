#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// typedef (int, int)> tuple;

int encodeTuples(int *tups){
int n  = sizeof(*tups);
int current_i=-1 ; // keep track of analyzed lines
int col[n]; // keep column values (size = len(tups))
int line[n]; // keep values' indexes (size = len(matrix lines))
// for (int t=0; t<sizeof(tups); t++)
//   {
//     while (tups[t][0]!=current_i)
//     {
//       current_i++ ;
//       line[current_i]=t;
//     }
//     col[t]=tups[t][1] ;
//
//   }

while ( current_i < sizeof(line) ) //finish filling line array with t maximal value
{
line[current_i]=sizeof(tups) ;
current_i++ ;
}    // for every remaining value of line array use
// return col ;
return 1;
}


int main(int argc, char *argv[])
{
    int age=23;
    int *premierpointeur = &age ;
    int tups[2]={2,2};
    printf("%d\n", *tups);
    printf("%d\n", *premierpointeur) ;
    printf("This is executed %s \n",  argv[0]);
}


/* static struct { char strVal[21]; int intVal; } tuple[10];
static int tupleCount = 0;

static void listTuples(void) {
    printf("==========\nTuple count is %d\n", tupleCount);
    for (int i = 0; i < tupleCount; ++i)
        printf("   [%s] -> %d\n", tuple[i].strVal, tuple[i].intVal);
    puts("==========");
}

static void addTuple(char *str, int val) {
    printf("Adding '%s', mapped to %d\n", str, val);
    strcpy(tuple[tupleCount].strVal, str);
    tuple[tupleCount++].intVal = val;
}

static void deleteTuple(char *str) {
    int index = 0;
    while (index < tupleCount) {
        if (strcmp(str, tuple[index].strVal) == 0) break;
        ++index;
    }
    if (index == tupleCount) return;

    printf("Deleting '%s', mapped to %d\n", str, tuple[index].intVal);
    if (index != tupleCount - 1) {
        strcpy(tuple[index].strVal, tuple[tupleCount - 1].strVal);
        tuple[index].intVal = tuple[tupleCount - 1].intVal;
    }
    --tupleCount;
} */
