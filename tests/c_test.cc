#include <iostream>
#include <map>

typedef std::map<int, int> Dict;
typedef Dict::const_iterator It;

int main()
{
   Dict d;
   for (int l=0; l<10; l ++ ){
   d[l] = l*2;
   }
   printf("Double dictionnary %d %d\n",12, d[12]);
   // for (It it(d.begin()); it != d.end(); ++it)
   // {
   //    int i(it->first.first);
   //    int j(it->first.second);
   //    std::cout <<it->second <<' '
   //              <<d[std::make_pair(j, i)] <<'\n';
   // }
}
