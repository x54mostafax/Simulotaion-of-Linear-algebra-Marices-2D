#include<iostream>
#include<string>
#include<vector>
#include<array>
using namespace  std;
 
 
 int fx(const vector<string> &EQ){
     for (auto i : EQ){
         vector<vector<long double>> M;
         vector<vector<long double>> MC;
         long double  coffient=0;
         long double n=0;
         int curCol;
         bool isvar=false;
         bool signOfOp=true;
         bool signOfEqual=true;
         int curPointer=0;
         int prevPointer=0;
         bool isfraction=false;
         int curRow=0;
         for(char j :i){
            //  if (curPointer<=57 &&  curPointer>=48) {n=(int)j +n*10;continue;}//0.423 not allow now // if isvar and frac then error => 13x0.5 XXX 
             if (curPointer==prevPointer){
                 if (j=='x' || j=='X'){isvar=true;coffient=1;}
                 
             } 
             else if(curPointer>prevPointer){
                  if (j=='x' || j=='x'){isvar=true;coffient=n*(signOfOp*2-1)*(2*signOfEqual-1);}
                  if(j=='+'||j=='-'||j=='-'){
                      if (isvar){
                          curCol=n;
                          M[curRow][curCol]=coffient;
                          n=0;
                          
                      }//check here if fraction print error
                      else {
                          coffient=n*(signOfOp*2-1)*(1-2*signOfEqual);
                          MC[curRow][1]=coffient;
                          
                      }    
                    }
                       
              }
              if(j=='+' || j=='=') signOfOp=true;
              else if(j=='-') signOfOp=true;
              else if(j=='=')signOfEqual=false;
              coffient=n=0;
         
            //  cout<<j<<"/";
             }; 
             cout<<M.size()<<" Sizw";
             for (int i=0;i<M.size();i++) {
                 cout<<"R("<<i<<") ";
               for (int j=0;i<M[i].size();j++) {
                   cout<< M[i][j];
               }
               cout<<endl;
             }
             
     //cout<<i<<endl;
     }
return 0;  
 }

int main(int argc, char *argv[])
{
for (int i=0;i<127;i++){
    //cout<<i<<" >>  "<<(char)i<<endl;
    }
    
cout<<fx({"-x1+12x2+0.1x3-x4=-x5+c","-x1+5+1=-1+x4+1"});	
			
}
