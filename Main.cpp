#include <iostream>
#include<vector>
#include<string>
#include<cmath>
#include <algorithm>
using namespace std;

/*
Matrix 

####Main Parts
--Det 
--dispaly +
--RREF ++

####attribiutes
--dim +
--Transpose +
--INVERSE => RREF + IDENTITY ++
--ADJOINT => Transpose + Cofactor
--COFACTOR ++

###Static items
--IDENTITY Matrix
--ZERO Matrix
--MultipleMatrices

###Opreations on Matrix ->RREF:
--addRow(R2,R1,C)
--multipleRow(R,C)
--swapRows(R1,R2)
--multipleMatrices(Matrix2D M1,Matrix2D M2)
--addMatrices(Matrix2D M1,Matrix2D M2)

###Sqaure Matrix =>
--Det
--Power
--Types(future)=>diag(elements)
1,0,1,0    1 0 0 x
0 1,1,9 => 0 1 0 9-x
5,8,6,0    0 0 1 x
*/
class LinearAlgebra{
    public:
    string Name;
    
};
class Matrix2D{
    private:
    vector<vector<double>>T;
    vector<vector<double>>_RREF;
    vector<vector<double>>_Cofactor;
    public:
    int dim[2];
    vector<vector<double>> Matrix;
    
    
    //initialize 
     Matrix2D(int dims[2]){
        copy(dims,dims+2,dim);
        Matrix.resize(dims[0], vector<double>(dims[1], 0));
        cout<<Matrix[5][5]<<endl;
    };
    
    
    Matrix2D(vector<vector<double>> &M){
        Matrix=M;
        dim[0]=M.size();
        dim[1]=M[0].size();
        
    };
    
    void display(){
        cout<<"X ";
         for (int i = 0; i <dim[1]; i++) {
             cout<<""<<i+1<<" ";
         }
         cout<<endl;
         for (int i = 0; i<dim[0];i++) {
              cout<<""<<i+1;
              
             for (int j = 0; j<dim[1];j++) {
                cout<<' '<<Matrix[i][j];
             }
             cout<<endl;
         }
        cout<<"Dimension :"<<dim[0]<<" Row "<<dim[1]<< " Column\n";
    }
    
    
    
    
    
    // Square type Matrices 
    static Matrix2D IDENTITY(int order){
        vector<vector<double>>I(order,vector<double>(order,0));
        for (int i = 0; i < order; i++) {
            I[i][i]=1;
        }
        return I;
    }
    
    static Matrix2D ZERO(int order){
         vector<vector<double>> Z (order,vector<double>(order,0));
         return Z;
    }
    static Matrix2D diag(double elements[] , int order ){
          vector<vector<double>>diag(order,vector<double>(order,0));
          for (int i = 0; i < order; i++) {
            diag[i][i]=elements[i];
        }
        return diag;
    }
    
    
    
    
    
    // Transpose 
    static Matrix2D Transpose(Matrix2D &M){
        if(M.T.empty()){
            cout<<"Iam empty\n";
            M.T.resize(M.dim[1],vector<double>(M.dim[0]));
            for (int i = 0; i<M.dim[0];i++) for (int j = 0; j<M.dim[1];j++) M.T[j][i]=M.Matrix[i][j];
             
         };
        return M.T;
    }
    static Matrix2D Cofactor(Matrix2D &M){
        if(!M._Cofactor.empty()) return M._Cofactor;
        vector<int> PassedCoulmn(M.dim[1],0);
        vector<vector<double>>CofactorMatrix(M.dim[0],vector<double>(M.dim[1],0));
        vector<int> PassedRow(M.dim[0],0) ;
        for (int r = 0; r < M.dim[0]; r++) {
            for (int c = 0; c < M.dim[1]; c++) {
                PassedRow[r]=PassedCoulmn[c]=1;
                CofactorMatrix[r][c]= Matrix2D::NDet(M,PassedCoulmn,PassedRow)*pow(-1,r-c);
                PassedCoulmn[c]=PassedRow[r]=0;
            }
        }
        return CofactorMatrix;
    }

    static Matrix2D Adjoint(Matrix2D &M){Matrix2D R= Matrix2D::Cofactor(M);return Matrix2D::Transpose(R);}
    
    static double Det(Matrix2D &M,vector<int> memo={},int cRow=0){
        
        if(memo.empty())memo.resize(M.dim[0],0);
        double Result=0;
        int sign=-1;
        for (int i = 0; i < M.dim[0]; i++) {
            if(memo[i]==0){
                sign*=-1;
                memo[i]=1;
                 if(cRow==M.dim[0]-1){
                    //cout <<  "End Det ("<<cRow<<","<<i<<") => "<<M.Matrix[cRow][i]<<" r " <<Result<< endl;
                    return M.Matrix[cRow][i];
                }
                
                
                Result+=M.Matrix[cRow][i]*Det(M,memo,cRow+1)*sign;
                //cout<<"I Det ("<<cRow<<","<<i<<") => "<< M.Matrix[cRow][i]*pow(-1,i)<<" r "<<Result<<endl;
                memo[i]=0;
            }
            
        }
        return Result;
    }
    
     static double NDet(Matrix2D &M,vector<int> PassedCoulmn={},vector<int> PassedRow={},int cRow=0){
        
        if(PassedCoulmn.empty())PassedCoulmn.resize(M.dim[0],0);
        if(PassedRow.empty())PassedRow.resize(M.dim[0],0);
        double Result=0;
        int sign=-1;
        if(PassedRow[cRow])++cRow;
        if(cRow==M.dim[0]){
                    // cout <<  "End Det ("<<cRow<<","<<") => "<<Result<< endl;
                    return 1;
        }
        for (int i = 0; i < M.dim[0]; i++) {
            if(PassedCoulmn[i]==0){
                sign*=-1;
                 
                PassedCoulmn[i]=PassedRow[cRow]=1;
                Result+=M.Matrix[cRow][i]*NDet(M,PassedCoulmn,PassedRow,cRow+1)*sign;
                // cout<<"I Det ("<<cRow<<","<<i<<") => "<< M.Matrix[cRow][i]*pow(-1,i)<<" r "<<Result<<endl;
                 PassedCoulmn[i]=0;
            }
            
        }
        return Result;
         
     }
     
   
   
    static Matrix2D MultipleMatrices(Matrix2D &M1,Matrix2D &M2){
        // m*n X n*p =. m*p
        //return M1;
        if(M1.dim[1]!=M2.dim[0]){ cout<<"NOT Allowed To multiple Them\n";throw "Multiple Error"; }
        vector<vector<double>>Result(M1.dim[0],vector<double>(M2.dim[1],0));
        for (int i = 0; i < M1.dim[0]; i++) 
        for (int j = 0; j < M2.dim[1]; j++) 
        for (int k = 0; k < M2.dim[0]; k++) 
        Result[i][j]+=M1.Matrix[i][k]*M2.Matrix[k][j];
        return Result;   
    }
    
    //5-OPERATION on Matrix ->RREF:
    static Matrix2D SwapRows(int R1,int R2,vector<vector<double>>&Matrix ){
         swap(Matrix[R1],Matrix[R2]);
         return Matrix;
    }
    static Matrix2D AddRow(int R2,int R1,double C,vector<vector<double>>&Matrix){
        for (int i = 0; i < Matrix[0].size(); i++) {
            Matrix[R2][i]+=Matrix[R1][i]*C;
        }
        return Matrix;
    }
    static Matrix2D MultipleRow(int R,double C,vector<vector<double>>&Matrix){
        for (int i = 0; i < Matrix[0].size(); i++) {
            //cout<<Matrix[R][i];
            if(Matrix[R][i]!=0)
            Matrix[R][i]*=C;
            //cout<<Matrix[R][i];
        }
         return Matrix;
    }
    
    static Matrix2D SegmentedMatrix(const Matrix2D &M ,const Matrix2D &M2){
    if(M.dim[0]!=M2.dim[0]){cout<<"unAllowed to Segmented them";throw "SegmentedMatrix Error";}
    vector<vector<double>> SegmentedM(M.dim[0],vector<double>(M.dim[1]+M2.dim[1]));
    
    for (int r = 0; r < M.dim[0]; r++) {
        for (int c = 0; c <M.dim[1]; c++) {
            SegmentedM[r][c]=M.Matrix[r][c];
        }
        for (int c = 0; c <M2.dim[1]; c++) {
            SegmentedM[r][c+M2.dim[1]]=M2.Matrix[r][c];
        }
    }
    return SegmentedM;
    }
    
    // 
    static void print(vector<vector<double>>M){
        for (auto i : M) {
            for (auto j : i) {
                cout<<" "<<j<<" ";
            }
            cout<<endl;
        }
    }
    
    
    Matrix2D RREF(){
        if(_RREF.empty()){
            cout<<"Iam RREF empty\n";
            _RREF=Matrix;
            for (int i = 0; i<*min_element(dim, dim + 2);i++){ 
                Matrix2D::print(_RREF);
                cout<<"\n Phase 1 Init\n";
                if(_RREF[i][i]==0)
                    for (int j = i+1; j < dim[0]; j++) 
                        if(_RREF[j][i]!=0){Matrix2D::SwapRows(i,j,_RREF);break;};
                        
                Matrix2D::print(_RREF);
                cout<<"\nPhase 2 swap \n";
                
                if(_RREF[i][i]!=0){
                    //cout<<"i :"<<i<<" ELE "<<1/_RREF[i][i]<<" "<<_RREF[i][i+1]*(1/_RREF[i][i])<<"\n ";
                    Matrix2D::MultipleRow(i,1/_RREF[i][i],_RREF);
                    Matrix2D::print(_RREF);
                    cout<<"\n Phase 3 Multiplication\n";
                    
                    for (int j = 0; j < dim[0]; j++) {
                        cout<<i<<" :"<<j<< "  "<< (_RREF[i][j]!=0 && i!=j) <<" "<<_RREF[i][j] ;
                        if(_RREF[j][i]!=0 && i!=j ){Matrix2D::AddRow(j,i,-1*(_RREF[j][i]/_RREF[i][i]),_RREF);cout<<" X ";}
                        cout<<endl;
                    }
                    Matrix2D::print(_RREF);
                    cout<<"\n Phase 4 Elimanation\n";
                }        
            }            
         };
        return _RREF;
        
    }
    
    
    
    
};




int main() 
{
    int arr[2]={8,8};
    vector<vector<double>>M={
        {1,5,0},
        {3,9,2},
        {9,4,5}
    };
     vector<vector<double>>M3={
        {1,1,1},
        {2,3,2},
        {1,2,2}
        
    };
    vector<vector<double>>M2={
        {2,0,0,0},
        {0,2,0,0},
        {0,0,2,0},
        {0,0,0,2}
    };
    
    Matrix2D m(M);
    Matrix2D m2(M2);
    Matrix2D m3(M3);
    // m.display();
    //m2.display();
   // m2.Transpose().display();
    //m2.Transpose();
    cout << "Multiple!\n";
   // Matrix2D::MultipleMatrices(m2,m);
    cout<<"Det\n";
    //cout<<Matrix2D::Det(m3);
    
    //Matrix2D::Cofactor(m).display();
    cout<<Matrix2D::Det(m);
    cout<<"\nRReF\n";
    //m.RREF().display();
    
    //m2.RREF();
    cout << "SegmentedMatrix" <<endl;
    Matrix2D::SegmentedMatrix(M3,Matrix2D::IDENTITY(3)).RREF().display();
    cout << "adjoint" <<endl;
    Matrix2D::Adjoint(m3).display();
    cout<<"End Program\n";
    return 0;
}
