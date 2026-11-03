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
    vector<vector<double>>_Inverse;
    public:
    int dim[2];
    vector<vector<double>> Matrix;
    
    
    //initialize 
     Matrix2D(int dims[2]){
        copy(dims,dims+2,dim);
        Matrix.resize(dims[0], vector<double>(dims[1], 0));
        cout<<Matrix[5][5]<<endl;
    };
    
    
    Matrix2D(const vector<vector<double>> &M){
        Matrix=M;
        dim[0]=M.size();
        dim[1]=M[0].size();
        
    };
    
    void display(){
        cout<<"   ";
         for (int i = 0; i <dim[1]; i++) {
             cout<<"C"<<i+1<<" ";
         }
         cout<<endl;
         for (int i = 0; i<dim[0];i++) {
              cout<<"R"<<i+1;
              
             for (int j = 0; j<dim[1];j++) {
                cout<<" "<<Matrix[i][j];
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
    Matrix2D TransposeMatix(){
        if(T.empty()){
            cout<<"Iam empty\n";
            T.resize(dim[1],vector<double>(dim[0]));
            for (int i = 0; i<dim[0];i++) for (int j = 0; j<dim[1];j++)  T[j][i]=Matrix[i][j];
             
         };
        return T;
    }
    static vector<vector<double>> Transpose(const vector<vector<double>> &M){
        vector<vector<double>>T(M[0].size(),vector<double>(M.size(),0));
        for (int i = 0; i<M.size();i++) for (int j = 0; j<M[0].size();j++)  T[j][i]=M[i][j];
        return T;
    }
    
    vector<vector<double>> Cofactor(){
        if(! _Cofactor.empty()) return _Cofactor;
        vector<int> PassedCoulmn(dim[1],0);
        vector<vector<double>>CofactorMatrix(dim[0],vector<double>(dim[1],0));
        vector<int> PassedRow(dim[0],0) ;
        for (int r = 0; r < dim[0]; r++) {
            for (int c = 0; c < dim[1]; c++) {
                PassedRow[r]=PassedCoulmn[c]=1;
                CofactorMatrix[r][c]= NDet(PassedCoulmn,PassedRow)*pow(-1,r-c);
                PassedCoulmn[c]=PassedRow[r]=0;
            }
        }
        return CofactorMatrix;
    }

    Matrix2D Adjoint(){return Matrix2D::Transpose(Cofactor());}
    
    vector<vector<double>> InveresMatrix(){
        if(!_Inverse.empty())return _Inverse;
        _Inverse.resize(dim[0],vector<double>(dim[0]));
        Matrix2D Segmented = Matrix2D::SegmentedMatrix(*this, Matrix2D::IDENTITY(dim[0])).RREF();
        for (int r = 0; r < dim[0]; r++) {
            for (int c = 0; c < dim[0]; c++) {
                _Inverse[r][c]=Segmented.Matrix[r][c+dim[0]];
            }
        }
        return _Inverse;
    }
    double Det(vector<int> memo={},int cRow=0){
        
        if(memo.empty())memo.resize(dim[0],0);
        double Result=0;
        int sign=-1;
        for (int i = 0; i < dim[0]; i++) {
            if(memo[i]==0){
                sign*=-1;
                memo[i]=1;
                 if(cRow==dim[0]-1){
                    //cout <<  "End Det ("<<cRow<<","<<i<<") => "<<Matrix[cRow][i]<<" r " <<Result<< endl;
                    return Matrix[cRow][i];
                }
                
                
                Result+=Matrix[cRow][i]*Det(memo,cRow+1)*sign;
                //cout<<"I Det ("<<cRow<<","<<i<<") => "<< Matrix[cRow][i]*pow(-1,i)<<" r "<<Result<<endl;
                memo[i]=0;
            }
            
        }
        return Result;
    }
    
    double NDet(vector<int>PassedCoulmn={},vector<int> PassedRow={},int cRow=0){
        if(PassedCoulmn.empty())PassedCoulmn.resize(dim[0],0);
        if(PassedRow.empty())PassedRow.resize(dim[0],0);
        double Result=0;
        int sign=-1;
        if(PassedRow[cRow])++cRow;
        if(cRow==dim[0]){
                    // cout <<  "End Det ("<<cRow<<","<<") => "<<Result<< endl;
                    return 1;
        }
        for (int i = 0; i < dim[0]; i++) {
            if(PassedCoulmn[i]==0){
                sign*=-1;
                 
                PassedCoulmn[i]=PassedRow[cRow]=1;
                Result+=Matrix[cRow][i]*NDet(PassedCoulmn,PassedRow,cRow+1)*sign;
                // cout<<"I Det ("<<cRow<<","<<i<<") => "<< Matrix[cRow][i]*pow(-1,i)<<" r "<<Result<<endl;
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
    static Matrix2D SwapRows(int R1,int R2,vector<vector<double>>&Matrix){
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
        for (int i = 0; i < Matrix.size(); i++) {
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
            SegmentedM[r][c+M.dim[1]]=M2.Matrix[r][c];
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
                        if(_RREF[j][i]!=0){ Matrix2D::SwapRows(i,j,_RREF);break;};
                        
                Matrix2D::print(_RREF);
                cout<<"\nPhase 2 swap \n";
                
                if(_RREF[i][i]!=0){
                    //cout<<"i :"<<i<<" ELE "<<1/_RREF[i][i]<<" "<<_RREF[i][i+1]*(1/_RREF[i][i])<<"\n ";
                    Matrix2D::MultipleRow(i,1/_RREF[i][i],_RREF);
                    Matrix2D::print(_RREF);
                    cout<<"\n Phase 3 Multiplication\n";
                    
                    for (int j = 0; j < dim[0]; j++) {
                        // cout<<i<<" :"<<j<< "  "<< (_RREF[i][j]!=0 && i!=j) <<" "<<_RREF[i][j] ;
                        if(_RREF[j][i]!=0 && i!=j ){ Matrix2D::AddRow(j,i,-1*(_RREF[j][i]/_RREF[i][i]),_RREF);}
                        
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
    vector<vector<double>>M2={
        {2,0,0,0},
        {0,2,0,0},
        {0,0,2,0},
        {0,0,0,2}
    };
     vector<vector<double>>M3={
        {1,1,1},
        {2,3,2},
        {1,2,2}
        
    };
    Matrix2D m(M);
    Matrix2D m2(M2);
    Matrix2D m3(M3);
    cout<<"M1 :\n";
    m.display();
    cout<<"M2 :\n";
    m2.display();
    cout<<"M3 :\n";
    m3.display();
    cout<<"Transpose M1 :\n";
    m.TransposeMatix().display();
    cout << "Multiple M3*M\n";
    Matrix2D::MultipleMatrices(m3,m).display();
    cout << "Multiple M*M3\n";
    Matrix2D::MultipleMatrices(m,m3).display();
    cout<<"Det M3 :";
    cout<< m3.Det()<<endl;
    cout<< "Cofactor M\n";
    Matrix2D::print(m.Cofactor());
    // cout<<Det(m);
    cout<<"\nRReF\n";
    m.RREF().display();
    cout << "SegmentedMatrix (M3|I)" <<endl;
    Matrix2D::SegmentedMatrix(M3,Matrix2D::IDENTITY(3)).RREF().display();
    cout << "Adjoint M3" <<endl;
     m3.Adjoint().display();
    cout<<"_Inverse M3\n";
    Matrix2D::print(m3.InveresMatrix());
    cout<<"End Program\n";
    return 0;
}
