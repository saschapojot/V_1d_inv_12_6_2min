//
// Created by polya on 8/16/24.
//
#include "potentialFunctionPrototype.hpp"

class V_1d_inv_12_6_2min:public potentialFunction
{
public:
    V_1d_inv_12_6_2min(const std::string &coefsStr):potentialFunction(){
        this->coefsInStr=coefsStr;
    }

    void json2Coefs(const std::string &coefsStr)override{
        std::stringstream iss;
        iss<<coefsStr;
        std::string temp;
        //read a
        if (std::getline(iss, temp, ',')){
            this->a=std::stod(temp);
        }

        //read b

        if (std::getline(iss, temp, ',')){
            this->b=std::stod(temp);
        }

        //read A1
        if (std::getline(iss, temp, ',')){
            this->A1=std::stod(temp);
        }

        //read A2
        if (std::getline(iss, temp, ',')){
            this->A2=std::stod(temp);
        }

        //read a1
        if (std::getline(iss, temp, ',')){
            this->a1=std::stod(temp);
        }

        //read a2
        if (std::getline(iss, temp, ',')){
            this->a2=std::stod(temp);
        }

        //read sigma1
        if (std::getline(iss, temp, ',')){
            this->sigma1=std::stod(temp);
        }

        //read sigma2
        if (std::getline(iss, temp, ',')){
            this->sigma2=std::stod(temp);
        }

        //read f

        if (std::getline(iss, temp, ',')){
            this->f=std::stod(temp);
        }

        //read d2

        if (std::getline(iss, temp, ',')){
            this->d2=std::stod(temp);
        }


        //read N

        if (std::getline(iss, temp, ',')){
            this->N=std::stoi(temp);
        }
    }

    void init() override{
        this->json2Coefs(coefsInStr);
        this->r1=2.0*d1;
        this->r2=1.1*d2;
        this->lm=(static_cast<double >(N)*(r1+r2))*2;
        this->eps=((r1+r2)/2.0)/8;
        //        pow_result_tmp=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

        std::cout << "a=" << a << ", b=" << b << ", A1=" << A1 << ", A2=" << A2
        <<", a1="<<a1<<", a2="<<a2
        <<", sigma1="<<sigma1<<", sigma2="<<sigma2
        <<", f="<<f<<", d2="<<d2<< std::endl;
        std::cout<<"r1="<<r1<<", r2="<<r2<<std::endl;
        std::cout<<"lm="<<lm<<std::endl;
        //        std::cout<<"eps="<<eps<<std::endl;


    }


    double operator() (const double * xVec, const int& j)override{
        //only computes the part of the potential that changes
        if(j<0 or j>=2*N){
            std::cerr<<"j="+std::to_string(j)+" out of range";
            std::exit(14);
        }
        if(j==0){
            return V1(xVec[j+1]-xVec[j]);
        }//end 0
        else if(j==2*N-1){
            return V1(xVec[j]-xVec[j-1]);
        }//end 2N-1
        else{
            if(j%2==1){
                //B, odd position
                double d1=xVec[j]-xVec[j-1];
                double d2=xVec[j+1]-xVec[j];
                return V1(d1)+V2(d2);
            }//end odd
            else{
                //A, even position
                double d2=xVec[j]-xVec[j-1];
                double d1=xVec[j+1]-xVec[j];
                return V2(d2)+V1(d1);
            }//end even



        }//end middle

    }//end () operator

    double potentialFull(const double * xVec)override
    {
        double val=0;

        for (int j=0;j<2*N;j+=2)
        {
            double d1=xVec[j+1]-xVec[j];
            val+=V1(d1);
        }
// std::cout<<"val1="<<val<<std::endl;
        for(int j=1;j<2*N-1;j+=2)
        {
            double d2=xVec[j+1]-xVec[j];
            val+=V2(d2);
        }

        return val;
    }
    double V1(const double &r){
        double val1=a*std::pow(r,-12.0)-b*std::pow(r,-6.0);

        double val2=-A1*std::exp(-1.0/std::pow(sigma1,2.0)*std::pow(r-a1,2.0))
                    -A2*std::exp(-1.0/std::pow(sigma2,2.0)*std::pow(r-a2,2.0));

        return val1+val2;
    }

    double V2(const double &r){
        double val=std::pow(r-d2,2.0);
        val*=f;
        return val;
    }

    double getLm() const override {
        return lm;
    }
    double get_eps() const override {
        return eps;
    }

public:
    double a;
    double b;
    double A1;
    double A2;
    double a1;
    double a2;
    double sigma1;
    double sigma2;
    double f;
    double d1;
    double d2;
    std::string coefsInStr;
    double r1;//min position of V1
    double r2;//min position of V2
    double lm;//range of distances
    double eps;//half interval length of uniform distribution
    int N;

};

std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &coefsJsonStr) {
    if (funcName == "V_1d_inv_12_6_2min") {

        return std::make_shared<V_1d_inv_12_6_2min>(coefsJsonStr);
    }

    else {
        throw std::invalid_argument("Unknown potential function type");
    }
}